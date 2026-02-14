#!/usr/bin/env python3

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
import anndata
import pandas as pd
from data import parse_anndata, parse_dataframe


def safe_truncate(X, lower_percentile=1, upper_percentile=99):
    '''
    Truncates the values in X to the specified percentiles.
    '''
    lower_bound = np.percentile(X, lower_percentile, axis=0)
    upper_bound = np.percentile(X, upper_percentile, axis=0)
    X_clipped = np.clip(X, lower_bound, upper_bound)
    return X_clipped

def calculate_bps(X, batches, default_preprocess=True, n_jobs=-1, compute_importance=False):
    '''
    Computes BPS-RF and BPS-LR metrics for the given dataset.

    Parameters:
    - X: Feature matrix (gene expression data).
    - batches: Labels (batch information).
    - default_preprocess: Whether to apply default preprocessing.
    - n_jobs: Number of parallel jobs.
    - compute_importance: If True, fit models on full data and return feature importances.

    Returns:
    - (BPS_RF, BPS_LR) when compute_importance=False
    - {"BPS_RF": float, "BPS_LR": float, "rf_importance": array, "lr_importance": array}
      when compute_importance=True
    '''

    if default_preprocess:
        print("Preprocessing data...")
        X = StandardScaler().fit_transform(X)


    rf = RandomForestClassifier(n_estimators=200, max_depth=100, class_weight="balanced_subsample",
                                n_jobs=n_jobs, random_state=42)
    lr = LogisticRegression(max_iter=10000, C=1.0, n_jobs=n_jobs, class_weight="balanced")

    # Compute BPS-RF (Random Forest)
    print("Computing BPS-RF...")
    rf_scores = cross_val_score(rf, X, batches, cv=5, scoring='roc_auc_ovr')  # Calcoliamo ROC AUC
    BPS_RF = (np.mean(rf_scores) - 0.5) / 0.5

    if default_preprocess:
        X = safe_truncate(X, lower_percentile=0.1, upper_percentile=99)
    # Compute BPS-LR (Logistic Regression)
    print("Computing BPS-LR...")
    lr_scores = cross_val_score(lr, X, batches, cv=5, scoring='roc_auc_ovr')
    BPS_LR = (np.mean(lr_scores) - 0.5) / 0.5

    if not compute_importance:
        return BPS_RF, BPS_LR

    # Fit models on full dataset to extract feature importances
    print("Computing feature importances...")
    rf.fit(X, batches)
    rf_importance = rf.feature_importances_

    lr.fit(X, batches)
    # Mean absolute coefficient across all classes (OVR), then normalize to sum to 1
    lr_abs_coef = np.mean(np.abs(lr.coef_), axis=0)
    lr_importance = lr_abs_coef / lr_abs_coef.sum()

    return {
        "BPS_RF": BPS_RF,
        "BPS_LR": BPS_LR,
        "rf_importance": rf_importance,
        "lr_importance": lr_importance,
    }




def bps(input_data, batch_key='batch', n_jobs=-1, default_preprocess=True, verbose=True, compute_importance=False):
    '''
    Computes BPS-RF and BPS-LR metrics for the given dataset.

    Parameters:
    - input_data: anndata object or a pandas dataframe.
    - batch_key: Key in input_data that contains the batch information.
    - n_jobs: Number of jobs to run in parallel.
    - default_preprocess: Whether to apply default preprocessing (standardization and truncation).
    - compute_importance: If True, also compute and return feature importances.

    Returns:
    - (BPS_RF, BPS_LR) when compute_importance=False
    - (BPS_RF, BPS_LR, importance_df) when compute_importance=True,
      where importance_df is a DataFrame with columns: gene, rf_importance, lr_coef_norm
    '''

    if isinstance(input_data, anndata.AnnData):
        # If the input data is an AnnData object
        adata = input_data
        X, batches, genes = parse_anndata(adata, batch_key=batch_key)
    elif isinstance(input_data, pd.DataFrame):
        # If the input data is a DataFrame
        df = input_data
        X, batches, genes = parse_dataframe(df, batch_key=batch_key)
    else:
        raise ValueError("Unsupported input data type. Please provide an AnnData object or a pandas DataFrame.")

    if verbose:
        print(f"Loaded data with shape {X.shape} and {len(np.unique(batches))} batches.")

    # Compute BPS scores
    result = calculate_bps(X, batches, default_preprocess=default_preprocess, n_jobs=n_jobs,
                           compute_importance=compute_importance)

    if not compute_importance:
        BPS_RF, BPS_LR = result
        if verbose:
            print(f"BPS-RF: {BPS_RF}")
            print(f"BPS-LR: {BPS_LR}")
        return BPS_RF, BPS_LR

    BPS_RF = result["BPS_RF"]
    BPS_LR = result["BPS_LR"]

    if verbose:
        print(f"BPS-RF: {BPS_RF}")
        print(f"BPS-LR: {BPS_LR}")

    # Build importance DataFrame
    gene_labels = genes if genes is not None else [f"feature_{i}" for i in range(X.shape[1])]
    importance_df = pd.DataFrame({
        "gene": gene_labels,
        "rf_importance": result["rf_importance"],
        "lr_coef_norm": result["lr_importance"],
    }).sort_values("rf_importance", ascending=False).reset_index(drop=True)

    return BPS_RF, BPS_LR, importance_df
