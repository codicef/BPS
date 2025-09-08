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

def calculate_bps(X, batches, default_preprocess=True, n_jobs=-1):
    '''
    Computes BPS-RF and BPS-LR metrics for the given dataset.

    Parameters:
    - X: Feature matrix (gene expression data).
    - batches: Labels (batch information).

    Returns:
    - BPS_RF: BPS-RF score.
    - BPS_LR: BPS-LR score.
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

    return BPS_RF, BPS_LR




def bps(input_data, batch_key='batch', n_jobs=-1, default_preprocess=True, verbose=True):
    '''
    Computes BPS-RF and BPS-LR metrics for the given dataset.

    Parameters:
    - input_data: anndata object or a pandas dataframe.
    - batch_key: Key in input_data that contains the batch information.
    - n_jobs: Number of jobs to run in parallel.
    - default_preprocess: Whether to apply default preprocessing (standardization and truncation).

    Returns:
    - BPS_RF: BPS-RF score.
    - BPS_LR: BPS-LR score.
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
    BPS_RF, BPS_LR = calculate_bps(X, batches, default_preprocess=default_preprocess, n_jobs=n_jobs)

    if verbose:
        print(f"BPS-RF: {BPS_RF}")
        print(f"BPS-LR: {BPS_LR}")
    return BPS_RF, BPS_LR
