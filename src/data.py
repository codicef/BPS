#!/usr/bin/env python3

import pandas as pd
import scanpy as sc
import numpy as np



def parse_anndata(adata, batch_key='batch'):
    """
    Parse an AnnData object to extract relevant information for downstream analysis.

    Parameters:
    - adata: AnnData object containing the dataset.
    - batch_key: Key in adata.obs that contains the batch information.

    Returns:
    - X: Feature matrix (gene expression data).
    - batches: Labels (batch information).
    - genes: List of gene names.
    """
    # Extract the feature matrix
    X = adata.X
    if X is None:
        print(f"adata.X is None, checking if adata.raw is available...")
        if adata.raw is not None:
            print(f"Using adata.raw.X...")
            X = adata.raw.X
        else:
            print(f"adata.raw is None, checking if layer is available...")
            if 'counts' in adata.layers:
                print(f"Using adata.layers['counts']...")
                X = adata.layers['counts']
            elif 'normalized' in adata.layers:
                print(f"Using adata.layers['normalized']...")
                X = adata.layers['normalized']
            else:
                raise ValueError("No suitable data found in adata. Please check the AnnData object.")

    # Check if X is a sparse matrix and convert to dense if necessary
    if hasattr(X, 'toarray'):
        X = X.toarray()
    elif hasattr(X, 'todense'):
        X = X.todense()
    elif not isinstance(X, np.ndarray):
        raise ValueError("X is not a valid matrix format. Please check the AnnData object.")



    # Check batch_key presence
    if batch_key not in adata.obs.columns and 'tech' not in adata.obs.columns:
        print(f"Available keys in adata.obs: {adata.obs.columns.tolist()}")
        raise KeyError(f"Batch key '{batch_key}' not found in adata.obs.")
    elif 'tech' in adata.obs.columns:
        print(f"Warning: 'tech' found in adata.obs, using it as batch_key instead of '{batch_key}'")
        batch_key = 'tech'

    # Extract the labels (batch information)
    batches = adata.obs[batch_key].values

    # Extract gene names
    genes = adata.var_names.tolist()
    return X, batches, genes


def parse_dataframe(df, batch_key='batch', gene_key='gene'):
    """
    Parse a DataFrame to extract relevant information for downstream analysis.

    Parameters:
    - df: DataFrame containing the dataset.
    - batch_key: Key in df that contains the batch information.
    - gene_key: Key in df that contains the gene names.

    Returns:
    - X: Feature matrix (gene expression data).
    - batches: Labels (batch information).
    - genes: List of gene names.
    """


    if gene_key in df.columns:
        # Extract gene names
        genes = df[gene_key].tolist()
        df = df.drop(columns=[gene_key])
    else:
        genes = None

    # Extract the feature matrix
    X = df.drop(columns=[batch_key]).values

    # Check batch_key presence
    if batch_key not in df.columns:
        raise KeyError(f"Batch key '{batch_key}' not found in DataFrame.")
    # Extract the labels (batch information)
    batches = df[batch_key].values

    return X, batches, genes


def load_data(file_path, batch_key='batch', gene_key='gene'):
    '''
    Load data from a file and return the feature matrix and batch labels.

    Parameters:
    - file_path: Path to the input data file.

    Returns:
    - X: Feature matrix (gene expression data).
    - batches: Labels (batch information).
    - genes: List of gene names (if available).

    '''
    if file_path.endswith('.h5ad'):
        # If the file is in .h5ad format
        print(f"Loading adata from {file_path}...")
        adata = sc.read(file_path)
        X, batches, genes = parse_anndata(adata)
    elif file_path.endswith('.csv') or file_path.endswith('.tsv'):
        print(f"Loading dataframe from {file_path}...")
        # If the file is in .csv or .tsv format
        sep = ',' if file_path.endswith('.csv') else '\t'
        df = pd.read_csv(file_path, sep=sep)
        X, batches, genes = parse_dataframe(df, batch_key=batch_key, gene_key=gene_key)
    else:
        raise ValueError("Unsupported file format. Please provide a .h5ad, .csv, or .tsv file.")
    print(f"Loaded {file_path} with shape {X.shape} and {len(np.unique(batches))} batches.")

    return X, batches, genes
