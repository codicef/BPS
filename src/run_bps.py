#!/usr/bin/env python3

import argparse
from data import load_data
from bps import calculate_bps

def main():
    parser = argparse.ArgumentParser(
        description="Compute BPS-RF and BPS-LR scores from single-cell gene expression data."
    )

    parser.add_argument(
        "input",
        help="Path to input file (.h5ad, .csv, or .tsv)"
    )
    parser.add_argument(
        "--batch_key", default="batch",
        help="Column name or obs key for batch labels (default: 'batch')"
    )
    parser.add_argument(
        "--n_jobs", type=int, default=1,
        help="Number of jobs to run in parallel for BPS calculation (default: 1)"
    )

    args = parser.parse_args()

    print("Loading data...")
    X, batches = load_data(args.input, batch_key=args.batch_key)

    print("Running BPS calculation...")
    BPS_RF, BPS_LR = calculate_bps(X, batches, n_jobs=args.n_jobs)

    print("\nResults:")
    print(f"BPS-RF: {BPS_RF:.4f}")
    print(f"BPS-LR: {BPS_LR:.4f} (linear batch signal)")

if __name__ == "__main__":
    main()
