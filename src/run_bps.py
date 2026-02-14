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
    parser.add_argument(
        "--feature_importance", action="store_true",
        help="Compute and display feature importance (RF) and normalized coefficients (LR)"
    )
    parser.add_argument(
        "--top_n", type=int, default=20,
        help="Number of top features to display (default: 20)"
    )
    parser.add_argument(
        "--output_csv", type=str, default=None,
        help="Path to save the full feature importance DataFrame as CSV"
    )

    args = parser.parse_args()

    print("Loading data...")
    X, batches, genes = load_data(args.input, batch_key=args.batch_key)

    print("Running BPS calculation...")
    result = calculate_bps(X, batches, n_jobs=args.n_jobs, compute_importance=args.feature_importance)

    if not args.feature_importance:
        BPS_RF, BPS_LR = result
        print("\nResults:")
        print(f"BPS-RF: {BPS_RF:.4f}")
        print(f"BPS-LR: {BPS_LR:.4f} (linear batch signal)")
    else:
        import pandas as pd

        BPS_RF = result["BPS_RF"]
        BPS_LR = result["BPS_LR"]

        print("\nResults:")
        print(f"BPS-RF: {BPS_RF:.4f}")
        print(f"BPS-LR: {BPS_LR:.4f} (linear batch signal)")

        # Build importance DataFrame
        gene_labels = genes if genes is not None else [f"feature_{i}" for i in range(X.shape[1])]
        importance_df = pd.DataFrame({
            "gene": gene_labels,
            "rf_importance": result["rf_importance"],
            "lr_coef_norm": result["lr_importance"],
        }).sort_values("rf_importance", ascending=False).reset_index(drop=True)

        # Display top N features
        top_n = min(args.top_n, len(importance_df))
        print(f"\nTop {top_n} features by RF importance:")
        print(importance_df.head(top_n).to_string(index=False))

        # Save full DataFrame to CSV if requested
        if args.output_csv:
            importance_df.to_csv(args.output_csv, index=False)
            print(f"\nFull feature importance saved to {args.output_csv}")

if __name__ == "__main__":
    main()
