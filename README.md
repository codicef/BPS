# BPS Calculator

**BPS Calculator** is a command-line tool for computing the **Batch Probing Score (BPS)**, a supervised metric designed to quantify residual batch effects in single-cell gene expression data. It provides two variants:

- `BPS-RF`: Batch signal detected via a Random Forest classifier  
- `BPS-LR`: Batch signal detected via a Logistic Regression classifier (linear decision boundary)

Both scores are based on the ROC AUC of a batch classifier, normalized between 0 and 1 as:

Higher BPS values indicate stronger residual batch signal after correction. A score near 0 implies minimal detectable batch signal, while values close to 1 indicate strong batch-specific structure.

## Input formats

The tool accepts:
- `.h5ad`: AnnData format (used in Scanpy workflows)
- `.csv` or `.tsv`: Tabular formats with one row per sample and one column indicating the batch

### Example `.csv` format:
| gene1 | gene2 | gene3 | ... | batch |
|-------|-------|-------|-----|--------|
|  2.3  |  1.0  |  0.7  | ... |   A    |
|  1.9  |  0.8  |  0.5  | ... |   B    |


## Usage

To compute the BPS-RF and BPS-LR scores, run the script with the following command:

```bash
python src/run_bps.py <input_file> [--batch_key <batch_column>] [--n_jobs <num_jobs>]
```

### Optional arguments

| Argument             | Description                                                         | Default   |
|----------------------|---------------------------------------------------------------------|-----------|
| `--batch_key <batch>`| Column name (CSV/TSV) or `.obs` key (h5ad) for batch labels         | `"batch"` |
| `--n_jobs <num_jobs>`| Number of jobs to run in parallel for the BPS calculation           | 1         |

### Requirements

To run this tool, you need the following Python libraries:

- `scanpy`
- `numpy`
- `scikit-learn`
- `pandas`

You can install the necessary dependencies using `pip`:

```bash
pip install -r requirements.txt
```


### Dataset example 
We can use the **human_pancreas_norm_complexBatch.h5ad** dataset as a test dataset; it contains five distinct batches (technologies/sources of cells) along with annotated cell types and raw count matrices, making it ideal for evaluating and benchmarking integration methods. Download it here:https://figshare.com/ndownloader/files/24539828 (linked from the LIGER benchmarking).

### Usage example

```bash
~/BPS (main) [1]> python3 src/run_bps.py ../sc_batch_correction/data/raw/dataset7/human_pancreas_5.h5ad
Loading data...
Loading adata from ../sc_batch_correction/data/raw/dataset7/human_pancreas_5.h5ad...
Warning: 'tech' found in adata.obs, using it as batch_key instead of 'batch'
Loaded ../sc_batch_correction/data/raw/dataset7/human_pancreas_5.h5ad with shape (16382, 19093) and 9 batches.
Running BPS calculation...
Preprocessing data...
Computing BPS-RF...
Computing BPS-LR...

Results:
BPS-RF: 0.9973
BPS-LR: 0.9999 (linear batch signal)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


