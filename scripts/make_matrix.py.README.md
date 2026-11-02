# make_matrix.py

## What it does
Creates a gene-presence matrix by HOG and haplotype. Unlike the count version, each matrix cell stores the gene IDs for a given HOG/genome combination as a semicolon-delimited string instead of an integer count.

## Inputs
- `--data`: input TSV file.
- `--output`: output CSV path.

Required data contract:
- The input needs columns `id`, `genome`, and `globHOG`.
- `globHOG` should contain a HOG ID followed by an orthogroup token; the script extracts only the first token.

## Outputs
CSV with header:
- `HOG_ID` followed by one column per genome in the input.

Each cell contains a semicolon-separated list of gene IDs present in that HOG/genome combination.

## Usage
`python make_matrix.py --data <input.tsv> --output <matrix.csv>`

## Dependencies
- Python standard library: `csv`, `argparse`, `collections.defaultdict`.

## Important processing steps / limitations
- Gene IDs are deduplicated within a single cell to avoid repeated insertion.
- Rows missing a `globHOG` are skipped.
- It prints a short summary of total haplotypes and HOG IDs at runtime.
- The output cells are string-based; downstream consumers must parse the semicolon-delimited gene ID lists.