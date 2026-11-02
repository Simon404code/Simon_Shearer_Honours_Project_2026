# make_matrix_count.py

## What it does
Builds a HOG-by-genome count matrix from a TSV file. For each HOG and genome, it counts the number of distinct genes in that HOG and writes a CSV matrix with rows as HOG IDs and columns as genomes.

## Inputs
- `--data`: tab-delimited TSV file.
- `--output`: output CSV path.

Required data contract:
- Input rows must include `id`, `genome`, and `globHOG`.
- `globHOG` is expected to look like `N0.HOG0008484 OG0000263`; only the first token is used as the HOG ID.

## Outputs
CSV with header:
- `HOG_ID` followed by one column per genome present in the input.

Each cell contains the count of distinct genes from that HOG in that genome.

## Usage
`python make_matrix_count.py --data <input.tsv> --output <matrix.csv>`

## Dependencies
- Python standard library: `csv`, `argparse`, `collections.defaultdict`.

## Important processing steps / limitations
- The script skips rows with blank `globHOG` values.
- Duplicate gene IDs within the same HOG/genome pair are deduplicated with a set, so counts reflect unique genes, not raw rows.
- Genome names and HOG IDs are sorted for consistent output.
- The script is matrix-oriented and does not carry gene IDs in the output values.