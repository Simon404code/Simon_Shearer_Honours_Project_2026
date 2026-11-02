# single_vs_multiple_HOGs.py

## What it does
Classifies tandem arrays into two categories based on how many distinct global HOG IDs they contain: `single-HOG` or `multi-HOG`. It groups rows by `(genome, arrayID)`, collects the HOG IDs from the `globHOG` field, and writes a summary table for each array.

## Inputs
- `--data`: input TSV file.
- Required columns in the data include at least `genome`, `arrayID`, `id`, and `globHOG`.
- `globHOG` is expected to look like `N0.HOG0008484 OG0000263`; the first token is treated as the HOG identifier.

## Outputs
A TSV file with columns:
- `genome`
- `arrayID`
- `gene_count`
- `HOG_count`
- `array_type`
- `globHOG_ids`
- `gene_ids`

Rows with a single gene in an array are skipped.

## Usage
`python single_vs_multiple_HOGs.py --data <input.tsv> --output <summary.tsv>`

## Dependencies
- Python standard library: `csv`, `argparse`, `collections.defaultdict`.

## Important processing steps / limitations
- The script uses only the first token of `globHOG` (`glob_hog.split()[0]`) to represent a HOG, ignoring the orthogroup component after the first whitespace.
- It does not check for malformed rows beyond field presence in the CSV reader.
- It preserves gene IDs in a semicolon-delimited string and may not be suitable for large files if memory usage is a concern because all arrays are collected in memory.