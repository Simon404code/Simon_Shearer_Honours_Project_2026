# classify_core_shell_HOGs.py

## What it does
Reads a HOG presence/absence matrix and classifies each HOG as `Core`, `Soft-core`, or `Shell` based on how many haplotypes are present.

## Inputs
- `--input`: a CSV matrix with the HOG ID in the first column and one column per haplotype.
- `--output`: output CSV path.

Expected matrix contract:
- The first column is the HOG ID.
- Remaining columns represent haplotypes.
- A value greater than zero indicates presence in that haplotype.

## Outputs
A CSV with columns:
- `HOG_ID`
- `Category`
- `Haplotypes_present`

## Usage
`python classify_core_shell_HOGs.py --input <matrix.csv> --output <classified.csv>`

## Dependencies
- Python standard library: `csv`, `argparse`

## Important processing steps / limitations
- If `present_count == haplotype_count`, the HOG is labeled `Core`.
- If `present_count == haplotype_count - 1`, it is labeled `Soft-core`.
- Otherwise, it is labeled `Shell`.
- The script assumes all matrix values can be converted to `float` via `float(value)`, so non-numeric placeholders may fail.