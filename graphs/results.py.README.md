# results.py

## What it does
This is not a reusable pipeline script; it is a one-off analysis snippet that loads a saved `combBed` file and two output tables, then prints counts of tandem arrays by category. It appears to have been used to inspect the distribution of tandem arrays in core, soft-core, and singleton classes.

## Inputs
The script hard-codes absolute file paths to the project data directory:
- `combBed_v5_complete-2.txt`
- `out_table2.5.csv`
- `out_table3.5.csv`

It expects pandas-readable CSV/TSV data with at least the columns used by the selected indices and filters.

## Outputs
Console output only, including counts for:
- total tandem arrays
- core tandem arrays
- soft-core tandem arrays
- core & soft-core tandem arrays
- singleton tandem array count

It also computes a `tandem_array_distribution` list but does not write it to disk.

## Usage
This script is effectively not CLI-driven. It is intended to be run from a Python environment where the hard-coded paths are valid.

## Dependencies
- `pandas`

## Important processing steps / limitations
- It uses hard-coded Windows paths and is therefore not portable without editing.
- It sets index columns (`HOG_ID`, `HOG_ID(array)`) and then inspects a specific known HOG (`N0.HOG0000001 OG0000000`) by label.
- It is not a general-purpose analysis tool and works only for the dataset layout assumed by the author.
- It does not validate file existence or column names before use; it expects the data files to be available in the exact paths shown.