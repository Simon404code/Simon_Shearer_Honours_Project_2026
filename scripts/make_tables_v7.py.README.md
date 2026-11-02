# make_tables_v7.py

## What it does
Generates summary tables for tandem-array HOGs from a gene-level CSV and a `combBed` file. It validates input columns, computes HOG-level and array-level summaries, and writes a set of CSV outputs with names derived from an `output_prefix`.

## Inputs
Positional arguments:
- `table1`: a CSV of all genes belonging to tandem-array HOGs.
- `combbed`: `combBed` TSV or TSV.GZ with gene-level array annotations.
- `output_prefix` (optional, default: `tandem`): filename prefix for generated outputs.
- `--n-genomes` (optional, default: `50`): expected genome count used for warnings.

Expected data contracts:
- `table1` should include at least: `HOG_ID`, `Gene_ID`, `Genome`, `array_ID`, `ofID`, `isArrayRep`, `chr`, `start`, `end`.
- `combBed` is expected to contain gene IDs and array metadata in a gene-centric form, with aliases recognized for `Gene_ID`, `HOG_ID`, `Genome`, `arrayID`, `chr`, `start`, and `end`.

## Outputs
Creates files named like:
- `<output_prefix>_table2.csv`
- `<output_prefix>_table3.csv`
- `<output_prefix>_table4.csv`

The script also writes a debug merge file:
- `<output_prefix>_merged_debug.csv`

These outputs include aggregated HOG summary statistics and tandem-array summaries.

## Usage
`python make_tables_v7.py <table1.csv> <combBed.tsv> [output_prefix] [--n-genomes 50]`

## Dependencies
- `pandas`
- Python standard library: `argparse`, `sys`

## Important processing steps / limitations
- Strict validation is performed: duplicate `Gene_ID`s, missing required columns, and invalid numeric coordinates raise errors.
- Table 2 classifies HOGs as `singleton`, `minority shell`, `majority shell`, `soft-core`, or `core` based on the number of genomes carrying the HOG.
- Table 3 classifies HOGs as `Core`, `Soft-core`, or `Shell` based on the number of genomes with tandem arrays.
- It assumes a maximum of 50 genomes per HOG and warns if the data size differs from the expected count.
- It is built around specific column aliases and may fail on schema deviations.