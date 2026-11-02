# group_by_tandem_array.py

## What it does
Identifies tandem-array HOGs from a `combBed`-style file and writes two outputs: a list of tandem-array HOG IDs and a detailed summary of the tandem arrays and representative genes for each array. The script is a direct workflow script rather than a library.

## Inputs
The script reads four positional arguments from `argv`:
1. `combbed_file`: a `combBed` TSV with gene-level rows.
2. `provenance_file`: a CSV file path, though the script also defaults to a local `Pangenome_provenaces_edited.csv` beside the script when it is used in the repository context.
3. `GENOME_MAP_FILE` is effectively part of the caller's environment but not exposed as an argparse option; the script is written to read it as `argv[3]`.
4. `OUTPUT_FILE` is the final CSV path, again provided as `argv[4]`.

Required columns in the combed data:
- `chr`, `start`, `end`, `id`, `ofID`, `genome`, `arrayID`, `globHOG`

Provenance file contract:
- CSV with `ProjectName`, `Location`, `Region`, `Subpopulation`

## Outputs
- `tandem_array_HOGs.txt`: one HOG ID per line.
- `tandem_array_HOG_summary.csv`: a detailed summary of tandem arrays.
- Final output file as provided in `argv[4]` (the script comments say the final output is `.csv.gz`, but the implemented code writes a CSV-like file using pandas and no compression logic is visible).

## Usage
This script is command-line driven by positional arguments, e.g.:
`python group_by_tandem_array.py <combBed.tsv> <provenance.csv> <genome_map.txt> <output.csv>`

The actual argument structure is inferred from the source rather than from an `argparse` parser, so the script should be treated as a hard-coded workflow script.

## Dependencies
- `pandas`
- Python standard library: `pathlib`, `sys`

## Important processing steps / limitations
- It groups rows by `HOG_ID`, `genome`, and `arrayID` and keeps only groups with at least two genes.
- It verifies that no tandem array spans multiple chromosomes and that every genome has provenance metadata.
- It writes arrays sorted by `HOG_ID`, `genome`, `arrayID`, and start coordinate.
- It uses hard-coded output filenames (`tandem_array_HOGs.txt` and `tandem_array_HOG_summary.csv`) in the current working directory.
- It is not robust to schema deviations or alternative naming conventions beyond the exact required columns.