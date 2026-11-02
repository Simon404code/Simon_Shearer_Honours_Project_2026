# all_genes_in_tandem_HOGs.py

## What it does
Scans a `combBed` TSV and keeps only rows belonging to HOGs that contain at least one tandem array of length 2 or more within a single genome. It writes all genes from those HOGs to a reduced output CSV.

## Inputs
- `combbed`: a `combBed` TSV file.
- `-o/--output` (optional): output CSV path, default `all_genes_in_tandem_HOGs.csv`.

Required columns in the input include at least:
- `globHOG`, `genome`, `arrayID`

## Outputs
A CSV with columns:
- `HOG_ID`
- `Gene_ID`
- `Genome`
- `array_ID`
- `ofID`
- `isArrayRep`
- `chr`
- `start`
- `end`

It only contains genes whose HOG appears in at least one `(HOG_ID, genome, arrayID)` group with two or more rows.

## Usage
`python all_genes_in_tandem_HOGs.py <combBed.tsv> -o <output.csv>`

## Dependencies
- Python standard library: `argparse`, `csv`, `collections.defaultdict`

## Important processing steps / limitations
- It defines a tandem array as a HOG/genome/array group with at least two gene rows.
- It uses `globHOG` directly as the HOG identifier and does not split off any orthogroup suffix; this is a potential schema dependency.
- It keeps all genes from every qualifying HOG, not just the tandem-array genes themselves.
- The script expects a full `combBed` file and does not handle missing fields gracefully beyond the default CSV reader behavior.