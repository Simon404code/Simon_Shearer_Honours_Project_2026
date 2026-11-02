# find_CNV.py

## What it does
Summarizes tandem-array copy number variation (CNV) per HOG. It counts genes grouped by `(Genome, array_ID, HOG_ID)`, ignores groups smaller than two genes, and reports minimum, maximum, mean array size, number of arrays, total genes captured, and the genomes involved.

## Inputs
Positional argument:
- `input`: a gene CSV containing at least `HOG_ID`, `Gene_ID`, `Genome`, and `array_ID`.
- `-o/--output` (optional): summary CSV path, default `CNV_summary.csv`.

## Outputs
A CSV with columns:
- `HOG_ID`
- `min_array_size`
- `max_array_size`
- `avg_CNV`
- `number_of_arrays`
- `number_of_genes_in_arrays`
- `genomes`

## Usage
`python find_CNV.py <input.csv> -o <output.csv>`

## Dependencies
- Python standard library: `argparse`, `csv`, `collections.Counter`, `collections.defaultdict`, `pathlib`

## Important processing steps / limitations
- It counts only arrays with at least two genes in the same genome and HOG.
- It ignores arrays with only one gene, which means it is not a full HOG copy-number table but a tandem-array summary.
- It does not compute statistics across chromosomes or array coordinates; it is purely count-based.
- It expects no blank required values in the input rows.