# blast_stats_per_hog.py

## What it does
Reads a BLAST results CSV and calculates per-HOG summary statistics for chromosome diversity, protein-length diversity, and similarity metrics. It aggregates by the first column, which is assumed to hold the HOG ID.

## Inputs
Positional argument:
- `blast_csv`: CSV containing BLAST results with a `HOG_ID`-like first column.
- `-o/--output` (optional): output CSV path, default `blast_stats_per_hog.csv`.

It tries to find chromosome and protein-length columns by alias names such as:
- `query_Chromosome`, `hit_Chromosome`, `query_Protein_Length`, `hit_Protein_Length`, `pident`

## Outputs
A CSV summarizing each HOG with columns such as:
- `HOG_ID`
- `blast_rows`
- `different_chromosomes`
- `different_protein_lengths`
- `protein_length_difference_min`
- `protein_length_difference_max`
- `protein_length_difference_mean`
- `similarity_min`
- `similarity_max`
- `similarity_mean`
- `similarity_range`

## Usage
`python blast_stats_per_hog.py <blast_results.csv> -o <summary.csv>`

## Dependencies
- Python standard library: `argparse`, `csv`, `math`, `re`, `sys`, `collections.defaultdict`

## Important processing steps / limitations
- The script is tolerant of multiple column naming variants but may warn if it does not find expected fields.
- It computes min/max/sum/count statistics without storing every value, which keeps memory usage lower.
- It only aggregates over rows encountered for a given HOG; any rows missing necessary metrics are treated as absent rather than fatal errors.
- The script assumes a generic BLAST output style and does not validate the full BLAST schema beyond the header names it recognizes.