# variation_between_arrays_within_HOGs_vs_gene_number.R

## Purpose
Generates a dot plot showing how tandem-array size varies across genomes for each HOG. The x-axis is the number of genomes contributing to a HOG, and the y-axis is the largest array-size difference within that HOG.

## Inputs and expected formats
This script expects a precomputed data frame named `tandem_array_HOG_summary` with at least these columns:
- `HOG_ID`
- `genome`
- `arrayID`
- `genes_in_array`

The script filters to `genes_in_array >= 2`, removes NA values, and keeps one record per `HOG_ID` + `genome` + `arrayID` with `distinct()`.

## Outputs
Writes these files to the current working directory:
- `tandem_array_variation_vs_genome_count.png`
- `tandem_array_variation_vs_genome_count.pdf`
- `tandem_array_variation_by_HOG.csv`

It also prints the plot.

## Usage
This is not a CLI script; it assumes the input table already exists in the session. Output files are saved to `getwd()`, so the working directory must be the intended output location.

## Dependencies
- `tidyverse`
- `ggplot2` (loaded via `tidyverse`)
- `scales`

## Processing steps and limitations
1. Keeps valid tandem arrays (genes >= 2) and deduplicates entries.
2. Computes, per HOG and genome, the largest tandem-array size.
3. Summarises each HOG's genome count and size range.
4. Filters out HOGs present in only one genome.
5. Produces a scatter plot with jitter and saves the graph and the aggregated table.

Limitations:
- It relies on an existing `tandem_array_HOG_summary` object rather than loading data from disk.
- It uses a simple `HOG_ID`-by-genome size-range summary and does not model complex distributional statistics.
