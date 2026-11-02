# isolated_runs.R

## Purpose
This script is a near-duplicate of the boxplot analysis for tandem-array counts by size, but it limits the plot to the non-zero counts used in the earlier visualisation and keeps the subpopulation annotation.

## Inputs and expected formats
Expected inputs are the same as in the other boxplot scripts:
- `tandem_array_HOG_summary` with `genome`, `arrayID`, `genes_in_array`, and `Subpopulation`
- a hard-coded `array_totals` table with size totals and percentages

It filters the data to valid array sizes and stores grouped counts in `plot_data`.

## Outputs
Writes these files to the current working directory:
- `tandem_array_boxplots_subpopulation_log10.png`
- `tandem_array_counts_by_haplotype_subpopulation_no_zeros.csv`

It also prints the plot.

## Usage
No command-line arguments are present. It assumes the relevant objects already exist and saves outputs relative to the session working directory.

## Dependencies
- `ggplot2`
- `tidyverse`
- `ggbeeswarm`

## Processing steps and limitations
1. Builds a fixed `array_totals` table and calculates percentages.
2. Filters tandem arrays to valid non-missing records and removes duplicates by `genome` and `arrayID`.
3. Warns if haplotypes map to multiple subpopulations.
4. Counts arrays by `genome`, `Subpopulation`, and `genes_in_array`.
5. Produces a boxplot with log10 y-axis, quasirandom points, and custom x-axis labels.
6. Saves the figure and a CSV summary table.

Limitations:
- Like the other boxplot scripts, it is tightly parameterised to a specific dataset.
- It does not read file inputs; all required numbers and labels are embedded in the script.
- Hard-coded array sizes mean it is not adaptable without editing the source.
