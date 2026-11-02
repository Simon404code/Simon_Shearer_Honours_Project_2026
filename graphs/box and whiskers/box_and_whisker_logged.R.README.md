# box_and_whisker_logged.R

## Purpose
Generates a log-scaled boxplot of tandem-array counts by array size, with subpopulation colour coding per haplotype, and includes the percentage and total count for each array size in the x-axis labels.

## Inputs and expected formats
The script assumes a `tandem_array_HOG_summary` table is already in the R environment. Relevant fields include:
- `genome`
- `arrayID`
- `genes_in_array`
- `Subpopulation`

It also defines a hard-coded `array_totals` tibble listing array sizes and their total counts, for example `2, 22871`, `3, 2104`, and so on. This is used to create size labels and to limit the plot to selected non-zero array sizes.

## Outputs
Writes these files to the current working directory:
- `tandem_array_boxplots_subpopulation_log10.png`
- `tandem_array_counts_by_haplotype_subpopulation_no_zeros.csv`

It also displays the plot in the graphics device.

## Usage
This script is not a CLI tool. It relies on precomputed data and writes plots to the current working directory. There are no command-line flags or external configuration files.

## Dependencies
- `tidyverse`
- `ggbeeswarm`
- `ggplot2`
- `scales`

## Processing steps and limitations
1. Predefines `array_totals` and converts them to percentages.
2. Filters `tandem_array_HOG_summary` to valid arrays and deduplicates by genome and array ID.
3. Checks whether a haplotype is assigned to multiple subpopulations.
4. Counts arrays by genome, subpopulation, and size.
5. Builds a log10 boxplot with quasirandom jittered points for each haplotype.
6. Labels each x-axis category with percentage, total count, and array size.
7. Saves the figure and the filtered summary table.

Limitations:
- Uses strongly hard-coded array-size totals and selected size ranges.
- Assumes the `Subpopulation` field is present and consistent.
- The output file locations are relative to the current working directory.
