# box_and_whisker.R

## Purpose
Produces a boxplot and quasirandom point display of tandem-array counts by size, with clear x-axis labels showing percentage and total counts for each array size.

## Inputs and expected formats
The script requires a `plot_data` object and a `tandem_array_HOG_summary` table in the R session. The main fields used are:
- `genome`
- `arrayID`
- `genes_in_array`
- `Subpopulation`
- `number_of_arrays`

It also defines a hard-coded `array_totals` tibble summarising total counts per size. The plot data is derived from counts by haplotype and array size.

## Outputs
Writes these files to the current working directory:
- `tandem_array_boxplots_subpopulation_with_totals.png`
- `tandem_array_boxplots_subpopulation_with_totals.pdf`

It also prints the plot.

## Usage
The script expects the analysis objects to already exist in memory and uses `ggsave()` to save outputs in the active working directory. No input paths are configured.

## Dependencies
- `tidyverse`
- `ggbeeswarm`

## Processing steps and limitations
1. Sets up the fixed `array_totals` table.
2. Converts totals to percentages for annotation.
3. Builds a plot-level data frame filtered to valid array sizes and zero-free totals.
4. Creates a labelled boxplot with `geom_quasirandom()` points coloured by subpopulation.
5. Saves a wide PNG and a PDF version.

Limitations:
- Relies on downstream data objects and hard-coded totals rather than reading from a data source.
- Most labelling is static and tuned to a specific dataset.
- Not a reusable generic analysis script for other assemblies or experiments.
