# violin_plots.R

## Purpose
This is an exploratory plotting script rather than a polished pipeline. It creates several ggplot figures from a prebuilt `out_table2.5` summary table, mostly to inspect tandem-array rates, HOG counts, gene counts, and CNV-related trends.

## Inputs and expected formats
The script assumes a data frame named `out_table2.5` is already in the R environment. Relevant columns appear to include:
- `classifier`
- `tandem_array_rate_percent`
- `number_of_Genomes`
- `number_of_Genes`
- `avg_CNV`
- `avg_gene_length`

These values are used in `ggplot()` calls without any file reading or validation.

## Outputs
No files are written. The script prints graphs to the active R plotting device (for example, RStudio or a graphics window).

## Usage
Not a standalone CLI script. The R session must already contain the prepared summary table, and the user likely runs the relevant ggplot block(s) individually. The script does not define a working-directory assumption or file paths.

## Dependencies
- `tidyverse`

## Processing steps and limitations
- Uses a series of ad hoc examples and duplicate plotting blocks.
- No data cleaning or validation is performed.
- Comments suggest intended analyses (box plots, violin plots, dot plots, CNV trends), but the script is a scratch notebook rather than a reproducible workflow.
- It does not save outputs or support command-line options.
