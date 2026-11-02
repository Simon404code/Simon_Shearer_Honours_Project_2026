# pie_charts.R

## Purpose
Builds three pie charts summarising the distribution of tandem-array carrying HOGs, HOGs with tandem arrays, and genes across the pangenome using fixed counts supplied in the script.

## Inputs and expected formats
No external files are read. The script defines three small data frames with columns:
- `category`
- `n`

Categories and counts are hard-coded for:
1. tandem HOG distribution (`shell`, `core`, `soft-core`, `singletons`)
2. all HOGs vs tandem HOGs
3. all genes vs tandem genes

## Outputs
No output files are saved. The pie charts are rendered directly in the graphics device.

## Usage
This script is effectively a static figure generator. No working-directory assumptions are needed because it does not read or save files. It is intended to be run as a whole in an R session.

## Dependencies
- `dplyr`
- `ggplot2`
- `ggrepel`

## Processing steps and limitations
- Builds a helper function `make_pie()` that converts counts to percentages and labels.
- Creates three separate plots using `geom_col()` with polar coordinates.
- Uses `geom_text_repel()` to place percentage labels.

Limitations:
- Counts are hard-coded rather than computed from input data.
- There is no export automation, so output is only visible in the active graphics device.
