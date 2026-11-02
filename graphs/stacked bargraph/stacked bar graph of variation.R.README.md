# stacked bar graph of variation.R

## Purpose
Creates a stacked bar chart showing, for each number of genomes containing a HOG, the split between zero-variation and non-zero-variation tandem arrays.

## Inputs and expected formats
This script expects a `plot_data` object already present in the workspace. The relevant fields are the count of genomes per HOG (`number_of_genomes`) and a variation metric (`variation`) used to classify each row as:
- `Zero variation`
- `Variation`

The data are transformed with `if_else()` and then expanded to include all genome counts from 2 to 50 using `complete()`.

## Outputs
Writes these files to the current working directory:
- `tandem_array_zero_vs_nonzero_variation_stacked.png`
- `tandem_array_zero_vs_nonzero_variation_stacked.pdf`

It also prints the graph.

## Usage
This script assumes the upstream analysis has already produced `plot_data` in the same R session. It saves output to the current working directory and does not accept arguments.

## Dependencies
- `ggplot2`
- `dplyr`

## Processing steps and limitations
1. Reclassifies HOGs as zero or non-zero variation.
2. Counts HOGs by `number_of_genomes` and variation category.
3. Completes missing categories for all genome counts from 2 to 50.
4. Adds total-HOG labels above each stacked bar.
5. Produces a grouped stacked bar chart and exports PNG/PDF versions.

Limitations:
- It depends on a preexisting `plot_data` structure and does not validate it.
- It hard-codes the genome-count range as `2:50` and uses a fixed chart size.
