# dot_plot_and_cumulative_line_var_of_HOG_ID_and_genome_involvment.R

## Purpose
Plots tandem-array size variation by HOG presence and overlays the cumulative frequency of HOG occurrence. The figure compares the number of genomes containing a HOG against the variation in tandem-array size within that HOG.

## Inputs and expected formats
The script requires two objects in the R environment:
- `combBed_v5_complete.2` with `globHOG` and `genome`
- `tandem_array_HOG_summary` with `HOG_ID`, `genome`, and `genes_in_array`

It extracts a HOG ID from `globHOG` strings such as `N0.HOG0000782 OG0000003`, then counts how many genomes contain each HOG. For tandem arrays, it keeps `genes_in_array >= 2` and uses the largest array per HOG per genome when multiple arrays exist.

## Outputs
Writes these files to the current working directory:
- `tandem_array_variation_vs_HOG_presence.png`
- `tandem_array_variation_vs_HOG_presence.pdf`
- `tandem_array_variation_by_HOG_presence.csv`
- `tandem_array_cumulative_frequency.csv`

It also prints the plot and summary diagnostics to the console.

## Usage
This script assumes the required tables are already loaded. Output is written relative to the current working directory, and there are no command-line arguments or config files.

## Dependencies
- `tidyverse`
- `ggplot2`
- `scales`

## Processing steps and limitations
1. Counts the number of genomes containing each HOG.
2. Determines each HOG's largest tandem array per genome and computes size variation.
3. Joins HOG presence and variation tables.
4. Calculates a cumulative frequency curve across genome counts.
5. Draws a dual-axis scatter/line plot with a 50% reference line.
6. Writes both plot images and supporting CSV tables.

Limitations:
- It is designed around a specific `combBed_v5_complete.2` table structure.
- It uses fixed plotting aesthetics and a hard-coded 50% reference line based on the data maximum.
- It does not validate data beyond required-column checks.
