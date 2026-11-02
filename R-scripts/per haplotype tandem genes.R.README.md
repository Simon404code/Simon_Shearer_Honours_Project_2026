# per haplotype tandem genes.R

## Purpose
Summarises gene counts and tandem-array counts per haplotype (or genome) using two preloaded data frames: the full gene table (`combBed_v5_complete.2`) and the tandem-array gene list (`all_genes_in_tandem_HOGs_2`).

## Inputs and expected formats
Expected objects in the R workspace:
- `combBed_v5_complete.2`: must contain `genome` and `id` columns.
- `all_genes_in_tandem_HOGs_2`: must contain `HOG_ID`, `Gene_ID`, `Genome`, and `array_ID`.

The script treats `genome` and `Genome` as haplotype identifiers and `id`/`Gene_ID` as gene identifiers. It trims whitespace and removes blank/NA values before counting.

## Outputs
Writes a CSV file to the current working directory:
- `gene_and_tandem_HOG_array_counts_per_haplotype.csv`

It also prints a summary to the console and displays the table with `View()`.

## Usage
This script expects the input objects to already be loaded into the R session. It uses `getwd()` for output placement, so the file is written relative to the current project directory. There are no arguments or CLI flags.

## Dependencies
- `tidyverse`
- `writexl` (imported but not used here)

## Processing steps and limitations
1. Validates presence of required columns.
2. Counts unique genes across the full gene table per haplotype.
3. Filters tandem-array rows to non-empty gene IDs and array IDs.
4. Counts unique genes in tandem arrays and distinct tandem arrays per haplotype.
5. Joins counts and calculates `percentage_genes_in_tandem_arrays`.

Limitations:
- It assumes the objects exist in memory before execution.
- It does not read any source files; all inputs are expected to be precomputed datasets.
- The file name is hard-coded and no output directory is created.
