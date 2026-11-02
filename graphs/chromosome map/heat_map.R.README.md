# heat_map.R

## Purpose
Creates a chromosome-level heatmap showing the number of tandem genes in 400 kb genomic bins for the 11 chromosomes in the study. The plot is karyotype-like and intended to visualise tandem-gene distribution across chromosomes.

## Inputs and expected formats
The script requires a prebuilt object named `tandem_array_HOG_summary` containing at least:
- `chr`
- `genome`
- `arrayID`
- `start`
- `end`
- `genes_in_array`

It also defines a hard-coded `chromosome_lengths` tibble with chromosome names and lengths for:
- `Chr01` to `Chr11`

The script filters out scaffolds/other chromosome labels not found in this table.

## Outputs
Writes the following files into a folder named `tandem_gene_chromosome_karyotype` under the current working directory:
- `tandem_gene_counts_per_400kb_bin.csv`
- `tandem_gene_karyotype.png`
- `tandem_gene_karyotype.pdf`

It also prints the heatmap.

## Usage
The script assumes the tandem-array summary table is already loaded in the R session. Output is written relative to `getwd()`, and a directory named `tandem_gene_chromosome_karyotype` is created automatically. No command-line arguments are supported.

## Dependencies
- `tidyverse`
- `viridis`

## Processing steps and limitations
1. Validates that required input columns are present.
2. Standardises chromosome names and drops invalid rows.
3. Calculates each tandem array midpoint and keeps only arrays that fall within the 11 chromosomal lengths.
4. Bins tandem genes into 400 kb windows and creates complete bins for each chromosome.
5. Joins counts to the full chromosome bin table and produces a tiled heatmap.
6. Saves a CSV and high-resolution PNG/PDF figures.

Limitations:
- It explicitly excludes chromosomes not listed in the hard-coded length table.
- It assumes the `chromosome_lengths` table is accurate for this genome assembly.
- It is tuned for the 11-chromosome Eucalyptus-like assembly and not generic across other assemblies.
