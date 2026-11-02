# Tandem gene array project

This repository contains scripts and analysis outputs for a study of tandem
gene arrays, gene-family (HOG) membership, and variation across haplotypes and
genomes. It includes Python analysis scripts, R plotting and GeneSpace
workflows, figures, and notes for running analyses on the UConn HPC system.

## Repository contents

| Path | Contents |
|---|---|
| [`isolate_tandem_genes.py`](isolate_tandem_genes.py) | Top-level Python tandem-gene isolation script. |
| [`scripts/`](scripts/) | Python workflows for grouping genes and HOGs, processing BLAST results, translation, CNV, matrices, and tables. |
| [`R-scripts/`](R-scripts/) | R analyses and plots, including per-haplotype tandem-gene and violin-plot workflows. |
| [`graphs/`](graphs/) | R plotting scripts and generated figures/tables, organized by plot type. |
| [`scripted run/`](scripted%20run/) | GeneSpace R workflow and its SLURM batch submission script. |
| [`TAG0014_HAP1vsHAP2/`](TAG0014_HAP1vsHAP2/) | Synteny, dot-plot, and riparian outputs for the TAG0014 haplotype comparison. |
| [`UCONN_HPC_Access/`](UCONN_HPC_Access/) | HPC access and interactive-session notes. |
| [`IDs of the haplotypes used.txt`](IDs%20of%20the%20haplotypes%20used.txt) | List of haplotype identifiers used in the analysis. |

Every Python and R script has a companion `<script filename>.README.md` in the
same directory. Consult that file for the script's purpose, expected inputs,
outputs, dependencies, and usage details.

## Workflow overview

The scripts cover several related stages rather than one turnkey pipeline:

1. Prepare or isolate tandem genes and associate genes with HOGs.
2. Process gene-family and BLAST-related results, including comparisons and
   summaries.
3. Generate matrices, tables, and statistics used in downstream analyses.
4. Analyze and visualize tandem-array patterns with R, including variation
   among HOGs, haplotypes, and genomes.
5. Run comparative genome/synteny analyses with GeneSpace where needed.

The precise inputs and ordering vary by script. Start with the companion
README for each script you plan to use; there is no single project-wide command
that runs every analysis.

## Requirements and reproducibility

Requirements depend on the selected workflow. Python modules and R packages
are listed in the relevant script READMEs. The GeneSpace workflow additionally
uses the `GENESPACE` R package and MCScanX; HPC setup notes mention other R
packages and environment-specific configuration.

Input genome annotations, sequence files, intermediate analysis data, and
cluster environments may be external to this repository. Before running a
script:

- Confirm that its documented inputs are available and match the expected
  format.
- Review and update any machine-specific paths, genome IDs, and output
  locations in the script.
- Run it from the documented working directory and preserve required input
  naming conventions.
- Keep generated files separate from source inputs and record the software
  versions and parameters used for a reproducible run.

Many research scripts are configured for a particular dataset or HPC
environment. Treat example commands and paths in the script documentation as
dataset-specific unless stated otherwise.

## Running the GeneSpace HPC workflow

[`scripted run/genespace_run_50.R`](<scripted run/genespace_run_50.R>) is
submitted by [`scripted run/scripted_genespace_run.pbs`](<scripted run/scripted_genespace_run.pbs>).
The batch file activates a Conda environment named `genespace1` and runs the R
script with `Rscript`. Before submission, verify the environment, input paths,
MCScanX location, and the scheduler resource settings for your cluster.
Additional setup and interactive-session examples are in
[`UCONN_HPC_Access/`](UCONN_HPC_Access/).

## Figures and analysis outputs

The `graphs/` and `TAG0014_HAP1vsHAP2/` directories contain generated plots,
tables, and comparative-genomics results. These are outputs and examples, not
inputs required by every script. Check each script's companion README to
identify which outputs it creates and where it writes them.
