# genespace_run_50.R

## Purpose
Runs a GENESPACE analysis over a set of 50 haplotype genomes and then performs additional riparian-window processing from FASTA files and `synBlkCoords` outputs.

## Inputs and expected formats
The script expects a genome repository structured under:
- `~/h_project/50_genespace_run/rss/`

The list `genomes2run` contains 50 genome identifiers, for example `RK05HAP1`, `RK25HAP2`, `RM02HAP1`, and `TAG0014HAP2`. The script also expects:
- `MCScanX` installed or available at `~/h_project/50_genespace_run/MCScanX/`
- a working directory at `/scratch/genespace2`
- a secondary FASTA repository at `50_genespace_run/other_rss/files_fa/`
- local files `clean_windows.R` and `riparian_windows.R` in the working directory

The script uses `parse_annotations()`, `init_genespace()`, `run_genespace()`, and `plot_riparian()`, then later calls `init_windows()` and `riparian_windows()`.

## Outputs
The script is likely intended to generate GENESPACE outputs in `/scratch/genespace2` and derived riparian-window products. The code does not explicitly save a final summary file, but it expects to populate the GENESPACE work directory and produce output from `synBlkCoords` and riparian-window functions.

## Usage
This script is hard-coded to specific local paths and a specific set of genomes. It is not a portable CLI script and should be run only in the corresponding project environment. It assumes the `GENESPACE` package and related dependencies are installed.

## Dependencies
- `GENESPACE`
- `BiocGenerics`
- `IRanges`
- `GenomeInfoDb`
- `GenomicRanges`
- `GenomeInfoDbData`
- `rtracklayer`
- `data.table`
- `XVector`
- `Biostrings`
- `ggplot2`
- `devtools`
- `parallel`
- `fpc`
- `dbscan`
- `dplyr`

## Processing steps and limitations
1. Defines a genome repository, working directory, and MCScanX path.
2. Calls `parse_annotations()` and `init_genespace()` for the 50-genome list.
3. Runs `run_genespace()` and then `plot_riparian()` with `refGenome = "TAGHAP1"`.
4. Finds FASTA files under a second repository and assigns them to genome IDs by stripping prefixes/suffixes.
5. Sources `clean_windows.R` and `riparian_windows.R` from the working directory.
6. Initialises windows and runs `clean_windows()` and `riparian_windows()`.

Limitations:
- Uses many absolute and user-specific paths (`~`, `/scratch`, hard-coded genome IDs).
- The script is not portable and may fail on another machine or in a different project layout.
- It assumes required helper files are present and in the current working directory.
- There may be a mismatch in the reference name (`TAGHAP1` vs `TAG0014HAP1` in the genome list), which should be checked before relying on the plot output.
