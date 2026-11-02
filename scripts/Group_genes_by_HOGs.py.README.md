# Group_genes_by_HOGs.py

## What it does
Builds a HOG-level summary table from a gene CSV and a combBed file. It joins in provenance metadata for each genome and computes summary statistics such as number of genes, average gene size, average gene position, number of genomes, tandem-array counts, and average array size.

## Inputs
Positional arguments:
- `genes`: CSV containing gene-level records for tandem-array HOGs.
- `combbed`: full `combBed` TSV or TSV.GZ.
- `-o/--output`: CSV output path (default: `Group_genes_by_HOGs.csv`).
- `--provenance`: semicolon-delimited genome metadata CSV; defaults to a file named `Pangenome_provenaces_edited.csv` next to the script.

Required columns:
- Gene CSV: `HOG_ID`, `Gene_ID`, `Genome`, `chr`, `start`, `end`
- combBed: `globHOG`, `genome`, `arrayID`, `id`
- Provenance CSV: `ProjectName`, `Location`, `Region`, `Subpopulation`

## Outputs
A CSV with one row per HOG containing columns such as:
- `HOG_ID`
- `number_of_genomes`
- `genomes`
- `locations`
- `regions`
- `subpopulations`
- `total_number_of_genes`
- `total_number_of_arrays`
- `average_array_size`
- `average_genes_per_genome`
- `average_gene_size`
- `average_gene_position`
- `chromosomes`

## Usage
`python Group_genes_by_HOGs.py <genes.csv> <combBed.tsv> -o <output.csv> [--provenance <metadata.csv>]`

## Dependencies
- `pandas`
- Python standard library: `argparse`, `pathlib`, `sys`

## Important processing steps / limitations
- The script validates that all gene IDs are unique and that provenance is present for every genome in the gene table.
- It groups by HOG and then by genome/array ID in combBed to compute tandem-array sizes.
- It raises errors if any genome lacks complete provenance metadata or if there are conflicting metadata values for the same genome.
- The script assumes a very specific metadata schema and exact column names.