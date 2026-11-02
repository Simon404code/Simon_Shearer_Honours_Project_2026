# find_hog_blasts_v3.py

## What it does
Finds BLAST hits for gene pairs belonging to the same HOG across the genomes represented in a tandem-array gene table. It normalizes genome names, maps them to BLAST genome numbers using a translation file, reads the corresponding BLAST files, and writes a CSV containing the matched gene pairs and BLAST metrics.

## Inputs
Positional arguments:
- `tandem_genes`: CSV file produced by the tandem-gene extraction workflow.
- `blast_directory`: directory containing BLAST files named like `Blast<query>_<subject>.txt` or `.txt.gz`.
- `translation_file`: file describing genome-to-BLAST-number mapping in a `number: genome` style.
- `-o/--output` (optional): output CSV path, default `tandem_hog_blast_results.csv`.

Expected data contract for `tandem_genes`:
- `HOG_ID`, `Gene_ID`, `Genome`, `ofID`, `isArrayRep`, `chr`, `start`, `end`

The BLAST translation file is expected to contain lines such as:
- `0: TAG0014Hap1.fa`
- `1: RM01Hap1.fa`

## Outputs
A CSV with columns including:
- `HOG_ID`
- `query_Gene_ID`, `query_Genome`, `query_ofID`
- `hit_Gene_ID`, `hit_Genome`, `hit_ofID`
- `query_blast_genome`, `hit_blast_genome`
- `pident`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, `bitscore`

Runtime output also prints counts for gene pairs examined, BLAST matches written, missing files, and processed files.

## Usage
`python find_hog_blasts_v3.py <tandem_genes.csv> <blast_dir> <translation_file> -o <output.csv>`

## Dependencies
- Python standard library: `argparse`, `csv`, `gzip`, `os`, `re`, `collections.defaultdict`

## Important processing steps / limitations
- Genome names are normalized by removing underscore, hyphen, and space characters and converting to uppercase, so the mapping must match the file naming conventions used by the BLAST outputs.
- It searches only for pairs of genes from the same HOG and only reads BLAST files for relevant genome pairs.
- If a genome has no BLAST number or a BLAST file is missing, it emits warnings and skips the affected comparisons.
- The script relies on the BLAST file format being tab-delimited and at least 12 columns wide.
- It stores only the last duplicate hit per query/subject pair and discards earlier duplicates.