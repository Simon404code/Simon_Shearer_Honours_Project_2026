# blast_results_between_HOG_families.py

## What it does
Filters all tandem-array HOG gene pairs against a BLAST directory and writes only those BLAST alignments that match valid HOG gene pairs. It expects a genome-number mapping file and a tandem-HOG table with gene IDs and `ofID`s, then reads compressed BLAST files to identify direction-specific matches between genes from the same HOG.

## Inputs
The script reads positional arguments in order:
1. `HOG_FILE`: a tandem-array HOG summary file with columns `HOG_ID`, `genome`, `id`, and `ofID`.
2. `BLAST_DIR`: directory containing BLAST files such as `Blast0_1.txt.gz`.
3. `GENOME_MAP_FILE`: whitespace-delimited mapping file with `genome_name genome_number`.
4. `OUTPUT_FILE`: destination CSV/CSV.GZ path.

## Outputs
A filtered set of BLAST results including the HOG and genome context of each matching pair. The code also prints runtime progress and counts of total BLAST lines, matching lines, and files with hits.

## Usage
`python blast_results_between_HOG_families.py <tandem_hog_summary.csv> <blast_dir> <genome_map.txt> <output.csv>`

## Dependencies
- `pandas`
- Python standard library: `gzip`, `re`, `pathlib`, `collections.defaultdict`

## Important processing steps / limitations
- It creates all pairwise combinations within each HOG and only retains those that appear in the relevant BLAST files.
- It requires exact genome numbers to match the BLAST filename conventions.
- It validates that every gene's genome is present in the mapping file and exits early if any are missing.
- The script is older and less modular than the other pipeline pieces: it uses global `argv` variables and a direct script-style workflow rather than an argparse parser.
- It assumes BLAST files are compressed `.txt.gz` and named in the expected `Blast<number>_<number>.txt.gz` pattern.