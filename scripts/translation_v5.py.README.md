# translation_v5.py

## What it does
Loads a directory of genome-specific gene ID dictionaries and uses them to translate gene IDs in a tab-delimited data file. It selects the correct dictionary from the genome name stored in column 8 of each row, replaces the gene ID in column 4, and writes a translated output plus a missing-ID report.

## Inputs
- `--data`: input TSV/TSV-like file.
- `--dictionaries`: directory containing dictionary files. Each file should contain two tab-delimited columns: original gene ID and translated gene ID. File names are used to infer genome names; common suffixes such as `_pre_to_post_mapping.tsv` are stripped.
- `--output`: translated output path.
- `--missing` (optional): missing-ID report, defaulting to `<output>.missing.tsv`.

Expected row format for the main data file:
- Tab-delimited rows.
- The gene ID to translate is in column 4 (`row[3]`), and the genome identifier used to select the dictionary is in column 8 (`row[7]`).
- The script skips rows with fewer than 8 columns and copies them unchanged.
- A header row is preserved as-is unless the first cell looks like `pre_id`, `geneid`, `gene_id`, or `original_id`.

## Outputs
- Translated data file written to `--output`.
- A TSV report of untranslated/missing genes at `--missing` with columns: `genome`, `geneID`, `reason`.
- Console logs showing counts of rows processed, translated IDs, missing dictionaries, and missing gene IDs.

## Usage
`python translation_v5.py --data <input.tsv> --dictionaries <dict_dir> --output <translated.tsv> [--missing <report.tsv>]`

Assumes the current working directory is the project directory or that paths are supplied explicitly. The script does not search for files relative to the script location.

## Dependencies
- Python standard library: `csv`, `argparse`, `os`, `collections` (the module is imported and used for dictionary grouping, though the script is mostly CSV-driven).

## Important processing steps / limitations
- Dictionary lookup is genome-based and uses substring matching as a fallback when a direct genome name match is unavailable.
- It can warn when the same original gene ID maps to multiple translated IDs in the same dictionary.
- Rows without a matching genome dictionary are copied unchanged and counted in the missing report.
- The script assumes the relevant genome identifier is in column 8 and the target gene is in column 4; it does not generalize to other schemas without code changes.
- It does not validate downstream file types beyond the dictionary directory and the input file path.