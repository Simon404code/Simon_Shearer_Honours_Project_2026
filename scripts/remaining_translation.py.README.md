# remaining_translation.py

## What it does
Translates residual gene IDs from a combined dictionary file. It targets rows whose gene ID starts with `PAV.` and strips the prefix before lookup, then restores it after a successful map. It also handles non-PAV IDs in the same way via a generic dictionary lookup.

## Inputs
- `--data`: input TSV data file.
- `--dictionaries`: a single TSV file with at least two columns: `original_id` and `translated_id`.
- `--output`: output TSV.

Expected data contract:
- The gene ID to replace is in column 4 (`row[3]`).
- The script reads the file row-by-row and writes all rows back out, replacing only column 4 as needed.

## Outputs
- A translated output TSV at `--output`.
- Console counts for total rows, replaced IDs, unchanged IDs, PAV counts, and PAV replacements.

## Usage
`python remaining_translation.py --data <input.tsv> --dictionaries <combined_dictionary.tsv> --output <translated.tsv>`

## Dependencies
- Python standard library: `argparse`, `csv`.

## Important processing steps / limitations
- It treats any `PAV.` prefix specially: `PAV.<id>` becomes `<id>` before dictionary lookup and is restored as `PAV.<translated_id>` if found.
- If a gene ID is not found, the original value is retained unchanged.
- The script is narrow: it assumes a single dictionary file and column 4 semantics. It does not infer genome-specific dictionaries.