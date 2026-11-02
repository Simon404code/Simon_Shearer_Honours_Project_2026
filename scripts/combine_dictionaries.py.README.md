# combine_dictionaries.py

## What it does
Combines every line from every regular file in a directory into a single output file. The script is a generic concatenation tool rather than a biology-specific pipeline step.

## Inputs
- `input_dir`: directory containing source files to concatenate.
- `output_file`: destination path for the merged content.

## Outputs
A single file containing the concatenated contents of all regular files in the input directory, in sorted filename order.

## Usage
`python combine_dictionaries.py <input_dir> <output_file>`

## Dependencies
- Python standard library: `argparse`, `pathlib`

## Important processing steps / limitations
- It skips the output file itself if it is inside the same directory and would otherwise be read recursively.
- It writes each source file sequentially and appends a newline if the source line does not end with one.
- It is not schema-aware and does not validate content structure; it simply concatenates lines.
- It assumes all files in the directory are plain text and suitable for concatenation.