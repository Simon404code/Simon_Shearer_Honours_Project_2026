"Combine every line from every file in a directory into one output file."

import argparse
from pathlib import Path

def combine_files(input_dir, output_file):
    """Write every line from every regular file in input_dir to output_file."""
    input_path = Path(input_dir)
    output_path = Path(output_file)

    if not input_path.is_dir():
        raise NotADirectoryError(f"Input directory does not exist: {input_dir}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines_written = 0

    with open(output_path, "w", encoding="utf-8", newline="") as outfile:
        for input_file in sorted(input_path.iterdir()):
            if not input_file.is_file() or input_file.resolve() == output_path.resolve():
                continue

            print(f"Reading {input_file.name}")
            with open(input_file, "r", encoding="utf-8", newline="") as infile:
                for line in infile:
                    outfile.write(line)
                    if not line.endswith(("\n", "\r")):
                        outfile.write("\n")
                    lines_written += 1

    return lines_written

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combine every line from every file in a directory into one file."
    )
    parser.add_argument(
        "input_dir",
        type=str,
        help="Directory containing the files to combine."
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="Path to the combined output file."
    )

    args = parser.parse_args()
    lines_written = combine_files(args.input_dir, args.output_file)
    print(f"Wrote {lines_written} lines to: {args.output_file}")
                    