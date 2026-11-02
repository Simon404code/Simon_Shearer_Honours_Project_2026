
#!/usr/bin/env python3

import csv
import argparse
from collections import defaultdict


def make_hog_matrix(data_file, output_file):

    # hog_matrix[HOG_ID][haplotype] = list of gene IDs
    hog_matrix = defaultdict(lambda: defaultdict(list))

    haplotypes = set()

    print(f"Reading: {data_file}")

    with open(data_file, "r", newline="", encoding="utf-8-sig") as f:

        # Read the .txt file as tab-separated
        reader = csv.DictReader(f, delimiter="\t")

        required_columns = {"id", "genome", "globHOG"}

        if not required_columns.issubset(reader.fieldnames):
            raise ValueError(
                f"Missing required columns. "
                f"Expected: {required_columns}. "
                f"Found: {reader.fieldnames}"
            )

        for row in reader:

            gene_id = row["id"]
            genome = row["genome"]
            glob_hog = row["globHOG"]

            # Skip rows without a HOG ID
            if not glob_hog:
                continue

            # Extract only the HOG ID
            # Example: N0.HOG0012215 OG0000517
            # Result:  N0.HOG0012215
            hog_id = glob_hog.split()[0]

            haplotypes.add(genome)

            # Avoid duplicate gene IDs in the same cell
            if gene_id not in hog_matrix[hog_id][genome]:
                hog_matrix[hog_id][genome].append(gene_id)

    # Sort for consistent output
    haplotypes = sorted(haplotypes)
    hog_ids = sorted(hog_matrix.keys())

    print()
    print(f"Total haplotypes: {len(haplotypes)}")
    print(f"Total HOG IDs:    {len(hog_ids)}")

    # Write Excel-compatible CSV
    with open(output_file, "w", newline="", encoding="utf-8-sig") as f:

        writer = csv.writer(
            f,
            delimiter=",",
            quoting=csv.QUOTE_MINIMAL,
            lineterminator="\n"
        )

        # Header
        writer.writerow(["HOG_ID"] + haplotypes)

        # Rows
        for hog_id in hog_ids:

            row = [hog_id]

            for haplotype in haplotypes:

                genes = hog_matrix[hog_id].get(haplotype, [])

                # Multiple gene IDs in one cell
                cell_value = ";".join(genes)

                row.append(cell_value)

            writer.writerow(row)

    print()
    print("Finished")
    print("------------------------------")
    print(f"HOG IDs:       {len(hog_ids)}")
    print(f"Haplotypes:    {len(haplotypes)}")
    print(f"Output file:   {output_file}")


def main():

    parser = argparse.ArgumentParser(
        description="Create a HOG-by-haplotype gene matrix from a TSV .txt file."
    )

    parser.add_argument(
        "--data",
        required=True,
        help="Input tab-separated .txt data file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV file"
    )

    args = parser.parse_args()

    make_hog_matrix(
        args.data,
        args.output
    )


if __name__ == "__main__":
    main()