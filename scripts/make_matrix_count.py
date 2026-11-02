import csv
import argparse
from collections import defaultdict


def make_hog_count_matrix(data_file, output_file):

    # Store gene IDs for each HOG and genome
    hog_matrix = defaultdict(lambda: defaultdict(set))
    haplotypes = set()

    # Read the tab-separated input file
    with open(data_file, "r", encoding="utf-8-sig", newline="") as infile:

        reader = csv.DictReader(infile, delimiter="\t")

        required_columns = {"id", "genome", "globHOG"}

        if not required_columns.issubset(reader.fieldnames):
            missing = required_columns - set(reader.fieldnames)
            raise ValueError(
                f"Missing required columns: {', '.join(missing)}"
            )

        for row in reader:

            gene_id = row["id"]
            genome = row["genome"]
            glob_hog = row["globHOG"]

            # Skip rows without a HOG
            if not glob_hog or glob_hog.strip() == "":
                continue

            # Example:
            # N0.HOG0008484 OG0000263
            #
            # We only use the HOG ID:
            # N0.HOG0008484
            hog_id = glob_hog.split()[0]

            haplotypes.add(genome)

            # Store the gene ID.
            # Using a set prevents the same gene being counted twice.
            hog_matrix[hog_id][genome].add(gene_id)

    # Sort genomes and HOGs
    haplotypes = sorted(haplotypes)
    hog_ids = sorted(hog_matrix.keys())

    # Write CSV output
    with open(
        output_file,
        "w",
        encoding="utf-8-sig",
        newline=""
    ) as outfile:

        writer = csv.writer(
            outfile,
            delimiter=",",
            quoting=csv.QUOTE_MINIMAL
        )

        # Header
        writer.writerow(["HOG_ID"] + haplotypes)

        # Rows
        for hog_id in hog_ids:

            row = [hog_id]

            for genome in haplotypes:

                # Number of genes belonging to this HOG
                count = len(hog_matrix[hog_id].get(genome, set()))

                row.append(count)

            writer.writerow(row)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Create a HOG x haplotype gene-count matrix."
    )

    parser.add_argument(
        "--data",
        required=True,
        help="Input tab-separated .txt file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV file"
    )

    args = parser.parse_args()

    make_hog_count_matrix(
        args.data,
        args.output
    )