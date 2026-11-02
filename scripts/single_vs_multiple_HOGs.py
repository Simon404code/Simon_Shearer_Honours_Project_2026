import csv
import argparse
from collections import defaultdict


def classify_arrays(input_file, output_file):

    arrays = defaultdict(lambda: {
        "genes": [],
        "hogs": set()
    })

    with open(input_file, "r", encoding="utf-8-sig", newline="") as infile:

        reader = csv.DictReader(infile, delimiter="\t")

        for row in reader:

            genome = row["genome"].strip()
            array_id = row["arrayID"].strip()
            gene_id = row["id"].strip()
            glob_hog = row["globHOG"].strip()

            if not array_id:
                continue

            key = (genome, array_id)

            arrays[key]["genes"].append(gene_id)

            if glob_hog:
                hog_id = glob_hog.split()[0]
                arrays[key]["hogs"].add(hog_id)

    with open(output_file, "w", encoding="utf-8", newline="") as outfile:

        writer = csv.writer(outfile, delimiter="\t")

        writer.writerow([
            "genome",
            "arrayID",
            "gene_count",
            "HOG_count",
            "array_type",
            "globHOG_ids",
            "gene_ids"
        ])

        for (genome, array_id), data in sorted(arrays.items()):

            gene_count = len(data["genes"])

            if gene_count < 2:
                continue

            hog_count = len(data["hogs"])

            if hog_count == 1:
                array_type = "single-HOG"
            else:
                array_type = "multi-HOG"

            writer.writerow([
                genome,
                array_id,
                gene_count,
                hog_count,
                array_type,
                ";".join(sorted(data["hogs"])),
                ";".join(data["genes"])
            ])


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Classify tandem arrays by their HOG composition."
    )

    parser.add_argument("--data", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    classify_arrays(args.data, args.output)
