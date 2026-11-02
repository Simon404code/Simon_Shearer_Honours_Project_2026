#!/usr/bin/env python3
# python3 all_genes_in_tandem_HOGs.py combBed_v5_complete-2.txt -o all_genes_in_tandem_HOGs.csv

#takes the combbed file as input and filters through the HOG_IDs to find those HOGs that contain 
#tandem arrays. This is defined as HOGs that have:
#at least 2 genes that share a HOG ID and an array ID in 1 genome. (there is a tandem array of at least length 2)


import argparse
import csv
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Extract genes from HOGs that contain at least one "
            "tandem array of length >= 2 within a single genome."
        )
    )

    parser.add_argument(
        "combbed",
        help="Input combBed TSV file"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="all_genes_in_tandem_HOGs.csv",
        help="Output CSV file"
    )

    args = parser.parse_args()

    # ------------------------------------------------------------
    # Read combBed
    # ------------------------------------------------------------

    with open(args.combbed, "r", newline="") as f:

        reader = csv.DictReader(f, delimiter="\t")

        rows = list(reader)

    print("Input genes: {:,}".format(len(rows)))

    # ------------------------------------------------------------
    # Group genes by:
    #
    #     HOG_ID
    #     Genome
    #     arrayID
    #
    # This is the definition of a potential tandem array.
    # ------------------------------------------------------------

    groups = defaultdict(list)

    for row in rows:

        # globHOG contains something like:
        #
        # N0.HOG0000763 OG0000003

        hog_id = row["globHOG"]
        genome = row["genome"]
        array_id = row["arrayID"]

        # The complete grouping key is:
        #
        # HOG + genome + array
        #
        key = (
            hog_id,
            genome,
            array_id
        )

        groups[key].append(row)

    # ------------------------------------------------------------
    # Identify groups containing >= 2 genes.
    #
    # These are tandem arrays of length >= 2.
    # ------------------------------------------------------------
    tandem_groups = {
        key: group
        for key, group in groups.items()
        if len(group) >= 2
    }

    print(
        f"Tandem HOG/genome/array groups: "
        f"{len(tandem_groups):,}"
    )

    # ------------------------------------------------------------
    # Determine which HOGs contain at least one tandem array.
    # ------------------------------------------------------------
    tandem_hogs = {
        key[0]
        for key in tandem_groups
    }

    print(
        f"HOGs containing tandem arrays: "
        f"{len(tandem_hogs):,}"
    )

    # ------------------------------------------------------------
    # Collect all genes that belong to any HOG with a tandem array.
    # ------------------------------------------------------------
    output_rows = []

    for row in rows:
        hog_id = row["globHOG"]
        if hog_id in tandem_hogs:
            output_rows.append(row)

    # ------------------------------------------------------------
    # Output
    # ------------------------------------------------------------
    output_columns = [
        "HOG_ID",
        "Gene_ID",
        "Genome",
        "array_ID",
        "ofID",
        "isArrayRep",
        "chr",
        "start",
        "end"
    ]

    with open(args.output, "w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=output_columns
        )

        writer.writeheader()

        for row in output_rows:

            writer.writerow({
                "HOG_ID": row["globHOG"],
                "Gene_ID": row["id"],
                "Genome": row["genome"],
                "array_ID": row["arrayID"],
                "ofID": row["ofID"],
                "isArrayRep": row["isArrayRep"],
                "chr": row["chr"],
                "start": row["start"],
                "end": row["end"]
            })

    print(
        f"Genes written from tandem-array HOGs: "
        f"{len(output_rows):,}"
    )

    print(
        f"Output: {args.output}"
    )


if __name__ == "__main__":
    main()