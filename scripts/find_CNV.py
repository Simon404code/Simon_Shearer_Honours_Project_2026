#!/usr/bin/env python3
#py -3 .\find_CNV.py .\all_genes_in_tandem_HOGs_2.csv -o .\CNV_summary.csv
"""Summarize tandem-array sizes for each HOG."""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


REQUIRED_COLUMNS = {"HOG_ID", "Gene_ID", "Genome", "array_ID"}
OUTPUT_COLUMNS = [
    "HOG_ID",
    "min_array_size",
    "max_array_size",
    "avg_CNV",
    "number_of_arrays",
    "number_of_genes_in_arrays",
    "genomes",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize gene counts in tandem arrays grouped by genome, "
            "array_ID, and HOG_ID."
        )
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="all_genes_in_tandem_HOGs_2.csv",
        help="Input gene CSV (default: all_genes_in_tandem_HOGs_2.csv)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="CNV_summary.csv",
        help="Output summary CSV (default: CNV_summary.csv)",
    )
    return parser.parse_args()


def summarize_arrays(input_path):
    group_sizes = Counter()
    with open(input_path, "r", newline="", encoding="utf-8-sig") as input_file:
        reader = csv.DictReader(input_file)
        if reader.fieldnames is None:
            raise ValueError("Input CSV is empty or has no header.")

        missing_columns = sorted(REQUIRED_COLUMNS - set(reader.fieldnames))
        if missing_columns:
            raise ValueError(
                f"Input CSV is missing required columns: {missing_columns}"
            )

        for line_number, row in enumerate(reader, start=2):
            values = {
                column: (row.get(column) or "").strip()
                for column in REQUIRED_COLUMNS
            }
            missing_values = [
                column for column, value in values.items() if not value
            ]
            if missing_values:
                raise ValueError(
                    f"Input CSV has blank required values at line "
                    f"{line_number}: {missing_values}"
                )

            group_key = (
                values["Genome"],
                values["array_ID"],
                values["HOG_ID"],
            )
            group_sizes[group_key] += 1

    hog_stats = defaultdict(
        lambda: {
            "array_sizes": [],
            "genomes": set(),
        }
    )
    for (genome, _array_id, hog_id), group_size in group_sizes.items():
        if group_size < 2:
            continue
        stats = hog_stats[hog_id]
        stats["array_sizes"].append(group_size)
        stats["genomes"].add(genome)

    return hog_stats


def write_summary(hog_stats, output_path):
    with open(output_path, "w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()

        for hog_id in sorted(hog_stats):
            stats = hog_stats[hog_id]
            array_sizes = stats["array_sizes"]
            writer.writerow(
                {
                    "HOG_ID": hog_id,
                    "min_array_size": min(array_sizes),
                    "max_array_size": max(array_sizes),
                    "avg_CNV": sum(array_sizes) / len(array_sizes),
                    "number_of_arrays": len(array_sizes),
                    "number_of_genes_in_arrays": sum(array_sizes),
                    "genomes": "; ".join(sorted(stats["genomes"])),
                }
            )


def main():
    args = parse_args()
    input_path = Path(args.input)
    if not input_path.is_file():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    hog_stats = summarize_arrays(input_path)
    write_summary(hog_stats, args.output)
    print(f"HOGs summarized: {len(hog_stats):,}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, csv.Error, ValueError) as error:
        raise SystemExit(f"ERROR: {error}")



    ### count groups
#     #Eg_RM02_H2.07G019300
# Eg_RM02_H2.07G018900
# Eg_RK86_H2.07G022400
# Eg_RK86_H2.07G022200
# Eg_RK86_H2.07G020800
# Eg_RK84_H1.07G023000
# Eg_RK792_H2.07G021700
# Eg_RK792_H2.07G021400
# Eg_RK792_H1.07G022000
# Eg_RK792_H1.07G021900
# Eg_RK715_H2.07G019800
# Eg_RK687_H2.07G024300
# Eg_RK687_H2.07G024200
# Eg_RK687_H2.07G024100
# Eg_RK687_H2.07G024000
# Eg_RK687_H2.07G023900
# Eg_RK687_H2.07G023800
# Eg_RK687_H2.07G023700
# Eg_RK687_H2.07G023600
# Eg_RK687_H2.07G023500
# Eg_RK687_H2.07G023400
# Eg_RK687_H2.07G023300
# Eg_RK687_H2.07G023200
# Eg_RK687_H2.07G023100
# Eg_RK687_H2.07G023000
# Eg_RK687_H2.07G022900
# Eg_RK687_H2.07G022800
# Eg_RK687_H2.07G022700
# Eg_RK687_H2.07G022600
# Eg_RK687_H2.07G022500
# Eg_RK687_H2.07G022400
# Eg_RK687_H2.07G022300
# Eg_RK687_H2.07G022200
# Eg_RK687_H2.07G022100
# Eg_RK687_H2.07G022000
# Eg_RK687_H2.07G021900
# Eg_RK687_H2.07G021800
# Eg_RK687_H2.07G021700
# Eg_RK687_H2.07G021600
# Eg_RK687_H2.07G021500
# Eg_RK687_H2.07G021400
# Eg_RK687_H2.07G021300
# Eg_RK687_H2.07G021200
# Eg_RK687_H2.07G021100
# Eg_RK687_H2.07G021000
# Eg_RK687_H2.07G020900
# Eg_RK687_H2.07G020700
# Eg_RK687_H2.07G020600
# Eg_RK687_H2.07G020500
# Eg_RK687_H2.07G020400
# Eg_RK687_H2.07G020300
# Eg_RK687_H2.07G020200
# Eg_RK687_H2.07G020100
# Eg_RK687_H2.07G020000
# Eg_RK687_H2.07G019900
# Eg_RK687_H2.07G019800
# Eg_RK687_H2.07G019700
# Eg_RK687_H2.07G019600
# Eg_RK687_H2.07G019500
# Eg_RK687_H2.07G019400
# Eg_RK687_H2.07G019300
# Eg_RK687_H2.07G019200
# Eg_RK687_H2.07G019100
# Eg_RK687_H2.07G019000
# Eg_RK687_H2.07G018900
# Eg_RK687_H2.07G018700
# Eg_RK687_H2.07G018600
# Eg_RK687_H2.07G018500
# Eg_RK687_H2.07G018400
# Eg_RK687_H1.07G021500
# Eg_RK687_H1.07G021200
# Eg_RK687_H1.07G020800
# Eg_RK687_H1.07G020500
# Eg_RK561_H1.07G020000
# Eg_RK561_H1.07G019800
# Eg_RK561_H1.07G019200
# Eg_RK561_H1.07G018700
# Eg_RK498_H2.07G021100
# Eg_RK498_H2.07G021000
# Eg_RK498_H2.07G020900
# Eg_RK498_H2.07G020800
# Eg_RK498_H2.07G019200
# Eg_RK498_H2.07G019100
# Eg_RK498_H2.07G019000
# Eg_RK498_H2.07G018900
# Eg_RK498_H2.07G018800
# Eg_RK498_H2.07G018500
# Eg_RK498_H2.07G018400
# Eg_RK498_H2.07G018300
# Eg_RK498_H2.07G018200
# Eg_RK498_H2.07G018100
# Eg_RK498_H1.07G024300
# Eg_RK498_H1.07G024200
# Eg_RK436_H2.07G029000
# Eg_RK436_H2.07G028200
# Eg_RK436_H2.07G028000
# Eg_RK436_H1.07G020300
# Eg_RK388_H1.07G022700
# Eg_RK388_H1.07G022500
# Eg_RK320_H2.07G024500
# Eg_RK320_H2.07G023500
# Eg_RK320_H1.07G017700
# Eg_RK25_H2.07G020700
# Eg_RK25_H2.07G019800
# Eg_RK25_H2.07G019400
# Eg_RK25_H1.07G022700
# Eg_RK25_H1.07G021900
# Eg_RK25_H1.07G020700
# Eg_RK25_H1.07G020600
# Eg_RK25_H1.07G020400
# Eg_RK25_H1.07G020300
# Eg_RK25_H1.07G019800
# Eg_RK25_H1.07G019400
# Eg_RK05_H1.07G019800
# Eg_RK05_H1.07G019700

    #