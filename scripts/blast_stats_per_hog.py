# This file takes the Blast-hog-tandem-results.csv as input, as well as the 
# HOG_IDs as input, then it 
# dervives stats about each HOG, and outputs a summary table of the 
# stats per HOG. 
# Stats inlcude, different chromosomes, different protein lengths, 
# max, mix and average similiartiy difference#
# #
# #python3 .\blast_stats_per_hog.py .\Blast-hog-tandem-results.csv -o .\blast_stats_per_hog.csv
import argparse
import csv
import math
import re
import sys
from collections import defaultdict


def normalize(name):
    return re.sub(r"[^a-z0-9]", "", name.lower())


def find_column(fieldnames, aliases):
    fields = {normalize(name): index for index, name in enumerate(fieldnames)}
    for alias in aliases:
        index = fields.get(normalize(alias))
        if index is not None:
            return index
    return None


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def format_number(value):
    return "" if value is None else f"{value:.6g}"


def update_range(item, prefix, value):
    """Update min, max, sum, and count without storing every value."""
    if value is None:
        return

    count_key = f"{prefix}_count"
    item[f"{prefix}_min"] = min(item.get(f"{prefix}_min", value), value)
    item[f"{prefix}_max"] = max(item.get(f"{prefix}_max", value), value)
    item[f"{prefix}_sum"] = item.get(f"{prefix}_sum", 0.0) + value
    item[count_key] = item.get(count_key, 0) + 1


def main():
    parser = argparse.ArgumentParser(
        description="Calculate chromosome, protein-length, and similarity "
                    "statistics for each HOG in the input CSV."
    )
    parser.add_argument("blast_csv", help="BLAST-results CSV; HOG_ID must be first")
    parser.add_argument(
        "-o", "--output",
        default="blast_stats_per_hog.csv",
        help="Output CSV path (default: blast_stats_per_hog.csv)",
    )
    args = parser.parse_args()

    stats = defaultdict(lambda: {
        "rows": 0,
        "chromosomes": set(),
        "protein_lengths": set(),
    })
    hog_order = []

    with open(args.blast_csv, "r", encoding="utf-8-sig", newline="") as file:
        reader = csv.reader(file)
        headers = next(reader, None)

        if not headers:
            parser.error("The BLAST-results file is empty or has no header.")

        headers = [header.strip() for header in headers]
        # HOG_ID is taken from the first column of each data row.
        query_chrom_col = find_column(
            headers, ["query_Chromosome", "query_Chr", "q_Chromosome", "q_Chr"]
        )
        hit_chrom_col = find_column(
            headers, ["hit_Chromosome", "hit_Chr", "subject_Chromosome", "s_Chr"]
        )
        query_length_col = find_column(
            headers,
            ["query_Protein_Length", "query_Prot_Length", "query_Protein_Len",
             "query_Prot_Len"],
        )
        hit_length_col = find_column(
            headers,
            ["hit_Protein_Length", "hit_Prot_Length", "hit_Protein_Len",
             "hit_Prot_Len", "subject_Protein_Length"],
        )
        similarity_col = find_column(headers, ["pident", "percent_identity"])

        for row in reader:
            if not row:
                continue

            hog_id = row[0].strip()
            if not hog_id or normalize(hog_id) in {"hogid", "hog"}:
                continue

            if hog_id not in stats:
                hog_order.append(hog_id)

            item = stats[hog_id]
            item["rows"] += 1

            for column in (query_chrom_col, hit_chrom_col):
                if column is not None and column < len(row):
                    value = row[column].strip()
                    if value:
                        item["chromosomes"].add(value)

            query_length = (
                number(row[query_length_col])
                if query_length_col is not None and query_length_col < len(row)
                else None
            )
            hit_length = (
                number(row[hit_length_col])
                if hit_length_col is not None and hit_length_col < len(row)
                else None
            )

            for value in (query_length, hit_length):
                if value is not None:
                    item["protein_lengths"].add(value)

            if query_length is not None and hit_length is not None:
                update_range(
                    item, "length_difference",
                    abs(query_length - hit_length),
                )

            similarity = (
                number(row[similarity_col])
                if similarity_col is not None and similarity_col < len(row)
                else None
            )
            update_range(item, "similarity", similarity)

    output_columns = [
        "HOG_ID",
        "blast_rows",
        "different_chromosomes",
        "different_protein_lengths",
        "protein_length_difference_min",
        "protein_length_difference_max",
        "protein_length_difference_mean",
        "similarity_min",
        "similarity_max",
        "similarity_mean",
        "similarity_range",
    ]

    with open(args.output, "w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=output_columns)
        writer.writeheader()

        for hog_id in hog_order:
            item = stats[hog_id]
            diff_count = item.get("length_difference_count", 0)
            sim_count = item.get("similarity_count", 0)

            diff_min = item.get("length_difference_min")
            diff_max = item.get("length_difference_max")
            sim_min = item.get("similarity_min")
            sim_max = item.get("similarity_max")

            writer.writerow({
                "HOG_ID": hog_id,
                "blast_rows": item["rows"],
                "different_chromosomes": len(item["chromosomes"]),
                "different_protein_lengths": len(item["protein_lengths"]),
                "protein_length_difference_min": format_number(diff_min),
                "protein_length_difference_max": format_number(diff_max),
                "protein_length_difference_mean": format_number(
                    item["length_difference_sum"] / diff_count
                    if diff_count else None
                ),
                "similarity_min": format_number(sim_min),
                "similarity_max": format_number(sim_max),
                "similarity_mean": format_number(
                    item["similarity_sum"] / sim_count if sim_count else None
                ),
                "similarity_range": format_number(
                    sim_max - sim_min if sim_min is not None else None
                ),
            })

    if not query_chrom_col or not hit_chrom_col:
        print("Warning: chromosome columns were not found.", file=sys.stderr)
    if not query_length_col or not hit_length_col:
        print("Warning: query/hit protein-length columns were not found.",
              file=sys.stderr)
    if similarity_col is None:
        print("Warning: pident/similarity column was not found.",
              file=sys.stderr)

    print(f"Summarized {len(hog_order)} HOG(s) to {args.output}")


if __name__ == "__main__":
    main()