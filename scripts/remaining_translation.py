"This file will isolate all the genes that haven't been translated, Ientified by their 'm.g'"
"in their ID. It will then find the corresponding translated ID in one of the dictionaries "
"It will then replace that geneID with the coresponding translated ID. If the geneID is not "
"found in any of the dictionaries, it will be left unchanged, stats about thr translation will be"
"kept, for skipped genes, translated genes, untranslated genes."
import argparse
import csv


def main():
    parser = argparse.ArgumentParser(
        description="Replace gene IDs in column 4 using one combined dictionary."
    )

    parser.add_argument(
        "--data",
        required=True,
        help="Input data TSV file"
    )

    parser.add_argument(
        "--dictionaries",
        required=True,
        help="TSV file containing the dictionaries"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file"
    )

    args = parser.parse_args()

    print("Loading dictionaries...")
    id_map = load_dictionaries(args.dictionaries)

    print()
    print("Processing data file...")

    replace_ids(
        args.data,
        args.output,
        id_map
    )

def replace_ids(data_file, output_file, id_map):
    """Replace gene IDs in column 4 of the data file using the provided ID map."""
    total = 0
    replaced = 0
    unchanged = 0
    pav_count = 0
    pav_replaced = 0

    with open(data_file, "r", encoding="utf-8") as infile, open(
        output_file, "w", encoding="utf-8", newline=""
    ) as outfile:

        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")

        for row in reader:
            total += 1
            original_id = row[3]

            if original_id.startswith("PAV."):
                pav_count += 1
                translated_id = id_map.get(original_id[4:])
                if translated_id:
                    row[3] = "PAV." + translated_id
                    pav_replaced += 1
                else:
                    row[3] = original_id  # Retain original if not found
            else:
                translated_id = id_map.get(original_id)
                if translated_id:
                    row[3] = translated_id
                    replaced += 1
                else:
                    unchanged += 1  # Retain original if not found

            writer.writerow(row)

def load_dictionaries(dictionaries_file):
    """Load multiple TSV dictionary files into a single dictionary."""
    id_map = {}
    with open(dictionaries_file, "r", encoding="utf-8") as dict_file:
        for line in dict_file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue  # Skip empty lines and comments
            parts = line.split("\t")
            if len(parts) < 2:
                continue  # Skip lines that don't have at least two columns
            original_id, translated_id = parts[0], parts[1]
            id_map[original_id] = translated_id


if __name__ == "__main__":

    main()
