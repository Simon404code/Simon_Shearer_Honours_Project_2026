#!/usr/bin/env python3

import csv
import argparse
import os
from collections import defaultdict


def load_dictionaries(dictionary_directory):
    """
    Load all dictionary files from a directory.

    Each dictionary must be tab-delimited with:
        column 1 = original gene ID
        column 2 = translated gene ID

    The dictionary filename must contain the genome name.

    Example:
        RK05_HAP1_dictionary.tsv
        RK05_HAP2_dictionary.tsv
    """

    dictionaries = {}

    for filename in sorted(os.listdir(dictionary_directory)):

        dictionary_path = os.path.join(
            dictionary_directory,
            filename
        )

        if not os.path.isfile(dictionary_path):
            continue

        # Ignore unrelated files
        if not filename.endswith((".txt", ".tsv", ".csv")):
            continue

        # The genome name is taken from the filename.
        # This removes common dictionary suffixes.
        genome = filename

        if genome.endswith("_pre_to_post_mapping.tsv"):
            genome = genome[:-len("_pre_to_post_mapping.tsv")]
            ##this adds the _ before the HAP1 or HAP2
            genome = genome[:-4] + "_" + genome[-4:] 
        id_map = {}

        with open(
            dictionary_path,
            "r",
            encoding="utf-8-sig",
            newline=""
        ) as infile:

            reader = csv.reader(infile, delimiter="\t")

            for line_number, row in enumerate(reader, start=1):

                if len(row) < 2:
                    continue

                original_id = row[0].strip()
                translated_id = row[1].strip()

                # Skip blank rows
                if not original_id or not translated_id:
                    continue

                # Skip a header if present
                if line_number == 1 and original_id.lower() in {
                    "pre_id",
                    "geneid",
                    "gene_id",
                    "original_id"
                }:
                    continue

                if original_id in id_map:
                    if id_map[original_id] != translated_id:
                        print(
                            f"WARNING: conflicting mapping in "
                            f"{filename}, line {line_number}: "
                            f"{original_id} -> "
                            f"{id_map[original_id]} / "
                            f"{translated_id}"
                        )

                id_map[original_id] = translated_id

        dictionaries[genome] = id_map

        print(
            f"Loaded {len(id_map):,} mappings for genome "
            f"{genome} from {filename}"
        )

    print()
    print(f"Total dictionaries loaded: {len(dictionaries)}")

    return dictionaries


def find_dictionary_for_genome(genome, dictionaries):
    """
    Find the dictionary corresponding to the genome.

    First tries an exact filename-derived genome match.
    Then tries matching the genome against dictionary names.
    """

    if genome in dictionaries:
        return dictionaries[genome]

    # Also allow minor naming differences
    for dictionary_genome, id_map in dictionaries.items():

        if genome == dictionary_genome:
            return id_map

        if genome in dictionary_genome:
            return id_map

        if dictionary_genome in genome:
            return id_map

    return None


def translate_data_file(
    data_file,
    output_file,
    missing_file,
    dictionaries
):
    """
    Translate the gene ID in column 4 using the dictionary
    selected from the genome in column 8.

    Both columns are counted from 1:
        column 4 = gene ID
        column 8 = genome

    The file is processed top to bottom.
    Each row is translated once only.
    """

    total_rows = 0
    translated_rows = 0
    unchanged_rows = 0
    missing_dictionary_rows = 0

    missing_genes = []
    missing_genes_seen = set()

    genomes_without_dictionary = set()

    with open(
        data_file,
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as infile, open(
        output_file,
        "w",
        encoding="utf-8",
        newline=""
    ) as outfile:

        reader = csv.reader(infile, delimiter="\t")

        writer = csv.writer(
            outfile,
            delimiter="\t",
            lineterminator="\n"
        )

        # Copy the header unchanged
        header = next(reader)
        writer.writerow(header)

        for line_number, row in enumerate(reader, start=2):

            if not row:
                continue

            total_rows += 1

            # Column numbers are converted to Python indexes:
            # column 4 -> index 3
            # column 8 -> index 7
            if len(row) < 8:

                print(
                    f"WARNING: line {line_number} has fewer than "
                    f"8 columns. Row copied unchanged."
                )

                writer.writerow(row)
                unchanged_rows += 1
                continue

            original_gene_id = row[3]
            genome = row[7].strip()

            # Select the dictionary using the genome column
            id_map = find_dictionary_for_genome(
                genome,
                dictionaries
            )

            if id_map is None:

                print(
                    f"NO DICTIONARY: line {line_number}, "
                    f"genome={genome}, gene={original_gene_id}"
                )

                genomes_without_dictionary.add(genome)
                missing_dictionary_rows += 1

                missing_key = (
                    genome,
                    original_gene_id,
                    "NO_DICTIONARY"
                )

                if missing_key not in missing_genes_seen:
                    missing_genes.append(missing_key)
                    missing_genes_seen.add(missing_key)

                writer.writerow(row)
                unchanged_rows += 1
                continue

            # Exact lookup.
            # PAV. prefixes and suffixes are retained.
            if original_gene_id in id_map:

                translated_gene_id = id_map[original_gene_id]

                # Replace column 4 only
                row[3] = translated_gene_id

                translated_rows += 1

            else:

                print(
                    f"NOT FOUND: line {line_number}, "
                    f"genome={genome}, gene={original_gene_id}"
                )

                missing_key = (
                    genome,
                    original_gene_id,
                    "ID_NOT_FOUND"
                )

                if missing_key not in missing_genes_seen:
                    missing_genes.append(missing_key)
                    missing_genes_seen.add(missing_key)

                unchanged_rows += 1

            writer.writerow(row)

    # Write untranslated genes to a report
    with open(
        missing_file,
        "w",
        encoding="utf-8",
        newline=""
    ) as outfile:

        writer = csv.writer(
            outfile,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow([
            "genome",
            "geneID",
            "reason"
        ])

        for genome, gene_id, reason in missing_genes:
            writer.writerow([
                genome,
                gene_id,
                reason
            ])

    print()
    print("Translation finished")
    print("------------------------------")
    print(f"Rows processed:              {total_rows:,}")
    print(f"IDs translated:              {translated_rows:,}")
    print(f"IDs not translated:          {unchanged_rows:,}")
    print(f"Rows without dictionary:     {missing_dictionary_rows:,}")
    print(f"Unique missing gene entries: {len(missing_genes):,}")

    if genomes_without_dictionary:
        print()
        print("Genomes without dictionaries:")

        for genome in sorted(genomes_without_dictionary):
            print(f"  {genome}")

    print()
    print(f"Translated output: {output_file}")
    print(f"Missing-ID report:  {missing_file}")


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Translate gene IDs using the genome-specific dictionary "
            "selected from column 8."
        )
    )

    parser.add_argument(
        "--data",
        required=True,
        help="Input tab-delimited data file"
    )

    parser.add_argument(
        "--dictionaries",
        required=True,
        help="Directory containing the genome-specific dictionaries"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Translated tab-delimited output file"
    )

    parser.add_argument(
        "--missing",
        default=None,
        help=(
            "Report of genes not translated. "
            "Defaults to OUTPUT.missing.tsv"
        )
    )

    args = parser.parse_args()

    if args.missing is None:
        missing_file = args.output + ".missing.tsv"
    else:
        missing_file = args.missing

    if not os.path.isfile(args.data):
        raise FileNotFoundError(
            f"Data file not found: {args.data}"
        )

    if not os.path.isdir(args.dictionaries):
        raise NotADirectoryError(
            f"Dictionary directory not found: {args.dictionaries}"
        )

    dictionaries = load_dictionaries(args.dictionaries)

    if not dictionaries:
        raise RuntimeError(
            "No dictionary files were found in the supplied directory."
        )

    translate_data_file(
        data_file=args.data,
        output_file=args.output,
        missing_file=missing_file,
        dictionaries=dictionaries
    )


if __name__ == "__main__":
    main()