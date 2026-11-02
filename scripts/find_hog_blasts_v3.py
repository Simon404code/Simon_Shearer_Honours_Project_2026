#!/usr/bin/env python3
#python find_hog_blasts.py all_genes_in_tandem_HOGs.csv /path/to/blast /path/to/translation.txt
# python3 find_hog_blasts.py all_genes_in_tandem_HOGs.csv ../../ SpeciesIDs_from_50_run.txt --output tandem_hog_blast_results.csv

import argparse
import csv
import gzip
import os
import re
from collections import defaultdict


def normalise_genome_name(name):
    """
    Make genome names easier to match.

    For example:

        RM01_HAP1
        RM01HAP1
        rm01hap1

    are converted to comparable uppercase strings.
    """

    return re.sub(r"[_\-\s]", "", str(name).upper())


def read_translation_file(filename):
    """
    Read a genome translation file such as:

        TAG0014Hap1 - 0
        RM01Hap1 - 1

    Returns:

        {
            "TAG0014HAP1": "0",
            "RM01HAP1": "1"
        }
    """

    genome_to_number = {}

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Accept colon separated entries like:
            # 0: TAG0014Hap1.fa
            # 1: RM01Hap1.fa

            match = re.match(
                r"^(.+?):(.+?)\s*$",
                line
            )

            if not match:
                print(
                    f"WARNING: Could not parse translation line: {line}"
                )
                continue

            number = match.group(1).strip()
            genome = match.group(2).strip()

            root, ext = os.path.splitext(genome)
            if ext.lower() in {".fa", ".fasta"}:
                genome = root

            genome_to_number[normalise_genome_name(genome)] = number

    return genome_to_number


def read_tandem_file(filename):
    """
    Read the CSV produced by extract_tandem_genes.py.
    Returns:
        genes_by_hog[HOG_ID] = list of gene records
    """

    genes_by_hog = defaultdict(list)

    with open(filename, "r", newline="") as f:

        reader = csv.DictReader(f)

        required_columns = {
            "HOG_ID",
            "Gene_ID",
            "Genome",
            "ofID",
            "isArrayRep",
            "chr",
            "start",
            "end"
        }

        missing = required_columns - set(reader.fieldnames)

        if missing:
            raise ValueError(
                f"Missing columns in tandem gene file: {', '.join(sorted(missing))}"
            )

        for row in reader:
            #takes the HOG_ID from the read in row and appends the Gene and the row to that HOG_ID
            genes_by_hog[row["HOG_ID"]].append(row)

    return genes_by_hog


def find_blast_file(blast_dir, query_number, subject_number):

    """
    Find a BLAST file corresponding to:

        BlastX_Y

    where X = query genome
    and Y = subject genome.

    The function allows common extensions such as:
        .txt
        .tsv
        .blast
        .out
    """

    possible_names = [
        f"Blast{query_number}_{subject_number}.txt",
        f"Blast{query_number}_{subject_number}.txt.gz"
    ]

    for name in possible_names:

        path = os.path.join(blast_dir, name)

        if os.path.isfile(path):
            return path

    return None


def read_blast_file(filename):
    """
    Read BLAST tabular output and index hits by query gene and target gene.

    Returns:

        {
            qseqid: {
                sseqid: [fields[:12]]
            }
        }

    so each gene pair can be looked up in constant time instead of scanning
    the entire BLAST file for every comparison.
    """

    blast_index = {}

    opener = gzip.open if filename.endswith(".gz") else open

    with opener(filename, "rt") as f:

        for line_number, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            fields = line.split("\t")

            if len(fields) < 12:
                print(
                    f"WARNING: {filename}, line {line_number} "
                    f"has only {len(fields)} columns"
                )
                continue

            blast_row = fields[:12]
            query_id = blast_row[0]
            hit_id = blast_row[1]

            blast_index.setdefault(query_id, {})
            blast_index[query_id][hit_id] = blast_row

    return blast_index


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Find BLAST results between every pair of genes "
            "belonging to the same HOG."
        )
    )

    parser.add_argument(
        "tandem_genes",
        help="CSV produced by extract_tandem_genes.py"
    )

    parser.add_argument(
        "blast_directory",
        help="Directory containing BlastX_Y files"
    )

    parser.add_argument(
        "translation_file",
        help="Genome-to-BLAST-number translation file"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="tandem_hog_blast_results.csv",
        help="Output CSV file"
    )

    args = parser.parse_args()

    # ------------------------------------------------------------
    # Read genome translation
    # Normalise names so that they don't have .fa or .fasta extensions and are uppercase
    # So that the genome names match the blast annoations.
    # ------------------------------------------------------------

    genome_to_number_raw = read_translation_file(
        args.translation_file
    )

    genome_to_number = {
        normalise_genome_name(genome): number
        for genome, number in genome_to_number_raw.items()
    }

    print(
        f"Genome translations loaded: {len(genome_to_number)}"
    )

    # ------------------------------------------------------------
    # Read tandem genes
    # ------------------------------------------------------------

    genes_by_hog = read_tandem_file(
        args.tandem_genes
    )

    print(
        f"HOGs containing tandem-array genes: "
        f"{len(genes_by_hog):,}"
    )

    # Index requested gene pairs by BLAST file and gene IDs.
    # ...existing code...

    requests = defaultdict(lambda: defaultdict(list))
    warned_genomes = set()

    for hog_id, genes in genes_by_hog.items():
        for query_gene in genes:
            query_key = normalise_genome_name(query_gene["Genome"])
            query_number = genome_to_number.get(query_key)

            if query_number is None:
                if query_key not in warned_genomes:
                    print(f"WARNING: No BLAST number for genome {query_gene['Genome']}")
                    warned_genomes.add(query_key)
                continue

            for hit_gene in genes:
                hit_key = normalise_genome_name(hit_gene["Genome"])
                hit_number = genome_to_number.get(hit_key)

                if hit_number is None:
                    if hit_key not in warned_genomes:
                        print(f"WARNING: No BLAST number for genome {hit_gene['Genome']}")
                        warned_genomes.add(hit_key)
                    continue

                pair = (query_gene["ofID"], hit_gene["ofID"])
                requests[(query_number, hit_number)][pair].append(
                    (hog_id, query_gene, hit_gene)
                )

# ...existing code...
    # ------------------------------------------------------------
    # Output columns
    # ------------------------------------------------------------

    output_columns = [
        "HOG_ID",

        "query_Gene_ID",
        "query_Genome",
        "query_ofID",

        "hit_Gene_ID",
        "hit_Genome",
        "hit_ofID",

        "query_blast_genome",
        "hit_blast_genome",

        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore"
    ]

    
    n_pairs = sum(
    len(records)
    for pair_requests in requests.values()
    for records in pair_requests.values()
    )
    n_blast_hits = 0
    n_blast_misses = 0
    n_missing_files = 0
    n_files_processed = 0

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=output_columns)
        writer.writeheader()

        for (query_number, hit_number), pair_requests in requests.items():
            blast_file = find_blast_file(
                args.blast_directory, query_number, hit_number
            )

            if blast_file is None:
                print(f"WARNING: BLAST file not found: Blast{query_number}_{hit_number}")
                n_missing_files += sum(
                    len(records) for records in pair_requests.values()
                )
                continue

            wanted = defaultdict(set)
            for query_id, hit_id in pair_requests:
                wanted[query_id].add(hit_id)

            print(f"Reading BLAST file: {blast_file}")
            opener = gzip.open if blast_file.endswith(".gz") else open
            found = {}

            # Stream the file; retain only requested hits.
            with opener(blast_file, "rt") as blast:
                for line_number, line in enumerate(blast, start=1):
                    if not line.strip() or line.startswith("#"):
                        continue

                    fields = line.rstrip("\r\n").split("\t")
                    if len(fields) < 12:
                        print(
                            f"WARNING: {blast_file}, line {line_number} "
                            f"has only {len(fields)} columns"
                        )
                        continue

                    query_id, hit_id = fields[0], fields[1]
                    if hit_id in wanted.get(query_id, ()):
                        # Preserve the original behavior: the last duplicate hit wins.
                        found[(query_id, hit_id)] = fields[:12]

            n_files_processed += 1

            for pair, records in pair_requests.items():
                blast_row = found.get(pair)
                if blast_row is None:
                    n_blast_misses += len(records)
                    continue

                for hog_id, query_gene, hit_gene in records:
                    writer.writerow({
                        "HOG_ID": hog_id,
                        "query_Gene_ID": query_gene["Gene_ID"],
                        "query_Genome": query_gene["Genome"],
                        "query_ofID": query_gene["ofID"],
                        "hit_Gene_ID": hit_gene["Gene_ID"],
                        "hit_Genome": hit_gene["Genome"],
                        "hit_ofID": hit_gene["ofID"],
                        "query_blast_genome": query_number,
                        "hit_blast_genome": hit_number,
                        "pident": blast_row[2],
                        "length": blast_row[3],
                        "mismatch": blast_row[4],
                        "gapopen": blast_row[5],
                        "qstart": blast_row[6],
                        "qend": blast_row[7],
                        "sstart": blast_row[8],
                        "send": blast_row[9],
                        "evalue": blast_row[10],
                        "bitscore": blast_row[11],
                    })
                    n_blast_hits += 1

# ...existing code...
    # ------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------
    print()
    print("Finished.")
    print(f"Gene pairs examined: {n_pairs:,}")
    print(f"BLAST matches written: {n_blast_hits:,}")
    print(f"BLAST searches with no hits: {n_blast_misses:,}")
    print(f"Missing BLAST file accesses: {n_missing_files:,}")
    print(f"BLAST files processed: {n_files_processed:,}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()