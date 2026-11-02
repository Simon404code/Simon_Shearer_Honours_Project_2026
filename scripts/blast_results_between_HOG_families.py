
#!/usr/bin/env python3

from sys import argv
import pandas as pd
import gzip
import re
from pathlib import Path
from collections import defaultdict


# ============================================================
# INPUTS
# ============================================================

# Tandem-array HOG summary produced by the previous script
HOG_FILE = argv[1]

# Directory containing Blast0_1.txt.gz, Blast0_2.txt.gz, etc.
BLAST_DIR = argv[2]

# Genome name <-> BLAST number mapping
GENOME_MAP_FILE = argv[3]

# Final output file type .csv.gz
OUTPUT_FILE = argv[4]


# ============================================================
# BLAST COLUMN NAMES
# ============================================================

BLAST_COLUMNS = [
    "Query_seq_ID",
    "Hit_ID",
    "percent_similarity",
    "Align_length",
    "Number_of_mismatches",
    "Gap_openings",
    "Query_start",
    "Query_end",
    "Subject_start",
    "Subject_end",
    "E-value",
    "Bit-score"
]


# ============================================================
# READ GENOME NUMBER MAPPING
# ============================================================

print("Reading genome-number mapping...")

genome_map = {}

with open(GENOME_MAP_FILE, "r") as f:

    for line in f:

        line = line.strip()

        if not line:
            continue

        parts = line.split()

        if len(parts) < 2:
            continue

        genome_name = parts[0]
        genome_number = int(parts[1])

        genome_map[genome_name] = genome_number


print(f"Genomes in mapping: {len(genome_map)}")


# ============================================================
# READ TANDEM HOG FILE
# ============================================================

print(f"\nReading {HOG_FILE}...")

hog_df = pd.read_csv(
    HOG_FILE,
    dtype=str
)

required_columns = [
    "HOG_ID",
    "genome",
    "id",
    "ofID"
]

missing = [
    col for col in required_columns
    if col not in hog_df.columns
]

if missing:

    print("\nERROR: Missing columns from HOG file:")

    for col in missing:
        print(f"  {col}")

    raise SystemExit(1)


# Remove whitespace
for col in ["HOG_ID", "genome", "id", "ofID"]:
    hog_df[col] = hog_df[col].astype(str).str.strip()


# ============================================================
# ADD BLAST GENOME NUMBER TO EACH GENE
# ============================================================

hog_df["genome_number"] = hog_df["genome"].map(genome_map)


missing_genomes = (
    hog_df.loc[
        hog_df["genome_number"].isna(),
        "genome"
    ]
    .drop_duplicates()
    .tolist()
)

if missing_genomes:

    print("\nERROR: These genomes are present in the HOG file")
    print("but are missing from the genome-number mapping:")

    for genome in missing_genomes:
        print(f"  {genome}")

    raise SystemExit(1)


hog_df["genome_number"] = (
    hog_df["genome_number"]
    .astype(int)
)


# ============================================================
# BUILD LOOKUP OF ALL OF_IDS BELONGING TO TANDEM HOGs
# ============================================================

# ofID -> information about the gene
#
# Example:
#
# 12_22397 ->
#     HOG_ID = N0.HOG0000001
#     genome = RK339_HAP1
#     genome_number = 12
#     gene_id = Eg_RK339_H1.07G016800

gene_lookup = {}

for _, row in hog_df.iterrows():

    ofid = row["ofID"]

    gene_lookup[ofid] = {
        "HOG_ID": row["HOG_ID"],
        "genome": row["genome"],
        "genome_number": row["genome_number"],
        "gene_id": row["id"]
    }


print(f"Unique of_IDs in tandem HOGs: {len(gene_lookup)}")


# ============================================================
# CREATE SET OF VALID HOG GENE PAIRS
# ============================================================

# For every HOG, create all pairs of genes.
#
# Example:
#
# HOG1:
#     gene A
#     gene B
#     gene C
#
# produces:
#
#     A -> B
#     A -> C
#     B -> A
#     B -> C
#     C -> A
#     C -> B
#
# The direction matters because BLAST files are directional.

print("Creating HOG gene-pair lookup...")

hog_pairs = set()

for hog, group in hog_df.groupby("HOG_ID"):

    genes = group["ofID"].unique().tolist()

    for query_id in genes:

        for subject_id in genes:

            if query_id == subject_id:
                continue

            hog_pairs.add(
                (query_id, subject_id)
            )


print(f"Gene pairs to search for: {len(hog_pairs)}")


# ============================================================
# FIND BLAST FILES
# ============================================================

blast_dir = Path(BLAST_DIR)

blast_files = sorted(
    blast_dir.glob("Blast*_*.txt.gz")
)

print(f"\nBLAST files found: {len(blast_files)}")

if len(blast_files) == 0:

    print(
        f"ERROR: No Blast*_*.txt.gz files found in "
        f"{blast_dir}"
    )

    raise SystemExit(1)


# ============================================================
# PROCESS BLAST FILES
# ============================================================

results = []

files_with_hits = 0
total_blast_lines = 0
matching_blast_lines = 0


for file_number, blast_file in enumerate(blast_files, start=1):

    # --------------------------------------------------------
    # Extract query and subject genome numbers
    #
    # Blast0_4.txt.gz
    #      ^ ^
    #      | |
    #      | subject = 4
    #      query = 0
    # --------------------------------------------------------

    match = re.match(
        r"Blast(\d+)_(\d+)\.txt\.gz$",
        blast_file.name
    )

    if not match:

        print(
            f"Skipping file with unexpected name: "
            f"{blast_file.name}"
        )

        continue

    query_genome_number = int(match.group(1))
    subject_genome_number = int(match.group(2))


    # --------------------------------------------------------
    # Determine whether this genome pair could contain
    # relevant HOG gene pairs.
    #
    # This prevents us from reading BLAST files where neither
    # genome contains genes from the tandem HOGs.
    # --------------------------------------------------------

    relevant_pairs = {
        pair
        for pair in hog_pairs
        if (
            gene_lookup[pair[0]]["genome_number"]
            == query_genome_number
            and
            gene_lookup[pair[1]]["genome_number"]
            == subject_genome_number
        )
    }

    if not relevant_pairs:
        continue


    print(
        f"[{file_number}/{len(blast_files)}] "
        f"{blast_file.name} "
        f"({len(relevant_pairs)} possible pairs)"
    )


    file_has_hit = False


    # --------------------------------------------------------
    # Read compressed BLAST file
    # --------------------------------------------------------

    with gzip.open(
        blast_file,
        "rt"
    ) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            parts = line.split()

            if len(parts) < 12:
                continue

            total_blast_lines += 1


            # ------------------------------------------------
            # Extract query and subject IDs
            # ------------------------------------------------

            query_id = parts[0]
            subject_id = parts[1]


            # ------------------------------------------------
            # Check whether this exact pair belongs to the
            # same HOG
            # ------------------------------------------------

            if (query_id, subject_id) not in relevant_pairs:
                continue


            matching_blast_lines += 1
            file_has_hit = True


            # ------------------------------------------------
            # Add HOG / genome information
            # ------------------------------------------------

            query_info = gene_lookup[query_id]
            subject_info = gene_lookup[subject_id]


            results.append({

                "HOG_ID":
                    query_info["HOG_ID"],

                "query_genome":
                    query_info["genome"],

                "query_genome_number":
                    query_genome_number,

                "query_gene":
                    query_info["gene_id"],

                "query_ofID":
                    query_id,

                "subject_genome":
                    subject_info["genome"],

                "subject_genome_number":
                    subject_genome_number,

                "subject_gene":
                    subject_info["gene_id"],

                "subject_ofID":
                    subject_id,

                "blast_file":
                    blast_file.name,

                "percent_similarity":
                    parts[2],

                "Align_length":
                    parts[3],

                "Number_of_mismatches":
                    parts[4],

                "Gap_openings":
                    parts[5],

                "Query_start":
                    parts[6],

                "Query_end":
                    parts[7],

                "Subject_start":
                    parts[8],

                "Subject_end":
                    parts[9],

                "E-value":
                    parts[10],

                "Bit-score":
                    parts[11]
            })


    if file_has_hit:
        files_with_hits += 1


# ============================================================
# WRITE OUTPUT
# ============================================================

print("\nFinished reading BLAST files.")

if len(results) == 0:

    print(
        "\nNo BLAST results were found between genes "
        "belonging to the same tandem-array HOG."
    )

    raise SystemExit(0)


result_df = pd.DataFrame(results)


# Sort in a useful order
result_df = result_df.sort_values(
    [
        "HOG_ID",
        "query_genome_number",
        "subject_genome_number",
        "query_ofID",
        "subject_ofID"
    ]
)


result_df.to_csv(
    OUTPUT_FILE,
    index=False
)


# ============================================================
# SUMMARY
# ============================================================

print("\n========================================")
print("SUMMARY")
print("========================================")

print(
    f"Tandem-array HOGs: "
    f"{hog_df['HOG_ID'].nunique()}"
)

print(
    f"Unique genes: "
    f"{hog_df['ofID'].nunique()}"
)

print(
    f"BLAST files found: "
    f"{len(blast_files)}"
)

print(
    f"BLAST files containing relevant hits: "
    f"{files_with_hits}"
)

print(
    f"BLAST lines examined: "
    f"{total_blast_lines}"
)

print(
    f"Matching BLAST lines: "
    f"{matching_blast_lines}"
)

print(
    f"Output rows: "
    f"{len(result_df)}"
)

print(
    f"\nOutput: {OUTPUT_FILE}"
)