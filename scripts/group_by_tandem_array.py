####program takes in a combbed file and provenance file as an argument and outputs a list of HOG IDs that contain at least one tandem array, as well as a detailed
#  summary of the tandem arrays, as well as the representative genes in each array per genome. The output files are written to the current
#  working directory.
from pathlib import Path
from sys import argv
import pandas as pd

combed_file = argv[1]
provenance_file = Path(__file__).with_name(
    "Pangenome_provenaces_edited.csv"
)

# HOG IDs that contain at least one tandem array
hog_output = "tandem_array_HOGs.txt"

# Detailed information about the tandem arrays
summary_output = "tandem_array_HOG_summary.csv"


# ============================================================
# READ COMBED FILE
# ============================================================

print("Reading combed file...")

df = pd.read_csv(
    combed_file,
    sep="\t",
    dtype=str
)

print(f"Genes read: {len(df)}")

provenance = pd.read_csv(
    provenance_file,
    sep=";",
    dtype=str,
    encoding="utf-8-sig",
)
provenance_columns = [
    "ProjectName",
    "Location",
    "Region",
    "Subpopulation",
]
missing_provenance_columns = [
    column for column in provenance_columns if column not in provenance.columns
]
if missing_provenance_columns:
    raise ValueError(
        "Provenance file is missing required columns: "
        f"{missing_provenance_columns}"
    )

for column in provenance_columns:
    provenance[column] = provenance[column].str.strip()

if provenance["ProjectName"].duplicated().any():
    duplicate_genomes = (
        provenance.loc[
            provenance["ProjectName"].duplicated(keep=False), "ProjectName"
        ]
        .drop_duplicates()
        .tolist()
    )
    raise ValueError(
        "Provenance file contains duplicate genome names: "
        f"{duplicate_genomes[:10]}"
    )

provenance_by_genome = provenance.set_index("ProjectName")[
    ["Location", "Region", "Subpopulation"]
]


# ============================================================
# CHECK COLUMNS
# ============================================================

required_columns = [
    "chr",
    "start",
    "end",
    "id",
    "ofID",
    "genome",
    "arrayID",
    "globHOG"
]

missing = [col for col in required_columns if col not in df.columns]

if missing:
    print("\nERROR: Missing required columns:")
    for col in missing:
        print(f"  {col}")
    raise SystemExit(1)


# ============================================================
# CLEAN COLUMNS
# ============================================================

for col in ["id", "genome", "arrayID", "globHOG"]:
    df[col] = df[col].astype(str).str.strip()

for column in ["start", "end"]:
    df[column] = pd.to_numeric(df[column], errors="coerce")
    if df[column].isna().any():
        raise ValueError(f"combBed contains missing or non-numeric {column} values.")


# ============================================================
# EXTRACT THE HOG ID
#
# globHOG looks like:
#
#     N0.HOG0000763 OG0000003
#
# We only want:
#
#     N0.HOG0000763
# ============================================================

df["HOG_ID"] = df["globHOG"].str.split().str[0]


# ============================================================
# COUNT GENES FOR EACH:
#
#     HOG + GENOME + ARRAY
#
# Example:
#
#     N0.HOG0000763 | RK05_HAP1 | Arr0000123
#
# ============================================================

array_counts = (
    df
    .groupby(
        ["HOG_ID", "genome", "arrayID"],
        as_index=False
    )
    .agg(
        genes_in_array=("id", "count")
    )
)


# ============================================================
# KEEP ONLY ARRAYS WITH >= 2 GENES
# ============================================================

tandem_arrays = array_counts[
    array_counts["genes_in_array"] >= 2
].copy()


print(
    f"Tandem arrays containing >=2 genes from the same HOG: "
    f"{len(tandem_arrays)}"
)


# ============================================================
# GET ALL HOG IDs WITH AT LEAST ONE TANDEM ARRAY
# ============================================================

tandem_hogs = (
    tandem_arrays["HOG_ID"]
    .drop_duplicates()
    .sort_values()
    .reset_index(drop=True)
)


# ============================================================
# WRITE HOG IDs
# ============================================================

tandem_hogs.to_csv(
    hog_output,
    index=False,
    header=False
)


# ============================================================
# ADD THE ACTUAL GENES TO THE SUMMARY
# ============================================================

# Join the tandem-array information back to the original
# combed file so we can see which genes caused each HOG
# to qualify.

tandem_gene_summary = df.merge(
    tandem_arrays[
        ["HOG_ID", "genome", "arrayID", "genes_in_array"]
    ],
    on=["HOG_ID", "genome", "arrayID"],
    how="inner"
)
tandem_gene_summary = tandem_gene_summary.join(
    provenance_by_genome,
    on="genome",
)

missing_metadata = tandem_gene_summary[
    ["Location", "Region", "Subpopulation"]
].isna().any(axis=1)
if missing_metadata.any():
    missing_genomes = (
        tandem_gene_summary.loc[missing_metadata, "genome"]
        .drop_duplicates()
        .tolist()
    )
    raise ValueError(
        "No complete provenance metadata found for genomes: "
        f"{missing_genomes[:10]}"
    )

array_keys = ["HOG_ID", "genome", "arrayID"]
chromosome_counts = tandem_gene_summary.groupby(array_keys)["chr"].nunique()
if chromosome_counts.gt(1).any():
    examples = chromosome_counts.index[chromosome_counts.gt(1)].tolist()[:10]
    raise ValueError(
        "A tandem array cannot span multiple chromosomes; "
        f"examples: {examples}"
    )


# ============================================================
# SORT OUTPUT
# ============================================================

tandem_gene_summary = tandem_gene_summary.sort_values(
    [
        "HOG_ID",
        "genome",
        "arrayID",
        "start"
    ]
)

tandem_gene_summary = (
    tandem_gene_summary.groupby(array_keys, sort=True, as_index=False)
    .agg(
        genes_in_array=("id", "nunique"),
        Location=("Location", "first"),
        Region=("Region", "first"),
        Subpopulation=("Subpopulation", "first"),
        chr=("chr", "first"),
        start=("start", "min"),
        end=("end", "max"),
        id=("id", lambda values: ", ".join(pd.unique(values))),
        ofID=("ofID", lambda values: ", ".join(pd.unique(values))),
        globOG=("globOG", "first"),
        globHOG=("globHOG", "first"),
    )
)

# ============================================================
# SELECT USEFUL COLUMNS
# ============================================================

output_columns = [
    "HOG_ID",
    "genome",
    "arrayID",
    "genes_in_array",
    "Location",
    "Region",
    "Subpopulation",
    "chr",
    "start",
    "end",
    "id",
    "ofID",
    "globOG",
    "globHOG"
]

tandem_gene_summary = tandem_gene_summary[
    output_columns
]


# ============================================================
# WRITE SUMMARY
# ============================================================

tandem_gene_summary.to_csv(
    summary_output,
    index=False
)


# ============================================================
# PRINT SUMMARY
# ============================================================

print("\nFinished.")
print(f"Unique tandem-array HOGs: {len(tandem_hogs)}")
print(f"HOG IDs written to: {hog_output}")
print(f"Detailed summary written to: {summary_output}")

print("\nFirst 20 tandem-array HOGs:")
print(tandem_hogs.head(20).to_string(index=False))

print("\nFirst 20 tandem arrays:")
print(
    tandem_arrays
    .head(20)
    .to_string(index=False)
)