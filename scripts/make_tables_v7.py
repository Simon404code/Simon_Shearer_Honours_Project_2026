#!/usr/bin/env python3
#$ py -3 ./make_tables_v7.py ./all_genes_in_tandem_HOGs_2.csv ./combBed_v5_complete-2.txt tandem --n-genomes 50

"""Generate summary tables from all genes in tandem HOGs and a combBed file."""

import argparse
import sys
import pandas as pd


TABLE1_COLUMNS = [
    "HOG_ID",
    "Gene_ID",
    "Genome",
    "array_ID",
    "ofID",
    "isArrayRep",
    "chr",
    "start",
    "end",
]

COMBBED_COLUMN_ALIASES = {
    "Gene_ID": ["id", "Gene_ID", "geneID", "gene_id"],
    "HOG_ID": ["globHOG", "HOG_ID", "hog", "HOG"],
    "Genome": ["Genome", "genome"],
    "arrayID": ["arrayID", "ArrayID", "array_id"],
    "chr": ["chr", "chromosome", "chrom"],
    "start": ["start", "Start"],
    "end": ["end", "End"],
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Create Tables 2 and 3 from all genes in tandem HOGs, and "
            "Table 4 from combBed."
        )
    )
    parser.add_argument(
        "table1",
        help="CSV of all genes belonging to HOGs with tandem arrays",
    )
    parser.add_argument("combbed", help="combBed TSV or TSV.GZ file")
    parser.add_argument(
        "output_prefix",
        nargs="?",
        default="tandem",
        help="Prefix for generated files (default: tandem)",
    )
    parser.add_argument(
        "--n-genomes",
        type=int,
        default=50,
        help="Expected number of genomes; used for a warning, not as the output denominator (default: 50)",
    )
    return parser.parse_args()


def find_column(df, possible_names, description):
    """Find a column by trying known source-file names in order."""
    for name in possible_names:
        if name in df.columns:
            return name

    raise ValueError(
        f"Could not find {description}. Available columns: {list(df.columns)}"
    )


def clean_string_columns(df, columns):
    """Trim text identifiers while preserving missing values as missing."""
    for column in columns:
        df[column] = df[column].astype("string").str.strip()
        df[column] = df[column].replace("", pd.NA)


def load_table1(path):
    table1 = pd.read_csv(path)
    missing = [column for column in TABLE1_COLUMNS if column not in table1.columns]
    if missing:
        raise ValueError(f"Table 1 is missing required columns: {missing}")

    clean_string_columns(
        table1,
        ["HOG_ID", "Gene_ID", "Genome", "array_ID", "ofID", "chr"],
    )

    for column in ["HOG_ID", "Gene_ID", "Genome", "array_ID"]:
        if table1[column].isna().any():
            raise ValueError(f"Table 1 contains blank or missing values in {column}.")

    duplicate_genes = table1["Gene_ID"].duplicated(keep=False)
    if duplicate_genes.any():
        examples = table1.loc[duplicate_genes, "Gene_ID"].head(10).tolist()
        raise ValueError(
            "Gene_ID must be unique in Table 1; duplicate examples: "
            f"{examples}"
        )

    table1["start"] = pd.to_numeric(table1["start"], errors="coerce")
    table1["end"] = pd.to_numeric(table1["end"], errors="coerce")
    table1["gene_length"] = (table1["end"] - table1["start"] + 1).abs()
    return table1


def make_table2(table1):
    table2 = (
        table1.groupby("HOG_ID", sort=True)
        .agg(
            list_of_gene_IDs=("Gene_ID", lambda values: ";".join(values)),
            number_of_Genes=("Gene_ID", "nunique"),
            number_of_Genomes=("Genome", "nunique"),
            unique_chromosomes=("chr", "nunique"),
            avg_gene_length=("gene_length", "mean")
        )
        .reset_index()
    )
    per_hog_genome_gene_counts = (
        table1.groupby(["HOG_ID", "Genome"])["Gene_ID"].nunique()
    )
    tandem_array_rates = (
        per_hog_genome_gene_counts.groupby("HOG_ID")
        .apply(lambda gene_counts: gene_counts.ge(2).mean() * 100)
        .rename("tandem_array_rate_percent")
        .reset_index()
    )
    table2 = table2.merge(tandem_array_rates, on="HOG_ID", how="left")

    if table2["number_of_Genomes"].gt(50).any():
        raise ValueError(
            "Table 2 classification expects at most 50 genomes per HOG."
        )

    def classify_hog(genome_count):
        if genome_count == 1:
            return "singleton"
        if genome_count <= 24:
            return "minority shell"
        if genome_count <= 48:
            return "majority shell"
        if genome_count == 49:
            return "soft-core"
        return "core"

    table2["classifier"] = table2["number_of_Genomes"].map(classify_hog)
    return table2


def load_combbed(path):
    combbed = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        low_memory=False,
    )

    rename_columns = {
        find_column(combbed, aliases, description): standard_name
        for standard_name, aliases in COMBBED_COLUMN_ALIASES.items()
        for description in [standard_name]
    }
    combbed = combbed.rename(columns=rename_columns)

    clean_string_columns(
        combbed,
        ["Gene_ID", "HOG_ID", "Genome", "arrayID", "chr"],
    )
    for column in ["Gene_ID", "HOG_ID", "Genome"]:
        if combbed[column].isna().any():
            # Rows without a usable identifier cannot participate in matching
            # or per-HOG/per-genome summaries.
            combbed = combbed.dropna(subset=[column])

    duplicate_genes = combbed["Gene_ID"].duplicated(keep=False)
    if duplicate_genes.any():
        examples = combbed.loc[duplicate_genes, "Gene_ID"].head(10).tolist()
        raise ValueError(
            "combBed Gene_ID values must be unique for the gene-based merge; "
            f"duplicate examples: {examples}"
        )

    return combbed


def make_table3(table1, expected_genomes):
    n_genomes = table1["Genome"].nunique()
    if n_genomes == 0:
        raise ValueError("The tandem-HOG CSV contains no usable genome identifiers.")
    if expected_genomes > 0 and n_genomes != expected_genomes:
        print(
            f"WARNING: the tandem-HOG CSV contains {n_genomes} genomes, "
            f"not the expected {expected_genomes}."
        )

    per_hog_array = (
        table1.groupby(["HOG_ID", "Genome", "array_ID"])["Gene_ID"]
        .nunique()
        .rename("array_copy_count")
        .reset_index()
    )
    tandem_arrays = per_hog_array.loc[
        per_hog_array["array_copy_count"] > 1
    ]

    frequency = (
        tandem_arrays.groupby("HOG_ID")["Genome"]
        .nunique()
        .rename("frequency_count")
    )
    cnv_stats = (
        tandem_arrays["array_copy_count"]
        .agg(min_cnv="min", max_cnv="max", avg_CNV="mean")
    )
    duplicate_genomes = frequency.rename("duplicate_Genomes")

    length_stats = (
        table1.groupby("HOG_ID")["gene_length"]
        .agg(
            min_gene_length="min",
            max_gene_length="max",
            avg_gene_length="mean",
        )
    )

    table3 = pd.DataFrame(index=sorted(table1["HOG_ID"].unique()))
    table3.index.name = "HOG_ID"
    table3 = table3.join(cnv_stats).join(frequency).join(duplicate_genomes).join(
        length_stats
    )
    table3["frequency_count"] = table3["frequency_count"].fillna(0).astype(int)
    table3["duplicate_Genomes"] = (
        table3["duplicate_Genomes"].fillna(0).astype(int)
    )
    table3["min_cnv"] = table3["min_cnv"].fillna(0)
    table3["max_cnv"] = table3["max_cnv"].fillna(0)
    table3["avg_CNV"] = table3["avg_CNV"].fillna(0)

    def classify_type(frequency_count):
        if frequency_count == n_genomes:
            return "Core"
        if frequency_count == n_genomes - 1:
            return "Soft-core"
        return "Shell"

    table3["Type"] = table3["frequency_count"].map(classify_type)
    table3["HOG_ID(array)"] = table3.index
    table3["CNV min"] = table3["min_cnv"].astype(int).astype(str)
    table3["CNV max"] = table3["max_cnv"].astype(int).astype(str)
    table3["frequency"] = (
        table3["frequency_count"].astype(str) + f"/{n_genomes}"
    )
    table3["duplicate_rate"] = (
        table3["duplicate_Genomes"].astype(str) + f"/{n_genomes}"
    )

    output = table3[
        [
            "HOG_ID(array)",
            "Type",
            "CNV min",
            "CNV max",
            "avg_CNV",
            "frequency",
            "duplicate_rate",
            "min_gene_length",
            "max_gene_length",
            "avg_gene_length",
        ]
    ].reset_index(drop=True)

    return output


def make_table4(table1, combbed):
    table1_genes = set(table1["Gene_ID"].dropna())
    matched = table1_genes.intersection(combbed["Gene_ID"])

    print("Genes in Table 1:", len(table1_genes))
    print("Genes found in combBed:", len(matched))
    print("Genes NOT found:", len(table1_genes - matched))
    if not matched:
        raise ValueError(
            "No Table 1 Gene_ID values matched combBed. "
            f"Table 1 examples: {table1['Gene_ID'].head(10).tolist()}; "
            f"combBed examples: {combbed['Gene_ID'].head(10).tolist()}"
        )

    combbed_for_merge = combbed.loc[
        combbed["Gene_ID"].isin(table1_genes),
        ["Gene_ID", "HOG_ID", "Genome", "arrayID"],
    ].rename(
        columns={
            "HOG_ID": "HOG_ID_combbed",
            "Genome": "Genome_combbed",
        }
    )
    merged = table1.merge(
        combbed_for_merge,
        on="Gene_ID",
        how="left",
        validate="one_to_one",
    )

    array_counts = (
        merged.dropna(
            subset=["HOG_ID_combbed", "Genome_combbed", "arrayID"]
        )
        .groupby(
            ["HOG_ID_combbed", "Genome_combbed", "arrayID"],
            sort=True,
            dropna=True,
        )["Gene_ID"]
        .nunique()
        .rename("array_gene_count")
        .reset_index()
    )
    tandem_arrays = array_counts.loc[
        array_counts["array_gene_count"] >= 2
    ].copy()

    array_keys = ["Genome_combbed", "arrayID", "HOG_ID_combbed"]
    array_members = merged.merge(
        tandem_arrays[array_keys],
        on=array_keys,
        how="inner",
        validate="many_to_one",
    )
    missing_location = array_members[["chr", "start", "end"]].isna().any(axis=1)
    if missing_location.any():
        examples = array_members.loc[
            missing_location, array_keys + ["Gene_ID"]
        ].head(10)
        raise ValueError(
            "Tandem-array genes must have chromosome, start, and end "
            f"coordinates; examples: {examples.to_dict(orient='records')}"
        )

    chromosome_counts = array_members.groupby(array_keys, sort=True)["chr"].nunique()
    if chromosome_counts.gt(1).any():
        examples = chromosome_counts.loc[chromosome_counts.gt(1)].head(10)
        raise ValueError(
            "A tandem array cannot span multiple chromosomes; "
            f"examples: {examples.index.tolist()}"
        )

    table4 = (
        array_members.groupby(array_keys, sort=True)
        .agg(
            gene_copy_number=("Gene_ID", "nunique"),
            start=("start", "min"),
            end=("end", "max"),
            chr=("chr", "first"),
        )
        .reset_index()
        .rename(
            columns={
                "Genome_combbed": "Genome",
                "HOG_ID_combbed": "HOG_ID",
                "arrayID": "Array_ID",
            }
        )
    )

    return tandem_arrays, table4


def main():
    args = parse_args()

    if args.n_genomes < 0:
        raise ValueError("--n-genomes must be zero or greater.")

    print("\n========================================")
    print("TANDEM TABLE GENERATION")
    print("========================================\n")
    print("Table 1:", args.table1)
    print("combBed:", args.combbed)
    print("Output prefix:", args.output_prefix)
    print()

    print("Reading Table 1...")
    table1 = load_table1(args.table1)
    print("Table 1 loaded.")
    print("Rows:", len(table1))
    print("Columns:", list(table1.columns))
    print("Unique HOGs in Table 1:", table1["HOG_ID"].nunique())
    print("Unique genes in Table 1:", table1["Gene_ID"].nunique())
    print("Unique Genomes in Table 1:", table1["Genome"].nunique())
    print()

    print("Reading combBed...")
    combbed = load_combbed(args.combbed)
    print("combBed loaded.")
    print("Rows:", len(combbed))
    print("Columns:", list(combbed.columns))
    print()

    table3 = make_table3(table1, args.n_genomes)
    tandem_arrays, table4 = make_table4(table1, combbed)
    table3_file = f"{args.output_prefix}_table3.csv"
    table3.to_csv(table3_file, index=False)
    table4_file = f"{args.output_prefix}_table4.csv"
    table4.to_csv(table4_file, index=False)
    print("Table 4 created:", table4_file)
    print("Rows:", len(table4))
    print()

    table2 = make_table2(table1).merge(
        table3[["HOG_ID(array)", "avg_CNV"]].rename(
            columns={"HOG_ID(array)": "HOG_ID"}
        ),
        on="HOG_ID",
        how="left",
        validate="one_to_one",
    )
    table2_file = f"{args.output_prefix}_table2.csv"
    table2.to_csv(table2_file, index=False)
    print("Table 2 created:", table2_file)
    print("Rows:", len(table2))
    print()

    merged_file = f"{args.output_prefix}_merged_debug.csv"
    debug_merge = table1.merge(
        combbed.loc[
            combbed["Gene_ID"].isin(set(table1["Gene_ID"])),
            ["Gene_ID", "HOG_ID", "Genome", "arrayID"],
        ].rename(
            columns={
                "HOG_ID": "HOG_ID_combbed",
                "Genome": "Genome_combbed",
            }
        ),
        on="Gene_ID",
        how="left",
        validate="one_to_one",
    )
    debug_merge.to_csv(merged_file, index=False)

    print("Table 3 created:", table3_file)
    print("Rows:", len(table3))
    print()
    print("========================================")
    print("DONE")
    print("========================================")
    print("Table 2:", table2_file)
    print("Table 3:", table3_file)
    print("Table 4:", table4_file)
    print("Debug:", merged_file)
    print("Table 2 HOGs:", len(table2))
    print("Table 3 HOGs:", len(table3))
    print("Qualifying tandem arrays:", len(tandem_arrays))
    print("Table 4 arrays:", len(table4))


if __name__ == "__main__":
    try:
        main()
    except (OSError, pd.errors.ParserError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
