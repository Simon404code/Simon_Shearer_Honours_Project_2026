#!/usr/bin/env python3
# python3 Group_genes_by_HOGs.py -o group_genes_by_HOGs.csv --provenance Pangenome_provenaces_edited.csv all_genes_in_tandem_HOGs.csv combBed_v5_complete-2.txt
"""Summarize genes, tandem arrays, and genome representation by HOG."""

import argparse
from pathlib import Path
import sys

import pandas as pd


DEFAULT_PROVENANCE = (
    Path(__file__).with_name("Pangenome_provenaces_edited.csv")
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize each HOG in all_genes_in_tandem_HOGs.csv, using "
            "combBed to calculate tandem-array statistics."
        )
    )
    parser.add_argument(
        "genes",
        help="CSV of all genes belonging to HOGs with tandem arrays",
    )
    parser.add_argument(
        "combbed",
        help="Full combBed TSV or TSV.GZ file",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="Group_genes_by_HOGs.csv",
        help="Output CSV path (default: Group_genes_by_HOGs.csv)",
    )
    parser.add_argument(
        "--provenance",
        default=str(DEFAULT_PROVENANCE),
        help=(
            "Semicolon-delimited genome metadata CSV with ProjectName, "
            "Location, Region, and Subpopulation columns "
            f"(default: {DEFAULT_PROVENANCE.name})"
        ),
    )
    return parser.parse_args()


def clean_strings(frame, columns, source_name):
    for column in columns:
        frame[column] = frame[column].astype("string").str.strip()
        frame[column] = frame[column].replace("", pd.NA)
        if frame[column].isna().any():
            raise ValueError(
                f"{source_name} contains blank or missing values "
                f"in {column}."
            )


def load_genes(path):
    genes = pd.read_csv(path, dtype={"HOG_ID": "string", "Gene_ID": "string"})
    required = {
        "HOG_ID",
        "Gene_ID",
        "Genome",
        "chr",
        "start",
        "end",
    }
    missing = sorted(required - set(genes.columns))
    if missing:
        raise ValueError(f"Gene CSV is missing required columns: {missing}")

    clean_strings(
        genes,
        ["HOG_ID", "Gene_ID", "Genome", "chr"],
        "Gene CSV",
    )
    if genes["Gene_ID"].duplicated().any():
        examples = genes.loc[
            genes["Gene_ID"].duplicated(keep=False), "Gene_ID"
        ].head(10).tolist()
        raise ValueError(
            f"Gene CSV contains duplicate Gene_ID values: {examples}"
        )

    for column in ("start", "end"):
        genes[column] = pd.to_numeric(genes[column], errors="coerce")
        if genes[column].isna().any():
            raise ValueError(
                f"Gene CSV contains missing or non-numeric {column} values."
            )

    genes["gene_size"] = (genes["end"] - genes["start"]).abs() + 1
    genes["gene_position"] = (genes["start"] + genes["end"]) / 2
    return genes


def load_combbed(path, hog_ids):
    combbed = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        low_memory=False,
        dtype={"globHOG": "string", "genome": "string", "arrayID": "string"},
    )
    required = {"globHOG", "genome", "arrayID", "id"}
    missing = sorted(required - set(combbed.columns))
    if missing:
        raise ValueError(f"combBed is missing required columns: {missing}")

    for column in ("globHOG", "genome", "arrayID", "id"):
        combbed[column] = combbed[column].astype("string").str.strip()
        combbed[column] = combbed[column].replace("", pd.NA)

    combbed = combbed.loc[
        combbed["globHOG"].isin(hog_ids)
        & combbed["genome"].notna()
        & combbed["arrayID"].notna()
        & combbed["id"].notna(),
        ["globHOG", "genome", "arrayID", "id"],
    ]
    return combbed


def load_provenance(path):
    metadata = pd.read_csv(path, sep=";", dtype="string")
    required = {"ProjectName", "Location", "Region", "Subpopulation"}
    missing = sorted(required - set(metadata.columns))
    if missing:
        raise ValueError(
            f"Provenance CSV is missing required columns: {missing}"
        )

    value_columns = ["Location", "Region", "Subpopulation"]
    clean_strings(
        metadata,
        ["ProjectName"] + value_columns,
        "Provenance CSV",
    )
    conflicts = metadata.groupby("ProjectName")[value_columns].nunique().gt(1)
    if conflicts.any().any():
        examples = conflicts.index[conflicts.any(axis=1)].tolist()[:10]
        raise ValueError(
            "Provenance CSV has conflicting location, region, or "
            "subpopulation values for genomes: "
            f"{examples}"
        )
    return metadata.drop_duplicates("ProjectName").set_index("ProjectName")[
        value_columns
    ]


def make_summary(genes, combbed, provenance):
    genome_metadata = genes[["HOG_ID", "Genome"]].drop_duplicates().copy()
    genome_metadata = genome_metadata.join(provenance, on="Genome")
    missing_genomes = sorted(
        genome_metadata.loc[
            genome_metadata[["Location", "Region", "Subpopulation"]]
            .isna()
            .any(axis=1),
            "Genome",
        ].unique()
    )
    if missing_genomes:
        raise ValueError(
            "No complete provenance metadata found for genomes: "
            f"{missing_genomes[:10]}"
        )

    hog_summary = (
        genes.groupby("HOG_ID", sort=True)
        .agg(
            total_number_of_genes=("Gene_ID", "nunique"),
            average_gene_size=("gene_size", "mean"),
            average_gene_position=("gene_position", "mean"),
            chromosomes=(
                "chr",
                lambda values: ", ".join(sorted(values.unique())),
            ),
            genomes=(
                "Genome",
                lambda values: ", ".join(sorted(values.unique())),
            ),
            number_of_genomes=("Genome", "nunique"),
        )
    )

    genome_gene_counts = (
        genes.groupby(["HOG_ID", "Genome"])["Gene_ID"]
        .nunique()
        .groupby("HOG_ID")
        .mean()
        .rename("average_genes_per_genome")
    )
    provenance_lists = (
        genome_metadata.groupby("HOG_ID")[
            ["Location", "Region", "Subpopulation"]
        ]
        .agg(lambda values: ", ".join(sorted(values.unique())))
        .rename(
            columns={
                "Location": "locations",
                "Region": "regions",
                "Subpopulation": "subpopulations",
            }
        )
    )

    array_sizes = (
        combbed.groupby(
            ["globHOG", "genome", "arrayID"],
            sort=True,
            dropna=True,
        )["id"]
        .nunique()
        .rename("array_size")
    )
    tandem_array_sizes = array_sizes.loc[array_sizes.ge(2)]
    array_summary = tandem_array_sizes.groupby("globHOG").agg(
        total_number_of_arrays="size",
        average_array_size="mean",
    )
    array_summary = array_summary.rename_axis("HOG_ID")

    summary = (
        hog_summary.join(genome_gene_counts)
        .join(provenance_lists)
        .join(array_summary)
    )
    summary["total_number_of_arrays"] = (
        summary["total_number_of_arrays"].fillna(0).astype(int)
    )
    summary["average_array_size"] = summary["average_array_size"].fillna(0)
    return summary.reset_index()[
        [
            "HOG_ID",
            "number_of_genomes",
            "genomes",
            "locations",
            "regions",
            "subpopulations",
            "total_number_of_genes",
            "total_number_of_arrays",
            "average_array_size",
            "average_genes_per_genome",
            "average_gene_size",
            "average_gene_position",
            "chromosomes",
        ]
    ]


def main():
    args = parse_args()
    genes = load_genes(args.genes)
    provenance = load_provenance(args.provenance)
    combbed = load_combbed(args.combbed, set(genes["HOG_ID"]))
    summary = make_summary(genes, combbed, provenance)
    summary.to_csv(args.output, index=False)
    print(f"HOGs summarized: {len(summary):,}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, pd.errors.ParserError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
