#!/usr/bin/env python3
"""
Remove reads from a FASTQ file using Kraken-classified taxon assignments.

This script extracts read IDs assigned to one or more target taxa from a
Kraken classifications file and removes those reads from an input FASTQ file.

It is intended for cases such as background clean-up, where a supposedly
depleted read set still contains reads classified to a taxonomic group of
interest.

Outputs
-------
All tabular outputs are written as tab-separated files.

1. A text file containing one read ID per line for all removed reads.
2. A tab-separated summary of removed read counts per taxon.
3. A filtered FASTQ file containing only reads not selected for removal.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Remove reads from a FASTQ file based on Kraken taxonomic "
            "assignments."
        )
    )
    parser.add_argument(
        "--kraken_classifications",
        required=True,
        help="Path to Kraken classifications TSV.",
    )
    parser.add_argument(
        "--kraken_report",
        required=True,
        help=(
            "Path to Kraken report TSV. This is used to map taxids to taxon "
            "names and to select taxa by pattern."
        ),
    )
    parser.add_argument(
        "--input_fastq",
        required=True,
        help="Path to the input FASTQ file to filter.",
    )
    parser.add_argument(
        "--output_fastq",
        required=True,
        help="Path to the filtered FASTQ file.",
    )
    parser.add_argument(
        "--read_ids_out",
        required=True,
        help="Path to write removed read IDs, one per line.",
    )
    parser.add_argument(
        "--summary_tsv",
        required=True,
        help="Path to write a tab-separated summary of removed reads.",
    )
    parser.add_argument(
        "--taxon_pattern",
        default="Plasmodium",
        help=(
            "Case-sensitive substring used to select taxa from the Kraken "
            "report. Default: Plasmodium"
        ),
    )
    return parser.parse_args()


def load_kraken_report(report_path: Path, taxon_pattern: str) -> pd.DataFrame:
    """
    Load the Kraken report and keep taxa matching the requested pattern.

    Parameters
    ----------
    report_path : Path
        Path to the Kraken report TSV.
    taxon_pattern : str
        Substring used to select taxa of interest.

    Returns
    -------
    pd.DataFrame
        DataFrame of matching taxa.
    """
    report_cols = [
        "percent",
        "clade_reads",
        "direct_reads",
        "rank_code",
        "taxid",
        "name",
    ]
    report_df = pd.read_csv(
        report_path,
        sep="\t",
        names=report_cols,
        header=None,
        dtype={"taxid": str},
    )
    report_df["taxid"] = report_df["taxid"].astype(str).str.strip()
    report_df["name"] = report_df["name"].astype(str).str.strip()
    report_df["rank_code"] = report_df["rank_code"].astype(str).str.strip()
    report_df["clade_reads"] = pd.to_numeric(
        report_df["clade_reads"],
        errors="coerce",
    )
    report_df["direct_reads"] = pd.to_numeric(
        report_df["direct_reads"],
        errors="coerce",
    )

    return report_df[
        report_df["name"].str.contains(taxon_pattern, na=False)
    ].copy()


def load_matching_read_ids(
    classifications_path: Path,
    selected_taxids: set[str],
) -> pd.DataFrame:
    """
    Load Kraken classifications and keep reads matching selected taxids.

    Parameters
    ----------
    classifications_path : Path
        Path to the Kraken classifications TSV.
    selected_taxids : set[str]
        Taxids to remove.

    Returns
    -------
    pd.DataFrame
        DataFrame of classified reads selected for removal.
    """
    class_cols = ["status", "read_id", "taxid", "length", "lca"]
    class_df = pd.read_csv(
        classifications_path,
        sep="\t",
        names=class_cols,
        header=None,
        dtype={"taxid": str, "read_id": str, "status": str},
    )
    return class_df[
        (class_df["status"] == "C")
        & (class_df["taxid"].astype(str).isin(selected_taxids))
    ].copy()


def write_filtered_fastq(
    input_fastq: Path,
    output_fastq: Path,
    read_ids_to_remove: set[str],
) -> tuple[int, int]:
    """
    Write a FASTQ file with selected read IDs removed.

    Parameters
    ----------
    input_fastq : Path
        Path to the input FASTQ file.
    output_fastq : Path
        Path to the filtered output FASTQ file.
    read_ids_to_remove : set[str]
        Read IDs to exclude.

    Returns
    -------
    tuple[int, int]
        Number of reads kept and removed.
    """
    kept = 0
    removed = 0

    with input_fastq.open("r") as in_handle, output_fastq.open("w") as out_handle:
        while True:
            header = in_handle.readline()
            if not header:
                break
            sequence = in_handle.readline()
            plus = in_handle.readline()
            quality = in_handle.readline()

            if not quality:
                raise ValueError("Input FASTQ appears truncated.")

            read_id = header.strip().split()[0].lstrip("@")
            if read_id in read_ids_to_remove:
                removed += 1
                continue

            out_handle.write(header)
            out_handle.write(sequence)
            out_handle.write(plus)
            out_handle.write(quality)
            kept += 1

    return kept, removed


def main() -> None:
    """
    Run the filtering workflow.
    """
    args = parse_args()

    selected_taxa = load_kraken_report(
        report_path=Path(args.kraken_report),
        taxon_pattern=args.taxon_pattern,
    )
    selected_taxids = set(selected_taxa["taxid"].dropna().astype(str))

    reads_df = load_matching_read_ids(
        classifications_path=Path(args.kraken_classifications),
        selected_taxids=selected_taxids,
    )

    read_ids = list(reads_df["read_id"].astype(str))
    read_id_set = set(read_ids)

    Path(args.read_ids_out).write_text("\n".join(read_ids) + "\n")

    summary_df = (
        reads_df["taxid"]
        .value_counts()
        .rename_axis("taxid")
        .reset_index(name="read_count")
        .merge(
            selected_taxa[["taxid", "name", "rank_code"]].drop_duplicates(),
            on="taxid",
            how="left",
        )
        .loc[:, ["taxid", "name", "rank_code", "read_count"]]
    )
    summary_df.to_csv(args.summary_tsv, sep="\t", index=False)

    kept, removed = write_filtered_fastq(
        input_fastq=Path(args.input_fastq),
        output_fastq=Path(args.output_fastq),
        read_ids_to_remove=read_id_set,
    )

    print(f"Selected taxa: {len(selected_taxids)}")
    print(f"Removed reads: {removed}")
    print(f"Kept reads: {kept}")


if __name__ == "__main__":
    main()
