#!/usr/bin/env python3
"""Summarise MetaMaps WIMP output for one or more requested target labels."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Optional


POSSIBLE_LEVEL_COLUMNS = ["AnalysisLevel", "Level", "level"]
POSSIBLE_LABEL_COLUMNS = ["Name", "TaxonName", "taxonName", "Label", "Taxon"]
POSSIBLE_COUNT_COLUMNS = ["Absolute", "absolute", "Reads", "Count", "n_reads"]


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Summarise MetaMaps WIMP output for target taxa."
    )
    parser.add_argument(
        "--wimp_tsv",
        required=True,
        help="MetaMaps .EM.WIMP file.",
    )
    parser.add_argument(
        "--target_label",
        action="append",
        default=[],
        help="Target label substring to match in the WIMP taxon label.",
    )
    parser.add_argument(
        "--analysis_level",
        default="definedGenomes",
        help="MetaMaps WIMP analysis level to query. Default: definedGenomes.",
    )
    parser.add_argument(
        "--out_tsv",
        required=True,
        help="Output summary TSV.",
    )
    return parser.parse_args()


def choose_column(fieldnames: list[str], candidates: list[str]) -> Optional[str]:
    """
    Return the first candidate column present in the table.

    Parameters
    ----------
    fieldnames : list[str]
        Available column names.
    candidates : list[str]
        Candidate column names.

    Returns
    -------
    Optional[str]
        Matching column name, or None if none are present.
    """
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    return None


def load_rows(path: str) -> tuple[list[dict[str, str]], list[str]]:
    """
    Load WIMP rows from a tab-separated file.

    Parameters
    ----------
    path : str
        WIMP TSV path.

    Returns
    -------
    tuple[list[dict[str, str]], list[str]]
        Parsed rows and field names.
    """
    with open(path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader) if fieldnames else []
    return rows, fieldnames


def to_int(value: str) -> int:
    """
    Convert a string count value to an integer robustly.

    Parameters
    ----------
    value : str
        Input value.

    Returns
    -------
    int
        Parsed integer count.
    """
    try:
        return int(round(float(value)))
    except Exception:
        return 0


def write_summary(
    out_tsv: str,
    wimp_tsv: str,
    analysis_level: str,
    target_labels: list[str],
    total_reads: int,
    target_counts: list[int],
    status: str,
) -> None:
    """
    Write summary TSV.

    Parameters
    ----------
    out_tsv : str
        Output TSV path.
    wimp_tsv : str
        Input WIMP TSV path.
    analysis_level : str
        Requested analysis level.
    target_labels : list[str]
        Requested target labels.
    total_reads : int
        Total reads at the selected level.
    target_counts : list[int]
        Target read counts.
    status : str
        Status string describing the parse outcome.
    """
    output_path = Path(out_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        header = [
            "wimp_tsv",
            "analysis_level",
            "total_reads_at_level",
            "status",
        ]
        header.extend([f"target_reads_{i + 1}" for i in range(len(target_labels))])
        header.extend([f"target_label_{i + 1}" for i in range(len(target_labels))])
        writer.writerow(header)

        row = [wimp_tsv, analysis_level, total_reads, status]
        row.extend(target_counts)
        row.extend(target_labels)
        writer.writerow(row)


def main() -> None:
    """
    Run WIMP summarisation.
    """
    args = parse_args()

    rows, fieldnames = load_rows(path=args.wimp_tsv)
    target_labels = args.target_label or []

    if not fieldnames:
        write_summary(
            out_tsv=args.out_tsv,
            wimp_tsv=args.wimp_tsv,
            analysis_level=args.analysis_level,
            target_labels=target_labels,
            total_reads=0,
            target_counts=[0] * len(target_labels),
            status="no_header",
        )
        return

    if not rows:
        write_summary(
            out_tsv=args.out_tsv,
            wimp_tsv=args.wimp_tsv,
            analysis_level=args.analysis_level,
            target_labels=target_labels,
            total_reads=0,
            target_counts=[0] * len(target_labels),
            status="no_data_rows",
        )
        return

    level_column = choose_column(fieldnames=fieldnames, candidates=POSSIBLE_LEVEL_COLUMNS)
    label_column = choose_column(fieldnames=fieldnames, candidates=POSSIBLE_LABEL_COLUMNS)
    count_column = choose_column(fieldnames=fieldnames, candidates=POSSIBLE_COUNT_COLUMNS)

    if level_column is None or label_column is None or count_column is None:
        write_summary(
            out_tsv=args.out_tsv,
            wimp_tsv=args.wimp_tsv,
            analysis_level=args.analysis_level,
            target_labels=target_labels,
            total_reads=0,
            target_counts=[0] * len(target_labels),
            status="missing_expected_columns",
        )
        return

    selected_rows = [
        row
        for row in rows
        if str(row.get(level_column, "")).strip() == args.analysis_level
    ]

    if not selected_rows:
        write_summary(
            out_tsv=args.out_tsv,
            wimp_tsv=args.wimp_tsv,
            analysis_level=args.analysis_level,
            target_labels=target_labels,
            total_reads=0,
            target_counts=[0] * len(target_labels),
            status="no_rows_at_requested_level",
        )
        return

    total_reads = sum(to_int(row.get(count_column, "0")) for row in selected_rows)

    target_counts = []
    for label in target_labels:
        label_cf = label.casefold()
        count = sum(
            to_int(row.get(count_column, "0"))
            for row in selected_rows
            if label_cf in str(row.get(label_column, "")).casefold()
        )
        target_counts.append(count)

    write_summary(
        out_tsv=args.out_tsv,
        wimp_tsv=args.wimp_tsv,
        analysis_level=args.analysis_level,
        target_labels=target_labels,
        total_reads=total_reads,
        target_counts=target_counts,
        status="ok",
    )


if __name__ == "__main__":
    main()