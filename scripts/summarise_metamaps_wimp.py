#!/usr/bin/env python3
"""Summarise MetaMaps WIMP output for one or more requested target labels."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


POSSIBLE_LEVEL_COLUMNS = ["AnalysisLevel", "Level", "level"]
POSSIBLE_LABEL_COLUMNS = ["Name", "TaxonName", "taxonName", "Label", "Taxon"]
POSSIBLE_COUNT_COLUMNS = ["Absolute", "absolute", "Reads", "Count", "n_reads"]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
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


def choose_column(fieldnames: list[str], candidates: list[str]) -> str:
    """Return the first candidate column present in the table."""
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    raise ValueError(
        f"Could not find any of the expected columns {candidates} in WIMP file."
    )


def load_rows(path: str) -> list[dict[str, str]]:
    """Load WIMP rows from a tab-separated file."""
    with open(path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No header found in WIMP file: {path}")
        rows = list(reader)
        if not rows:
            raise ValueError(f"No data rows found in WIMP file: {path}")
        return rows


def to_int(value: str) -> int:
    """Convert a string count value to an integer robustly."""
    return int(round(float(value)))


def main() -> None:
    """Run WIMP summarisation."""
    args = parse_args()
    rows = load_rows(path=args.wimp_tsv)
    fieldnames = list(rows[0].keys())
    level_column = choose_column(fieldnames=fieldnames, candidates=POSSIBLE_LEVEL_COLUMNS)
    label_column = choose_column(fieldnames=fieldnames, candidates=POSSIBLE_LABEL_COLUMNS)
    count_column = choose_column(fieldnames=fieldnames, candidates=POSSIBLE_COUNT_COLUMNS)

    selected_rows = [
        row for row in rows if str(row.get(level_column, "")).strip() == args.analysis_level
    ]
    if not selected_rows:
        raise ValueError(
            f"No rows found for analysis level '{args.analysis_level}' in {args.wimp_tsv}"
        )

    total_reads = sum(to_int(row[count_column]) for row in selected_rows)
    target_labels = args.target_label
    target_counts = []
    for label in target_labels:
        label_cf = label.casefold()
        count = sum(
            to_int(row[count_column])
            for row in selected_rows
            if label_cf in str(row.get(label_column, "")).casefold()
        )
        target_counts.append(count)

    output_path = Path(args.out_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        header = ["wimp_tsv", "analysis_level", "total_reads_at_level"]
        header.extend([f"target_reads_{i + 1}" for i in range(len(target_labels))])
        header.extend([f"target_label_{i + 1}" for i in range(len(target_labels))])
        writer.writerow(header)

        row = [args.wimp_tsv, args.analysis_level, total_reads]
        row.extend(target_counts)
        row.extend(target_labels)
        writer.writerow(row)


if __name__ == "__main__":
    main()
