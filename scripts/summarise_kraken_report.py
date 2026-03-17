#!/usr/bin/env python3
"""Summarise Kraken2 report files and count target-associated assignments."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass
class KrakenRow:
    """Container for one Kraken2 report row."""

    percent: float
    clade_reads: int
    direct_reads: int
    rank_code: str
    taxid: str
    name: str


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarise a Kraken2 report file."
    )
    parser.add_argument(
        "--kraken_report_tsv",
        required=True,
        help="Input Kraken2 report TSV.",
    )
    parser.add_argument(
        "--target_label",
        action="append",
        default=[],
        help=(
            "Target label substring to match in taxon names. "
            "May be provided multiple times."
        ),
    )
    parser.add_argument(
        "--out_tsv",
        required=True,
        help="Output summary TSV.",
    )
    return parser.parse_args()


def read_kraken_report(path: str) -> list[KrakenRow]:
    """Parse Kraken2 report rows."""
    rows: list[KrakenRow] = []
    with open(path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for raw in reader:
            if len(raw) < 6:
                continue
            rows.append(
                KrakenRow(
                    percent=float(raw[0].strip()),
                    clade_reads=int(raw[1].strip()),
                    direct_reads=int(raw[2].strip()),
                    rank_code=raw[3].strip(),
                    taxid=raw[4].strip(),
                    name=raw[5].strip(),
                )
            )
    if not rows:
        raise ValueError(f"No rows parsed from Kraken report: {path}")
    return rows


def get_unclassified_reads(rows: Iterable[KrakenRow]) -> int:
    """Return the number of unclassified reads from a Kraken report."""
    for row in rows:
        if row.rank_code == "U":
            return row.clade_reads
    return 0


def get_total_reads(rows: list[KrakenRow]) -> int:
    """Return the total number of sequences reported by Kraken2."""
    unclassified_reads = get_unclassified_reads(rows)
    classified_reads = 0
    for row in rows:
        if row.taxid == "1" and row.name == "root":
            classified_reads = row.clade_reads
            break
    if classified_reads == 0:
        classified_reads = max(
            (row.clade_reads for row in rows if row.rank_code != "U"),
            default=0,
        )
    return classified_reads + unclassified_reads


def get_classified_reads(rows: list[KrakenRow]) -> int:
    """Return the number of classified reads from a Kraken report."""
    total_reads = get_total_reads(rows)
    unclassified_reads = get_unclassified_reads(rows)
    classified_reads = total_reads - unclassified_reads
    return max(classified_reads, 0)


def count_target_reads(rows: Iterable[KrakenRow], target_label: str) -> int:
    """Count Kraken clade reads for rows whose names contain the target label."""
    label = target_label.casefold()
    return sum(
        row.clade_reads for row in rows if label in row.name.casefold()
    )


def write_summary(
    kraken_report: str,
    classified_reads: int,
    target_labels: list[str],
    target_counts: list[int],
    out_tsv: str,
) -> None:
    """Write summary TSV for Kraken report parsing."""
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        header = ["kraken_report", "classified_reads"]
        header.extend([f"target_reads_{idx + 1}" for idx in range(len(target_labels))])
        header.extend([f"target_label_{idx + 1}" for idx in range(len(target_labels))])
        writer.writerow(header)

        row = [kraken_report, classified_reads]
        row.extend(target_counts)
        row.extend(target_labels)
        writer.writerow(row)


def main() -> None:
    """Run the Kraken report summarisation workflow."""
    args = parse_args()
    rows = read_kraken_report(path=args.kraken_report_tsv)
    classified_reads = get_classified_reads(rows=rows)
    target_labels = args.target_label
    target_counts = [
        count_target_reads(rows=rows, target_label=label)
        for label in target_labels
    ]
    write_summary(
        kraken_report=args.kraken_report_tsv,
        classified_reads=classified_reads,
        target_labels=target_labels,
        target_counts=target_counts,
        out_tsv=args.out_tsv,
    )


if __name__ == "__main__":
    main()
