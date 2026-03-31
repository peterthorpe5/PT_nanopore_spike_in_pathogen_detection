#!/usr/bin/env python3
"""
Summarise Kraken2 report files and additionally write a reported-taxa long table.

This wrapper preserves the original target-centric summary TSV but also writes
a second long-format TSV containing all report rows that match a set of broad
taxon-name filters, so off-target and cross-hit behaviour can be analysed
downstream.

Default reported-taxa filter behaviour
--------------------------------------
If ``--reported_name_prefix`` is not supplied, the script derives one or more
prefixes from the first token of each ``--target_label`` value. For example,
targets such as "Plasmodium vivax" and "Plasmodium falciparum" will lead to a
single broad filter of "Plasmodium".

Outputs
-------
1. Standard target-centric summary TSV via ``--out_tsv``.
2. A sidecar reported-taxa long TSV. By default this is written next to the
   original Kraken report and is named:
   ``<kraken_report_stem>.reported_taxa_long.tsv``

The sidecar contains one row per report taxon per metric, with metrics:
- kraken_reported_clade_reads
- kraken_reported_direct_reads
- kraken_reported_found
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable, Optional

import summarise_kraken_report as base


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Summarise a Kraken2 report file and also write a reported-taxa "
            "long table for off-target analysis."
        )
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
        help="Output target-centric summary TSV.",
    )
    parser.add_argument(
        "--reported_taxa_tsv",
        default=None,
        help=(
            "Optional output path for the reported-taxa long TSV. "
            "Defaults to <kraken_report_stem>.reported_taxa_long.tsv."
        ),
    )
    parser.add_argument(
        "--reported_name_prefix",
        action="append",
        default=[],
        help=(
            "Optional broad taxon-name prefix to retain in the reported-taxa "
            "long table. May be supplied multiple times."
        ),
    )
    return parser.parse_args()


def derive_reported_prefixes(
    target_labels: list[str],
    explicit_prefixes: list[str],
) -> list[str]:
    """
    Derive broad taxon-name prefixes for reported-taxa export.

    Parameters
    ----------
    target_labels : list[str]
        Requested target labels.
    explicit_prefixes : list[str]
        User-supplied prefixes.

    Returns
    -------
    list[str]
        Ordered unique list of prefixes.
    """
    prefixes: list[str] = []

    for prefix in explicit_prefixes:
        prefix = str(prefix).strip()
        if prefix and prefix not in prefixes:
            prefixes.append(prefix)

    if prefixes:
        return prefixes

    for label in target_labels:
        label = str(label).strip()
        if not label:
            continue
        first_token = label.split()[0]
        if first_token and first_token not in prefixes:
            prefixes.append(first_token)

    if not prefixes:
        prefixes.append("Plasmodium")

    return prefixes


def default_reported_taxa_path(kraken_report_tsv: str) -> Path:
    """
    Build the default sidecar reported-taxa TSV path.

    Parameters
    ----------
    kraken_report_tsv : str
        Input Kraken report TSV path.

    Returns
    -------
    Path
        Default reported-taxa TSV path.
    """
    report_path = Path(kraken_report_tsv)
    return report_path.with_name(
        f"{report_path.stem}.reported_taxa_long.tsv"
    )


def row_matches_any_prefix(name: str, prefixes: list[str]) -> tuple[bool, Optional[str]]:
    """
    Test whether a Kraken taxon name matches any requested broad prefix.

    Parameters
    ----------
    name : str
        Kraken report taxon name.
    prefixes : list[str]
        Prefixes to test.

    Returns
    -------
    tuple[bool, Optional[str]]
        Match flag and the first matching prefix.
    """
    name_cf = name.casefold()
    for prefix in prefixes:
        prefix_cf = prefix.casefold()
        if prefix_cf in name_cf:
            return True, prefix
    return False, None


def matched_target_labels(
    name: str,
    target_labels: list[str],
) -> list[str]:
    """
    Return requested target labels that match a Kraken taxon name.

    Parameters
    ----------
    name : str
        Kraken report taxon name.
    target_labels : list[str]
        Requested target labels.

    Returns
    -------
    list[str]
        Matching target labels.
    """
    name_cf = name.casefold()
    matches: list[str] = []
    for label in target_labels:
        label_cf = str(label).strip().casefold()
        if label_cf and label_cf in name_cf and label not in matches:
            matches.append(label)
    return matches


def write_reported_taxa_long(
    rows: Iterable[base.KrakenRow],
    kraken_report_tsv: str,
    out_tsv: Path,
    target_labels: list[str],
    prefixes: list[str],
) -> None:
    """
    Write the reported-taxa long TSV.

    Parameters
    ----------
    rows : Iterable[base.KrakenRow]
        Parsed Kraken report rows.
    kraken_report_tsv : str
        Source Kraken report path.
    out_tsv : Path
        Output reported-taxa TSV path.
    target_labels : list[str]
        Requested target labels.
    prefixes : list[str]
        Broad prefixes retained in the sidecar output.
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "classifier",
        "report_tsv",
        "filter_prefix",
        "matched_target_labels",
        "taxid",
        "taxon_name",
        "rank_code",
        "percent",
        "metric_name",
        "metric_value",
    ]

    with out_tsv.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)

        for row in rows:
            is_match, matched_prefix = row_matches_any_prefix(
                name=row.name,
                prefixes=prefixes,
            )
            if not is_match:
                continue

            matched_labels = matched_target_labels(
                name=row.name,
                target_labels=target_labels,
            )
            matched_label_text = "; ".join(matched_labels)

            metric_rows = [
                ("kraken_reported_clade_reads", row.clade_reads),
                ("kraken_reported_direct_reads", row.direct_reads),
                ("kraken_reported_found", 1 if row.clade_reads > 0 else 0),
            ]

            for metric_name, metric_value in metric_rows:
                writer.writerow(
                    [
                        "kraken",
                        kraken_report_tsv,
                        matched_prefix,
                        matched_label_text,
                        row.taxid,
                        row.name,
                        row.rank_code,
                        row.percent,
                        metric_name,
                        metric_value,
                    ]
                )


def main() -> None:
    """Run the Kraken report summarisation workflow."""
    args = parse_args()

    rows = base.read_kraken_report(path=args.kraken_report_tsv)
    classified_reads = base.get_classified_reads(rows=rows)
    target_labels = args.target_label
    target_counts = [
        base.count_target_reads(rows=rows, target_label=label)
        for label in target_labels
    ]
    base.write_summary(
        kraken_report=args.kraken_report_tsv,
        classified_reads=classified_reads,
        target_labels=target_labels,
        target_counts=target_counts,
        out_tsv=args.out_tsv,
    )

    prefixes = derive_reported_prefixes(
        target_labels=target_labels,
        explicit_prefixes=args.reported_name_prefix,
    )
    reported_taxa_tsv = (
        Path(args.reported_taxa_tsv)
        if args.reported_taxa_tsv is not None
        else default_reported_taxa_path(args.kraken_report_tsv)
    )

    write_reported_taxa_long(
        rows=rows,
        kraken_report_tsv=args.kraken_report_tsv,
        out_tsv=reported_taxa_tsv,
        target_labels=target_labels,
        prefixes=prefixes,
    )


if __name__ == "__main__":
    main()
