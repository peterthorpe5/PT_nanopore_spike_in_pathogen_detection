#!/usr/bin/env python3
"""
Summarise Metabuli report output and additionally write a reported-taxa long table.

This wrapper preserves the original target-centric summary TSV but also writes
a second long-format TSV containing all Metabuli report rows that match a set
of broad taxon-name filters, so off-target and cross-hit behaviour can be
analysed downstream.

Default reported-taxa filter behaviour
--------------------------------------
If ``--reported_name_prefix`` is not supplied, the script derives one or more
prefixes from the first token of each ``--target-label`` value. For example,
targets such as "Plasmodium vivax" and "Plasmodium falciparum" will lead to a
single broad filter of "Plasmodium".

Outputs
-------
1. Standard target-centric summary TSV via ``--out_tsv``.
2. A sidecar reported-taxa long TSV. By default this is written next to the
   original Metabuli report and is named:
   ``<report_stem>.reported_taxa_long.tsv``

The sidecar contains one row per report taxon per metric, with metrics:
- metabuli_reported_clade_reads
- metabuli_reported_direct_reads
- metabuli_reported_found
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Optional

import summarise_metabuli_report as base


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
            "Summarise Metabuli report output for requested targets and also "
            "write a reported-taxa long table for off-target analysis."
        )
    )
    parser.add_argument(
        "--report_tsv",
        required=True,
        help="Metabuli JobID_report.tsv file.",
    )
    parser.add_argument(
        "--out_tsv",
        required=True,
        help="Output target-centric summary TSV.",
    )
    parser.add_argument(
        "--target-label",
        action="append",
        dest="target_labels",
        default=[],
        help="Target species label to match exactly after whitespace stripping.",
    )
    parser.add_argument(
        "--target-taxid",
        action="append",
        dest="target_taxids",
        default=[],
        help=(
            "Optional target taxon ID. If provided, taxon ID matching is used "
            "for the corresponding target."
        ),
    )
    parser.add_argument(
        "--reported_taxa_tsv",
        default=None,
        help=(
            "Optional output path for the reported-taxa long TSV. "
            "Defaults to <report_stem>.reported_taxa_long.tsv."
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


def default_reported_taxa_path(report_tsv: str) -> Path:
    """
    Build the default sidecar reported-taxa TSV path.

    Parameters
    ----------
    report_tsv : str
        Input Metabuli report TSV path.

    Returns
    -------
    Path
        Default reported-taxa TSV path.
    """
    report_path = Path(report_tsv)
    return report_path.with_name(
        f"{report_path.stem}.reported_taxa_long.tsv"
    )


def row_matches_any_prefix(name: str, prefixes: list[str]) -> tuple[bool, Optional[str]]:
    """
    Test whether a Metabuli taxon name matches any requested broad prefix.

    Parameters
    ----------
    name : str
        Report taxon name.
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
    Return requested target labels that match a report taxon name.

    Parameters
    ----------
    name : str
        Report taxon name.
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
    rows: list[dict[str, str]],
    report_tsv: str,
    out_tsv: Path,
    target_labels: list[str],
    prefixes: list[str],
) -> None:
    """
    Write the reported-taxa long TSV.

    Parameters
    ----------
    rows : list[dict[str, str]]
        Parsed Metabuli report rows.
    report_tsv : str
        Source report TSV path.
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
            taxon_name = row["name"]
            is_match, matched_prefix = row_matches_any_prefix(
                name=taxon_name,
                prefixes=prefixes,
            )
            if not is_match:
                continue

            matched_labels = matched_target_labels(
                name=taxon_name,
                target_labels=target_labels,
            )
            matched_label_text = "; ".join(matched_labels)

            clade_reads = base.to_int(row["clade_reads"])
            direct_reads = base.to_int(row["direct_reads"])

            metric_rows = [
                ("metabuli_reported_clade_reads", clade_reads),
                ("metabuli_reported_direct_reads", direct_reads),
                ("metabuli_reported_found", 1 if clade_reads > 0 else 0),
            ]

            for metric_name, metric_value in metric_rows:
                writer.writerow(
                    [
                        "metabuli",
                        report_tsv,
                        matched_prefix,
                        matched_label_text,
                        row["taxid"],
                        taxon_name,
                        row["rank_code"],
                        row["percent"],
                        metric_name,
                        metric_value,
                    ]
                )


def main() -> None:
    """
    Run report summarisation.
    """
    args = parse_args()

    target_labels = [str(x).strip() for x in args.target_labels]
    target_taxids = base.normalise_taxid_list(
        target_labels=target_labels,
        target_taxids=args.target_taxids,
    )

    rows = base.load_report_rows(report_tsv=args.report_tsv)
    total_reads_in_report = 0
    if rows:
        total_reads_in_report = max(base.to_int(row["clade_reads"]) for row in rows)

    clade_counts: list[int] = []
    direct_counts: list[int] = []
    found_flags: list[int] = []

    for target_label, target_taxid in zip(target_labels, target_taxids):
        row = base.find_target_row(
            rows=rows,
            target_label=target_label,
            target_taxid=target_taxid,
        )
        if row is None:
            clade_counts.append(0)
            direct_counts.append(0)
            found_flags.append(0)
        else:
            clade_value = base.to_int(row["clade_reads"])
            direct_value = base.to_int(row["direct_reads"])
            clade_counts.append(clade_value)
            direct_counts.append(direct_value)
            found_flags.append(1 if clade_value > 0 else 0)

    base.write_summary(
        out_tsv=args.out_tsv,
        total_reads_in_report=total_reads_in_report,
        clade_counts=clade_counts,
        direct_counts=direct_counts,
        found_flags=found_flags,
    )

    prefixes = derive_reported_prefixes(
        target_labels=target_labels,
        explicit_prefixes=args.reported_name_prefix,
    )
    reported_taxa_tsv = (
        Path(args.reported_taxa_tsv)
        if args.reported_taxa_tsv is not None
        else default_reported_taxa_path(args.report_tsv)
    )

    write_reported_taxa_long(
        rows=rows,
        report_tsv=args.report_tsv,
        out_tsv=reported_taxa_tsv,
        target_labels=target_labels,
        prefixes=prefixes,
    )


if __name__ == "__main__":
    main()
