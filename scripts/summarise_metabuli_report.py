#!/usr/bin/env python3
"""
Summarise Metabuli Kraken-format report output for one or more targets.

This script parses the `JobID_report.tsv` file produced by Metabuli and
extracts clade and direct counts for one or more requested targets.

Targets can be matched by taxon name, taxon ID, or both. If a taxon ID is
provided, it is preferred over the name for matching. The output is designed
to be shell-friendly for integration into a larger spike-in benchmarking
workflow.

Outputs
-------
A two-line TSV file:
- header
- one summary row

Columns
-------
- total_reads_in_report
- for each target:
  - target_clade_reads_N
  - target_direct_reads_N
  - target_found_N

Notes
-----
Metabuli produces a Kraken-format report file. This helper uses the clade
count as the main per-target detection count because that most closely matches
the usual taxon-level summary used in metagenomic benchmarking.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Optional


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Summarise Metabuli report output for requested targets."
    )
    parser.add_argument(
        "--report_tsv",
        required=True,
        help="Metabuli JobID_report.tsv file.",
    )
    parser.add_argument(
        "--out_tsv",
        required=True,
        help="Output summary TSV.",
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
    return parser.parse_args()


def normalise_taxid_list(
    target_labels: list[str],
    target_taxids: list[str],
) -> list[Optional[str]]:
    """
    Normalise the target taxon ID list to the same length as target labels.

    Parameters
    ----------
    target_labels : list[str]
        Target labels.
    target_taxids : list[str]
        Target taxon IDs.

    Returns
    -------
    list[Optional[str]]
        Taxon IDs aligned to the target label list.
    """
    taxids: list[Optional[str]] = [None] * len(target_labels)
    for idx, value in enumerate(target_taxids[: len(target_labels)]):
        if value is None:
            continue
        value = str(value).strip()
        taxids[idx] = value if value else None
    return taxids


def load_report_rows(report_tsv: str) -> list[dict[str, str]]:
    """
    Load a Kraken-format report file.

    Parameters
    ----------
    report_tsv : str
        Path to report TSV.

    Returns
    -------
    list[dict[str, str]]
        Parsed row dictionaries.
    """
    rows: list[dict[str, str]] = []

    with open(report_tsv, mode="rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for fields in reader:
            if len(fields) < 6:
                continue

            rows.append(
                {
                    "percent": fields[0].strip(),
                    "clade_reads": fields[1].strip(),
                    "direct_reads": fields[2].strip(),
                    "rank_code": fields[3].strip(),
                    "taxid": fields[4].strip(),
                    "name": fields[5].strip(),
                }
            )

    return rows


def to_int(value: str) -> int:
    """
    Convert a string to an integer count.

    Parameters
    ----------
    value : str
        Input value.

    Returns
    -------
    int
        Parsed integer, or zero if parsing fails.
    """
    try:
        return int(value)
    except Exception:
        return 0


def find_target_row(
    rows: list[dict[str, str]],
    target_label: str,
    target_taxid: Optional[str],
) -> Optional[dict[str, str]]:
    """
    Find a target row in the Kraken-format report.

    Parameters
    ----------
    rows : list[dict[str, str]]
        Parsed report rows.
    target_label : str
        Target label.
    target_taxid : Optional[str]
        Optional target taxon ID.

    Returns
    -------
    Optional[dict[str, str]]
        Matching row, or None if not found.
    """
    if target_taxid is not None:
        for row in rows:
            if row["taxid"] == target_taxid:
                return row

    label_norm = target_label.strip().casefold()
    for row in rows:
        if row["name"].strip().casefold() == label_norm:
            return row

    return None


def write_summary(
    out_tsv: str,
    total_reads_in_report: int,
    clade_counts: list[int],
    direct_counts: list[int],
    found_flags: list[int],
) -> None:
    """
    Write the summary TSV.

    Parameters
    ----------
    out_tsv : str
        Output TSV path.
    total_reads_in_report : int
        Total reads represented by the report.
    clade_counts : list[int]
        Per-target clade read counts.
    direct_counts : list[int]
        Per-target direct read counts.
    found_flags : list[int]
        Per-target found flags.
    """
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    header = ["total_reads_in_report"]
    row = [str(total_reads_in_report)]

    for idx in range(len(clade_counts)):
        header.extend(
            [
                f"target_clade_reads_{idx + 1}",
                f"target_direct_reads_{idx + 1}",
                f"target_found_{idx + 1}",
            ]
        )
        row.extend(
            [
                str(clade_counts[idx]),
                str(direct_counts[idx]),
                str(found_flags[idx]),
            ]
        )

    with open(out_tsv, mode="wt", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        handle.write("\t".join(row) + "\n")


def main() -> None:
    """
    Run report summarisation.
    """
    args = parse_args()

    target_labels = [str(x).strip() for x in args.target_labels]
    target_taxids = normalise_taxid_list(
        target_labels=target_labels,
        target_taxids=args.target_taxids,
    )

    rows = load_report_rows(report_tsv=args.report_tsv)
    total_reads_in_report = 0
    if rows:
        total_reads_in_report = max(to_int(row["clade_reads"]) for row in rows)

    clade_counts: list[int] = []
    direct_counts: list[int] = []
    found_flags: list[int] = []

    for target_label, target_taxid in zip(target_labels, target_taxids):
        row = find_target_row(
            rows=rows,
            target_label=target_label,
            target_taxid=target_taxid,
        )
        if row is None:
            clade_counts.append(0)
            direct_counts.append(0)
            found_flags.append(0)
        else:
            clade_value = to_int(row["clade_reads"])
            direct_value = to_int(row["direct_reads"])
            clade_counts.append(clade_value)
            direct_counts.append(direct_value)
            found_flags.append(1 if clade_value > 0 else 0)

    write_summary(
        out_tsv=args.out_tsv,
        total_reads_in_report=total_reads_in_report,
        clade_counts=clade_counts,
        direct_counts=direct_counts,
        found_flags=found_flags,
    )


if __name__ == "__main__":
    main()
