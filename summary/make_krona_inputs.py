#!/usr/bin/env python3
"""Prepare simple Krona input tables from Kraken report files.

This script scans one or more directories for Kraken report TSV files and
writes Krona-compatible text tables. Each output line contains a count followed
by a taxonomic path split on semicolons.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


KRANKEN_REPORT_NAME = "kraken.report.tsv"



def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for Krona input preparation."""
    parser = argparse.ArgumentParser(
        description="Create Krona input tables from Kraken report TSV files."
    )
    parser.add_argument(
        "--input_dirs",
        nargs="+",
        required=True,
        help="One or more directories to scan recursively for kraken.report.tsv files.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for Krona input text files.",
    )
    parser.add_argument(
        "--min_count",
        type=int,
        default=1,
        help="Minimum clade count required to include a report row.",
    )
    return parser.parse_args()



def discover_reports(input_dirs: list[str]) -> list[Path]:
    """Discover Kraken report files recursively."""
    found: list[Path] = []
    for input_dir in input_dirs:
        base = Path(input_dir).expanduser().resolve()
        if not base.exists():
            continue
        found.extend(sorted(base.rglob(KRANKEN_REPORT_NAME)))
    return found



def parse_kraken_report(path: Path) -> pd.DataFrame:
    """Parse a Kraken report into a pandas DataFrame."""
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            "percent",
            "clade_count",
            "taxon_count",
            "rank_code",
            "taxid",
            "name",
        ],
    )



def write_krona_input(dataframe: pd.DataFrame, out_path: Path, min_count: int) -> None:
    """Write a simple Krona input file from a Kraken report DataFrame."""
    with open(out_path, "w", encoding="utf-8") as handle:
        for _, row in dataframe.iterrows():
            clade_count = int(row["clade_count"])
            if clade_count < min_count:
                continue
            raw_name = str(row["name"]).strip()
            clean_name = " ".join(raw_name.split())
            if not clean_name:
                continue
            lineage_parts = [part.strip() for part in clean_name.split(";") if part.strip()]
            if not lineage_parts:
                lineage_parts = [clean_name]
            handle.write("\t".join([str(clade_count)] + lineage_parts) + "\n")



def main() -> None:
    """Run the Krona input preparation workflow."""
    args = parse_args()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    reports = discover_reports(input_dirs=args.input_dirs)
    manifest_records: list[dict[str, str]] = []

    for report_path in reports:
        try:
            dataframe = parse_kraken_report(path=report_path)
            out_path = out_dir / f"{report_path.parent.name}.krona_input.txt"
            write_krona_input(dataframe=dataframe, out_path=out_path, min_count=args.min_count)
            manifest_records.append(
                {
                    "kraken_report": str(report_path),
                    "krona_input": str(out_path),
                }
            )
        except Exception as exc:  # noqa: BLE001
            manifest_records.append(
                {
                    "kraken_report": str(report_path),
                    "krona_input": f"ERROR: {type(exc).__name__}: {exc}",
                }
            )

    manifest_df = pd.DataFrame.from_records(manifest_records)
    manifest_df.to_csv(out_dir / "krona_manifest.tsv", sep="\t", index=False)
    print(f"[INFO] Wrote Krona inputs to: {out_dir}")


if __name__ == "__main__":
    main()
