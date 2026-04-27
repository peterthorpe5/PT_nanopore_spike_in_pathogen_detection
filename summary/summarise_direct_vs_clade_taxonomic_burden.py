#!/usr/bin/env python3
"""Summarise expected versus off-target taxonomic signal.

This script scans reported-taxa long-format files produced by the spike-in
benchmark and summarises the balance between expected target signal and
non-expected Plasmodium signal.

The main purpose is to separate two related but different interpretations of
Kraken-style taxonomic reports:

1. Clade-level burden: uses hierarchical clade counts such as
   ``kraken_reported_clade_reads``. These are useful for understanding the
   broader reporting burden visible in a taxonomic report, but they can count
   the same evidence at several taxonomic levels.
2. Direct-assignment burden: uses direct counts such as
   ``kraken_reported_direct_reads``. These provide a less double-counted view
   of where reads or contigs are assigned directly.

The script writes TSV, Excel, and HTML outputs. All delimited output files are
written as tab-separated text.
"""

from __future__ import annotations

import argparse
import html
import math
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


DEFAULT_METRIC_GROUPS = {
    "clade": [
        "kraken_reported_clade_reads",
        "metabuli_reported_clade_reads",
        "reported_clade_reads",
    ],
    "direct": [
        "kraken_reported_direct_reads",
        "metabuli_reported_direct_reads",
        "reported_direct_reads",
    ],
    "found": [
        "kraken_reported_found",
        "metabuli_reported_found",
        "reported_found",
    ],
}


REQUIRED_COLUMNS = {
    "classifier",
    "report_tsv",
    "matched_target_labels",
    "taxon_name",
    "metric_name",
    "metric_value",
}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Summarise expected versus off-target taxonomic signal using "
            "reported_taxa_long-style files."
        )
    )
    parser.add_argument(
        "--input_roots",
        nargs="+",
        required=True,
        help=(
            "Run directories or summary directories to scan recursively for "
            "reported_taxa_long TSV files."
        ),
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for TSV, Excel, and HTML summaries.",
    )
    parser.add_argument(
        "--file_patterns",
        nargs="+",
        default=[
            "*reported_taxa_long.tsv",
            "*.reported_taxa_long.tsv",
            "relative_taxonomic_burden_input_classified.tsv",
        ],
        help="Filename patterns to scan under each input root.",
    )
    parser.add_argument(
        "--taxon_regex",
        default="Plasmodium",
        help="Regex used to define taxa considered relevant off-targets.",
    )
    parser.add_argument(
        "--min_metric_value",
        type=float,
        default=0.0,
        help="Minimum metric value to include in burden calculations.",
    )
    parser.add_argument(
        "--summary_name",
        default="direct_vs_clade_taxonomic_burden",
        help="Prefix for output filenames.",
    )
    return parser.parse_args()


def find_input_files(
    *,
    input_roots: Iterable[str],
    file_patterns: Iterable[str],
) -> list[Path]:
    """Find unique input files matching one or more patterns."""
    files: list[Path] = []
    seen: set[Path] = set()
    for root_text in input_roots:
        root = Path(root_text).expanduser().resolve()
        if root.is_file():
            candidates = [root]
        else:
            candidates = []
            for pattern in file_patterns:
                candidates.extend(root.rglob(pattern))
        for candidate in candidates:
            path = candidate.resolve()
            if path in seen or not path.is_file() or path.stat().st_size == 0:
                continue
            seen.add(path)
            files.append(path)
    return sorted(files)


def safe_read_tsv(path: Path) -> pd.DataFrame:
    """Read a TSV file, returning an empty frame if it is unusable."""
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return pd.DataFrame()
    if df.empty:
        return pd.DataFrame()
    df.columns = [str(col).strip() for col in df.columns]
    return df


def infer_run_fields(report_tsv: str, source_file: str) -> dict[str, object]:
    """Infer run name, replicate, and spike level from paths."""
    path_text = report_tsv or source_file
    parts = Path(path_text).parts

    run_name = "unknown_run"
    for part in parts:
        if part.startswith("spikein_") or part.startswith("runs_minimap"):
            run_name = part

    mix_match = re.search(r"mix_rep(\d+)_n(\d+)", path_text)
    replicate = int(mix_match.group(1)) if mix_match else math.nan
    spike_n = int(mix_match.group(2)) if mix_match else math.nan

    if "multi" in run_name:
        workflow = "multi_assembly" if "flye" in run_name else "multi_read"
    elif "single" in run_name:
        workflow = "single_assembly" if "flye" in run_name else "single_read"
    elif "shuffle" in run_name:
        workflow = "shuffled_control"
    else:
        workflow = "unknown"

    is_shuffled = "shuffle" in run_name.lower() or "shuffled" in run_name.lower()

    return {
        "run_name": run_name,
        "workflow": workflow,
        "replicate": replicate,
        "spike_n": spike_n,
        "is_shuffled_control": is_shuffled,
    }


def normalise_input_frame(*, path: Path, taxon_regex: str) -> pd.DataFrame:
    """Read and normalise one reported-taxa file."""
    df = safe_read_tsv(path)
    if df.empty or not REQUIRED_COLUMNS.issubset(df.columns):
        return pd.DataFrame()

    df = df.copy()
    df["source_file"] = str(path)
    df["metric_value"] = pd.to_numeric(df["metric_value"], errors="coerce").fillna(0)
    df["matched_target_labels"] = df["matched_target_labels"].fillna("").astype(str)
    df["taxon_name"] = df["taxon_name"].fillna("").astype(str)
    df["metric_name"] = df["metric_name"].fillna("").astype(str)
    df["classifier"] = df["classifier"].fillna("unknown").astype(str)
    df["report_tsv"] = df["report_tsv"].fillna("").astype(str)

    inferred = df.apply(
        lambda row: infer_run_fields(
            report_tsv=row.get("report_tsv", ""),
            source_file=row.get("source_file", ""),
        ),
        axis=1,
        result_type="expand",
    )
    df = pd.concat([df, inferred], axis=1)

    taxon_pattern = re.compile(taxon_regex, flags=re.IGNORECASE)
    df["is_expected"] = df["matched_target_labels"].str.strip().ne("")
    df["is_relevant_taxon"] = df["taxon_name"].map(
        lambda value: bool(taxon_pattern.search(value))
    )
    df["signal_class"] = "ignored"
    df.loc[df["is_expected"], "signal_class"] = "expected"
    df.loc[
        (~df["is_expected"]) & df["is_relevant_taxon"],
        "signal_class",
    ] = "off_target"

    return df


def assign_metric_group(metric_name: str) -> str:
    """Assign a metric name to a broad burden group."""
    for group, names in DEFAULT_METRIC_GROUPS.items():
        if metric_name in names:
            return group
    if "clade" in metric_name:
        return "clade"
    if "direct" in metric_name:
        return "direct"
    if "found" in metric_name:
        return "found"
    return "other"


def summarise_observations(
    *,
    df: pd.DataFrame,
    min_metric_value: float,
) -> pd.DataFrame:
    """Create one burden row per observation and metric group."""
    if df.empty:
        return pd.DataFrame()

    work = df.copy()
    work = work[work["metric_value"] >= min_metric_value]
    work["metric_group"] = work["metric_name"].map(assign_metric_group)
    work = work[work["metric_group"].isin(["clade", "direct", "found"])]

    group_cols = [
        "workflow",
        "classifier",
        "run_name",
        "replicate",
        "spike_n",
        "is_shuffled_control",
        "metric_group",
        "report_tsv",
        "source_file",
    ]

    rows = []
    for keys, sub in work.groupby(group_cols, dropna=False):
        key_dict = dict(zip(group_cols, keys))
        expected = sub[sub["signal_class"] == "expected"]
        off = sub[sub["signal_class"] == "off_target"]
        expected_signal = expected["metric_value"].sum()
        off_signal = off["metric_value"].sum()
        total_signal = expected_signal + off_signal
        ratio = off_signal / expected_signal if expected_signal > 0 else math.nan
        expected_fraction = expected_signal / total_signal if total_signal > 0 else math.nan
        off_fraction = off_signal / total_signal if total_signal > 0 else math.nan

        if off.empty:
            dominant_taxon = ""
            dominant_signal = 0.0
            dominant_fraction = math.nan
        else:
            off_by_taxon = (
                off.groupby("taxon_name", dropna=False)["metric_value"]
                .sum()
                .sort_values(ascending=False)
            )
            dominant_taxon = str(off_by_taxon.index[0])
            dominant_signal = float(off_by_taxon.iloc[0])
            dominant_fraction = (
                dominant_signal / off_signal if off_signal > 0 else math.nan
            )

        rows.append(
            {
                **key_dict,
                "n_expected_taxa": expected["taxon_name"].nunique(),
                "n_off_target_taxa": off["taxon_name"].nunique(),
                "expected_signal": expected_signal,
                "off_target_signal": off_signal,
                "total_expected_plus_off_target_signal": total_signal,
                "off_target_to_expected_ratio": ratio,
                "expected_signal_fraction": expected_fraction,
                "off_target_signal_fraction": off_fraction,
                "dominant_off_target_taxon": dominant_taxon,
                "dominant_off_target_signal": dominant_signal,
                "dominant_off_target_fraction_of_offtarget": dominant_fraction,
                "expected_taxa": "; ".join(sorted(expected["taxon_name"].unique())),
                "off_target_taxa": "; ".join(sorted(off["taxon_name"].unique())),
            }
        )
    return pd.DataFrame(rows)


def summarise_group(
    *,
    observation_df: pd.DataFrame,
    group_cols: list[str],
) -> pd.DataFrame:
    """Summarise observation-level burden across grouping columns."""
    if observation_df.empty:
        return pd.DataFrame()

    summaries = []
    for keys, sub in observation_df.groupby(group_cols, dropna=False):
        key_dict = dict(zip(group_cols, keys))
        summaries.append(
            {
                **key_dict,
                "n_observations": len(sub),
                "median_expected_signal": sub["expected_signal"].median(),
                "median_off_target_signal": sub["off_target_signal"].median(),
                "median_expected_signal_fraction": sub[
                    "expected_signal_fraction"
                ].median(),
                "median_off_target_signal_fraction": sub[
                    "off_target_signal_fraction"
                ].median(),
                "median_off_target_to_expected_ratio": sub[
                    "off_target_to_expected_ratio"
                ].median(),
                "median_n_off_target_taxa": sub["n_off_target_taxa"].median(),
                "mean_expected_signal": sub["expected_signal"].mean(),
                "mean_off_target_signal": sub["off_target_signal"].mean(),
                "mean_expected_signal_fraction": sub[
                    "expected_signal_fraction"
                ].mean(),
                "mean_off_target_signal_fraction": sub[
                    "off_target_signal_fraction"
                ].mean(),
                "mean_off_target_to_expected_ratio": sub[
                    "off_target_to_expected_ratio"
                ].mean(),
                "mean_n_off_target_taxa": sub["n_off_target_taxa"].mean(),
            }
        )
    return pd.DataFrame(summaries)


def write_tsv(*, df: pd.DataFrame, path: Path) -> None:
    """Write a DataFrame as a tab-separated file."""
    df.to_csv(path, sep="\t", index=False)


def write_excel(*, tables: dict[str, pd.DataFrame], path: Path) -> None:
    """Write all result tables to a styled Excel workbook."""
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, df in tables.items():
            safe_name = sheet_name[:31]
            df.to_excel(writer, sheet_name=safe_name, index=False)
            ws = writer.book[safe_name]
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions
            for cell in ws[1]:
                cell.fill = __import__("openpyxl").styles.PatternFill(
                    start_color="1F4E79",
                    end_color="1F4E79",
                    fill_type="solid",
                )
                cell.font = __import__("openpyxl").styles.Font(
                    bold=True,
                    color="FFFFFF",
                )
                cell.alignment = __import__("openpyxl").styles.Alignment(
                    wrap_text=True,
                    vertical="center",
                )
            for col in ws.columns:
                col_letter = col[0].column_letter
                max_len = max(
                    len(str(cell.value)) if cell.value is not None else 0
                    for cell in col[:200]
                )
                ws.column_dimensions[col_letter].width = min(max(max_len + 2, 10), 45)


def dataframe_to_html_table(df: pd.DataFrame, *, max_rows: int = 100) -> str:
    """Render a compact HTML table."""
    if df.empty:
        return "<p>No rows available.</p>"
    display_df = df.head(max_rows).copy()
    return display_df.to_html(index=False, escape=True, classes="data-table")


def write_html_report(*, tables: dict[str, pd.DataFrame], path: Path) -> None:
    """Write an HTML report summarising the burden outputs."""
    css = """
    body { font-family: Arial, Helvetica, sans-serif; margin: 28px; color: #1a1a1a; }
    h1, h2 { color: #1f4e79; }
    .note { max-width: 1100px; line-height: 1.45; }
    .table-wrap { overflow-x: auto; border: 1px solid #d9e2ef; margin: 18px 0; }
    table { border-collapse: collapse; font-size: 13px; width: 100%; }
    th { background: #1f4e79; color: white; padding: 7px; text-align: left; white-space: nowrap; }
    td { border-bottom: 1px solid #e6edf5; padding: 6px; vertical-align: top; }
    tr:nth-child(even) { background: #fbfdff; }
    .small { color: #555; font-size: 13px; }
    """
    sections = []
    for title, df in tables.items():
        escaped = html.escape(title.replace("_", " ").title())
        sections.append(
            f"<h2>{escaped}</h2>\n<div class='table-wrap'>"
            f"{dataframe_to_html_table(df=df, max_rows=200)}</div>"
        )
    body = "\n".join(sections)
    path.write_text(
        f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Direct versus clade taxonomic burden</title>
<style>{css}</style>
</head>
<body>
<h1>Direct versus clade taxonomic burden</h1>
<p class="note">
This report separates hierarchical clade-level burden from direct-assignment
burden. Clade-level counts describe the broader taxonomic report seen by a
user, but can count the same evidence at multiple taxonomic levels. Direct
counts provide a stricter signal-to-noise view of taxa assigned directly to the
expected species versus non-expected Plasmodium taxa.
</p>
{body}
</body>
</html>
""",
        encoding="utf-8",
    )


def main() -> None:
    """Run the direct-versus-clade burden summary."""
    args = parse_args()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    input_files = find_input_files(
        input_roots=args.input_roots,
        file_patterns=args.file_patterns,
    )
    manifest = pd.DataFrame({"input_file": [str(path) for path in input_files]})

    frames = [
        normalise_input_frame(path=path, taxon_regex=args.taxon_regex)
        for path in input_files
    ]
    frames = [frame for frame in frames if not frame.empty]
    all_rows = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    observations = summarise_observations(
        df=all_rows,
        min_metric_value=args.min_metric_value,
    )
    workflow_summary = summarise_group(
        observation_df=observations,
        group_cols=["workflow", "classifier", "metric_group"],
    )
    run_summary = summarise_group(
        observation_df=observations,
        group_cols=["workflow", "classifier", "run_name", "metric_group"],
    )
    spike_summary = summarise_group(
        observation_df=observations,
        group_cols=["workflow", "classifier", "metric_group", "spike_n"],
    )

    tables = {
        "manifest": manifest,
        "input_rows": all_rows,
        "observation_burden": observations,
        "workflow_summary": workflow_summary,
        "run_summary": run_summary,
        "spike_summary": spike_summary,
    }

    prefix = args.summary_name
    for name, table in tables.items():
        write_tsv(df=table, path=out_dir / f"{prefix}_{name}.tsv")

    write_excel(tables=tables, path=out_dir / f"{prefix}.xlsx")
    write_html_report(
        tables={
            "workflow_summary": workflow_summary,
            "run_summary": run_summary,
            "spike_summary": spike_summary,
            "observation_burden": observations,
            "manifest": manifest,
        },
        path=out_dir / f"{prefix}.html",
    )

    print(f"Wrote outputs to: {out_dir}")


if __name__ == "__main__":
    main()
