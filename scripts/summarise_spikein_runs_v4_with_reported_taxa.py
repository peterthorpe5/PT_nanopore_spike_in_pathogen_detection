#!/usr/bin/env python3
"""
Summarise ONT spike-in runs and additionally collate reported-taxa long tables.

This wrapper keeps the existing harmonised summary outputs from
``summarise_spikein_runs_v4.py`` unchanged and adds one extra file:

- reported_taxa_long.tsv

That additional table is intended for off-target and cross-hit analysis and is
constructed by collating sidecar ``*.reported_taxa_long.tsv`` files written by
the Kraken and Metabuli helper scripts.

Design
------
- ``combined_long.tsv`` remains target-centric for benchmarking.
- ``reported_taxa_long.tsv`` is taxon-centric for off-target review.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import pandas as pd

import summarise_spikein_runs_v4 as base


REPORTED_TAXA_PATTERN = "*.reported_taxa_long.tsv"
MIX_DIR_PATTERN = re.compile(r"^mix_rep(?P<replicate>\d+)_n(?P<spike>\d+(?:\.\d+)?)$")


def discover_reported_taxa_files(input_dirs: list[str]) -> list[Path]:
    """
    Recursively discover reported-taxa sidecar files.

    Parameters
    ----------
    input_dirs : list[str]
        Input directories to scan.

    Returns
    -------
    list[Path]
        Sorted list of discovered reported-taxa TSV files.
    """
    found: list[Path] = []
    seen: set[Path] = set()

    for input_dir in input_dirs:
        base_dir = Path(input_dir).expanduser().resolve()
        if not base_dir.exists():
            continue
        for path in base_dir.rglob(REPORTED_TAXA_PATTERN):
            if path.is_file() and path not in seen:
                found.append(path)
                seen.add(path)

    return sorted(found)


def parse_mix_metadata_from_path(path: Path) -> tuple[Optional[int], Optional[float]]:
    """
    Parse replicate and per-genome spike level from a mix directory name.

    Parameters
    ----------
    path : Path
        Reported-taxa TSV path.

    Returns
    -------
    tuple[Optional[int], Optional[float]]
        Replicate and spike_n_per_genome, if parsable.
    """
    match = MIX_DIR_PATTERN.match(path.parent.name)
    if match is None:
        return None, None

    replicate = int(match.group("replicate"))
    spike_text = match.group("spike")
    spike_value = float(spike_text)
    return replicate, spike_value


def find_run_root_for_sidecar(
    sidecar_path: Path,
    known_run_dirs: set[str],
) -> Optional[Path]:
    """
    Find the nearest ancestor directory that matches a known run root.

    Parameters
    ----------
    sidecar_path : Path
        Reported-taxa sidecar path.
    known_run_dirs : set[str]
        Known summary run-root directory paths.

    Returns
    -------
    Optional[Path]
        Matching run-root directory, if found.
    """
    for ancestor in sidecar_path.parents:
        if str(ancestor) in known_run_dirs:
            return ancestor
    return None


def lookup_run_row(
    combined_wide: pd.DataFrame,
    run_root: Path,
    replicate: Optional[int],
    spike_n_per_genome: Optional[float],
) -> Optional[pd.Series]:
    """
    Look up the matching wide-table row for a reported-taxa sidecar.

    Parameters
    ----------
    combined_wide : pd.DataFrame
        Combined wide summary table.
    run_root : Path
        Run-root directory path.
    replicate : Optional[int]
        Parsed replicate number.
    spike_n_per_genome : Optional[float]
        Parsed spike level per genome.

    Returns
    -------
    Optional[pd.Series]
        Matching row, or None if no suitable row is found.
    """
    if combined_wide.empty:
        return None

    subset = combined_wide.loc[combined_wide["run_dir"].astype(str) == str(run_root)].copy()
    if subset.empty:
        return None

    if replicate is not None and "replicate" in subset.columns:
        subset = subset.loc[pd.to_numeric(subset["replicate"], errors="coerce") == float(replicate)]

    if spike_n_per_genome is not None and "spike_n_per_genome" in subset.columns:
        subset = subset.loc[
            pd.to_numeric(subset["spike_n_per_genome"], errors="coerce") == float(spike_n_per_genome)
        ]

    if subset.empty:
        return None

    return subset.iloc[0]


def parse_reported_taxa_file(
    path: Path,
    combined_wide: pd.DataFrame,
    issues: list[base.IssueRecord],
) -> Optional[pd.DataFrame]:
    """
    Parse one reported-taxa sidecar file and attach run-level metadata.

    Parameters
    ----------
    path : Path
        Reported-taxa sidecar TSV path.
    combined_wide : pd.DataFrame
        Combined wide summary table.
    issues : list[base.IssueRecord]
        Mutable list of issue records.

    Returns
    -------
    Optional[pd.DataFrame]
        Harmonised reported-taxa long table for this sidecar.
    """
    dataframe = base.safe_read_tsv(path=path, issues=issues)
    if dataframe is None or dataframe.empty:
        return None

    known_run_dirs = set(combined_wide["run_dir"].astype(str)) if not combined_wide.empty else set()
    run_root = find_run_root_for_sidecar(path=path, known_run_dirs=known_run_dirs)

    replicate, spike_n_per_genome = parse_mix_metadata_from_path(path=path)
    matched_row = None
    if run_root is not None:
        matched_row = lookup_run_row(
            combined_wide=combined_wide,
            run_root=run_root,
            replicate=replicate,
            spike_n_per_genome=spike_n_per_genome,
        )

    dataframe = dataframe.copy()
    dataframe["reported_taxa_tsv"] = str(path)
    dataframe["run_root_dir"] = str(run_root) if run_root is not None else pd.NA
    dataframe["run_name"] = (
        str(matched_row.get("run_name"))
        if matched_row is not None
        else (run_root.name if run_root is not None else path.parent.name)
    )
    dataframe["run_dir"] = (
        str(matched_row.get("run_dir"))
        if matched_row is not None
        else (str(run_root) if run_root is not None else pd.NA)
    )
    dataframe["workflow_type"] = (
        str(matched_row.get("workflow_type"))
        if matched_row is not None
        else pd.NA
    )
    dataframe["summary_path"] = (
        str(matched_row.get("summary_path"))
        if matched_row is not None
        else pd.NA
    )
    dataframe["replicate"] = (
        matched_row.get("replicate")
        if matched_row is not None
        else replicate
    )
    dataframe["spike_n_per_genome"] = (
        matched_row.get("spike_n_per_genome")
        if matched_row is not None
        else spike_n_per_genome
    )
    dataframe["total_spike_n"] = (
        matched_row.get("total_spike_n")
        if matched_row is not None
        else pd.NA
    )
    dataframe["n_genomes"] = (
        matched_row.get("n_genomes")
        if matched_row is not None
        else pd.NA
    )
    dataframe["is_shuffled_control"] = (
        matched_row.get("is_shuffled_control")
        if matched_row is not None
        else False
    )

    numeric_columns = [
        "metric_value",
        "replicate",
        "spike_n_per_genome",
        "total_spike_n",
        "n_genomes",
        "percent",
    ]
    dataframe = base.coerce_numeric_columns(dataframe=dataframe, columns=numeric_columns)
    return dataframe


def main() -> None:
    """
    Run the discovery, parsing, collation, and plotting workflow.
    """
    args = base.parse_args()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    issues: list[base.IssueRecord] = []

    summary_paths = base.discover_summary_files(
        input_dirs=args.input_dirs,
        glob_patterns=args.glob_patterns,
        verbose=args.verbose,
    )

    if not summary_paths:
        issues.append(
            base.IssueRecord(
                "no_summary_files_found",
                str(out_dir),
                "No recognised summary TSV files were found under the input directories.",
            )
        )

    parsed_tables: list[pd.DataFrame] = []
    for summary_path in summary_paths:
        dataframe = base.parse_summary_file(
            path=summary_path,
            issues=issues,
            verbose=args.verbose,
        )
        if dataframe is not None and not dataframe.empty:
            parsed_tables.append(dataframe)

    run_manifest = base.build_run_manifest(dataframes=parsed_tables)

    if parsed_tables:
        combined_wide = pd.concat(parsed_tables, ignore_index=True, sort=False)
    else:
        combined_wide = pd.DataFrame()

    combined_long = (
        base.wide_to_long(dataframe=combined_wide)
        if not combined_wide.empty
        else pd.DataFrame()
    )

    plot_records = base.generate_plots(
        long_df=combined_long,
        out_dir=out_dir,
        plot_min_points=args.plot_min_points,
        issues=issues,
        verbose=args.verbose,
    )

    reported_taxa_paths = discover_reported_taxa_files(input_dirs=args.input_dirs)
    reported_taxa_tables: list[pd.DataFrame] = []
    for reported_taxa_path in reported_taxa_paths:
        dataframe = parse_reported_taxa_file(
            path=reported_taxa_path,
            combined_wide=combined_wide,
            issues=issues,
        )
        if dataframe is not None and not dataframe.empty:
            reported_taxa_tables.append(dataframe)

    if reported_taxa_tables:
        reported_taxa_long = pd.concat(
            reported_taxa_tables,
            ignore_index=True,
            sort=False,
        )
    else:
        reported_taxa_long = pd.DataFrame(
            columns=[
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
                "reported_taxa_tsv",
                "run_root_dir",
                "run_name",
                "run_dir",
                "workflow_type",
                "summary_path",
                "replicate",
                "spike_n_per_genome",
                "total_spike_n",
                "n_genomes",
                "is_shuffled_control",
            ]
        )

    issues_df = pd.DataFrame.from_records(
        [
            {
                "category": issue.category,
                "path": issue.path,
                "detail": issue.detail,
            }
            for issue in issues
        ]
    )
    if issues_df.empty:
        issues_df = pd.DataFrame(columns=["category", "path", "detail"])

    plot_manifest_df = pd.DataFrame.from_records(
        [
            {
                "workflow_type": record.workflow_type,
                "metric_name": record.metric_name,
                "target_label": record.target_label,
                "png_path": record.png_path,
                "pdf_path": record.pdf_path,
                "n_points": record.n_points,
            }
            for record in plot_records
        ]
    )
    if plot_manifest_df.empty:
        plot_manifest_df = pd.DataFrame(
            columns=[
                "workflow_type",
                "metric_name",
                "target_label",
                "png_path",
                "pdf_path",
                "n_points",
            ]
        )

    base.write_tsv(run_manifest, out_dir / "run_manifest.tsv")
    base.write_tsv(combined_wide, out_dir / "combined_wide.tsv")
    base.write_tsv(combined_long, out_dir / "combined_long.tsv")
    base.write_tsv(reported_taxa_long, out_dir / "reported_taxa_long.tsv")
    base.write_tsv(issues_df, out_dir / "missing_or_problematic_files.tsv")
    base.write_tsv(plot_manifest_df, out_dir / "plot_manifest.tsv")

    print(f"[INFO] Summary complete. Output directory: {out_dir}")
    print(f"[INFO] Run manifest: {out_dir / 'run_manifest.tsv'}")
    print(f"[INFO] Combined wide table: {out_dir / 'combined_wide.tsv'}")
    print(f"[INFO] Combined long table: {out_dir / 'combined_long.tsv'}")
    print(f"[INFO] Reported taxa long table: {out_dir / 'reported_taxa_long.tsv'}")
    print(f"[INFO] Issues table: {out_dir / 'missing_or_problematic_files.tsv'}")
    print(f"[INFO] Plot manifest: {out_dir / 'plot_manifest.tsv'}")
    print(f"[INFO] N summary files found: {len(summary_paths)}")
    print(f"[INFO] N parsed tables: {len(parsed_tables)}")
    print(f"[INFO] N reported-taxa files found: {len(reported_taxa_paths)}")
    print(f"[INFO] N generated plots: {len(plot_records)}")


if __name__ == "__main__":
    main()
