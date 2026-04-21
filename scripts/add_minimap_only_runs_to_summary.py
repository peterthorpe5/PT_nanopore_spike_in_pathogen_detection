#!/usr/bin/env python3
"""Add minimap-only run summaries into the existing spike-in summary tables.

This script reads one or more minimap-only run directories, converts their
per-run target and taxon sidecar TSVs into the ``combined_long.tsv`` and
``reported_taxa_long.tsv`` formats expected by the downstream reporting
workflow, and optionally merges them into an existing summary directory.
"""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path
from typing import Iterable

import pandas as pd


SUMMARY_PATTERNS = (
    "spikein_minimap_only_summary.tsv",
    "spikein_multi_minimap_only_summary.tsv",
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Fold minimap-only spike-in runs into combined_long.tsv and "
            "reported_taxa_long.tsv for the HTML reporting workflow."
        )
    )
    parser.add_argument(
        "--input_dirs",
        nargs="+",
        required=True,
        help=(
            "One or more run directories or parent directories to search "
            "recursively for minimap-only summary TSV files."
        ),
    )
    parser.add_argument(
        "--existing_summary_dir",
        default=None,
        help=(
            "Optional existing summary directory containing combined_long.tsv "
            "and reported_taxa_long.tsv. If supplied, its contents are copied "
            "into --out_dir before merged tables are written."
        ),
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for the merged summary tables.",
    )
    parser.add_argument(
        "--workflow_type",
        default="minimap2_readlevel_only",
        help="Workflow label to assign to the imported minimap-only rows.",
    )
    parser.add_argument(
        "--metric_name",
        default="minimap2_target_alignments",
        help="Metric name to use for target-level rows in combined_long.tsv.",
    )
    parser.add_argument(
        "--reported_metric_name",
        default="minimap2_reported_clade_reads",
        help=(
            "Metric name to use for taxon-level rows in reported_taxa_long.tsv. "
            "It should end with _reported_clade_reads for the current wrapper."
        ),
    )
    return parser.parse_args()


def discover_summary_files(input_dirs: Iterable[str]) -> list[Path]:
    """Discover minimap-only summary TSV files under the supplied paths.

    Parameters
    ----------
    input_dirs : Iterable[str]
        Run directories or parent directories.

    Returns
    -------
    list[Path]
        Sorted unique summary TSV paths.
    """
    discovered: set[Path] = set()
    for input_dir in input_dirs:
        base_path = Path(input_dir).expanduser().resolve()
        if base_path.is_file() and base_path.name in SUMMARY_PATTERNS:
            discovered.add(base_path)
            continue
        if not base_path.exists():
            continue
        for pattern in SUMMARY_PATTERNS:
            discovered.update(base_path.rglob(pattern))
    return sorted(discovered)


def infer_is_shuffled(summary_row: pd.Series, run_name: str) -> bool:
    """Infer whether a run represents a shuffled control.

    Parameters
    ----------
    summary_row : pd.Series
        One row from the minimap-only summary TSV.
    run_name : str
        Name of the run directory.

    Returns
    -------
    bool
        True if the run appears to be a shuffled control.
    """
    if "is_shuffled_control" in summary_row.index:
        value = str(summary_row["is_shuffled_control"]).strip().lower()
        return value in {"true", "1", "yes"}
    return "shuffle" in run_name.lower() or "shuffled" in run_name.lower()


def load_target_rows(
    target_summary_tsv: Path,
    summary_row: pd.Series,
    run_name: str,
    workflow_type: str,
    metric_name: str,
) -> list[dict[str, object]]:
    """Convert one target summary TSV into combined_long rows.

    Parameters
    ----------
    target_summary_tsv : Path
        Per-run target summary TSV.
    summary_row : pd.Series
        One row from the minimap-only run summary.
    run_name : str
        Name of the run directory.
    workflow_type : str
        Workflow label for the imported rows.
    metric_name : str
        Metric name for the imported rows.

    Returns
    -------
    list[dict[str, object]]
        Combined-long rows for each target label.
    """
    target_df = pd.read_csv(target_summary_tsv, sep="\t")
    rows: list[dict[str, object]] = []
    is_shuffled_control = infer_is_shuffled(summary_row=summary_row, run_name=run_name)
    spike_n = int(summary_row["spike_n"])
    total_spike_n = int(summary_row.get("total_spiked_reads", spike_n))
    n_genomes = int(summary_row.get("n_genomes", 1))
    replicate = int(summary_row["replicate"])

    for _, target_row in target_df.iterrows():
        rows.append(
            {
                "workflow_type": workflow_type,
                "classifier": "minimap2",
                "run_name": run_name,
                "target_label": str(target_row["target_label"]),
                "metric_name": metric_name,
                "metric_value": int(target_row.get("target_alignments", 0)),
                "replicate": replicate,
                "spike_n_per_genome": spike_n,
                "total_spike_n": total_spike_n,
                "n_genomes": n_genomes,
                "report_tsv": str(target_summary_tsv),
                "is_shuffled_control": is_shuffled_control,
            }
        )
    return rows


def load_reported_taxa_rows(
    reported_taxa_tsv: Path,
    summary_row: pd.Series,
    run_name: str,
    workflow_type: str,
    reported_metric_name: str,
) -> list[dict[str, object]]:
    """Convert one reported-taxa TSV into reported_taxa_long rows.

    Parameters
    ----------
    reported_taxa_tsv : Path
        Per-run taxon summary TSV.
    summary_row : pd.Series
        One row from the minimap-only run summary.
    run_name : str
        Name of the run directory.
    workflow_type : str
        Workflow label for the imported rows.
    reported_metric_name : str
        Metric name for the imported rows.

    Returns
    -------
    list[dict[str, object]]
        Reported-taxa-long rows.
    """
    taxa_df = pd.read_csv(reported_taxa_tsv, sep="\t")
    rows: list[dict[str, object]] = []
    is_shuffled_control = infer_is_shuffled(summary_row=summary_row, run_name=run_name)
    spike_n = int(summary_row["spike_n"])
    total_spike_n = int(summary_row.get("total_spiked_reads", spike_n))
    n_genomes = int(summary_row.get("n_genomes", 1))
    replicate = int(summary_row["replicate"])

    for _, taxon_row in taxa_df.iterrows():
        rows.append(
            {
                "workflow_type": workflow_type,
                "classifier": "minimap2",
                "run_name": run_name,
                "taxon_name": str(taxon_row["taxon_name"]),
                "taxid": str(taxon_row.get("taxid", "unknown_taxid")),
                "metric_name": reported_metric_name,
                "metric_value": int(taxon_row.get("alignments", 0)),
                "replicate": replicate,
                "spike_n_per_genome": spike_n,
                "total_spike_n": total_spike_n,
                "n_genomes": n_genomes,
                "matched_target_labels": "" if pd.isna(taxon_row.get("matched_target_labels", "")) else str(
                    taxon_row.get("matched_target_labels", "")
                ),
                "report_tsv": str(reported_taxa_tsv),
                "is_shuffled_control": is_shuffled_control,
            }
        )
    return rows


def convert_minimap_runs(
    summary_paths: Iterable[Path],
    workflow_type: str,
    metric_name: str,
    reported_metric_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convert all discovered minimap-only run summaries.

    Parameters
    ----------
    summary_paths : Iterable[Path]
        Summary TSV paths.
    workflow_type : str
        Workflow label for target and taxon rows.
    metric_name : str
        Target metric name.
    reported_metric_name : str
        Reported-taxa metric name.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Combined-long rows and reported-taxa-long rows.
    """
    combined_rows: list[dict[str, object]] = []
    reported_rows: list[dict[str, object]] = []

    for summary_path in summary_paths:
        summary_df = pd.read_csv(summary_path, sep="\t")
        run_name = summary_path.parent.name
        for _, summary_row in summary_df.iterrows():
            target_summary_tsv = Path(summary_row["minimap_target_summary_tsv"]).expanduser().resolve()
            reported_taxa_tsv = Path(summary_row["minimap_reported_taxa_tsv"]).expanduser().resolve()
            combined_rows.extend(
                load_target_rows(
                    target_summary_tsv=target_summary_tsv,
                    summary_row=summary_row,
                    run_name=run_name,
                    workflow_type=workflow_type,
                    metric_name=metric_name,
                )
            )
            reported_rows.extend(
                load_reported_taxa_rows(
                    reported_taxa_tsv=reported_taxa_tsv,
                    summary_row=summary_row,
                    run_name=run_name,
                    workflow_type=workflow_type,
                    reported_metric_name=reported_metric_name,
                )
            )

    combined_df = pd.DataFrame(combined_rows)
    reported_df = pd.DataFrame(reported_rows)
    return combined_df, reported_df


def copy_existing_summary(existing_summary_dir: Path, out_dir: Path) -> None:
    """Copy an existing summary directory into the output directory.

    Parameters
    ----------
    existing_summary_dir : Path
        Existing summary directory.
    out_dir : Path
        Output directory.
    """
    if out_dir.exists() and out_dir.samefile(existing_summary_dir):
        return
    shutil.copytree(existing_summary_dir, out_dir, dirs_exist_ok=True)


def read_existing_table(path: Path) -> pd.DataFrame:
    """Read an existing TSV if present, otherwise return an empty frame.

    Parameters
    ----------
    path : Path
        TSV path.

    Returns
    -------
    pd.DataFrame
        Existing table or an empty data frame.
    """
    if path.exists() and path.stat().st_size > 0:
        return pd.read_csv(path, sep="\t")
    return pd.DataFrame()


def merge_tables(existing_df: pd.DataFrame, new_df: pd.DataFrame) -> pd.DataFrame:
    """Merge existing and new rows while dropping exact duplicates.

    Parameters
    ----------
    existing_df : pd.DataFrame
        Existing table.
    new_df : pd.DataFrame
        Newly generated rows.

    Returns
    -------
    pd.DataFrame
        Merged table.
    """
    if existing_df.empty:
        merged_df = new_df.copy()
    elif new_df.empty:
        merged_df = existing_df.copy()
    else:
        merged_df = pd.concat([existing_df, new_df], ignore_index=True, sort=False)
    if merged_df.empty:
        return merged_df
    return merged_df.drop_duplicates().reset_index(drop=True)


def write_tsv(dataframe: pd.DataFrame, out_path: Path) -> None:
    """Write a TSV file.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Data to write.
    out_path : Path
        Output TSV path.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    dataframe.to_csv(out_path, sep="\t", index=False)


def write_manifest(summary_paths: list[Path], out_path: Path) -> None:
    """Write a manifest of minimap-only summary files that were imported.

    Parameters
    ----------
    summary_paths : list[Path]
        Imported summary paths.
    out_path : Path
        Manifest TSV path.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["summary_tsv"])
        for summary_path in summary_paths:
            writer.writerow([str(summary_path)])


def main() -> None:
    """Run the minimap-only summary merge workflow."""
    args = parse_args()

    summary_paths = discover_summary_files(args.input_dirs)
    if not summary_paths:
        raise FileNotFoundError(
            "No minimap-only summary TSV files were found under the supplied input directories."
        )

    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.existing_summary_dir:
        existing_summary_dir = Path(args.existing_summary_dir).expanduser().resolve()
        if not existing_summary_dir.exists():
            raise FileNotFoundError(
                f"Existing summary directory not found: {existing_summary_dir}"
            )
        copy_existing_summary(existing_summary_dir=existing_summary_dir, out_dir=out_dir)

    combined_new_df, reported_new_df = convert_minimap_runs(
        summary_paths=summary_paths,
        workflow_type=args.workflow_type,
        metric_name=args.metric_name,
        reported_metric_name=args.reported_metric_name,
    )

    combined_path = out_dir / "combined_long.tsv"
    reported_path = out_dir / "reported_taxa_long.tsv"

    combined_existing_df = read_existing_table(combined_path)
    reported_existing_df = read_existing_table(reported_path)

    combined_merged_df = merge_tables(
        existing_df=combined_existing_df,
        new_df=combined_new_df,
    )
    reported_merged_df = merge_tables(
        existing_df=reported_existing_df,
        new_df=reported_new_df,
    )

    write_tsv(combined_merged_df, combined_path)
    write_tsv(reported_merged_df, reported_path)
    write_manifest(summary_paths=summary_paths, out_path=out_dir / "minimap_only_import_manifest.tsv")


if __name__ == "__main__":
    main()
