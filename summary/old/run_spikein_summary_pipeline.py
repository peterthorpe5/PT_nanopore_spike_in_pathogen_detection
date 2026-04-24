#!/usr/bin/env python3
"""Run the full spike-in summary and reporting pipeline from one command.

This script consolidates the current multi-step summary workflow into a single
entry point that starts from a runs directory and produces the standard summary
outputs, method-performance tables, HTML reports, replicate-resolved report,
threshold-calibration report, combined real-world report, and optional Krona
inputs.

The underlying analytical scripts remain separate, but the user-facing
interface is reduced to one command. This keeps the existing functionality
while removing the need to remember or maintain a long chain of manual calls.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd


SCRIPT_NAMES = {
    "summariser": "summarise_spikein_runs_v4_reported_taxa_full.py",
    "method_performance": "build_method_performance_table_with_reported_taxa.py",
    "main_report": "make_spikein_report_v5.py",
    "replicate_report": "make_spikein_replicate_report_v2.py",
    "threshold_report": "make_spikein_threshold_calibration_report_v3.py",
    "real_world_report": "build_combined_real_world_report.py",
    "krona_inputs": "make_krona_inputs.py",
}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Run the complete ONT spike-in summary and report pipeline from a "
            "runs directory."
        )
    )
    parser.add_argument(
        "--runs_dir",
        required=True,
        help="Directory containing run output folders.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for the combined summary and reports.",
    )
    parser.add_argument(
        "--scripts_dir",
        default=None,
        help=(
            "Directory containing the existing summary/report scripts. "
            "Defaults to the directory containing this script."
        ),
    )
    parser.add_argument(
        "--title",
        default="ONT spike-in summary report",
        help="Title used in the HTML reports.",
    )
    parser.add_argument(
        "--threshold_mode",
        choices=["fixed", "baseline_mean", "baseline_max", "baseline_mean_plus_sd"],
        default="fixed",
        help="Thresholding rule used in downstream reports.",
    )
    parser.add_argument(
        "--min_detect_value",
        type=float,
        default=1.0,
        help="Minimum detection value used in downstream reports.",
    )
    parser.add_argument(
        "--sd_multiplier",
        type=float,
        default=2.0,
        help="Standard-deviation multiplier for thresholding rules that use it.",
    )
    parser.add_argument(
        "--target_fpr",
        type=float,
        default=0.05,
        help="Target false-positive rate for the threshold-calibration report.",
    )
    parser.add_argument(
        "--make_krona_inputs",
        action="store_true",
        help="Also generate Krona input tables from Kraken reports.",
    )
    parser.add_argument(
        "--krona_out_dir",
        default=None,
        help=(
            "Output directory for Krona input tables. Defaults to "
            "<out_dir>/krona_inputs."
        ),
    )
    parser.add_argument(
        "--krona_min_count",
        type=int,
        default=1,
        help="Minimum clade count for Krona input generation.",
    )
    parser.add_argument(
        "--copy_runs_into_out_dir",
        action="store_true",
        help=(
            "Copy the runs directory into the output directory before running. "
            "Useful if you want a self-contained report bundle."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print each command before running it.",
    )
    return parser.parse_args()


def resolve_scripts_dir(*, scripts_dir: str | None) -> Path:
    """Resolve the directory containing the pipeline scripts.

    Parameters
    ----------
    scripts_dir : str | None
        User-supplied scripts directory.

    Returns
    -------
    Path
        Resolved scripts directory.
    """
    if scripts_dir is not None:
        return Path(scripts_dir).expanduser().resolve()
    return Path(__file__).resolve().parent


def require_existing_paths(*, paths: Iterable[Path]) -> None:
    """Require that all supplied paths exist.

    Parameters
    ----------
    paths : Iterable[Path]
        Paths that must exist.

    Raises
    ------
    FileNotFoundError
        If any path is missing.
    """
    missing = [path for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Required path(s) not found: " + ", ".join(str(path) for path in missing)
        )


def run_command(*, command: list[str], verbose: bool) -> None:
    """Run a subprocess command and fail clearly if it exits non-zero.

    Parameters
    ----------
    command : list[str]
        Command to run.
    verbose : bool
        Whether to print the command first.

    Raises
    ------
    subprocess.CalledProcessError
        If the command exits non-zero.
    """
    if verbose:
        print("[INFO] Running:")
        print(" ".join(command))
    subprocess.run(command, check=True)


def write_pipeline_manifest(*, out_dir: Path, rows: list[dict[str, str]]) -> None:
    """Write a tab-separated pipeline manifest.

    Parameters
    ----------
    out_dir : Path
        Output directory.
    rows : list[dict[str, str]]
        Manifest rows.
    """
    manifest_df = pd.DataFrame(rows)
    manifest_df.to_csv(out_dir / "pipeline_manifest.tsv", sep="\t", index=False)


def maybe_copy_runs_dir(*, runs_dir: Path, out_dir: Path) -> Path:
    """Optionally copy the runs directory into the output bundle.

    Parameters
    ----------
    runs_dir : Path
        Original runs directory.
    out_dir : Path
        Output directory.

    Returns
    -------
    Path
        Runs directory to use downstream.
    """
    copied_runs_dir = out_dir / "runs"
    if copied_runs_dir.exists():
        shutil.rmtree(copied_runs_dir)
    shutil.copytree(src=runs_dir, dst=copied_runs_dir)
    return copied_runs_dir


def main() -> None:
    """Run the full spike-in summary pipeline."""
    args = parse_args()

    runs_dir = Path(args.runs_dir).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    scripts_dir = resolve_scripts_dir(scripts_dir=args.scripts_dir)

    out_dir.mkdir(parents=True, exist_ok=True)

    required_script_paths = {
        name: scripts_dir / script_name
        for name, script_name in SCRIPT_NAMES.items()
    }
    require_existing_paths(paths=[runs_dir, *required_script_paths.values()])

    active_runs_dir = runs_dir
    if args.copy_runs_into_out_dir:
        active_runs_dir = maybe_copy_runs_dir(runs_dir=runs_dir, out_dir=out_dir)

    summary_dir = out_dir
    replicate_out_dir = summary_dir / "replicate_resolved_report"
    threshold_out_dir = summary_dir / "threshold_calibration_report_v3"
    real_world_out_dir = summary_dir / "combined_real_world_report"
    krona_out_dir = (
        Path(args.krona_out_dir).expanduser().resolve()
        if args.krona_out_dir is not None
        else summary_dir / "krona_inputs"
    )

    manifest_rows: list[dict[str, str]] = []

    summariser_cmd = [
        sys.executable,
        str(required_script_paths["summariser"]),
        "--input_dirs",
        str(active_runs_dir),
        "--out_dir",
        str(summary_dir),
        "--verbose",
    ]
    run_command(command=summariser_cmd, verbose=args.verbose)
    manifest_rows.append(
        {
            "step": "summarise_runs",
            "output": str(summary_dir),
            "key_file": str(summary_dir / "combined_long.tsv"),
        }
    )

    combined_long_tsv = summary_dir / "combined_long.tsv"
    reported_taxa_long_tsv = summary_dir / "reported_taxa_long.tsv"
    require_existing_paths(paths=[combined_long_tsv])

    method_cmd = [
        sys.executable,
        str(required_script_paths["method_performance"]),
        "--combined_long_tsv",
        str(combined_long_tsv),
        "--out_dir",
        str(summary_dir),
        "--threshold_mode",
        args.threshold_mode,
        "--min_detect_value",
        str(args.min_detect_value),
        "--sd_multiplier",
        str(args.sd_multiplier),
    ]
    if reported_taxa_long_tsv.exists() and reported_taxa_long_tsv.stat().st_size > 0:
        method_cmd.extend([
            "--reported_taxa_long_tsv",
            str(reported_taxa_long_tsv),
        ])
    run_command(command=method_cmd, verbose=args.verbose)
    manifest_rows.append(
        {
            "step": "method_performance",
            "output": str(summary_dir),
            "key_file": str(summary_dir / "method_performance.xlsx"),
        }
    )

    main_report_cmd = [
        sys.executable,
        str(required_script_paths["main_report"]),
        "--summary_dir",
        str(summary_dir),
        "--title",
        args.title,
    ]
    run_command(command=main_report_cmd, verbose=args.verbose)
    manifest_rows.append(
        {
            "step": "main_html_report",
            "output": str(summary_dir),
            "key_file": str(summary_dir / "spikein_report_v5.html"),
        }
    )

    replicate_cmd = [
        sys.executable,
        str(required_script_paths["replicate_report"]),
        "--summary_dir",
        str(summary_dir),
        "--out_dir",
        str(replicate_out_dir),
        "--title",
        "ONT spike-in replicate-resolved report",
        "--threshold_mode",
        args.threshold_mode,
        "--min_detect_value",
        str(args.min_detect_value),
        "--sd_multiplier",
        str(args.sd_multiplier),
    ]
    run_command(command=replicate_cmd, verbose=args.verbose)
    manifest_rows.append(
        {
            "step": "replicate_report",
            "output": str(replicate_out_dir),
            "key_file": str(replicate_out_dir / "replicate_resolved_report.xlsx"),
        }
    )

    threshold_cmd = [
        sys.executable,
        str(required_script_paths["threshold_report"]),
        "--summary_dir",
        str(summary_dir),
        "--out_dir",
        str(threshold_out_dir),
        "--title",
        "ONT spike-in threshold calibration report",
        "--min_detect_value",
        str(args.min_detect_value),
        "--target_fpr",
        str(args.target_fpr),
        "--sd_multiplier",
        str(args.sd_multiplier),
    ]
    run_command(command=threshold_cmd, verbose=args.verbose)
    manifest_rows.append(
        {
            "step": "threshold_report",
            "output": str(threshold_out_dir),
            "key_file": str(threshold_out_dir / "threshold_calibration_report.xlsx"),
        }
    )

    method_performance_xlsx = summary_dir / "method_performance.xlsx"
    replicate_report_xlsx = replicate_out_dir / "replicate_resolved_report.xlsx"
    threshold_report_xlsx = threshold_out_dir / "threshold_calibration_report.xlsx"
    require_existing_paths(
        paths=[
            method_performance_xlsx,
            replicate_report_xlsx,
            threshold_report_xlsx,
        ]
    )

    real_world_cmd = [
        sys.executable,
        str(required_script_paths["real_world_report"]),
        "--method_performance_xlsx",
        str(method_performance_xlsx),
        "--replicate_report_xlsx",
        str(replicate_report_xlsx),
        "--threshold_report_xlsx",
        str(threshold_report_xlsx),
        "--out_dir",
        str(real_world_out_dir),
        "--report_title",
        f"{args.title} with real-world taxonomic burden",
    ]
    run_command(command=real_world_cmd, verbose=args.verbose)
    manifest_rows.append(
        {
            "step": "combined_real_world_report",
            "output": str(real_world_out_dir),
            "key_file": str(real_world_out_dir / "combined_real_world_report.xlsx"),
        }
    )

    if args.make_krona_inputs:
        krona_cmd = [
            sys.executable,
            str(required_script_paths["krona_inputs"]),
            "--input_dirs",
            str(active_runs_dir),
            "--out_dir",
            str(krona_out_dir),
            "--min_count",
            str(args.krona_min_count),
        ]
        run_command(command=krona_cmd, verbose=args.verbose)
        manifest_rows.append(
            {
                "step": "krona_inputs",
                "output": str(krona_out_dir),
                "key_file": str(krona_out_dir / "krona_manifest.tsv"),
            }
        )

    write_pipeline_manifest(out_dir=summary_dir, rows=manifest_rows)

    print(f"[INFO] Pipeline complete. Summary directory: {summary_dir}")
    print(f"[INFO] Main report: {summary_dir / 'spikein_report_v5.html'}")
    print(f"[INFO] Method workbook: {summary_dir / 'method_performance.xlsx'}")
    print(f"[INFO] Replicate report: {replicate_out_dir / 'replicate_resolved_report.html'}")
    print(f"[INFO] Threshold report: {threshold_out_dir / 'threshold_calibration_report.html'}")
    print(f"[INFO] Real-world report: {real_world_out_dir / 'combined_real_world_report.html'}")
    print(f"[INFO] Pipeline manifest: {summary_dir / 'pipeline_manifest.tsv'}")


if __name__ == "__main__":
    main()
