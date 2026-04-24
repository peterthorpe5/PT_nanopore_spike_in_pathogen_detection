#!/usr/bin/env python3
"""Run the full spike-in summary and reporting pipeline from one command.

This updated version optionally performs a minimap-specific rescue step before
building method-performance tables and HTML reports. The rescue step runs the
standalone minimap-specific summariser, then replaces generic minimap rows in
``combined_long.tsv`` and ``reported_taxa_long.tsv`` with the minimap-specific
rows so that all downstream reports use the corrected minimap stream.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd


SCRIPT_NAME_CHOICES = {
    "summariser": [
        "summarise_spikein_runs_v4_reported_taxa_full_minimap_patched.py",
        "summarise_spikein_runs_v4_reported_taxa_full.py",
    ],
    "method_performance": ["build_method_performance_table_with_reported_taxa.py"],
    "main_report": ["make_spikein_report_v5.py"],
    "replicate_report": ["make_spikein_replicate_report_v2.py"],
    "threshold_report": ["make_spikein_threshold_calibration_report_v3.py"],
    "real_world_report": ["build_combined_real_world_report.py"],
    "krona_inputs": ["make_krona_inputs.py"],
    "minimap_summary": ["minimap_specific_summary.py"],
}

BED_FILENAMES = (
    "minimap_sorted.filtered.bed",
    "minimap_sorted.MASKED.bed",
)


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
            "Directory containing the existing summary and report scripts. "
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
        "--include_minimap_rescue",
        action="store_true",
        help="Run the minimap-specific rescue summary and merge it upstream of all reports.",
    )
    parser.add_argument(
        "--minimap_input_roots",
        nargs="*",
        default=None,
        help=(
            "Optional minimap run roots. If omitted, minimap roots are discovered "
            "automatically from the runs directory."
        ),
    )
    parser.add_argument(
        "--minimap_out_dir",
        default=None,
        help=(
            "Output directory for minimap-specific tables. Defaults to "
            "<out_dir>/minimap_specific_summary."
        ),
    )
    parser.add_argument(
        "--single_target_label",
        default="Plasmodium vivax",
        help="Fallback single-target label for minimap rescue runs.",
    )
    parser.add_argument(
        "--panel2_tsv",
        default=None,
        help="Optional pathogen panel 2 TSV for minimap rescue.",
    )
    parser.add_argument(
        "--panel3_tsv",
        default=None,
        help="Optional pathogen panel 3 TSV for minimap rescue.",
    )
    parser.add_argument(
        "--target_threshold_alignments",
        type=int,
        default=1,
        help="Minimap rescue threshold for target alignments.",
    )
    parser.add_argument(
        "--target_threshold_unique_reads",
        type=int,
        default=1,
        help="Minimap rescue threshold for target unique reads.",
    )
    parser.add_argument(
        "--real_world_threshold_alignments",
        type=int,
        default=1,
        help="Minimap rescue threshold for real-world alignment presence.",
    )
    parser.add_argument(
        "--real_world_threshold_unique_reads",
        type=int,
        default=1,
        help="Minimap rescue threshold for real-world unique-read presence.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print each command before running.",
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


def resolve_script_paths(*, scripts_dir: Path) -> dict[str, Path]:
    """Resolve script paths, preferring newer script names when available.

    Parameters
    ----------
    scripts_dir : Path
        Directory containing pipeline scripts.

    Returns
    -------
    dict[str, Path]
        Resolved script paths.

    Raises
    ------
    FileNotFoundError
        If any required script cannot be resolved.
    """
    resolved: dict[str, Path] = {}
    missing: list[str] = []
    for role, candidates in SCRIPT_NAME_CHOICES.items():
        chosen = None
        for candidate in candidates:
            candidate_path = scripts_dir / candidate
            if candidate_path.exists():
                chosen = candidate_path
                break
        if chosen is None:
            missing.append(f"{role} ({', '.join(candidates)})")
        else:
            resolved[role] = chosen
    if missing:
        raise FileNotFoundError(
            "Required script(s) not found: " + "; ".join(missing)
        )
    return resolved


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


def discover_minimap_run_roots(*, runs_dir: Path) -> list[Path]:
    """Discover minimap run roots from BED files.

    Parameters
    ----------
    runs_dir : Path
        Runs directory to scan recursively.

    Returns
    -------
    list[Path]
        Sorted unique minimap run roots.
    """
    roots: set[Path] = set()
    for bed_filename in BED_FILENAMES:
        for bed_path in runs_dir.rglob(bed_filename):
            mix_dir = bed_path.parent
            run_root = mix_dir.parent
            if run_root.name == "main_run":
                run_root = run_root.parent
            lower_path = str(run_root).lower()
            if "broken" in lower_path or "superseded" in lower_path:
                continue
            roots.add(run_root)
    return sorted(roots)


def resolve_panel_path(*, user_path: str | None, scripts_dir: Path, panel_name: str) -> Path:
    """Resolve a panel TSV path.

    Parameters
    ----------
    user_path : str | None
        User-supplied panel path.
    scripts_dir : Path
        Scripts directory.
    panel_name : str
        File name of the panel TSV.

    Returns
    -------
    Path
        Resolved panel TSV path.
    """
    if user_path:
        return Path(user_path).expanduser().resolve()
    candidate = scripts_dir.parent / "configs" / panel_name
    return candidate.resolve()


def remove_generic_minimap_rows_from_combined_long(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Remove generic minimap rows from the main combined-long table.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Existing combined-long table.

    Returns
    -------
    pd.DataFrame
        Combined-long table with generic minimap rows removed.
    """
    out_df = dataframe.copy()
    metric_column = "metric" if "metric" in out_df.columns else "metric_name"
    if metric_column not in out_df.columns:
        return out_df
    metric_series = out_df[metric_column].fillna("").astype(str)
    keep_mask = ~metric_series.str.startswith("minimap_")
    return out_df.loc[keep_mask].copy()


def remove_generic_minimap_rows_from_reported_taxa(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Remove generic minimap rows from the main reported-taxa table.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Existing reported-taxa long table.

    Returns
    -------
    pd.DataFrame
        Reported-taxa table with generic minimap rows removed.
    """
    out_df = dataframe.copy()
    classifier_mask = pd.Series(False, index=out_df.index)
    metric_mask = pd.Series(False, index=out_df.index)
    if "classifier" in out_df.columns:
        classifier_mask = out_df["classifier"].fillna("").astype(str).str.lower().eq(
            "minimap"
        )
    if "metric_name" in out_df.columns:
        metric_mask = out_df["metric_name"].fillna("").astype(str).str.startswith(
            "minimap_"
        )
    keep_mask = ~(classifier_mask | metric_mask)
    return out_df.loc[keep_mask].copy()


def merge_minimap_outputs_into_summary(
    *,
    summary_dir: Path,
    minimap_out_dir: Path,
) -> dict[str, Path]:
    """Merge minimap-specific outputs into the main summary tables.

    Parameters
    ----------
    summary_dir : Path
        Main summary directory.
    minimap_out_dir : Path
        Minimap-specific summary directory.

    Returns
    -------
    dict[str, Path]
        Paths to the merged combined-long and reported-taxa tables.
    """
    combined_long_path = summary_dir / "combined_long.tsv"
    reported_taxa_path = summary_dir / "reported_taxa_long.tsv"
    minimap_combined_path = minimap_out_dir / "minimap_combined_long.tsv"
    minimap_reported_taxa_path = minimap_out_dir / "minimap_reported_taxa_long.tsv"

    require_existing_paths(paths=[combined_long_path, minimap_combined_path])

    combined_long_df = pd.read_csv(combined_long_path, sep="\t")
    minimap_combined_df = pd.read_csv(minimap_combined_path, sep="\t")
    merged_combined_df = pd.concat(
        [
            remove_generic_minimap_rows_from_combined_long(dataframe=combined_long_df),
            minimap_combined_df,
        ],
        ignore_index=True,
        sort=False,
    )

    combined_backup_path = summary_dir / "combined_long.pre_minimap.tsv"
    shutil.copyfile(combined_long_path, combined_backup_path)
    merged_combined_path = summary_dir / "combined_long.tsv"
    merged_combined_df.to_csv(merged_combined_path, sep="\t", index=False)

    if reported_taxa_path.exists() and reported_taxa_path.stat().st_size > 0:
        reported_taxa_df = pd.read_csv(reported_taxa_path, sep="\t")
    else:
        reported_taxa_df = pd.DataFrame()

    if minimap_reported_taxa_path.exists() and minimap_reported_taxa_path.stat().st_size > 0:
        minimap_reported_taxa_df = pd.read_csv(minimap_reported_taxa_path, sep="\t")
    else:
        minimap_reported_taxa_df = pd.DataFrame()

    if reported_taxa_path.exists():
        reported_backup_path = summary_dir / "reported_taxa_long.pre_minimap.tsv"
        shutil.copyfile(reported_taxa_path, reported_backup_path)

    if reported_taxa_df.empty:
        merged_reported_df = minimap_reported_taxa_df.copy()
    else:
        merged_reported_df = pd.concat(
            [
                remove_generic_minimap_rows_from_reported_taxa(
                    dataframe=reported_taxa_df,
                ),
                minimap_reported_taxa_df,
            ],
            ignore_index=True,
            sort=False,
        )

    merged_reported_path = summary_dir / "reported_taxa_long.tsv"
    merged_reported_df.to_csv(merged_reported_path, sep="\t", index=False)

    return {
        "combined_long": merged_combined_path,
        "reported_taxa_long": merged_reported_path,
        "combined_long_backup": combined_backup_path,
    }


def main() -> None:
    """Run the full spike-in summary pipeline."""
    args = parse_args()

    runs_dir = Path(args.runs_dir).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    scripts_dir = resolve_scripts_dir(scripts_dir=args.scripts_dir)

    out_dir.mkdir(parents=True, exist_ok=True)

    required_script_paths = resolve_script_paths(scripts_dir=scripts_dir)
    require_existing_paths(paths=[runs_dir])

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
    minimap_out_dir = (
        Path(args.minimap_out_dir).expanduser().resolve()
        if args.minimap_out_dir is not None
        else summary_dir / "minimap_specific_summary"
    )

    panel2_tsv = resolve_panel_path(
        user_path=args.panel2_tsv,
        scripts_dir=scripts_dir,
        panel_name="pathogen_panel_2.tsv",
    )
    panel3_tsv = resolve_panel_path(
        user_path=args.panel3_tsv,
        scripts_dir=scripts_dir,
        panel_name="pathogen_panel_3.tsv",
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

    minimap_roots: list[Path] = []
    if args.include_minimap_rescue:
        if args.minimap_input_roots:
            minimap_roots = [
                Path(path).expanduser().resolve() for path in args.minimap_input_roots
            ]
        else:
            minimap_roots = discover_minimap_run_roots(runs_dir=active_runs_dir)

        if minimap_roots:
            minimap_cmd = [
                sys.executable,
                str(required_script_paths["minimap_summary"]),
                "--input_roots",
                *[str(path) for path in minimap_roots],
                "--out_dir",
                str(minimap_out_dir),
                "--single_target_label",
                args.single_target_label,
                "--panel2_tsv",
                str(panel2_tsv),
                "--panel3_tsv",
                str(panel3_tsv),
                "--target_threshold_alignments",
                str(args.target_threshold_alignments),
                "--target_threshold_unique_reads",
                str(args.target_threshold_unique_reads),
                "--real_world_threshold_alignments",
                str(args.real_world_threshold_alignments),
                "--real_world_threshold_unique_reads",
                str(args.real_world_threshold_unique_reads),
            ]
            run_command(command=minimap_cmd, verbose=args.verbose)
            manifest_rows.append(
                {
                    "step": "minimap_specific_summary",
                    "output": str(minimap_out_dir),
                    "key_file": str(minimap_out_dir / "minimap_summary_run_info.tsv"),
                }
            )
            merge_paths = merge_minimap_outputs_into_summary(
                summary_dir=summary_dir,
                minimap_out_dir=minimap_out_dir,
            )
            manifest_rows.append(
                {
                    "step": "merge_minimap_into_summary",
                    "output": str(summary_dir),
                    "key_file": str(merge_paths["combined_long"]),
                }
            )
        else:
            manifest_rows.append(
                {
                    "step": "minimap_specific_summary",
                    "output": str(minimap_out_dir),
                    "key_file": "no_minimap_roots_found",
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
    print(
        f"[INFO] Replicate report: "
        f"{replicate_out_dir / 'replicate_resolved_report.html'}"
    )
    print(
        f"[INFO] Threshold report: "
        f"{threshold_out_dir / 'threshold_calibration_report.html'}"
    )
    print(
        f"[INFO] Real-world report: "
        f"{real_world_out_dir / 'combined_real_world_report.html'}"
    )
    if args.include_minimap_rescue:
        print(f"[INFO] Minimap summary: {minimap_out_dir}")
    print(f"[INFO] Pipeline manifest: {summary_dir / 'pipeline_manifest.tsv'}")


if __name__ == "__main__":
    main()
