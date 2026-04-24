#!/usr/bin/env python3
"""Create publication-focused summary plots for the ONT spike-in benchmark.

This script is intentionally selective. Rather than reproducing every plot from
all reports, it generates a compact set of figures that carry the main story of
the benchmark:

1. Sensitivity versus taxonomic burden across workflow families.
2. Species-by-complexity sensitivity heatmap.
3. Earliest reliable detection heatmap.
4. Flye versus Flye+Medaka comparison.
5. Threshold-band summary for a small set of representative method-target pairs.

The script uses the existing summary workbooks produced by the reporting
pipeline and can optionally add the focused minimap reanalysis as a separate
layer.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


SPECIES_ORDER = [
    "Plasmodium falciparum",
    "Plasmodium vivax",
    "Plasmodium knowlesi",
]

COMPLEXITY_ORDER = ["single", "panel2", "panel3"]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Create publication-focused summary plots from the spike-in report "
            "workbooks and optional minimap-specific summary tables."
        )
    )
    parser.add_argument(
        "--method_performance_xlsx",
        required=True,
        help="Path to method_performance.xlsx.",
    )
    parser.add_argument(
        "--real_world_xlsx",
        required=True,
        help="Path to combined_real_world_report.xlsx.",
    )
    parser.add_argument(
        "--replicate_report_xlsx",
        required=True,
        help="Path to replicate_resolved_report.xlsx.",
    )
    parser.add_argument(
        "--threshold_report_xlsx",
        required=True,
        help="Path to threshold_calibration_report.xlsx.",
    )
    parser.add_argument(
        "--minimap_tracked_tsv",
        required=False,
        default=None,
        help="Optional minimap_tracked_target_performance.tsv.",
    )
    parser.add_argument(
        "--minimap_real_world_tsv",
        required=False,
        default=None,
        help="Optional minimap_real_world_reference_summary.tsv.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for figures and the manifest TSV.",
    )
    return parser.parse_args()


def infer_complexity(*, run_name: str) -> Optional[str]:
    """Infer the panel complexity from a run name.

    Parameters
    ----------
    run_name : str
        Run name.

    Returns
    -------
    Optional[str]
        One of ``single``, ``panel2``, ``panel3``, or ``None`` if the run
        cannot be classified.
    """
    value = str(run_name).lower()
    if "panel2" in value:
        return "panel2"
    if "panel3" in value:
        return "panel3"
    if "single" in value or "panel1" in value or value == "main_run":
        return "single"
    return None


def normalise_target(*, target_label: str) -> Optional[str]:
    """Normalise target labels to the three principal species.

    Parameters
    ----------
    target_label : str
        Raw target label.

    Returns
    -------
    Optional[str]
        Normalised species label, or ``None`` if the target should be skipped.
    """
    value = str(target_label)
    if value in SPECIES_ORDER:
        return value
    if value == "single_target":
        return "Plasmodium vivax"
    if value.endswith(" clade"):
        base = value.replace(" clade", "")
        if base in SPECIES_ORDER:
            return base
    return None


def workflow_display_label(
    *, workflow: str, metric: str, run_name: str, classifier: Optional[str] = None
) -> str:
    """Create a readable workflow label.

    Parameters
    ----------
    workflow : str
        Workflow family.
    metric : str
        Metric name.
    run_name : str
        Run name.
    classifier : Optional[str]
        Optional classifier name from real-world summaries.

    Returns
    -------
    str
        Readable workflow label.
    """
    workflow = str(workflow)
    metric = str(metric)
    run_name_lower = str(run_name).lower()
    classifier = None if classifier is None else str(classifier)

    if classifier == "minimap" or metric == "minimap_target_alignments":
        if workflow == "single_read":
            return "Minimap focused single-read"
        return "Minimap focused multi-read"

    if classifier == "metabuli" or metric == "metabuli_target_found":
        return "Metabuli"

    if metric == "kraken_target_reads":
        if workflow == "single_read":
            return "Kraken2 single-read"
        return "Kraken2 multi-read"

    if metric == "kraken_target_contigs":
        if "medaka" in run_name_lower:
            return "Flye+Medaka"
        return "Flye"

    return f"{workflow} | {metric}"


def load_workbook_sheet(*, path: Path, sheet_name: str) -> pd.DataFrame:
    """Load one Excel worksheet.

    Parameters
    ----------
    path : Path
        Path to workbook.
    sheet_name : str
        Worksheet name.

    Returns
    -------
    pd.DataFrame
        Loaded worksheet.
    """
    return pd.read_excel(path, sheet_name=sheet_name)


def build_tradeoff_dataframe(
    *,
    real_world_df: pd.DataFrame,
    minimap_real_world_df: Optional[pd.DataFrame],
) -> pd.DataFrame:
    """Build the workflow-level trade-off summary table.

    Parameters
    ----------
    real_world_df : pd.DataFrame
        ``workflow_real_world_summary`` sheet.
    minimap_real_world_df : Optional[pd.DataFrame]
        Optional minimap-specific real-world summary.

    Returns
    -------
    pd.DataFrame
        Table used for the trade-off scatter.
    """
    rows: list[dict[str, object]] = []

    for _, row in real_world_df.iterrows():
        label = workflow_display_label(
            workflow=row["workflow"],
            metric="generic",
            run_name="",
            classifier=row.get("classifier"),
        )
        rows.append(
            {
                "label": label,
                "tracked_or_expected_sensitivity": row["all_expected_sensitivity"],
                "mean_off_target_taxa_positive": row[
                    "mean_off_target_taxa_positive"
                ],
                "mean_reported_taxa_positive": row["mean_reported_taxa_positive"],
                "positive_off_target_rate": row["positive_off_target_rate"],
                "source": "main",
            }
        )

    if minimap_real_world_df is not None and not minimap_real_world_df.empty:
        subset = minimap_real_world_df.loc[
            minimap_real_world_df["presence_metric"].astype(str) == "alignments"
        ].copy()
        for _, row in subset.iterrows():
            label = workflow_display_label(
                workflow=row["workflow"],
                metric="minimap_target_alignments",
                run_name="focused",
                classifier=row.get("classifier"),
            )
            rows.append(
                {
                    "label": label,
                    "tracked_or_expected_sensitivity": row[
                        "all_expected_sensitivity"
                    ],
                    "mean_off_target_taxa_positive": row[
                        "mean_off_target_taxa_positive"
                    ],
                    "mean_reported_taxa_positive": row[
                        "mean_reported_taxa_positive"
                    ],
                    "positive_off_target_rate": row["positive_off_target_rate"],
                    "source": "minimap",
                }
            )

    out_df = pd.DataFrame(rows)
    out_df = out_df.dropna(
        subset=[
            "tracked_or_expected_sensitivity",
            "mean_off_target_taxa_positive",
        ]
    ).copy()
    out_df = out_df.sort_values(
        by=["tracked_or_expected_sensitivity", "mean_off_target_taxa_positive"],
        ascending=[False, True],
    )
    return out_df


def make_tradeoff_plot(*, plot_df: pd.DataFrame, out_path: Path) -> None:
    """Plot sensitivity versus off-target burden.

    Parameters
    ----------
    plot_df : pd.DataFrame
        Trade-off summary table.
    out_path : Path
        Output SVG path.
    """
    fig, ax = plt.subplots(figsize=(8.5, 6.5))

    size_scale = 50
    sources = plot_df["source"].unique().tolist()
    markers = {"main": "o", "minimap": "s"}

    for source in sources:
        subset = plot_df.loc[plot_df["source"] == source].copy()
        ax.scatter(
            subset["tracked_or_expected_sensitivity"],
            subset["mean_off_target_taxa_positive"],
            s=subset["mean_reported_taxa_positive"] * size_scale,
            marker=markers.get(source, "o"),
            alpha=0.85,
            label="Focused minimap" if source == "minimap" else "Main workflows",
        )
        for _, row in subset.iterrows():
            ax.annotate(
                row["label"],
                xy=(
                    row["tracked_or_expected_sensitivity"],
                    row["mean_off_target_taxa_positive"],
                ),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
            )

    ax.set_xlabel("Recovered expected species set (all-expected sensitivity)")
    ax.set_ylabel("Mean off-target taxa per positive sample")
    ax.set_title("Sensitivity versus taxonomic burden")
    ax.grid(True, alpha=0.25)
    ax.set_xlim(-0.02, 1.05)
    ax.set_ylim(bottom=0)

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="", label="Main workflows"),
        Line2D([0], [0], marker="s", linestyle="", label="Focused minimap"),
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def build_species_complexity_dataframe(
    *,
    run_performance_df: pd.DataFrame,
    minimap_tracked_df: Optional[pd.DataFrame],
) -> pd.DataFrame:
    """Build a species-by-complexity sensitivity table.

    Parameters
    ----------
    run_performance_df : pd.DataFrame
        ``run_performance`` sheet.
    minimap_tracked_df : Optional[pd.DataFrame]
        Optional minimap tracked-target table.

    Returns
    -------
    pd.DataFrame
        Long-form table with one row per method, complexity, and species.
    """
    rows: list[dict[str, object]] = []

    df = run_performance_df.copy()
    for _, row in df.iterrows():
        target = normalise_target(target_label=row["target_label"])
        if target is None:
            continue
        if str(row["metric"]) == "metabuli_target_found":
            continue
        complexity = infer_complexity(run_name=row["run_name"])
        if complexity is None:
            continue
        label = workflow_display_label(
            workflow=row["workflow"], metric=row["metric"], run_name=row["run_name"]
        )
        rows.append(
            {
                "method": label,
                "complexity": complexity,
                "species": target,
                "sensitivity": row["sensitivity"],
                "lod50": row["lod50_spike_n"],
            }
        )

    if minimap_tracked_df is not None and not minimap_tracked_df.empty:
        subset = minimap_tracked_df.loc[
            (minimap_tracked_df["metric"].astype(str) == "minimap_target_alignments")
        ].copy()
        for _, row in subset.iterrows():
            target = normalise_target(target_label=row["target_label"])
            if target is None:
                continue
            complexity = infer_complexity(run_name=row["run_name"])
            if complexity is None:
                continue
            rows.append(
                {
                    "method": workflow_display_label(
                        workflow=row["workflow"],
                        metric=row["metric"],
                        run_name=row["run_name"],
                        classifier="minimap",
                    ),
                    "complexity": complexity,
                    "species": target,
                    "sensitivity": row["sensitivity"],
                    "lod50": row["lod50_spike_n"],
                }
            )

    out_df = pd.DataFrame(rows)
    if out_df.empty:
        return out_df

    out_df = (
        out_df.groupby(["method", "complexity", "species"], as_index=False)
        .agg({"sensitivity": "max", "lod50": "min"})
    )

    preferred_order = [
        "Kraken2 single-read",
        "Kraken2 multi-read",
        "Flye",
        "Flye+Medaka",
        "Minimap focused single-read",
        "Minimap focused multi-read",
    ]
    out_df["method"] = pd.Categorical(
        out_df["method"], categories=preferred_order, ordered=True
    )
    out_df["complexity"] = pd.Categorical(
        out_df["complexity"], categories=COMPLEXITY_ORDER, ordered=True
    )
    out_df["species"] = pd.Categorical(
        out_df["species"], categories=SPECIES_ORDER, ordered=True
    )
    out_df = out_df.sort_values(by=["method", "complexity", "species"])
    return out_df


def make_species_complexity_plot(*, plot_df: pd.DataFrame, out_path: Path) -> None:
    """Create a heatmap of species sensitivity across panel complexity.

    Parameters
    ----------
    plot_df : pd.DataFrame
        Species-by-complexity long table.
    out_path : Path
        Output SVG path.
    """
    methods = [x for x in plot_df["method"].cat.categories if x in plot_df["method"].astype(str).unique()]
    n_methods = len(methods)
    fig, axes = plt.subplots(
        nrows=n_methods,
        ncols=1,
        figsize=(8.2, max(3.0, 1.35 * n_methods + 1.5)),
        squeeze=False,
    )

    vmin = 0.0
    vmax = 1.0

    for idx, method in enumerate(methods):
        ax = axes[idx, 0]
        subset = plot_df.loc[plot_df["method"].astype(str) == method].copy()
        pivot = (
            subset.pivot(index="species", columns="complexity", values="sensitivity")
            .reindex(index=SPECIES_ORDER, columns=COMPLEXITY_ORDER)
        )
        data = pivot.to_numpy(dtype=float)
        im = ax.imshow(data, aspect="auto", vmin=vmin, vmax=vmax)
        ax.set_xticks(range(len(COMPLEXITY_ORDER)))
        ax.set_xticklabels(COMPLEXITY_ORDER)
        ax.set_yticks(range(len(SPECIES_ORDER)))
        ax.set_yticklabels(SPECIES_ORDER)
        ax.set_title(str(method), fontsize=10, loc="left")

        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                value = data[i, j]
                text = "NA" if math.isnan(value) else f"{value:.2f}"
                ax.text(j, i, text, ha="center", va="center", fontsize=8)

    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85)
    cbar.set_label("Sensitivity")
    fig.suptitle("Species sensitivity across panel complexity", y=0.995)
    fig.tight_layout()
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def build_first_detection_dataframe(*, run_thresholds_df: pd.DataFrame) -> pd.DataFrame:
    """Build a compact first-detection table for selected representative runs.

    Parameters
    ----------
    run_thresholds_df : pd.DataFrame
        ``run_thresholds`` sheet.

    Returns
    -------
    pd.DataFrame
        Compact first-detection table.
    """
    df = run_thresholds_df.copy()
    rows: list[dict[str, object]] = []
    for _, row in df.iterrows():
        target = normalise_target(target_label=row["target_label"])
        if target is None:
            continue
        label = workflow_display_label(
            workflow=row["workflow"], metric=row["metric"], run_name=row["run_name"]
        )
        complexity = infer_complexity(run_name=row["run_name"])
        if complexity is None:
            continue
        if label not in {"Kraken2 single-read", "Kraken2 multi-read", "Flye", "Flye+Medaka"}:
            continue
        rows.append(
            {
                "label": f"{label} | {complexity} | {target}",
                "first_any": row["first_detected_spike_any"],
                "first_50pct": row["first_detected_spike_50pct"],
                "first_95pct": row["first_detected_spike_95pct"],
                "first_all": row["first_detected_spike_all"],
            }
        )
    out_df = pd.DataFrame(rows).drop_duplicates()
    return out_df.sort_values(by="label")


def make_first_detection_plot(*, plot_df: pd.DataFrame, out_path: Path) -> None:
    """Create a heatmap of earliest reproducible detection.

    Parameters
    ----------
    plot_df : pd.DataFrame
        First-detection table.
    out_path : Path
        Output SVG path.
    """
    columns = ["first_any", "first_50pct", "first_95pct", "first_all"]
    labels = ["Any", "50%", "95%", "All"]
    data = plot_df[columns].to_numpy(dtype=float)

    finite_values = data[np.isfinite(data)]
    vmax = np.max(finite_values) if finite_values.size else 1.0
    fig, ax = plt.subplots(figsize=(8.5, max(4.0, 0.28 * len(plot_df) + 1.8)))
    im = ax.imshow(data, aspect="auto", vmin=0, vmax=vmax)
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels(labels)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["label"].tolist(), fontsize=8)
    ax.set_title("Earliest reproducible detection")

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            value = data[i, j]
            text = "NA" if math.isnan(value) else f"{int(value)}"
            ax.text(j, i, text, ha="center", va="center", fontsize=8)

    cbar = fig.colorbar(im, ax=ax, shrink=0.9)
    cbar.set_label("First detected spike")
    fig.tight_layout()
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def build_polishing_dataframe(*, polishing_df: pd.DataFrame) -> pd.DataFrame:
    """Reshape the polishing comparison into a compact plotting table.

    Parameters
    ----------
    polishing_df : pd.DataFrame
        ``polishing_comparison`` sheet.

    Returns
    -------
    pd.DataFrame
        Long-form polishing comparison table.
    """
    records: list[dict[str, object]] = []
    for _, row in polishing_df.iterrows():
        panel = row["panel"]
        for species in SPECIES_ORDER:
            safe = species.replace(" ", "_")
            unpolished_sens = row.get(f"unpolished_{safe}_sensitivity")
            polished_sens = row.get(f"polished_{safe}_sensitivity")
            if pd.notna(unpolished_sens):
                records.append(
                    {
                        "panel": panel,
                        "species": species,
                        "workflow": "Flye",
                        "sensitivity": unpolished_sens,
                        "mean_off_target_taxa_positive": row[
                            "unpolished_mean_off_target_taxa_positive"
                        ],
                    }
                )
            if pd.notna(polished_sens):
                records.append(
                    {
                        "panel": panel,
                        "species": species,
                        "workflow": "Flye+Medaka",
                        "sensitivity": polished_sens,
                        "mean_off_target_taxa_positive": row[
                            "polished_mean_off_target_taxa_positive"
                        ],
                    }
                )
    return pd.DataFrame(records)


def make_polishing_plot(*, plot_df: pd.DataFrame, out_path: Path) -> None:
    """Plot Flye versus Flye+Medaka sensitivities with taxonomic burden.

    Parameters
    ----------
    plot_df : pd.DataFrame
        Long-form polishing comparison data.
    out_path : Path
        Output SVG path.
    """
    panels = ["panel2", "panel3"]
    workflows = ["Flye", "Flye+Medaka"]
    x_positions = np.arange(len(SPECIES_ORDER))
    width = 0.35

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex="col")

    for col_idx, panel in enumerate(panels):
        subset = plot_df.loc[plot_df["panel"] == panel].copy()
        sens_ax = axes[0, col_idx]
        burden_ax = axes[1, col_idx]

        for shift_idx, workflow in enumerate(workflows):
            sub = subset.loc[subset["workflow"] == workflow].set_index("species")
            sens_values = [sub.loc[s, "sensitivity"] if s in sub.index else np.nan for s in SPECIES_ORDER]
            burden_values = [
                sub.loc[s, "mean_off_target_taxa_positive"] if s in sub.index else np.nan
                for s in SPECIES_ORDER
            ]
            shift = (shift_idx - 0.5) * width
            sens_ax.bar(x_positions + shift, sens_values, width=width, label=workflow)
            burden_ax.bar(x_positions + shift, burden_values, width=width, label=workflow)

        sens_ax.set_title(panel)
        sens_ax.set_ylabel("Sensitivity")
        sens_ax.set_ylim(0, 1.05)
        sens_ax.grid(True, axis="y", alpha=0.25)

        burden_ax.set_ylabel("Mean off-target taxa")
        burden_ax.set_xticks(x_positions)
        burden_ax.set_xticklabels([x.replace("Plasmodium ", "") for x in SPECIES_ORDER])
        burden_ax.grid(True, axis="y", alpha=0.25)

    axes[0, 0].legend(frameon=False, loc="upper left")
    fig.suptitle("Flye versus Flye+Medaka", y=0.995)
    fig.tight_layout()
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def build_threshold_dataframe(*, interpretive_df: pd.DataFrame) -> pd.DataFrame:
    """Select a small set of representative threshold-band rows.

    Parameters
    ----------
    interpretive_df : pd.DataFrame
        ``interpretive_bands`` sheet.

    Returns
    -------
    pd.DataFrame
        Selected threshold-band table.
    """
    df = interpretive_df.copy()
    keep_rows = []
    for _, row in df.iterrows():
        workflow = str(row["workflow"])
        metric = str(row["metric"])
        target = normalise_target(target_label=row["target_label"])
        if target is None:
            continue

        label = workflow_display_label(
            workflow=workflow,
            metric=metric,
            run_name=f"{workflow}_{target}",
        )
        if label == "Kraken2 multi-read" and target in {"Plasmodium vivax", "Plasmodium knowlesi"}:
            keep_rows.append((label, target, row))
        elif label == "Kraken2 single-read" and target == "Plasmodium vivax":
            keep_rows.append((label, target, row))
        elif label == "Flye" and target == "Plasmodium falciparum":
            keep_rows.append((label, target, row))

    records = []
    for label, target, row in keep_rows:
        records.append(
            {
                "label": f"{label} | {target}",
                "usually_background_at_or_below": row[
                    "usually_background_at_or_below"
                ],
                "grey_zone_from": row["grey_zone_from"],
                "grey_zone_to": row["grey_zone_to"],
                "likely_true_signal_at_or_above": row[
                    "likely_true_signal_at_or_above"
                ],
            }
        )
    return pd.DataFrame(records)


def make_threshold_plot(*, plot_df: pd.DataFrame, out_path: Path) -> None:
    """Plot simplified threshold bands for selected representative methods.

    Parameters
    ----------
    plot_df : pd.DataFrame
        Threshold-band table.
    out_path : Path
        Output SVG path.
    """
    fig, ax = plt.subplots(figsize=(9, max(3.5, 0.55 * len(plot_df) + 1.8)))

    y_positions = np.arange(len(plot_df))
    for idx, row in plot_df.iterrows():
        bg = row["usually_background_at_or_below"]
        likely = row["likely_true_signal_at_or_above"]
        grey_from = row["grey_zone_from"]
        grey_to = row["grey_zone_to"]

        ax.hlines(y=idx, xmin=0, xmax=bg, linewidth=6, alpha=0.7)
        if pd.notna(grey_from) and pd.notna(grey_to):
            ax.hlines(y=idx, xmin=grey_from, xmax=grey_to, linewidth=6, alpha=0.5)
        ax.hlines(y=idx, xmin=likely, xmax=max(likely, likely + 1), linewidth=6, alpha=0.9)
        ax.text(bg, idx + 0.12, f"≤{int(bg)}", fontsize=8)
        ax.text(likely, idx - 0.18, f"≥{int(likely)}", fontsize=8)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(plot_df["label"].tolist(), fontsize=8)
    ax.set_xlabel("Raw signal value")
    ax.set_title("Representative threshold bands")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def write_manifest(*, manifest_rows: list[dict[str, str]], out_path: Path) -> None:
    """Write a tab-separated plot manifest.

    Parameters
    ----------
    manifest_rows : list[dict[str, str]]
        Manifest entries.
    out_path : Path
        Output path.
    """
    pd.DataFrame(manifest_rows).to_csv(out_path, sep="\t", index=False)


def main() -> None:
    """Run the plotting workflow."""
    args = parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    method_path = Path(args.method_performance_xlsx)
    real_world_path = Path(args.real_world_xlsx)
    replicate_path = Path(args.replicate_report_xlsx)
    threshold_path = Path(args.threshold_report_xlsx)

    real_world_df = load_workbook_sheet(
        path=real_world_path, sheet_name="workflow_real_world_summary"
    )
    run_performance_df = load_workbook_sheet(
        path=replicate_path, sheet_name="run_performance"
    )
    run_thresholds_df = load_workbook_sheet(
        path=replicate_path, sheet_name="run_thresholds"
    )
    polishing_df = load_workbook_sheet(
        path=real_world_path, sheet_name="polishing_comparison"
    )
    interpretive_df = load_workbook_sheet(
        path=threshold_path, sheet_name="interpretive_bands"
    )

    minimap_tracked_df = None
    if args.minimap_tracked_tsv:
        minimap_tracked_df = pd.read_csv(args.minimap_tracked_tsv, sep="\t")

    minimap_real_world_df = None
    if args.minimap_real_world_tsv:
        minimap_real_world_df = pd.read_csv(args.minimap_real_world_tsv, sep="\t")

    manifest_rows = []

    tradeoff_df = build_tradeoff_dataframe(
        real_world_df=real_world_df,
        minimap_real_world_df=minimap_real_world_df,
    )
    tradeoff_table_path = out_dir / "plot01_tradeoff_table.tsv"
    tradeoff_df.to_csv(tradeoff_table_path, sep="\t", index=False)
    tradeoff_plot_path = out_dir / "plot01_sensitivity_vs_taxonomic_burden.svg"
    make_tradeoff_plot(plot_df=tradeoff_df, out_path=tradeoff_plot_path)
    manifest_rows.append(
        {
            "plot_id": "plot01",
            "title": "Sensitivity versus taxonomic burden",
            "figure_path": str(tradeoff_plot_path),
            "table_path": str(tradeoff_table_path),
            "why_useful": (
                "Best single summary of the main paper story: early recovery "
                "versus taxonomic clutter."
            ),
        }
    )

    species_df = build_species_complexity_dataframe(
        run_performance_df=run_performance_df,
        minimap_tracked_df=minimap_tracked_df,
    )
    species_table_path = out_dir / "plot02_species_complexity_table.tsv"
    species_df.to_csv(species_table_path, sep="\t", index=False)
    species_plot_path = out_dir / "plot02_species_by_complexity_heatmap.svg"
    make_species_complexity_plot(plot_df=species_df, out_path=species_plot_path)
    manifest_rows.append(
        {
            "plot_id": "plot02",
            "title": "Species sensitivity across panel complexity",
            "figure_path": str(species_plot_path),
            "table_path": str(species_table_path),
            "why_useful": (
                "Shows the biologically important collapse of P. falciparum "
                "as complexity rises, while vivax and knowlesi stay stronger."
            ),
        }
    )

    first_df = build_first_detection_dataframe(run_thresholds_df=run_thresholds_df)
    first_table_path = out_dir / "plot03_first_detection_table.tsv"
    first_df.to_csv(first_table_path, sep="\t", index=False)
    first_plot_path = out_dir / "plot03_first_detection_heatmap.svg"
    make_first_detection_plot(plot_df=first_df, out_path=first_plot_path)
    manifest_rows.append(
        {
            "plot_id": "plot03",
            "title": "Earliest reproducible detection",
            "figure_path": str(first_plot_path),
            "table_path": str(first_table_path),
            "why_useful": (
                "Summarises operational usefulness by distinguishing first "
                "signal from dependable replicate-level detection."
            ),
        }
    )

    polish_df = build_polishing_dataframe(polishing_df=polishing_df)
    polish_table_path = out_dir / "plot04_polishing_comparison_table.tsv"
    polish_df.to_csv(polish_table_path, sep="\t", index=False)
    polish_plot_path = out_dir / "plot04_flye_vs_medaka.svg"
    make_polishing_plot(plot_df=polish_df, out_path=polish_plot_path)
    manifest_rows.append(
        {
            "plot_id": "plot04",
            "title": "Flye versus Flye+Medaka",
            "figure_path": str(polish_plot_path),
            "table_path": str(polish_table_path),
            "why_useful": (
                "Best compact summary of the polishing question: narrower "
                "taxonomy in places, but mixed sensitivity effects."
            ),
        }
    )

    threshold_df = build_threshold_dataframe(interpretive_df=interpretive_df)
    threshold_table_path = out_dir / "plot05_threshold_bands_table.tsv"
    threshold_df.to_csv(threshold_table_path, sep="\t", index=False)
    threshold_plot_path = out_dir / "plot05_threshold_bands.svg"
    make_threshold_plot(plot_df=threshold_df, out_path=threshold_plot_path)
    manifest_rows.append(
        {
            "plot_id": "plot05",
            "title": "Representative threshold bands",
            "figure_path": str(threshold_plot_path),
            "table_path": str(threshold_table_path),
            "why_useful": (
                "Shows where low-count calls are usually background, grey-zone, "
                "or more likely to be real."
            ),
        }
    )

    write_manifest(
        manifest_rows=manifest_rows,
        out_path=out_dir / "plot_manifest.tsv",
    )


if __name__ == "__main__":
    main()
