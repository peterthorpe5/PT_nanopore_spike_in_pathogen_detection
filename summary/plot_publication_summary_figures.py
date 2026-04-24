#!/usr/bin/env python3
"""Create publication-quality summary figures for the ONT spike-in benchmark.

This script produces a compact set of cleaner, more manuscript-ready plots from
existing summary workbooks and optional minimap-specific summary tables.

Compared with the earlier plotting helper, this version focuses on:
- improved typography and spacing
- clearer grouping of methods and species
- log-scaled first-detection heatmaps for wide dynamic ranges
- a more interpretable Flye versus Flye+Medaka comparison
- segmented threshold-band plots with explicit background/grey/likely bands
- export to SVG, PDF, and PNG for each figure
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D


SPECIES_ORDER = [
    "Plasmodium falciparum",
    "Plasmodium vivax",
    "Plasmodium knowlesi",
]
COMPLEXITY_ORDER = ["single", "panel2", "panel3"]
METHOD_ORDER = [
    "Kraken2 single-read",
    "Kraken2 multi-read",
    "Flye",
    "Flye+Medaka",
    "Minimap focused single-read",
    "Minimap focused multi-read",
]
FAMILY_COLOURS = {
    "Kraken2 single-read": "#35608D",
    "Kraken2 multi-read": "#4F81BD",
    "Kraken2 single-assembly": "#7F9E44",
    "Kraken2 multi-assembly": "#9BBB59",
    "Flye": "#7F9E44",
    "Flye+Medaka": "#C0504D",
    "Metabuli": "#D28C28",
    "Minimap focused single-read": "#7A68A6",
    "Minimap focused multi-read": "#5D4B8A",
}
BAND_COLOURS = {
    "background": "#d9d9d9",
    "grey": "#f0c987",
    "likely": "#8fbc8f",
}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Create publication-quality summary plots from the spike-in report "
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
        help="Output directory for figures and tables.",
    )
    parser.add_argument(
        "--png_dpi",
        type=int,
        default=300,
        help="Resolution for PNG export. Default: 300.",
    )
    return parser.parse_args()


def set_publication_style() -> None:
    """Apply a restrained, publication-oriented Matplotlib style."""
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "figure.titlesize": 12,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def short_species_name(*, species: str) -> str:
    """Shorten a Plasmodium species label for plotting.

    Parameters
    ----------
    species : str
        Full species label.

    Returns
    -------
    str
        Shortened species label.
    """
    mapping = {
        "Plasmodium falciparum": "P. falciparum",
        "Plasmodium vivax": "P. vivax",
        "Plasmodium knowlesi": "P. knowlesi",
    }
    return mapping.get(str(species), str(species))


def clean_axis(*, ax: plt.Axes, grid_axis: Optional[str] = "y") -> None:
    """Apply consistent axis formatting.

    Parameters
    ----------
    ax : plt.Axes
        Axis to format.
    grid_axis : Optional[str]
        Axis along which the grid should be drawn.
    """
    ax.set_facecolor("white")
    if grid_axis is not None:
        ax.grid(True, axis=grid_axis, color="#d9d9d9", linewidth=0.6, alpha=0.7)
    ax.tick_params(direction="out")


def save_figure_outputs(
    *, fig: plt.Figure, out_path: Path, png_dpi: int
) -> tuple[Path, Path, Path]:
    """Save a figure as SVG, PDF, and PNG.

    Parameters
    ----------
    fig : plt.Figure
        Figure object.
    out_path : Path
        Base output path, usually ending in ``.svg``.
    png_dpi : int
        DPI for PNG export.

    Returns
    -------
    tuple[Path, Path, Path]
        Paths to the SVG, PDF, and PNG outputs.
    """
    svg_path = out_path.with_suffix(".svg")
    pdf_path = out_path.with_suffix(".pdf")
    png_path = out_path.with_suffix(".png")
    fig.savefig(svg_path)
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=png_dpi)
    return svg_path, pdf_path, png_path


def load_workbook_sheet(*, path: Path, candidate_names: list[str]) -> pd.DataFrame:
    """Load the first matching worksheet from an Excel workbook.

    Parameters
    ----------
    path : Path
        Workbook path.
    candidate_names : list[str]
        Sheet names to try in order.

    Returns
    -------
    pd.DataFrame
        Loaded worksheet.

    Raises
    ------
    ValueError
        If none of the candidate sheet names are found.
    """
    workbook = pd.ExcelFile(path)
    for sheet_name in candidate_names:
        if sheet_name in workbook.sheet_names:
            return pd.read_excel(path, sheet_name=sheet_name)
    raise ValueError(
        f"None of the requested sheets were found in {path}: {candidate_names}"
    )


def infer_complexity(*, run_name: str) -> Optional[str]:
    """Infer panel complexity from a run name.

    Parameters
    ----------
    run_name : str
        Run name.

    Returns
    -------
    Optional[str]
        ``single``, ``panel2``, ``panel3``, or ``None``.
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
    if value in {"target", "total", "total unique"}:
        return None
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
        Optional classifier name.

    Returns
    -------
    str
        Readable workflow label.
    """
    workflow = str(workflow)
    metric = str(metric)
    run_name_lower = str(run_name).lower()
    classifier = None if classifier is None else str(classifier)

    if classifier == "minimap" or metric.startswith("minimap_"):
        if workflow == "single_read":
            return "Minimap focused single-read"
        return "Minimap focused multi-read"

    if classifier == "metabuli" or metric == "metabuli_target_found":
        return "Metabuli"

    if workflow == "single_assembly":
        return "Kraken2 single-assembly"
    if workflow == "multi_assembly":
        return "Kraken2 multi-assembly"

    if metric == "kraken_target_contigs":
        if "medaka" in run_name_lower:
            return "Flye+Medaka"
        return "Flye"

    if metric == "kraken_target_reads":
        if workflow == "single_read":
            return "Kraken2 single-read"
        return "Kraken2 multi-read"

    return f"{workflow} | {metric}"


def add_direct_labels(
    *,
    ax: plt.Axes,
    dataframe: pd.DataFrame,
    x_col: str,
    y_col: str,
    label_col: str,
    min_separation_fraction: float = 0.05,
) -> None:
    """Add direct labels to scatter points with simple overlap reduction.

    Parameters
    ----------
    ax : plt.Axes
        Axis to annotate.
    dataframe : pd.DataFrame
        Data containing positions and labels.
    x_col : str
        Column used for x positions.
    y_col : str
        Column used for y positions.
    label_col : str
        Column used for labels.
    min_separation_fraction : float
        Minimum separation between final label y-positions as a fraction of the
        current y-range.
    """
    if dataframe.empty:
        return

    data = dataframe.sort_values(by=y_col).copy()
    y_values = data[y_col].astype(float).to_numpy()
    y_min = float(np.nanmin(y_values))
    y_max = float(np.nanmax(y_values))
    y_range = max(y_max - y_min, 1.0)
    min_sep = y_range * min_separation_fraction

    adjusted_y: list[float] = []
    for value in y_values:
        if not adjusted_y:
            adjusted_y.append(value)
            continue
        adjusted_value = max(value, adjusted_y[-1] + min_sep)
        adjusted_y.append(adjusted_value)

    if adjusted_y:
        overflow = adjusted_y[-1] - (y_max + y_range * 0.12)
        if overflow > 0:
            adjusted_y = [value - overflow for value in adjusted_y]

    for (_, row), y_text in zip(data.iterrows(), adjusted_y):
        ax.annotate(
            row[label_col],
            xy=(row[x_col], row[y_col]),
            xytext=(row[x_col] + 0.02, y_text),
            textcoords="data",
            ha="left",
            va="center",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "none", "alpha": 0.85},
            arrowprops={"arrowstyle": "-", "lw": 0.5, "color": "#666666", "shrinkA": 0, "shrinkB": 0},
        )


def build_tradeoff_dataframe(
    *,
    real_world_df: pd.DataFrame,
    minimap_real_world_df: Optional[pd.DataFrame],
) -> pd.DataFrame:
    """Build the workflow-level trade-off summary table.

    Parameters
    ----------
    real_world_df : pd.DataFrame
        Workflow-level real-world summary table.
    minimap_real_world_df : Optional[pd.DataFrame]
        Optional minimap-specific real-world summary.

    Returns
    -------
    pd.DataFrame
        Trade-off summary table.
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
                "all_expected_sensitivity": row["all_expected_sensitivity"],
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
                    "all_expected_sensitivity": row["all_expected_sensitivity"],
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

    out_df = pd.DataFrame(rows).dropna(
        subset=["all_expected_sensitivity", "mean_off_target_taxa_positive"]
    )
    return out_df.sort_values(
        by=["all_expected_sensitivity", "mean_off_target_taxa_positive"],
        ascending=[False, True],
    ).reset_index(drop=True)


def make_tradeoff_plot(
    *, plot_df: pd.DataFrame, out_path: Path, png_dpi: int
) -> None:
    """Plot all-expected sensitivity versus mean off-target burden.

    Parameters
    ----------
    plot_df : pd.DataFrame
        Trade-off summary table.
    out_path : Path
        Output base path.
    png_dpi : int
        PNG resolution.
    """
    fig, ax = plt.subplots(figsize=(7.6, 5.4))

    for _, row in plot_df.iterrows():
        colour = FAMILY_COLOURS.get(row["label"], "#666666")
        marker = "s" if row["source"] == "minimap" else "o"
        ax.scatter(
            row["all_expected_sensitivity"],
            row["mean_off_target_taxa_positive"],
            s=70,
            marker=marker,
            color=colour,
            edgecolor="black",
            linewidth=0.5,
            zorder=3,
        )

    add_direct_labels(
        ax=ax,
        dataframe=plot_df,
        x_col="all_expected_sensitivity",
        y_col="mean_off_target_taxa_positive",
        label_col="label",
        min_separation_fraction=0.045,
    )

    ax.axvline(0.5, color="#bfbfbf", linestyle="--", linewidth=0.8, zorder=1)
    ax.axhline(2.0, color="#bfbfbf", linestyle="--", linewidth=0.8, zorder=1)
    ax.set_xlim(-0.02, 1.05)
    ymax = max(1.0, float(plot_df["mean_off_target_taxa_positive"].max()) * 1.15)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("Recovered expected species set (all-expected sensitivity)")
    ax.set_ylabel("Mean off-target taxa per positive sample")
    ax.set_title("Sensitivity versus taxonomic burden", pad=10)
    clean_axis(ax=ax, grid_axis="both")

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor="#666666", markeredgecolor="black", label="Main workflows"),
        Line2D([0], [0], marker="s", linestyle="", markerfacecolor="#666666", markeredgecolor="black", label="Focused minimap"),
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="upper right")

    fig.tight_layout()
    save_figure_outputs(fig=fig, out_path=out_path, png_dpi=png_dpi)
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
        Run-level tracked-target performance table.
    minimap_tracked_df : Optional[pd.DataFrame]
        Optional minimap tracked-target summary table.

    Returns
    -------
    pd.DataFrame
        Long-form table for the grouped heatmap.
    """
    rows: list[dict[str, object]] = []

    for _, row in run_performance_df.iterrows():
        target = normalise_target(target_label=row["target_label"])
        if target is None:
            continue
        if "minimap" in str(row["run_name"]).lower() or str(row["metric"]).startswith("minimap_"):
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
            }
        )

    if minimap_tracked_df is not None and not minimap_tracked_df.empty:
        subset = minimap_tracked_df.loc[
            minimap_tracked_df["metric"].astype(str) == "minimap_target_alignments"
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
                }
            )

    out_df = pd.DataFrame(rows)
    if out_df.empty:
        return out_df

    out_df = (
        out_df.groupby(["method", "species", "complexity"], as_index=False)
        .agg({"sensitivity": "max"})
    )
    out_df["method"] = pd.Categorical(
        out_df["method"], categories=METHOD_ORDER, ordered=True
    )
    out_df["species"] = pd.Categorical(
        out_df["species"], categories=SPECIES_ORDER, ordered=True
    )
    out_df["complexity"] = pd.Categorical(
        out_df["complexity"], categories=COMPLEXITY_ORDER, ordered=True
    )
    return out_df.sort_values(by=["method", "species", "complexity"]).reset_index(drop=True)


def make_species_complexity_plot(
    *, plot_df: pd.DataFrame, out_path: Path, png_dpi: int
) -> None:
    """Create a grouped heatmap of species sensitivity by panel complexity.

    Parameters
    ----------
    plot_df : pd.DataFrame
        Species-by-complexity summary table.
    out_path : Path
        Output base path.
    png_dpi : int
        PNG resolution.
    """
    methods = [
        method
        for method in METHOD_ORDER
        if method in plot_df["method"].astype(str).unique().tolist()
    ]

    row_meta: list[tuple[str, str]] = []
    data_rows: list[list[float]] = []
    for method in methods:
        for species in SPECIES_ORDER:
            subset = plot_df.loc[
                (plot_df["method"].astype(str) == method)
                & (plot_df["species"].astype(str) == species)
            ].copy()
            if subset.empty:
                continue
            pivot = (
                subset.pivot(index="species", columns="complexity", values="sensitivity")
                .reindex(index=[species], columns=COMPLEXITY_ORDER)
            )
            data_rows.append(pivot.to_numpy(dtype=float).flatten().tolist())
            row_meta.append((method, species))

    matrix = np.asarray(data_rows, dtype=float)
    cmap = mpl.colormaps["viridis"].copy()
    cmap.set_bad("#eeeeee")

    fig_height = max(4.8, 0.42 * len(row_meta) + 1.6)
    fig, ax = plt.subplots(figsize=(7.3, fig_height))
    im = ax.imshow(matrix, aspect="auto", vmin=0.0, vmax=1.0, cmap=cmap)

    ax.set_xticks(np.arange(len(COMPLEXITY_ORDER)))
    ax.set_xticklabels(["Single", "Panel 2", "Panel 3"])
    ax.set_yticks(np.arange(len(row_meta)))
    ax.set_yticklabels([short_species_name(species=species) for _, species in row_meta])
    ax.set_title("Species sensitivity across panel complexity", pad=10)
    clean_axis(ax=ax, grid_axis=None)

    method_start = 0
    for idx, method in enumerate(methods):
        n_rows = sum(1 for name, _ in row_meta if name == method)
        if n_rows == 0:
            continue
        y_mid = method_start + (n_rows - 1) / 2
        ax.text(
            -0.85,
            y_mid,
            method,
            ha="right",
            va="center",
            fontsize=8,
            fontweight="bold",
            transform=ax.transData,
        )
        if method_start > 0:
            ax.axhline(method_start - 0.5, color="#c7c7c7", linewidth=0.8)
        method_start += n_rows

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            text = "NA" if math.isnan(value) else f"{value:.2f}"
            ax.text(
                j,
                i,
                text,
                ha="center",
                va="center",
                fontsize=7.5,
                color="white" if (not math.isnan(value) and value < 0.55) else "black",
            )

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Sensitivity")
    fig.subplots_adjust(left=0.34, right=0.90)
    save_figure_outputs(fig=fig, out_path=out_path, png_dpi=png_dpi)
    plt.close(fig)


def build_first_detection_dataframe(*, run_thresholds_df: pd.DataFrame) -> pd.DataFrame:
    """Build a compact first-detection summary table.

    Parameters
    ----------
    run_thresholds_df : pd.DataFrame
        Run-level threshold summary table.

    Returns
    -------
    pd.DataFrame
        First-detection summary table.
    """
    rows: list[dict[str, object]] = []
    for _, row in run_thresholds_df.iterrows():
        target = normalise_target(target_label=row["target_label"])
        if target is None:
            continue
        if "minimap" in str(row["run_name"]).lower() or str(row["metric"]).startswith("minimap_"):
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
                "Any": row.get("first_detected_spike_any", math.nan),
                "50%": row.get("first_detected_spike_50pct", math.nan),
                "95%": row.get("first_detected_spike_95pct", math.nan),
                "All": row.get("first_detected_spike_all", math.nan),
            }
        )

    out_df = pd.DataFrame(rows)
    if out_df.empty:
        return out_df

    out_df = (
        out_df.groupby(["method", "complexity", "species"], as_index=False)
        .agg({"Any": "min", "50%": "min", "95%": "min", "All": "min"})
    )

    keep_methods = [
        "Kraken2 single-read",
        "Kraken2 multi-read",
        "Flye",
        "Flye+Medaka",
    ]
    out_df = out_df.loc[out_df["method"].isin(keep_methods)].copy()
    out_df["method"] = pd.Categorical(
        out_df["method"], categories=keep_methods, ordered=True
    )
    out_df["complexity"] = pd.Categorical(
        out_df["complexity"], categories=COMPLEXITY_ORDER, ordered=True
    )
    out_df["species"] = pd.Categorical(
        out_df["species"], categories=SPECIES_ORDER, ordered=True
    )
    return out_df.sort_values(by=["method", "complexity", "species"]).reset_index(drop=True)


def make_first_detection_plot(
    *, plot_df: pd.DataFrame, out_path: Path, png_dpi: int
) -> None:
    """Create a log-scaled heatmap of earliest reproducible detection.

    Parameters
    ----------
    plot_df : pd.DataFrame
        First-detection summary table.
    out_path : Path
        Output base path.
    png_dpi : int
        PNG resolution.
    """
    columns = ["Any", "50%", "95%", "All"]
    row_labels = [
        f"{row['complexity']} | {short_species_name(species=row['species'])}"
        for _, row in plot_df.iterrows()
    ]
    matrix = plot_df[columns].to_numpy(dtype=float)

    finite_values = matrix[np.isfinite(matrix)]
    vmin = max(1.0, float(np.nanmin(finite_values))) if finite_values.size else 1.0
    vmax = max(vmin, float(np.nanmax(finite_values))) if finite_values.size else 1.0
    cmap = mpl.colormaps["viridis"].copy()
    cmap.set_bad("#eeeeee")

    fig_height = max(4.8, 0.34 * len(plot_df) + 1.8)
    fig, ax = plt.subplots(figsize=(7.8, fig_height))
    im = ax.imshow(
        matrix,
        aspect="auto",
        cmap=cmap,
        norm=LogNorm(vmin=vmin, vmax=vmax),
    )
    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels(columns)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels)
    ax.set_title("Earliest reproducible detection", pad=10)
    clean_axis(ax=ax, grid_axis=None)

    method_start = 0
    for method in plot_df["method"].cat.categories:
        n_rows = int((plot_df["method"].astype(str) == str(method)).sum())
        if n_rows == 0:
            continue
        y_mid = method_start + (n_rows - 1) / 2
        ax.text(
            -1.05,
            y_mid,
            str(method),
            ha="right",
            va="center",
            fontsize=8,
            fontweight="bold",
            transform=ax.transData,
        )
        if method_start > 0:
            ax.axhline(method_start - 0.5, color="#c7c7c7", linewidth=0.8)
        method_start += n_rows

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            text = "NA" if math.isnan(value) else f"{int(value)}"
            colour = "black"
            if not math.isnan(value):
                norm_val = (math.log10(value) - math.log10(vmin)) / max(
                    1e-6, math.log10(vmax) - math.log10(vmin)
                )
                colour = "white" if norm_val < 0.45 else "black"
            ax.text(j, i, text, ha="center", va="center", fontsize=7.2, color=colour)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("First detected spike (log scale)")
    fig.subplots_adjust(left=0.38, right=0.90)
    save_figure_outputs(fig=fig, out_path=out_path, png_dpi=png_dpi)
    plt.close(fig)


def build_polishing_dataframe(*, polishing_df: pd.DataFrame) -> pd.DataFrame:
    """Reshape the polishing comparison sheet into a long plotting table.

    Parameters
    ----------
    polishing_df : pd.DataFrame
        Polishing comparison table.

    Returns
    -------
    pd.DataFrame
        Long-form polishing comparison data.
    """
    records: list[dict[str, object]] = []
    for _, row in polishing_df.iterrows():
        panel = str(row["panel"])
        for species in SPECIES_ORDER:
            safe = species.replace(" ", "_")
            for prefix, workflow in [("unpolished", "Flye"), ("polished", "Flye+Medaka")]:
                sensitivity = row.get(f"{prefix}_{safe}_sensitivity")
                burden = row.get(f"{prefix}_mean_off_target_taxa_positive")
                if pd.notna(sensitivity):
                    records.append(
                        {
                            "panel": panel,
                            "workflow": workflow,
                            "species": species,
                            "sensitivity": sensitivity,
                            "mean_off_target_taxa_positive": burden,
                        }
                    )
    return pd.DataFrame(records)


def make_polishing_plot(
    *, plot_df: pd.DataFrame, out_path: Path, png_dpi: int
) -> None:
    """Create a paired Flye versus Flye+Medaka comparison figure.

    Parameters
    ----------
    plot_df : pd.DataFrame
        Long-form polishing comparison data.
    out_path : Path
        Output base path.
    png_dpi : int
        PNG resolution.
    """
    panels = ["panel2", "panel3"]
    workflows = ["Flye", "Flye+Medaka"]
    colours = {"Flye": FAMILY_COLOURS["Flye"], "Flye+Medaka": FAMILY_COLOURS["Flye+Medaka"]}

    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.6), sharey="row")
    y_positions = np.arange(len(SPECIES_ORDER))

    for row_idx, panel in enumerate(panels):
        subset = plot_df.loc[plot_df["panel"] == panel].copy()
        for col_idx, metric in enumerate(["sensitivity", "mean_off_target_taxa_positive"]):
            ax = axes[row_idx, col_idx]
            clean_axis(ax=ax, grid_axis="x")
            for species_idx, species in enumerate(SPECIES_ORDER):
                species_df = subset.loc[subset["species"] == species].set_index("workflow")
                if set(workflows).issubset(species_df.index):
                    x1 = float(species_df.loc["Flye", metric])
                    x2 = float(species_df.loc["Flye+Medaka", metric])
                    ax.plot([x1, x2], [species_idx, species_idx], color="#b0b0b0", linewidth=1.0, zorder=1)
                for workflow in workflows:
                    if workflow not in species_df.index:
                        continue
                    value = float(species_df.loc[workflow, metric])
                    ax.scatter(
                        value,
                        species_idx,
                        s=46,
                        color=colours[workflow],
                        edgecolor="black",
                        linewidth=0.5,
                        zorder=3,
                        marker="o" if workflow == "Flye" else "s",
                    )
            ax.set_yticks(y_positions)
            ax.set_yticklabels([short_species_name(species=species) for species in SPECIES_ORDER])
            if metric == "sensitivity":
                ax.set_xlim(0, max(0.2, float(plot_df[metric].max()) * 1.15))
                ax.set_xlabel("Sensitivity")
            else:
                ax.set_xlim(0, float(plot_df[metric].max()) * 1.15)
                ax.set_xlabel("Mean off-target taxa")
            if col_idx == 0:
                ax.set_title(f"{panel}: tracked-target sensitivity", loc="left", fontsize=10)
            else:
                ax.set_title(f"{panel}: taxonomic burden", loc="left", fontsize=10)

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=colours["Flye"], markeredgecolor="black", label="Flye"),
        Line2D([0], [0], marker="s", linestyle="", markerfacecolor=colours["Flye+Medaka"], markeredgecolor="black", label="Flye+Medaka"),
    ]
    fig.legend(handles=legend_handles, frameon=False, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.01))
    fig.suptitle("Flye versus Flye+Medaka", y=1.04)
    fig.tight_layout()
    save_figure_outputs(fig=fig, out_path=out_path, png_dpi=png_dpi)
    plt.close(fig)


def build_threshold_dataframe(*, interpretive_df: pd.DataFrame) -> pd.DataFrame:
    """Select representative threshold-band rows.

    Parameters
    ----------
    interpretive_df : pd.DataFrame
        Interpretive-band table.

    Returns
    -------
    pd.DataFrame
        Representative threshold-band rows.
    """
    selected_rows: list[dict[str, object]] = []
    desired = [
        ("single_read", "kraken_target_reads", "Plasmodium vivax", "Kraken2 single-read | P. vivax"),
        ("multi_read", "kraken_target_reads", "target", "Kraken2 multi-read | target"),
        ("multi_read", "metabuli_target_found", "Plasmodium vivax clade", "Metabuli | P. vivax clade"),
        ("multi_read", "metabuli_target_found", "target", "Metabuli | target"),
        ("multi_assembly", "kraken_target_contigs", "target", "Flye/Flye+Medaka | target contigs"),
    ]

    for workflow, metric, target, display_label in desired:
        subset = interpretive_df.loc[
            (interpretive_df["workflow"].astype(str) == workflow)
            & (interpretive_df["metric"].astype(str) == metric)
            & (interpretive_df["target_label"].astype(str) == target)
        ].copy()
        if subset.empty:
            continue
        row = subset.iloc[0]
        selected_rows.append(
            {
                "label": display_label,
                "usually_background_at_or_below": row["usually_background_at_or_below"],
                "grey_zone_from": row["grey_zone_from"],
                "grey_zone_to": row["grey_zone_to"],
                "likely_true_signal_at_or_above": row["likely_true_signal_at_or_above"],
                "plain_language_recommendation": row.get("plain_language_recommendation", ""),
            }
        )

    return pd.DataFrame(selected_rows)


def make_threshold_plot(
    *, plot_df: pd.DataFrame, out_path: Path, png_dpi: int
) -> None:
    """Create a segmented threshold-band plot.

    Parameters
    ----------
    plot_df : pd.DataFrame
        Threshold-band summary table.
    out_path : Path
        Output base path.
    png_dpi : int
        PNG resolution.
    """
    plot_df = plot_df.copy()
    max_value = float(
        np.nanmax(
            plot_df[
                [
                    "usually_background_at_or_below",
                    "grey_zone_to",
                    "likely_true_signal_at_or_above",
                ]
            ].to_numpy(dtype=float)
        )
    )
    x_max = max(max_value + 2, 5)

    fig_height = max(4.2, 0.65 * len(plot_df) + 1.8)
    fig, ax = plt.subplots(figsize=(9.2, fig_height))
    y_positions = np.arange(len(plot_df))

    for idx, row in plot_df.iterrows():
        background_end = float(row["usually_background_at_or_below"])
        likely_start = float(row["likely_true_signal_at_or_above"])
        grey_from = row["grey_zone_from"]
        grey_to = row["grey_zone_to"]

        ax.broken_barh([(0, max(background_end, 0.05))], (idx - 0.32, 0.64), facecolors=BAND_COLOURS["background"], edgecolors="none")
        if pd.notna(grey_from) and pd.notna(grey_to) and float(grey_to) > float(grey_from):
            ax.broken_barh(
                [(float(grey_from), float(grey_to) - float(grey_from))],
                (idx - 0.32, 0.64),
                facecolors=BAND_COLOURS["grey"],
                edgecolors="none",
            )
        ax.broken_barh(
            [(likely_start, max(x_max - likely_start, 0.1))],
            (idx - 0.32, 0.64),
            facecolors=BAND_COLOURS["likely"],
            edgecolors="none",
            alpha=0.9,
        )

        ax.text(background_end, idx + 0.20, f"≤{int(background_end)}", fontsize=7.5, ha="center")
        ax.text(likely_start, idx - 0.24, f"≥{int(likely_start)}", fontsize=7.5, ha="center")

    clean_axis(ax=ax, grid_axis="x")
    ax.set_xlim(0, x_max)
    ax.set_ylim(-0.8, len(plot_df) - 0.2)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(plot_df["label"].tolist())
    ax.invert_yaxis()
    ax.set_xlabel("Raw signal value")
    ax.set_title("Representative threshold bands", pad=10)

    legend_handles = [
        Line2D([0], [0], color=BAND_COLOURS["background"], linewidth=7, label="Usually background"),
        Line2D([0], [0], color=BAND_COLOURS["grey"], linewidth=7, label="Grey zone"),
        Line2D([0], [0], color=BAND_COLOURS["likely"], linewidth=7, label="Likely true signal"),
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="upper right")

    fig.tight_layout()
    save_figure_outputs(fig=fig, out_path=out_path, png_dpi=png_dpi)
    plt.close(fig)


def write_manifest(*, manifest_rows: list[dict[str, str]], out_path: Path) -> None:
    """Write a tab-separated manifest of produced plots.

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
    set_publication_style()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Keep the method-performance path as an explicit required argument so the
    # command stays aligned with the benchmark pipeline, even though the current
    # plots rely mainly on the real-world, replicate, and threshold workbooks.
    _ = Path(args.method_performance_xlsx)
    real_world_path = Path(args.real_world_xlsx)
    replicate_path = Path(args.replicate_report_xlsx)
    threshold_path = Path(args.threshold_report_xlsx)

    real_world_df = load_workbook_sheet(
        path=real_world_path,
        candidate_names=["workflow_real_world_summary"],
    )
    polishing_df = load_workbook_sheet(
        path=real_world_path,
        candidate_names=["polishing_comparison"],
    )
    run_performance_df = load_workbook_sheet(
        path=replicate_path,
        candidate_names=["run_performance"],
    )
    run_thresholds_df = load_workbook_sheet(
        path=replicate_path,
        candidate_names=["run_thresholds"],
    )
    interpretive_df = load_workbook_sheet(
        path=threshold_path,
        candidate_names=["interpretive_bands", "threshold_interpretive_bands"],
    )

    minimap_tracked_df = None
    if args.minimap_tracked_tsv:
        minimap_tracked_df = pd.read_csv(args.minimap_tracked_tsv, sep="\t")

    minimap_real_world_df = None
    if args.minimap_real_world_tsv:
        minimap_real_world_df = pd.read_csv(args.minimap_real_world_tsv, sep="\t")

    manifest_rows: list[dict[str, str]] = []

    tradeoff_df = build_tradeoff_dataframe(
        real_world_df=real_world_df,
        minimap_real_world_df=minimap_real_world_df,
    )
    tradeoff_table_path = out_dir / "plot01_tradeoff_table.tsv"
    tradeoff_df.to_csv(tradeoff_table_path, sep="\t", index=False)
    tradeoff_plot_path = out_dir / "plot01_sensitivity_vs_taxonomic_burden.svg"
    make_tradeoff_plot(
        plot_df=tradeoff_df,
        out_path=tradeoff_plot_path,
        png_dpi=args.png_dpi,
    )
    manifest_rows.append(
        {
            "plot_id": "plot01",
            "title": "Sensitivity versus taxonomic burden",
            "figure_svg": str(tradeoff_plot_path),
            "figure_pdf": str(tradeoff_plot_path.with_suffix(".pdf")),
            "figure_png": str(tradeoff_plot_path.with_suffix(".png")),
            "table_tsv": str(tradeoff_table_path),
            "why_useful": "Single-figure summary of recovered expected signal versus off-target burden.",
        }
    )

    species_df = build_species_complexity_dataframe(
        run_performance_df=run_performance_df,
        minimap_tracked_df=minimap_tracked_df,
    )
    species_table_path = out_dir / "plot02_species_complexity_table.tsv"
    species_df.to_csv(species_table_path, sep="\t", index=False)
    species_plot_path = out_dir / "plot02_species_by_complexity_heatmap.svg"
    make_species_complexity_plot(
        plot_df=species_df,
        out_path=species_plot_path,
        png_dpi=args.png_dpi,
    )
    manifest_rows.append(
        {
            "plot_id": "plot02",
            "title": "Species sensitivity across panel complexity",
            "figure_svg": str(species_plot_path),
            "figure_pdf": str(species_plot_path.with_suffix(".pdf")),
            "figure_png": str(species_plot_path.with_suffix(".png")),
            "table_tsv": str(species_table_path),
            "why_useful": "Shows which species remain robust and which weaken as mixed-panel complexity rises.",
        }
    )

    first_df = build_first_detection_dataframe(run_thresholds_df=run_thresholds_df)
    first_table_path = out_dir / "plot03_first_detection_table.tsv"
    first_df.to_csv(first_table_path, sep="\t", index=False)
    first_plot_path = out_dir / "plot03_first_detection_heatmap.svg"
    make_first_detection_plot(
        plot_df=first_df,
        out_path=first_plot_path,
        png_dpi=args.png_dpi,
    )
    manifest_rows.append(
        {
            "plot_id": "plot03",
            "title": "Earliest reproducible detection",
            "figure_svg": str(first_plot_path),
            "figure_pdf": str(first_plot_path.with_suffix(".pdf")),
            "figure_png": str(first_plot_path.with_suffix(".png")),
            "table_tsv": str(first_table_path),
            "why_useful": "Distinguishes first signal from dependable detection across replicates.",
        }
    )

    polish_df = build_polishing_dataframe(polishing_df=polishing_df)
    polish_table_path = out_dir / "plot04_polishing_comparison_table.tsv"
    polish_df.to_csv(polish_table_path, sep="\t", index=False)
    polish_plot_path = out_dir / "plot04_flye_vs_medaka.svg"
    make_polishing_plot(
        plot_df=polish_df,
        out_path=polish_plot_path,
        png_dpi=args.png_dpi,
    )
    manifest_rows.append(
        {
            "plot_id": "plot04",
            "title": "Flye versus Flye+Medaka",
            "figure_svg": str(polish_plot_path),
            "figure_pdf": str(polish_plot_path.with_suffix(".pdf")),
            "figure_png": str(polish_plot_path.with_suffix(".png")),
            "table_tsv": str(polish_table_path),
            "why_useful": "Shows whether polishing changes tracked-target sensitivity and off-target burden.",
        }
    )

    threshold_df = build_threshold_dataframe(interpretive_df=interpretive_df)
    threshold_table_path = out_dir / "plot05_threshold_bands_table.tsv"
    threshold_df.to_csv(threshold_table_path, sep="\t", index=False)
    threshold_plot_path = out_dir / "plot05_threshold_bands.svg"
    make_threshold_plot(
        plot_df=threshold_df,
        out_path=threshold_plot_path,
        png_dpi=args.png_dpi,
    )
    manifest_rows.append(
        {
            "plot_id": "plot05",
            "title": "Representative threshold bands",
            "figure_svg": str(threshold_plot_path),
            "figure_pdf": str(threshold_plot_path.with_suffix(".pdf")),
            "figure_png": str(threshold_plot_path.with_suffix(".png")),
            "table_tsv": str(threshold_table_path),
            "why_useful": "Illustrates which low-count values are usually background, ambiguous, or more likely to be real.",
        }
    )

    write_manifest(manifest_rows=manifest_rows, out_path=out_dir / "plot_manifest.tsv")


if __name__ == "__main__":
    main()
