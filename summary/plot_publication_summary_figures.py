#!/usr/bin/env python3
"""Create publication-ready summary plots for the ONT spike-in benchmark.

This script replaces the earlier quick-look plotting layer with a more
selective and publication-focused set of figures. It keeps the same command
line interface, but the figures are redesigned to be clearer in print and more
useful in the manuscript.

The output set is:

1. Sensitivity versus taxonomic burden:
   a labelled trade-off scatter highlighting the main workflows.
2. Species sensitivity across panel complexity:
   a compact heatmap with grouped rows and clearer labels.
3. Earliest reproducible detection:
   an interval plot on a log-scale x-axis, showing when signal first appears
   and when it becomes dependable across replicates.
4. Flye versus Flye+Medaka:
   dumbbell plots comparing tracked-target sensitivity and panel-level
   off-target burden.
5. Representative threshold bands:
   horizontal interpretive bands for a small set of method-target examples.

All plot source tables are written as tab-separated files.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable, Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

SPECIES_ORDER = [
    "Plasmodium falciparum",
    "Plasmodium vivax",
    "Plasmodium knowlesi",
]

SPECIES_SHORT = {
    "Plasmodium falciparum": "P. falciparum",
    "Plasmodium vivax": "P. vivax",
    "Plasmodium knowlesi": "P. knowlesi",
}

COMPLEXITY_ORDER = ["single", "panel2", "panel3"]
COMPLEXITY_DISPLAY = {"single": "Single", "panel2": "Panel 2", "panel3": "Panel 3"}

MAIN_METHOD_ORDER = [
    "Kraken2 single-read",
    "Kraken2 multi-read",
    "Minimap focused single-read",
    "Minimap focused multi-read",
    "Metabuli",
    "Kraken2 single-assembly",
    "Kraken2 multi-assembly",
]

METHOD_DISPLAY = {
    "Kraken2 single-read": "Kraken2\nsingle-read",
    "Kraken2 multi-read": "Kraken2\nmulti-read",
    "Minimap focused single-read": "Focused minimap\nsingle-read",
    "Minimap focused multi-read": "Focused minimap\nmulti-read",
    "Metabuli": "Metabuli",
    "Kraken2 single-assembly": "Assembly-first\nsingle",
    "Kraken2 multi-assembly": "Assembly-first\nmulti",
    "Flye": "Flye",
    "Flye+Medaka": "Flye+Medaka",
}

METHOD_COLOURS = {
    "Kraken2 single-read": "#1f77b4",
    "Kraken2 multi-read": "#003f5c",
    "Minimap focused single-read": "#7a5195",
    "Minimap focused multi-read": "#955196",
    "Metabuli": "#ef5675",
    "Kraken2 single-assembly": "#2f4b7c",
    "Kraken2 multi-assembly": "#665191",
    "Flye": "#4c956c",
    "Flye+Medaka": "#b56576",
}

METHOD_MARKERS = {
    "Kraken2 single-read": "o",
    "Kraken2 multi-read": "o",
    "Minimap focused single-read": "D",
    "Minimap focused multi-read": "D",
    "Metabuli": "^",
    "Kraken2 single-assembly": "s",
    "Kraken2 multi-assembly": "s",
    "Flye": "o",
    "Flye+Medaka": "s",
}

HEATMAP_CMAP = mpl.colormaps["viridis"]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Create publication-ready summary plots from the spike-in report "
            "workbooks and optional minimap-specific summary tables."
        )
    )
    parser.add_argument("--method_performance_xlsx", required=True)
    parser.add_argument("--real_world_xlsx", required=True)
    parser.add_argument("--replicate_report_xlsx", required=True)
    parser.add_argument("--threshold_report_xlsx", required=True)
    parser.add_argument("--minimap_tracked_tsv", required=False, default=None)
    parser.add_argument("--minimap_real_world_tsv", required=False, default=None)
    parser.add_argument("--out_dir", required=True)
    return parser.parse_args()


def apply_plot_style() -> None:
    """Apply a clean, publication-focused Matplotlib style."""
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
        }
    )


def load_sheet(*, path: Path, sheet_name: str) -> pd.DataFrame:
    """Load one worksheet from an Excel workbook.

    Parameters
    ----------
    path : Path
        Workbook path.
    sheet_name : str
        Sheet name.

    Returns
    -------
    pd.DataFrame
        Loaded worksheet.
    """
    return pd.read_excel(path, sheet_name=sheet_name)


def save_figure_set(*, fig: plt.Figure, out_path: Path) -> tuple[Path, Path, Path]:
    """Save a figure as SVG, PDF, and PNG.

    Parameters
    ----------
    fig : plt.Figure
        Figure to save.
    out_path : Path
        Base output path.

    Returns
    -------
    tuple[Path, Path, Path]
        SVG, PDF, and PNG output paths.
    """
    svg_path = out_path.with_suffix(".svg")
    pdf_path = out_path.with_suffix(".pdf")
    png_path = out_path.with_suffix(".png")
    fig.savefig(svg_path)
    fig.savefig(pdf_path)
    fig.savefig(png_path)
    return svg_path, pdf_path, png_path


def infer_complexity(*, run_name: str) -> Optional[str]:
    """Infer panel complexity from a run name.

    Parameters
    ----------
    run_name : str
        Run name.

    Returns
    -------
    Optional[str]
        One of ``single``, ``panel2``, ``panel3``, or ``None``.
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
    """Normalise a raw target label.

    Parameters
    ----------
    target_label : str
        Raw target label.

    Returns
    -------
    Optional[str]
        One of the principal species, or ``None`` if the row should be skipped.
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
        Optional classifier string.

    Returns
    -------
    str
        Display label.
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
        if workflow == "single_assembly":
            return "Kraken2 single-assembly"
        return "Kraken2 multi-assembly"

    return f"{workflow} | {metric}"


def shorten_species(*, species: str) -> str:
    """Return a short display name for a species.

    Parameters
    ----------
    species : str
        Full species name.

    Returns
    -------
    str
        Short species name.
    """
    return SPECIES_SHORT.get(str(species), str(species))


def lighten_colour(*, hex_colour: str, fraction: float) -> tuple[float, float, float]:
    """Lighten a hex colour by mixing it with white.

    Parameters
    ----------
    hex_colour : str
        Base colour.
    fraction : float
        Mixing fraction towards white.

    Returns
    -------
    tuple[float, float, float]
        RGB tuple.
    """
    rgb = np.array(mpl.colors.to_rgb(hex_colour))
    return tuple(rgb + (1.0 - rgb) * fraction)


def write_tsv(*, dataframe: pd.DataFrame, path: Path) -> None:
    """Write a table as TSV.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Table to write.
    path : Path
        Output path.
    """
    dataframe.to_csv(path, sep="\t", index=False)


def build_tradeoff_dataframe(
    *,
    real_world_df: pd.DataFrame,
    minimap_real_world_df: Optional[pd.DataFrame],
) -> pd.DataFrame:
    """Build the trade-off scatter table.

    Parameters
    ----------
    real_world_df : pd.DataFrame
        Workflow-level real-world summary.
    minimap_real_world_df : Optional[pd.DataFrame]
        Optional focused minimap real-world summary.

    Returns
    -------
    pd.DataFrame
        Trade-off table.
    """
    rows: list[dict[str, object]] = []

    for _, row in real_world_df.iterrows():
        label = workflow_display_label(
            workflow=row["workflow"],
            metric="generic",
            run_name="",
            classifier=row.get("classifier"),
        )
        if label not in MAIN_METHOD_ORDER:
            continue
        rows.append(
            {
                "label": label,
                "source": "main",
                "all_expected_sensitivity": row["all_expected_sensitivity"],
                "clean_sensitivity": row.get("clean_sensitivity", np.nan),
                "mean_off_target_taxa_positive": row["mean_off_target_taxa_positive"],
                "mean_reported_taxa_positive": row["mean_reported_taxa_positive"],
                "positive_off_target_rate": row["positive_off_target_rate"],
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
                    "source": "minimap",
                    "all_expected_sensitivity": row["all_expected_sensitivity"],
                    "clean_sensitivity": row.get("clean_sensitivity", np.nan),
                    "mean_off_target_taxa_positive": row["mean_off_target_taxa_positive"],
                    "mean_reported_taxa_positive": row["mean_reported_taxa_positive"],
                    "positive_off_target_rate": row["positive_off_target_rate"],
                }
            )

    out_df = pd.DataFrame(rows).dropna(
        subset=["all_expected_sensitivity", "mean_off_target_taxa_positive"]
    )
    out_df["label"] = pd.Categorical(
        out_df["label"], categories=MAIN_METHOD_ORDER, ordered=True
    )
    out_df = out_df.sort_values(by="label").reset_index(drop=True)
    return out_df


def make_tradeoff_plot(*, plot_df: pd.DataFrame, out_path: Path) -> None:
    """Plot sensitivity versus off-target burden.

    Parameters
    ----------
    plot_df : pd.DataFrame
        Trade-off table.
    out_path : Path
        Output path stem.
    """
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.grid(True, axis="both", alpha=0.18, linewidth=0.6)
    ax.axvline(0.5, color="#999999", linewidth=0.8, linestyle="--", alpha=0.6)
    ax.axhline(2.0, color="#999999", linewidth=0.8, linestyle="--", alpha=0.6)

    label_offsets = {
        "Kraken2 single-read": (0.010, 0.15),
        "Kraken2 multi-read": (0.010, -0.25),
        "Minimap focused single-read": (0.010, 0.10),
        "Minimap focused multi-read": (0.010, -0.20),
        "Metabuli": (0.010, 0.15),
        "Kraken2 single-assembly": (0.010, 0.10),
        "Kraken2 multi-assembly": (0.010, -0.15),
    }

    for _, row in plot_df.iterrows():
        label = str(row["label"])
        x_value = float(row["all_expected_sensitivity"])
        y_value = float(row["mean_off_target_taxa_positive"])
        size = 55 + 12 * float(row["mean_reported_taxa_positive"])
        colour = METHOD_COLOURS[label]
        marker = METHOD_MARKERS[label]
        edge_colour = "#222222"

        ax.scatter(
            x_value,
            y_value,
            s=size,
            marker=marker,
            color=colour,
            edgecolor=edge_colour,
            linewidth=0.7,
            zorder=3,
        )

        dx, dy = label_offsets.get(label, (0.010, 0.10))
        display = METHOD_DISPLAY.get(label, label)
        ax.annotate(
            display,
            xy=(x_value, y_value),
            xytext=(x_value + dx, y_value + dy),
            fontsize=8,
            ha="left",
            va="center",
            arrowprops={
                "arrowstyle": "-",
                "linewidth": 0.6,
                "color": "#666666",
                "shrinkA": 0,
                "shrinkB": 0,
            },
        )

    ax.set_xlabel("Recovered expected species set (all-expected sensitivity)")
    ax.set_ylabel("Mean off-target taxa per positive sample")
    ax.set_xlim(-0.02, 1.03)
    ax.set_ylim(0, max(8.2, plot_df["mean_off_target_taxa_positive"].max() + 0.8))
    ax.set_title("Sensitivity versus taxonomic burden", loc="left")

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=METHOD_COLOURS["Kraken2 multi-read"],
            markeredgecolor="#222222",
            label="Read-level classifier",
            markersize=7,
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            linestyle="",
            markerfacecolor=METHOD_COLOURS["Minimap focused multi-read"],
            markeredgecolor="#222222",
            label="Focused minimap",
            markersize=7,
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            markerfacecolor=METHOD_COLOURS["Kraken2 multi-assembly"],
            markeredgecolor="#222222",
            label="Assembly-first",
            markersize=7,
        ),
        Line2D(
            [0],
            [0],
            marker="^",
            linestyle="",
            markerfacecolor=METHOD_COLOURS["Metabuli"],
            markeredgecolor="#222222",
            label="Metabuli",
            markersize=7,
        ),
    ]
    ax.legend(
        handles=legend_handles,
        frameon=False,
        loc="upper right",
        title="Workflow class",
        title_fontsize=9,
    )
    fig.tight_layout()
    save_figure_set(fig=fig, out_path=out_path)
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
        Run-level performance sheet.
    minimap_tracked_df : Optional[pd.DataFrame]
        Optional minimap tracked-target table.

    Returns
    -------
    pd.DataFrame
        Long-form sensitivity table.
    """
    rows: list[dict[str, object]] = []

    for _, row in run_performance_df.iterrows():
        target = normalise_target(target_label=row["target_label"])
        if target is None:
            continue
        label = workflow_display_label(
            workflow=row["workflow"],
            metric=row["metric"],
            run_name=row["run_name"],
        )
        if label not in {
            "Kraken2 single-read",
            "Kraken2 multi-read",
            "Kraken2 single-assembly",
            "Kraken2 multi-assembly",
        }:
            continue
        complexity = infer_complexity(run_name=row["run_name"])
        if complexity is None:
            continue
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

    order = [
        "Kraken2 single-read",
        "Kraken2 multi-read",
        "Minimap focused single-read",
        "Minimap focused multi-read",
        "Kraken2 single-assembly",
        "Kraken2 multi-assembly",
    ]
    out_df["method"] = pd.Categorical(out_df["method"], categories=order, ordered=True)
    out_df["species"] = pd.Categorical(out_df["species"], categories=SPECIES_ORDER, ordered=True)
    out_df["complexity"] = pd.Categorical(
        out_df["complexity"], categories=COMPLEXITY_ORDER, ordered=True
    )
    out_df = out_df.sort_values(by=["method", "species", "complexity"]).reset_index(drop=True)
    return out_df


def make_species_complexity_plot(*, plot_df: pd.DataFrame, out_path: Path) -> None:
    """Create a grouped heatmap of species sensitivity across panel complexity.

    Parameters
    ----------
    plot_df : pd.DataFrame
        Long-form species sensitivity table.
    out_path : Path
        Output path stem.
    """
    methods = [
        method
        for method in plot_df["method"].cat.categories
        if method in set(plot_df["method"].astype(str))
    ]
    rows: list[str] = []
    row_labels: list[str] = []
    matrix_rows: list[list[float]] = []

    for method in methods:
        method_df = plot_df.loc[plot_df["method"].astype(str) == method].copy()
        for species in SPECIES_ORDER:
            species_df = method_df.loc[method_df["species"].astype(str) == species].copy()
            pivot = (
                species_df.set_index("complexity")["sensitivity"]
                .reindex(COMPLEXITY_ORDER)
                .tolist()
            )
            rows.append(f"{method} | {species}")
            row_labels.append(shorten_species(species=species))
            matrix_rows.append([np.nan if pd.isna(x) else float(x) for x in pivot])

    data = np.array(matrix_rows, dtype=float)
    fig_height = max(5.6, 0.52 * len(rows) + 1.5)
    fig, ax = plt.subplots(figsize=(6.8, fig_height))
    cmap = HEATMAP_CMAP.copy()
    cmap.set_bad("#efefef")
    im = ax.imshow(data, aspect="auto", vmin=0.0, vmax=1.0, cmap=cmap)

    ax.set_xticks(range(len(COMPLEXITY_ORDER)))
    ax.set_xticklabels([COMPLEXITY_DISPLAY[x] for x in COMPLEXITY_ORDER])
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels)

    for method_index, method in enumerate(methods):
        first_row = method_index * len(SPECIES_ORDER)
        ax.text(
            -0.72,
            first_row + 1.0,
            METHOD_DISPLAY.get(method, method),
            ha="right",
            va="center",
            fontsize=9,
            fontweight="bold",
            transform=ax.transData,
        )
        if method_index > 0:
            ax.axhline(first_row - 0.5, color="white", linewidth=2.0)

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            value = data[i, j]
            text = "NA" if math.isnan(value) else f"{value:.2f}"
            text_colour = "#111111" if math.isnan(value) or value > 0.55 else "white"
            ax.text(j, i, text, ha="center", va="center", fontsize=8, color=text_colour)

    ax.set_title("Species sensitivity across panel complexity", loc="left")
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Sensitivity")
    fig.tight_layout()
    save_figure_set(fig=fig, out_path=out_path)
    plt.close(fig)


def build_first_detection_dataframe(*, run_thresholds_df: pd.DataFrame) -> pd.DataFrame:
    """Build a clean interval table for earliest detection.

    Parameters
    ----------
    run_thresholds_df : pd.DataFrame
        Run thresholds sheet.

    Returns
    -------
    pd.DataFrame
        Selected earliest-detection table.
    """
    rows: list[dict[str, object]] = []

    for _, row in run_thresholds_df.iterrows():
        target = normalise_target(target_label=row["target_label"])
        if target is None:
            continue
        label = workflow_display_label(
            workflow=row["workflow"],
            metric=row["metric"],
            run_name=row["run_name"],
        )
        if label not in {"Kraken2 single-read", "Kraken2 multi-read"}:
            continue
        complexity = infer_complexity(run_name=row["run_name"])
        if complexity is None:
            continue
        rows.append(
            {
                "method": label,
                "complexity": complexity,
                "species": target,
                "first_any": row["first_detected_spike_any"],
                "first_50pct": row["first_detected_spike_50pct"],
                "first_95pct": row["first_detected_spike_95pct"],
                "first_all": row["first_detected_spike_all"],
            }
        )

    out_df = pd.DataFrame(rows)
    out_df = out_df.drop_duplicates().copy()

    sort_key = {
        ("Kraken2 single-read", "single", "Plasmodium vivax"): 1,
        ("Kraken2 multi-read", "panel2", "Plasmodium falciparum"): 2,
        ("Kraken2 multi-read", "panel2", "Plasmodium vivax"): 3,
        ("Kraken2 multi-read", "panel3", "Plasmodium falciparum"): 4,
        ("Kraken2 multi-read", "panel3", "Plasmodium vivax"): 5,
        ("Kraken2 multi-read", "panel3", "Plasmodium knowlesi"): 6,
    }
    out_df["sort_order"] = [
        sort_key.get((r["method"], r["complexity"], r["species"]), 999)
        for _, r in out_df.iterrows()
    ]
    out_df = out_df.sort_values(by="sort_order").reset_index(drop=True)
    return out_df.drop(columns="sort_order")


def make_first_detection_plot(*, plot_df: pd.DataFrame, out_path: Path) -> None:
    """Create a log-scale interval plot for earliest reproducible detection.

    Parameters
    ----------
    plot_df : pd.DataFrame
        Earliest-detection table.
    out_path : Path
        Output path stem.
    """
    fig_height = max(4.0, 0.55 * len(plot_df) + 1.2)
    fig, ax = plt.subplots(figsize=(7.8, fig_height))

    x_map = {
        "Any": "first_any",
        "50%": "first_50pct",
        "95%": "first_95pct",
        "All": "first_all",
    }
    stage_colours = {
        "Any": "#7b3294",
        "50%": "#4575b4",
        "95%": "#1a9850",
        "All": "#d73027",
    }
    stage_markers = {"Any": "o", "50%": "s", "95%": "^", "All": "D"}

    y_positions = np.arange(len(plot_df))[::-1]
    ax.grid(True, axis="x", which="both", alpha=0.18)

    for y_pos, (_, row) in zip(y_positions, plot_df.iterrows()):
        finite_values = [
            row[col]
            for col in x_map.values()
            if pd.notna(row[col]) and float(row[col]) > 0
        ]
        if finite_values:
            ax.hlines(
                y=y_pos,
                xmin=min(finite_values),
                xmax=max(finite_values),
                color="#bbbbbb",
                linewidth=1.2,
                zorder=1,
            )
        for stage, column in x_map.items():
            value = row[column]
            if pd.isna(value) or float(value) <= 0:
                continue
            ax.scatter(
                float(value),
                y_pos,
                s=55,
                color=stage_colours[stage],
                marker=stage_markers[stage],
                edgecolor="#222222",
                linewidth=0.5,
                zorder=3,
            )
            ax.text(
                float(value),
                y_pos + 0.13,
                f"{int(value)}",
                fontsize=8,
                ha="center",
                va="bottom",
            )

    labels = [
        f"{METHOD_DISPLAY.get(row['method'], row['method']).replace(chr(10), ' ')} | "
        f"{COMPLEXITY_DISPLAY.get(row['complexity'], row['complexity'])} | "
        f"{shorten_species(species=row['species'])}"
        for _, row in plot_df.iterrows()
    ]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels)
    ax.set_xscale("log")
    ax.set_xlabel("First detected spike")
    ax.set_title("Earliest reproducible detection", loc="left")

    legend_handles = [
        Line2D([0], [0], marker=stage_markers[stage], linestyle="", color=stage_colours[stage],
               markeredgecolor="#222222", label=stage, markersize=7)
        for stage in x_map
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="lower right", title="Detection stage")
    fig.tight_layout()
    save_figure_set(fig=fig, out_path=out_path)
    plt.close(fig)


def build_polishing_dataframe(*, polishing_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build tables for Flye versus Medaka comparison.

    Parameters
    ----------
    polishing_df : pd.DataFrame
        Polishing comparison sheet.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Species-level sensitivity table and panel-level burden table.
    """
    sens_rows: list[dict[str, object]] = []
    burden_rows: list[dict[str, object]] = []

    for _, row in polishing_df.iterrows():
        panel = str(row["panel"])
        burden_rows.extend(
            [
                {
                    "panel": panel,
                    "workflow": "Flye",
                    "mean_off_target_taxa_positive": row["unpolished_mean_off_target_taxa_positive"],
                },
                {
                    "panel": panel,
                    "workflow": "Flye+Medaka",
                    "mean_off_target_taxa_positive": row["polished_mean_off_target_taxa_positive"],
                },
            ]
        )
        for species in SPECIES_ORDER:
            safe_species = species.replace(" ", "_")
            sens_rows.extend(
                [
                    {
                        "panel": panel,
                        "species": species,
                        "workflow": "Flye",
                        "sensitivity": row.get(f"unpolished_{safe_species}_sensitivity"),
                    },
                    {
                        "panel": panel,
                        "species": species,
                        "workflow": "Flye+Medaka",
                        "sensitivity": row.get(f"polished_{safe_species}_sensitivity"),
                    },
                ]
            )

    sens_df = pd.DataFrame(sens_rows)
    burden_df = pd.DataFrame(burden_rows)
    return sens_df, burden_df


def make_polishing_plot(
    *, sens_df: pd.DataFrame, burden_df: pd.DataFrame, out_path: Path
) -> None:
    """Create a cleaner Flye versus Flye+Medaka comparison plot.

    Parameters
    ----------
    sens_df : pd.DataFrame
        Species-level sensitivity table.
    burden_df : pd.DataFrame
        Panel-level burden table.
    out_path : Path
        Output path stem.
    """
    fig = plt.figure(figsize=(8.8, 5.6))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 0.85], wspace=0.28)

    axes = {
        "panel2": fig.add_subplot(grid[0, 0]),
        "panel3": fig.add_subplot(grid[0, 1]),
        "burden": fig.add_subplot(grid[0, 2]),
    }

    workflows = ["Flye", "Flye+Medaka"]
    panel_order = ["panel2", "panel3"]

    for panel in panel_order:
        ax = axes[panel]
        subset = sens_df.loc[sens_df["panel"] == panel].copy()
        y_positions = np.arange(len(SPECIES_ORDER))[::-1]

        for idx, species in enumerate(SPECIES_ORDER):
            species_subset = subset.loc[subset["species"] == species].copy()
            values = {
                workflow: species_subset.loc[
                    species_subset["workflow"] == workflow, "sensitivity"
                ].iloc[0]
                for workflow in workflows
            }
            y_val = y_positions[idx]
            x1 = float(values["Flye"]) if pd.notna(values["Flye"]) else np.nan
            x2 = float(values["Flye+Medaka"]) if pd.notna(values["Flye+Medaka"]) else np.nan

            if pd.notna(x1) and pd.notna(x2):
                ax.plot([x1, x2], [y_val, y_val], color="#aaaaaa", linewidth=1.1, zorder=1)

            if pd.notna(x1):
                ax.scatter(
                    x1,
                    y_val,
                    s=55,
                    color=METHOD_COLOURS["Flye"],
                    edgecolor="#222222",
                    linewidth=0.6,
                    zorder=3,
                )
            if pd.notna(x2):
                ax.scatter(
                    x2,
                    y_val,
                    s=55,
                    color=METHOD_COLOURS["Flye+Medaka"],
                    marker="s",
                    edgecolor="#222222",
                    linewidth=0.6,
                    zorder=3,
                )

        ax.set_yticks(y_positions)
        if panel == "panel2":
            ax.set_yticklabels([shorten_species(species=s) for s in SPECIES_ORDER])
        else:
            ax.set_yticklabels([])
        ax.set_xlim(0, max(0.2, sens_df["sensitivity"].max() + 0.03))
        ax.grid(True, axis="x", alpha=0.18)
        ax.set_xlabel("Sensitivity")
        ax.set_title(f"{COMPLEXITY_DISPLAY[panel]}: tracked-target sensitivity", fontsize=10)

    burden_ax = axes["burden"]
    burden_ax.grid(True, axis="x", alpha=0.18)
    burden_y = np.array([1, 0])

    for idx, panel in enumerate(panel_order):
        subset = burden_df.loc[burden_df["panel"] == panel].copy()
        flye_val = float(
            subset.loc[subset["workflow"] == "Flye", "mean_off_target_taxa_positive"].iloc[0]
        )
        medaka_val = float(
            subset.loc[
                subset["workflow"] == "Flye+Medaka", "mean_off_target_taxa_positive"
            ].iloc[0]
        )
        y_val = burden_y[idx]
        burden_ax.plot([flye_val, medaka_val], [y_val, y_val], color="#aaaaaa", linewidth=1.1)
        burden_ax.scatter(
            flye_val,
            y_val,
            s=55,
            color=METHOD_COLOURS["Flye"],
            edgecolor="#222222",
            linewidth=0.6,
            zorder=3,
        )
        burden_ax.scatter(
            medaka_val,
            y_val,
            s=55,
            marker="s",
            color=METHOD_COLOURS["Flye+Medaka"],
            edgecolor="#222222",
            linewidth=0.6,
            zorder=3,
        )
        burden_ax.text(max(flye_val, medaka_val) + 0.05, y_val, COMPLEXITY_DISPLAY[panel], va="center")

    burden_ax.set_yticks([])
    burden_ax.set_xlabel("Mean off-target taxa")
    burden_ax.set_title("Panel-level taxonomic burden", fontsize=10)
    burden_ax.set_xlim(0, max(3.0, burden_df["mean_off_target_taxa_positive"].max() + 0.45))

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="", color=METHOD_COLOURS["Flye"],
               markeredgecolor="#222222", label="Flye", markersize=7),
        Line2D([0], [0], marker="s", linestyle="", color=METHOD_COLOURS["Flye+Medaka"],
               markeredgecolor="#222222", label="Flye+Medaka", markersize=7),
    ]
    fig.legend(handles=legend_handles, frameon=False, loc="upper centre", ncol=2)
    fig.suptitle("Flye versus Flye+Medaka", y=0.99, x=0.06, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure_set(fig=fig, out_path=out_path)
    plt.close(fig)


def build_threshold_dataframe(*, interpretive_df: pd.DataFrame) -> pd.DataFrame:
    """Select a representative threshold-band table.

    Parameters
    ----------
    interpretive_df : pd.DataFrame
        Interpretive bands sheet.

    Returns
    -------
    pd.DataFrame
        Selected threshold-band table.
    """
    rows: list[dict[str, object]] = []
    for _, row in interpretive_df.iterrows():
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
            rows.append(
                {
                    "label": f"{METHOD_DISPLAY[label].replace(chr(10), ' ')} | {shorten_species(species=target)}",
                    "usually_background_at_or_below": row["usually_background_at_or_below"],
                    "grey_zone_from": row["grey_zone_from"],
                    "grey_zone_to": row["grey_zone_to"],
                    "likely_true_signal_at_or_above": row["likely_true_signal_at_or_above"],
                }
            )
        elif label == "Kraken2 single-read" and target == "Plasmodium vivax":
            rows.append(
                {
                    "label": f"{METHOD_DISPLAY[label].replace(chr(10), ' ')} | {shorten_species(species=target)}",
                    "usually_background_at_or_below": row["usually_background_at_or_below"],
                    "grey_zone_from": row["grey_zone_from"],
                    "grey_zone_to": row["grey_zone_to"],
                    "likely_true_signal_at_or_above": row["likely_true_signal_at_or_above"],
                }
            )
        elif label == "Metabuli" and target in {"Plasmodium falciparum", "Plasmodium knowlesi"}:
            rows.append(
                {
                    "label": f"{METHOD_DISPLAY[label]} | {shorten_species(species=target)}",
                    "usually_background_at_or_below": row["usually_background_at_or_below"],
                    "grey_zone_from": row["grey_zone_from"],
                    "grey_zone_to": row["grey_zone_to"],
                    "likely_true_signal_at_or_above": row["likely_true_signal_at_or_above"],
                }
            )

    out_df = pd.DataFrame(rows)
    return out_df.drop_duplicates().reset_index(drop=True)


def make_threshold_plot(*, plot_df: pd.DataFrame, out_path: Path) -> None:
    """Create a segmented threshold-band plot.

    Parameters
    ----------
    plot_df : pd.DataFrame
        Threshold-band table.
    out_path : Path
        Output path stem.
    """
    max_signal = float(
        np.nanmax(
            [
                plot_df["usually_background_at_or_below"].max(),
                plot_df["grey_zone_to"].max(),
                plot_df["likely_true_signal_at_or_above"].max(),
            ]
        )
    )
    x_max = max(6.0, max_signal + 1.0)

    fig, ax = plt.subplots(figsize=(8.0, max(3.2, 0.65 * len(plot_df) + 1.2)))
    y_positions = np.arange(len(plot_df))[::-1]
    band_height = 0.58

    colours = {
        "background": "#c7dcef",
        "grey": "#d9d9d9",
        "likely": "#b8e3b1",
    }

    for y_pos, (_, row) in zip(y_positions, plot_df.iterrows()):
        bg_end = float(row["usually_background_at_or_below"])
        likely_start = float(row["likely_true_signal_at_or_above"])
        grey_from = row["grey_zone_from"]
        grey_to = row["grey_zone_to"]

        ax.add_patch(
            Rectangle(
                (0, y_pos - band_height / 2),
                bg_end,
                band_height,
                facecolor=colours["background"],
                edgecolor="white",
                linewidth=1.0,
            )
        )
        if pd.notna(grey_from) and pd.notna(grey_to) and float(grey_to) > float(grey_from):
            ax.add_patch(
                Rectangle(
                    (float(grey_from), y_pos - band_height / 2),
                    float(grey_to) - float(grey_from),
                    band_height,
                    facecolor=colours["grey"],
                    edgecolor="white",
                    linewidth=1.0,
                )
            )
        ax.add_patch(
            Rectangle(
                (likely_start, y_pos - band_height / 2),
                x_max - likely_start,
                band_height,
                facecolor=colours["likely"],
                edgecolor="white",
                linewidth=1.0,
            )
        )

        ax.text(bg_end, y_pos + 0.28, f"≤{int(bg_end)}", fontsize=8, ha="right")
        ax.text(likely_start, y_pos + 0.28, f"≥{int(likely_start)}", fontsize=8, ha="left")

    ax.set_yticks(y_positions)
    ax.set_yticklabels(plot_df["label"].tolist())
    ax.set_xlim(0, x_max)
    ax.set_xlabel("Raw signal value")
    ax.set_title("Representative threshold bands", loc="left")
    ax.grid(True, axis="x", alpha=0.18)

    legend_handles = [
        Patch(facecolor=colours["background"], edgecolor="none", label="Usually background"),
        Patch(facecolor=colours["grey"], edgecolor="none", label="Grey zone"),
        Patch(facecolor=colours["likely"], edgecolor="none", label="Likely true signal"),
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="upper right")
    fig.tight_layout()
    save_figure_set(fig=fig, out_path=out_path)
    plt.close(fig)


def write_manifest(*, rows: Iterable[dict[str, str]], out_path: Path) -> None:
    """Write a plot manifest as TSV.

    Parameters
    ----------
    rows : Iterable[dict[str, str]]
        Manifest rows.
    out_path : Path
        Output path.
    """
    pd.DataFrame(list(rows)).to_csv(out_path, sep="\t", index=False)


def main() -> None:
    """Run the publication plotting workflow."""
    apply_plot_style()
    args = parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    real_world_df = load_sheet(
        path=Path(args.real_world_xlsx),
        sheet_name="workflow_real_world_summary",
    )
    run_performance_df = load_sheet(
        path=Path(args.replicate_report_xlsx),
        sheet_name="run_performance",
    )
    run_thresholds_df = load_sheet(
        path=Path(args.replicate_report_xlsx),
        sheet_name="run_thresholds",
    )
    polishing_df = load_sheet(
        path=Path(args.real_world_xlsx),
        sheet_name="polishing_comparison",
    )
    interpretive_df = load_sheet(
        path=Path(args.threshold_report_xlsx),
        sheet_name="interpretive_bands",
    )

    minimap_tracked_df = (
        pd.read_csv(args.minimap_tracked_tsv, sep="\t")
        if args.minimap_tracked_tsv
        else None
    )
    minimap_real_world_df = (
        pd.read_csv(args.minimap_real_world_tsv, sep="\t")
        if args.minimap_real_world_tsv
        else None
    )

    manifest_rows = []

    tradeoff_df = build_tradeoff_dataframe(
        real_world_df=real_world_df,
        minimap_real_world_df=minimap_real_world_df,
    )
    tradeoff_table = out_dir / "plot01_tradeoff_table.tsv"
    write_tsv(dataframe=tradeoff_df, path=tradeoff_table)
    tradeoff_plot = out_dir / "plot01_sensitivity_vs_taxonomic_burden"
    save_paths = save_figure_set  # noqa: F841
    make_tradeoff_plot(plot_df=tradeoff_df, out_path=tradeoff_plot)
    manifest_rows.append(
        {
            "plot_id": "plot01",
            "title": "Sensitivity versus taxonomic burden",
            "figure_svg": str(tradeoff_plot.with_suffix(".svg")),
            "figure_pdf": str(tradeoff_plot.with_suffix(".pdf")),
            "figure_png": str(tradeoff_plot.with_suffix(".png")),
            "table_tsv": str(tradeoff_table),
        }
    )

    species_df = build_species_complexity_dataframe(
        run_performance_df=run_performance_df,
        minimap_tracked_df=minimap_tracked_df,
    )
    species_table = out_dir / "plot02_species_complexity_table.tsv"
    write_tsv(dataframe=species_df, path=species_table)
    species_plot = out_dir / "plot02_species_by_complexity_heatmap"
    make_species_complexity_plot(plot_df=species_df, out_path=species_plot)
    manifest_rows.append(
        {
            "plot_id": "plot02",
            "title": "Species sensitivity across panel complexity",
            "figure_svg": str(species_plot.with_suffix(".svg")),
            "figure_pdf": str(species_plot.with_suffix(".pdf")),
            "figure_png": str(species_plot.with_suffix(".png")),
            "table_tsv": str(species_table),
        }
    )

    first_df = build_first_detection_dataframe(run_thresholds_df=run_thresholds_df)
    first_table = out_dir / "plot03_first_detection_table.tsv"
    write_tsv(dataframe=first_df, path=first_table)
    first_plot = out_dir / "plot03_first_detection_heatmap"
    make_first_detection_plot(plot_df=first_df, out_path=first_plot)
    manifest_rows.append(
        {
            "plot_id": "plot03",
            "title": "Earliest reproducible detection",
            "figure_svg": str(first_plot.with_suffix(".svg")),
            "figure_pdf": str(first_plot.with_suffix(".pdf")),
            "figure_png": str(first_plot.with_suffix(".png")),
            "table_tsv": str(first_table),
        }
    )

    sens_df, burden_df = build_polishing_dataframe(polishing_df=polishing_df)
    polishing_table = out_dir / "plot04_polishing_comparison_table.tsv"
    joined_df = sens_df.merge(
        burden_df, on="panel", how="left", suffixes=("", "_panel")
    )
    write_tsv(dataframe=joined_df, path=polishing_table)
    polishing_plot = out_dir / "plot04_flye_vs_medaka"
    make_polishing_plot(sens_df=sens_df, burden_df=burden_df, out_path=polishing_plot)
    manifest_rows.append(
        {
            "plot_id": "plot04",
            "title": "Flye versus Flye+Medaka",
            "figure_svg": str(polishing_plot.with_suffix(".svg")),
            "figure_pdf": str(polishing_plot.with_suffix(".pdf")),
            "figure_png": str(polishing_plot.with_suffix(".png")),
            "table_tsv": str(polishing_table),
        }
    )

    threshold_df = build_threshold_dataframe(interpretive_df=interpretive_df)
    threshold_table = out_dir / "plot05_threshold_bands_table.tsv"
    write_tsv(dataframe=threshold_df, path=threshold_table)
    threshold_plot = out_dir / "plot05_threshold_bands"
    make_threshold_plot(plot_df=threshold_df, out_path=threshold_plot)
    manifest_rows.append(
        {
            "plot_id": "plot05",
            "title": "Representative threshold bands",
            "figure_svg": str(threshold_plot.with_suffix(".svg")),
            "figure_pdf": str(threshold_plot.with_suffix(".pdf")),
            "figure_png": str(threshold_plot.with_suffix(".png")),
            "table_tsv": str(threshold_table),
        }
    )

    write_manifest(rows=manifest_rows, out_path=out_dir / "plot_manifest.tsv")


if __name__ == "__main__":
    main()
