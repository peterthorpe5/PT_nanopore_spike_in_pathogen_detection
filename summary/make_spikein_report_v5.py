#!/usr/bin/env python3
"""Generate an enhanced HTML report from spike-in summary outputs.

This script keeps the stronger layout and tables from the later HTML report
versions, while adding experiment-level plots that give a rapid visual sense
of whether each spike-in behaved as expected.

The report:
- tolerates old and new column naming schemes
- recomputes shuffled-control status rather than trusting existing flags
- keeps the existing useful tables
- adds per-experiment detection plots
- optionally includes discovered Krona outputs if present
- does not fail if optional files are absent

Inputs expected inside ``--summary_dir``:
    - run_manifest.tsv
    - combined_wide.tsv (optional but preferred)
    - combined_long.tsv
    - summary_statistics.tsv (optional)
    - missing_or_problematic_files.tsv (optional)
    - plot_manifest.tsv (optional)
    - plots/ (optional)
    - krona_* directories or files (optional)

Outputs:
    - spikein_report_v5.html
    - report_plots/*.png
"""

from __future__ import annotations

import argparse
import html
import math
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Create an enhanced HTML spike-in report."
    )
    parser.add_argument(
        "--summary_dir",
        required=True,
        help="Directory containing summariser outputs.",
    )
    parser.add_argument(
        "--title",
        default="Spike-in analysis summary report",
        help="Title to display in the report.",
    )
    parser.add_argument(
        "--output_html",
        default=None,
        help="Optional output HTML path. Defaults to "
             "<summary_dir>/spikein_report_v5.html.",
    )
    parser.add_argument(
        "--max_plots_per_section",
        type=int,
        default=8,
        help="Maximum number of discovered plots to show per section.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    """Read a TSV file into a DataFrame.

    Args:
        path: Path to the TSV file.

    Returns:
        DataFrame loaded from the file.
    """
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def coerce_numeric_columns(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Convert numeric-looking columns to numeric dtype where possible.

    Args:
        dataframe: Input DataFrame.

    Returns:
        DataFrame with numeric-looking columns coerced where possible.
    """
    df = dataframe.copy()
    skip_columns = {
        "run_path",
        "run_dir",
        "run_directory",
        "run_root",
        "source_dir",
        "summary_path",
        "source_file",
        "plot_path",
        "path",
        "workflow",
        "workflow_type",
        "run_name",
        "target_label",
        "plot_kind",
        "metric",
        "metric_name",
        "issue",
        "file_path",
        "run_id",
    }

    for column in df.columns:
        if column in skip_columns:
            continue

        try:
            df[column] = pd.to_numeric(df[column])
        except Exception:
            pass

    return df


def get_first_existing_column(
    dataframe: pd.DataFrame,
    candidates: list[str],
) -> str | None:
    """Return the first matching column name from a candidate list.

    Args:
        dataframe: Input DataFrame.
        candidates: Ordered candidate column names.

    Returns:
        The first matching column name, or None if no match is found.
    """
    for candidate in candidates:
        if candidate in dataframe.columns:
            return candidate
    return None


def read_metadata_key_values(metadata_path: Path) -> dict[str, str]:
    """Read a two-column metadata TSV into a dictionary.

    Args:
        metadata_path: Path to run_metadata.tsv.

    Returns:
        Dictionary of metadata values. Empty if the file is absent or malformed.
    """
    metadata: dict[str, str] = {}
    if not metadata_path.exists():
        return metadata

    try:
        metadata_df = pd.read_csv(
            metadata_path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )
    except Exception:
        return metadata

    if "parameter" not in metadata_df.columns or "value" not in metadata_df.columns:
        return metadata

    for _, row in metadata_df.iterrows():
        metadata[str(row["parameter"])] = str(row["value"])
    return metadata


def looks_like_shuffle_path(text: str) -> bool:
    """Determine whether a path or run name indicates a shuffled control.

    Args:
        text: Run path or run name.

    Returns:
        True if the text strongly suggests a shuffled control run.
    """
    text_lower = text.lower()
    tokens = [
        "spikein_shuffle_",
        "shuffle_read",
        "shuffle_flye",
        "shuffled_control",
        "pathogen.shuffled.fasta",
        "/main_run",
    ]
    return any(token in text_lower for token in tokens)


def recompute_shuffled_flag(run_path: str) -> bool:
    """Recompute shuffled-control status for a run.

    Args:
        run_path: Path to the run directory.

    Returns:
        True if the run appears to be a shuffled control.
    """
    run_dir = Path(run_path)
    candidates = [run_dir, run_dir.parent]

    for candidate in candidates:
        if looks_like_shuffle_path(str(candidate)):
            return True
        if (candidate / "shuffle_metadata.tsv").exists():
            return True
        metadata = read_metadata_key_values(candidate / "run_metadata.tsv")
        metadata_text = " ".join(metadata.values()).lower()
        if "shuffle" in metadata_text or "shuffled" in metadata_text:
            return True

    return False


def compact_path(path_text: str) -> str:
    """Compact a long path for display.

    Args:
        path_text: Full path string.

    Returns:
        Shortened path string suitable for display.
    """
    if not path_text:
        return ""
    path = Path(path_text)
    parts = path.parts
    if len(parts) <= 4:
        return path_text
    return str(Path("…") / Path(*parts[-4:]))


def html_escape(value: object) -> str:
    """Escape a value for HTML display.

    Args:
        value: Value to escape.

    Returns:
        HTML-safe string.
    """
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return html.escape(str(value))


def format_number(value: object) -> str:
    """Format a number or return an escaped string.

    Args:
        value: Input value.

    Returns:
        Nicely formatted number or escaped string.
    """
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        if value.is_integer():
            return f"{int(value):,}"
        return f"{value:,.3f}"
    return html_escape(value)


def infer_target_label_from_run(run_path: str, fallback: str = "") -> str:
    """Infer a readable target label from run metadata or path.

    Args:
        run_path: Run directory path.
        fallback: Fallback label if inference fails.

    Returns:
        Target label string.
    """
    run_dir = Path(run_path)
    candidates = [run_dir, run_dir.parent]
    for candidate in candidates:
        metadata = read_metadata_key_values(candidate / "run_metadata.tsv")
        if metadata.get("target_label"):
            return metadata["target_label"]
    if fallback:
        return fallback
    return "Unspecified target"


def normalise_manifest(run_manifest: pd.DataFrame) -> pd.DataFrame:
    """Normalise manifest fields for display and report logic.

    Args:
        run_manifest: Raw run manifest.

    Returns:
        Normalised manifest.
    """
    df = run_manifest.copy()

    run_path_column = get_first_existing_column(
        dataframe=df,
        candidates=["run_path", "run_dir", "run_directory", "run_root", "source_dir"],
    )
    workflow_column = get_first_existing_column(
        dataframe=df,
        candidates=["workflow", "workflow_type"],
    )
    target_column = get_first_existing_column(
        dataframe=df,
        candidates=["target_label"],
    )

    if run_path_column is None:
        raise ValueError(
            "run_manifest.tsv must contain a run path column. "
            f"Observed columns: {', '.join(df.columns)}"
        )

    if run_path_column != "run_path":
        df = df.rename(columns={run_path_column: "run_path"})
    if workflow_column and workflow_column != "workflow":
        df = df.rename(columns={workflow_column: "workflow"})
    elif "workflow" not in df.columns:
        df["workflow"] = "unknown"
    if target_column and target_column != "target_label":
        df = df.rename(columns={target_column: "target_label"})
    elif "target_label" not in df.columns:
        df["target_label"] = ""

    if "run_name" not in df.columns:
        df["run_name"] = df["run_path"].map(lambda x: Path(str(x)).name)

    df["is_shuffled_control_recomputed"] = df["run_path"].map(
        recompute_shuffled_flag
    )

    df["target_label_display"] = [
        infer_target_label_from_run(
            run_path=str(run_path),
            fallback=str(target_label),
        )
        for run_path, target_label in zip(
            df["run_path"],
            df["target_label"],
        )
    ]

    df["workflow"] = df["workflow"].replace("", "unknown")
    df["run_path_compact"] = df["run_path"].map(compact_path)
    return df


def normalise_long_table(combined_long: pd.DataFrame) -> pd.DataFrame:
    """Normalise the combined long table.

    Args:
        combined_long: Raw combined long table.

    Returns:
        Normalised long-format table.
    """
    df = combined_long.copy()

    run_path_column = get_first_existing_column(
        dataframe=df,
        candidates=["run_path", "run_dir", "run_directory", "run_root", "source_dir"],
    )
    workflow_column = get_first_existing_column(
        dataframe=df,
        candidates=["workflow", "workflow_type"],
    )
    metric_column = get_first_existing_column(
        dataframe=df,
        candidates=["metric", "metric_name"],
    )
    value_column = get_first_existing_column(
        dataframe=df,
        candidates=["value", "metric_value", "observed_value", "count"],
    )
    spike_column = get_first_existing_column(
        dataframe=df,
        candidates=["spike_n", "total_spike_n", "spike_level", "n_reads_spiked"],
    )

    rename_map = {}
    if run_path_column and run_path_column != "run_path":
        rename_map[run_path_column] = "run_path"
    if workflow_column and workflow_column != "workflow":
        rename_map[workflow_column] = "workflow"
    if metric_column and metric_column != "metric":
        rename_map[metric_column] = "metric"
    if value_column and value_column != "value":
        rename_map[value_column] = "value"
    if spike_column and spike_column != "spike_n":
        rename_map[spike_column] = "spike_n"

    if rename_map:
        df = df.rename(columns=rename_map)

    if "run_path" in df.columns:
        df["is_shuffled_control_recomputed"] = df["run_path"].map(
            recompute_shuffled_flag
        )
    else:
        df["is_shuffled_control_recomputed"] = False

    if "workflow" not in df.columns:
        df["workflow"] = "unknown"
    if "target_label" not in df.columns:
        df["target_label"] = "Unspecified target"
    if "metric" not in df.columns:
        df["metric"] = "Unspecified metric"

    df["target_label"] = df["target_label"].replace("", "Unspecified target")
    df["workflow"] = df["workflow"].replace("", "unknown")
    df["metric"] = df["metric"].replace("", "Unspecified metric")
    return df



def summarise_key_metrics(combined_long: pd.DataFrame) -> pd.DataFrame:
    """Create a compact summary table of key metrics by workflow and target.

    Args:
        combined_long: Normalised combined long table.

    Returns:
        Summary DataFrame.
    """
    if combined_long.empty:
        return pd.DataFrame(
            columns=[
                "workflow",
                "target_label",
                "metric",
                "n_points",
                "min_spike_n",
                "max_spike_n",
                "mean_value",
                "max_value",
                "baseline_value",
            ]
        )

    required_cols = {"workflow", "target_label", "metric", "spike_n", "value"}
    if not required_cols.issubset(combined_long.columns):
        return pd.DataFrame()

    df = combined_long.copy()
    df["spike_n"] = pd.to_numeric(df["spike_n"], errors="coerce")
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["spike_n", "value"])

    rows: list[dict[str, object]] = []
    grouped = df.groupby(["workflow", "target_label", "metric"], dropna=False)

    for (workflow, target_label, metric), group in grouped:
        baseline_values = group.loc[group["spike_n"] == 0, "value"]
        baseline_value = baseline_values.mean() if not baseline_values.empty else math.nan
        rows.append(
            {
                "workflow": workflow,
                "target_label": target_label,
                "metric": metric,
                "n_points": int(group.shape[0]),
                "min_spike_n": int(group["spike_n"].min()),
                "max_spike_n": int(group["spike_n"].max()),
                "mean_value": float(group["value"].mean()),
                "max_value": float(group["value"].max()),
                "baseline_value": float(baseline_value) if not math.isnan(baseline_value) else math.nan,
            }
        )

    summary_df = pd.DataFrame(rows)
    if not summary_df.empty:
        summary_df = summary_df.sort_values(
            by=["workflow", "target_label", "metric"]
        ).reset_index(drop=True)
    return summary_df


def build_detection_overview(combined_long: pd.DataFrame) -> pd.DataFrame:
    """Build a compact detection overview by workflow, target, and metric.

    Args:
        combined_long: Normalised long-format table.

    Returns:
        Detection overview table.
    """
    required_cols = {"workflow", "target_label", "metric", "spike_n", "value"}
    if combined_long.empty or not required_cols.issubset(combined_long.columns):
        return pd.DataFrame()

    df = combined_long.copy()
    df["spike_n"] = pd.to_numeric(df["spike_n"], errors="coerce")
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["spike_n", "value"])

    rows: list[dict[str, object]] = []

    for (workflow, target_label, metric), group in df.groupby(
        ["workflow", "target_label", "metric"],
        dropna=False,
    ):
        baseline_values = group.loc[group["spike_n"] == 0, "value"]
        baseline_mean = baseline_values.mean() if not baseline_values.empty else 0.0

        non_zero = group.loc[group["spike_n"] > 0].copy()
        non_zero["delta_from_baseline"] = non_zero["value"] - baseline_mean

        detected = non_zero.loc[non_zero["delta_from_baseline"] > 0]
        first_detected_spike = (
            detected["spike_n"].min() if not detected.empty else math.nan
        )

        rows.append(
            {
                "workflow": workflow,
                "target_label": target_label,
                "metric": metric,
                "baseline_mean": float(baseline_mean),
                "max_observed_value": float(group["value"].max()),
                "max_delta_from_baseline": float((group["value"] - baseline_mean).max()),
                "first_spike_above_baseline": (
                    float(first_detected_spike)
                    if not math.isnan(first_detected_spike)
                    else math.nan
                ),
                "detected_above_baseline": bool(not detected.empty),
            }
        )

    overview = pd.DataFrame(rows)
    if not overview.empty:
        overview = overview.sort_values(
            by=["workflow", "target_label", "metric"]
        ).reset_index(drop=True)
    return overview


def top_signal_rows(combined_long: pd.DataFrame, n_rows: int = 20) -> pd.DataFrame:
    """Return the top rows by observed metric value.

    Args:
        combined_long: Normalised long-format table.
        n_rows: Number of rows to return.

    Returns:
        Top signal rows sorted by value.
    """
    if combined_long.empty or "value" not in combined_long.columns:
        return pd.DataFrame()

    df = combined_long.copy()
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["value"])

    columns_to_keep = [
        column for column in [
            "workflow",
            "target_label",
            "metric",
            "spike_n",
            "replicate",
            "value",
            "run_path",
        ] if column in df.columns
    ]

    if not columns_to_keep:
        return pd.DataFrame()

    out_df = (
        df[columns_to_keep]
        .sort_values(by="value", ascending=False)
        .head(n_rows)
        .reset_index(drop=True)
    )

    if "run_path" in out_df.columns:
        out_df["run_path"] = out_df["run_path"].map(compact_path)
    return out_df


def make_html_table(
    dataframe: pd.DataFrame,
    max_rows: int | None = None,
    class_name: str = "data-table",
) -> str:
    """Render a DataFrame as a simple HTML table.

    Args:
        dataframe: Table to render.
        max_rows: Optional maximum number of rows to display.
        class_name: HTML class name for the table.

    Returns:
        HTML string.
    """
    if dataframe.empty:
        return '<p class="muted">No data available.</p>'

    df = dataframe.copy()
    if max_rows is not None:
        df = df.head(max_rows)

    header_html = "".join(
        f"<th>{html_escape(column)}</th>" for column in df.columns
    )
    row_html_parts: list[str] = []

    for _, row in df.iterrows():
        row_html_parts.append(
            "<tr>" + "".join(
                f"<td>{format_number(value)}</td>" for value in row
            ) + "</tr>"
        )

    return (
        f'<table class="{class_name}">'
        f"<thead><tr>{header_html}</tr></thead>"
        f"<tbody>{''.join(row_html_parts)}</tbody>"
        f"</table>"
    )


def categorise_plot(plot_path: str, metric: str, plot_kind: str) -> str:
    """Assign a plot to a high-level report section.

    Args:
        plot_path: Plot file path.
        metric: Plot metric.
        plot_kind: Plot kind.

    Returns:
        Category label.
    """
    text = " ".join([plot_path or "", metric or "", plot_kind or ""]).lower()

    if "baseline" in text:
        return "Baseline-adjusted plots"
    if "efficiency" in text:
        return "Efficiency plots"
    if "assembly" in text or "contig" in text:
        return "Assembly-based results"
    if "minimap" in text:
        return "Read-level minimap2 results"
    if "kraken" in text:
        return "Read-level Kraken2 results"
    if "qc" in text or "completion" in text:
        return "Quality control"
    return "Other plots"


def select_plot_groups(
    plot_manifest: pd.DataFrame,
    max_plots_per_section: int,
) -> dict[str, list[dict[str, str]]]:
    """Group discovered plots for display in the report.

    Args:
        plot_manifest: Plot manifest table.
        max_plots_per_section: Maximum plots to include in each section.

    Returns:
        Dictionary mapping section titles to lists of plot records.
    """
    if plot_manifest.empty:
        return {}

    df = plot_manifest.copy()

    plot_path_column = get_first_existing_column(
        dataframe=df,
        candidates=["plot_path", "png_path", "path"],
    )
    workflow_column = get_first_existing_column(
        dataframe=df,
        candidates=["workflow", "workflow_type"],
    )
    target_column = get_first_existing_column(
        dataframe=df,
        candidates=["target_label"],
    )
    metric_column = get_first_existing_column(
        dataframe=df,
        candidates=["metric", "metric_name"],
    )
    plot_kind_column = get_first_existing_column(
        dataframe=df,
        candidates=["plot_kind", "plot_group"],
    )

    if plot_path_column is None:
        return {}

    rename_map = {}
    if plot_path_column != "plot_path":
        rename_map[plot_path_column] = "plot_path"
    if workflow_column and workflow_column != "workflow":
        rename_map[workflow_column] = "workflow"
    if metric_column and metric_column != "metric":
        rename_map[metric_column] = "metric"
    if target_column and target_column != "target_label":
        rename_map[target_column] = "target_label"
    if plot_kind_column and plot_kind_column != "plot_kind":
        rename_map[plot_kind_column] = "plot_kind"

    if rename_map:
        df = df.rename(columns=rename_map)

    for column in ["workflow", "target_label", "metric", "plot_kind"]:
        if column not in df.columns:
            df[column] = ""

    df = df[df["plot_path"].notna()].copy()
    df = df[df["plot_path"].astype(str).str.strip() != ""].copy()
    df = df[df["plot_path"].astype(str).str.lower() != "nan"].copy()
    df = df.drop_duplicates(subset=["plot_path"])

    section_map = {
        "raw": "Read-level and assembly raw plots",
        "baseline_adjusted": "Baseline-adjusted plots",
        "efficiency": "Efficiency plots",
        "qc": "Quality control",
    }

    df["section"] = [
        section_map.get(str(plot_kind), categorise_plot(
            plot_path=str(plot_path),
            metric=str(metric),
            plot_kind=str(plot_kind),
        ))
        for plot_path, metric, plot_kind in zip(
            df["plot_path"], df["metric"], df["plot_kind"]
        )
    ]

    sections: dict[str, list[dict[str, str]]] = {}
    for section, group in df.groupby("section", dropna=False):
        group = group.head(max_plots_per_section)
        sections[str(section)] = group.to_dict(orient="records")

    return sections



def relative_href(target_path: Path, from_dir: Path) -> str:
    """Create a relative href path for HTML links.

    Args:
        target_path: Destination file path.
        from_dir: Directory containing the HTML report.

    Returns:
        Relative path string for HTML.
    """
    try:
        return os.path.relpath(
            path=str(target_path.resolve()),
            start=str(from_dir.resolve()),
        )
    except Exception:
        return str(target_path)


def find_krona_outputs(summary_dir: Path) -> pd.DataFrame:
    """Discover optional Krona outputs in or near the summary directory.

    Args:
        summary_dir: Summary directory to inspect.

    Returns:
        DataFrame describing discovered Krona outputs.
    """
    search_roots = [summary_dir, summary_dir.parent]
    records: list[dict[str, str]] = []

    for root in search_roots:
        if not root.exists():
            continue

        patterns = [
            "krona_*",
            "**/krona_*",
            "**/*krona*.html",
            "**/*krona*.txt",
        ]
        seen: set[str] = set()

        for pattern in patterns:
            for path in root.glob(pattern):
                resolved = str(path.resolve())
                if resolved in seen:
                    continue
                seen.add(resolved)

                records.append(
                    {
                        "name": path.name,
                        "path": resolved,
                        "path_compact": compact_path(resolved),
                        "type": "directory" if path.is_dir() else path.suffix.lstrip(".") or "file",
                    }
                )

    if not records:
        return pd.DataFrame(columns=["name", "path", "path_compact", "type"])

    return pd.DataFrame(records).drop_duplicates(subset=["path"]).reset_index(drop=True)


def create_experiment_plots(
    combined_long: pd.DataFrame,
    output_dir: Path,
) -> pd.DataFrame:
    """Create per-experiment plots that summarise detection performance.

    Args:
        combined_long: Normalised combined long table.
        output_dir: Output directory for generated plots.

    Returns:
        Manifest DataFrame for generated report-specific plots.
    """
    if combined_long.empty:
        return pd.DataFrame(
            columns=["plot_path", "workflow", "target_label", "metric", "plot_kind"]
        )

    df = combined_long.copy()

    workflow_column = get_first_existing_column(
        dataframe=df,
        candidates=["workflow", "workflow_type"],
    )
    target_column = get_first_existing_column(
        dataframe=df,
        candidates=["target_label"],
    )
    metric_column = get_first_existing_column(
        dataframe=df,
        candidates=["metric", "metric_name"],
    )
    spike_column = get_first_existing_column(
        dataframe=df,
        candidates=["spike_n", "total_spike_n", "spike_level", "n_reads_spiked"],
    )
    value_column = get_first_existing_column(
        dataframe=df,
        candidates=["value", "metric_value", "observed_value", "count"],
    )
    replicate_column = get_first_existing_column(
        dataframe=df,
        candidates=["replicate", "rep"],
    )

    rename_map = {}
    if workflow_column and workflow_column != "workflow":
        rename_map[workflow_column] = "workflow"
    if target_column and target_column != "target_label":
        rename_map[target_column] = "target_label"
    if metric_column and metric_column != "metric":
        rename_map[metric_column] = "metric"
    if spike_column and spike_column != "spike_n":
        rename_map[spike_column] = "spike_n"
    if value_column and value_column != "value":
        rename_map[value_column] = "value"
    if replicate_column and replicate_column != "replicate":
        rename_map[replicate_column] = "replicate"

    if rename_map:
        df = df.rename(columns=rename_map)

    required_cols = {"workflow", "target_label", "metric", "spike_n", "value"}
    if not required_cols.issubset(df.columns):
        print(
            "[WARN] create_experiment_plots(): missing required columns. "
            f"Observed columns: {', '.join(df.columns)}"
        )
        return pd.DataFrame(
            columns=["plot_path", "workflow", "target_label", "metric", "plot_kind"]
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    df["spike_n"] = pd.to_numeric(df["spike_n"], errors="coerce")
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    if "replicate" in df.columns:
        df["replicate"] = pd.to_numeric(df["replicate"], errors="coerce")
    else:
        df["replicate"] = math.nan

    df = df.dropna(subset=["spike_n", "value"])

    manifests: list[dict[str, str]] = []

    grouped = df.groupby(["workflow", "target_label", "metric"], dropna=False)
    for (workflow, target_label, metric), group in grouped:
        summary_df = (
            group.groupby("spike_n", dropna=False)["value"]
            .agg(["mean", "min", "max", "count"])
            .reset_index()
            .sort_values("spike_n")
        )

        if summary_df.empty:
            continue

        baseline_rows = summary_df.loc[summary_df["spike_n"] == 0, "mean"]
        baseline_value = (
            float(baseline_rows.iloc[0]) if not baseline_rows.empty else 0.0
        )
        summary_df["delta_from_baseline"] = summary_df["mean"] - baseline_value
        summary_df["detected_above_baseline"] = (
            summary_df["delta_from_baseline"] > 0
        )

        safe_stub = re.sub(
            pattern=r"[^A-Za-z0-9._-]+",
            repl="_",
            string=f"{workflow}__{target_label}__{metric}",
        ).strip("_")

        fig, ax = plt.subplots(figsize=(7.5, 4.5))
        ax.plot(summary_df["spike_n"], summary_df["mean"], marker="o")
        ax.fill_between(
            summary_df["spike_n"],
            summary_df["min"],
            summary_df["max"],
            alpha=0.2,
        )
        ax.axhline(baseline_value, linestyle="--")
        ax.set_xlabel("Spiked reads")
        ax.set_ylabel("Observed signal")
        ax.set_title(f"{workflow} | {target_label} | {metric}")
        observed_path = output_dir / f"{safe_stub}.observed_vs_spike.png"
        fig.tight_layout()
        fig.savefig(observed_path, dpi=180)
        plt.close(fig)

        manifests.append(
            {
                "plot_path": str(observed_path.resolve()),
                "workflow": str(workflow),
                "target_label": str(target_label),
                "metric": str(metric),
                "plot_kind": "report_observed_vs_spike",
            }
        )

        fig, ax = plt.subplots(figsize=(7.5, 4.5))
        ax.plot(summary_df["spike_n"], summary_df["delta_from_baseline"], marker="o")
        ax.axhline(0, linestyle="--")
        ax.set_xlabel("Spiked reads")
        ax.set_ylabel("Delta from baseline")
        ax.set_title(f"{workflow} | {target_label} | {metric}")
        delta_path = output_dir / f"{safe_stub}.delta_from_baseline.png"
        fig.tight_layout()
        fig.savefig(delta_path, dpi=180)
        plt.close(fig)

        manifests.append(
            {
                "plot_path": str(delta_path.resolve()),
                "workflow": str(workflow),
                "target_label": str(target_label),
                "metric": str(metric),
                "plot_kind": "report_delta_from_baseline",
            }
        )

        fig, ax = plt.subplots(figsize=(7.5, 2.8))
        y_vals = [
            1 if value else 0
            for value in summary_df["detected_above_baseline"]
        ]
        ax.plot(summary_df["spike_n"], y_vals, marker="o")
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["No", "Yes"])
        ax.set_xlabel("Spiked reads")
        ax.set_ylabel("Detected")
        ax.set_title(f"{workflow} | {target_label} | {metric}")
        detection_path = output_dir / f"{safe_stub}.detected_boolean.png"
        fig.tight_layout()
        fig.savefig(detection_path, dpi=180)
        plt.close(fig)

        manifests.append(
            {
                "plot_path": str(detection_path.resolve()),
                "workflow": str(workflow),
                "target_label": str(target_label),
                "metric": str(metric),
                "plot_kind": "report_detected_boolean",
            }
        )

    if not manifests:
        return pd.DataFrame(
            columns=["plot_path", "workflow", "target_label", "metric", "plot_kind"]
        )

    return pd.DataFrame(manifests)


def build_executive_summary(
    run_manifest: pd.DataFrame,
    detection_overview: pd.DataFrame,
    issues_df: pd.DataFrame,
) -> str:
    """Build an executive summary block.

    Args:
        run_manifest: Normalised run manifest.
        detection_overview: Detection overview table.
        issues_df: Issues table.

    Returns:
        HTML fragment.
    """
    n_runs = int(run_manifest.shape[0]) if not run_manifest.empty else 0
    n_shuffled = (
        int(run_manifest["is_shuffled_control_recomputed"].sum())
        if not run_manifest.empty
        else 0
    )
    n_non_shuffled = n_runs - n_shuffled
    n_workflows = (
        int(run_manifest["workflow"].replace("", pd.NA).dropna().nunique())
        if not run_manifest.empty
        else 0
    )
    n_targets = (
        int(run_manifest["target_label_display"].replace("", pd.NA).dropna().nunique())
        if not run_manifest.empty
        else 0
    )
    n_issues = int(issues_df.shape[0]) if not issues_df.empty else 0
    n_detected = (
        int(detection_overview["detected_above_baseline"].fillna(False).sum())
        if not detection_overview.empty and "detected_above_baseline" in detection_overview.columns
        else 0
    )

    detection_text = (
        "The report summarises read-level, assembly-based, multi-genome, and "
        "shuffled-control benchmarking runs. The detection overview below is "
        "intended to provide an immediate view of whether observed signal rose "
        "above the baseline present at zero spiked reads."
    )

    return f"""
    <div class="summary-grid">
      <div class="summary-card">
        <div class="summary-label">Runs discovered</div>
        <div class="summary-value">{n_runs}</div>
      </div>
      <div class="summary-card">
        <div class="summary-label">Non-shuffled runs</div>
        <div class="summary-value">{n_non_shuffled}</div>
      </div>
      <div class="summary-card">
        <div class="summary-label">Shuffled controls</div>
        <div class="summary-value">{n_shuffled}</div>
      </div>
      <div class="summary-card">
        <div class="summary-label">Workflow categories</div>
        <div class="summary-value">{n_workflows}</div>
      </div>
      <div class="summary-card">
        <div class="summary-label">Distinct targets</div>
        <div class="summary-value">{n_targets}</div>
      </div>
      <div class="summary-card">
        <div class="summary-label">Detected metric series</div>
        <div class="summary-value">{n_detected}</div>
      </div>
      <div class="summary-card">
        <div class="summary-label">Logged issues</div>
        <div class="summary-value">{n_issues}</div>
      </div>
    </div>
    <p>{html_escape(detection_text)}</p>
    """


def make_plot_sections_html(
    plot_sections: dict[str, list[dict[str, str]]],
    summary_dir: Path,
) -> str:
    """Render grouped plots into HTML.

    Args:
        plot_sections: Plot groups to render.
        summary_dir: Summary directory containing the report.

    Returns:
        HTML fragment for plot sections.
    """
    if not plot_sections:
        return '<p class="muted">No plots were found for display.</p>'

    ordered_sections = [
        "Experiment-level detection plots",
        "Read-level Kraken2 results",
        "Read-level minimap2 results",
        "Assembly-based results",
        "Baseline-adjusted plots",
        "Efficiency plots",
        "Quality control",
        "Other plots",
    ]

    sections_html: list[str] = []
    for section_name in ordered_sections:
        if section_name not in plot_sections:
            continue

        plots = plot_sections[section_name]
        plots_html: list[str] = []
        for record in plots:
            raw_plot_path = record.get("plot_path", "")

            if raw_plot_path is None:
                continue
            if isinstance(raw_plot_path, float) and math.isnan(raw_plot_path):
                continue

            raw_plot_path = str(raw_plot_path).strip()
            if raw_plot_path == "" or raw_plot_path.lower() == "nan":
                continue

            plot_path = Path(raw_plot_path)
            href = relative_href(target_path=plot_path, from_dir=summary_dir)

            caption_bits = [
                bit for bit in [
                    record.get("workflow", ""),
                    record.get("target_label", ""),
                    record.get("metric", ""),
                    record.get("plot_kind", ""),
                ] if bit
            ]
            caption = " | ".join(caption_bits)

            plots_html.append(
                f"""
                <figure class="plot-card">
                <a href="{html_escape(href)}" target="_blank">
                    <img src="{html_escape(href)}" alt="{html_escape(caption)}">
                </a>
                <figcaption>{html_escape(caption)}</figcaption>
                </figure>
                """
            )

        sections_html.append(
            f"""
            <section class="report-section">
              <h2>{html_escape(section_name)}</h2>
              <div class="plot-grid">
                {''.join(plots_html)}
              </div>
            </section>
            """
        )

    for section_name, plots in plot_sections.items():
        if section_name in ordered_sections:
            continue

        plots_html: list[str] = []
        for record in plots:
            raw_plot_path = record.get("plot_path", "")

            if raw_plot_path is None:
                continue
            if isinstance(raw_plot_path, float) and math.isnan(raw_plot_path):
                continue

            raw_plot_path = str(raw_plot_path).strip()
            if raw_plot_path == "" or raw_plot_path.lower() == "nan":
                continue

            plot_path = Path(raw_plot_path)
            href = relative_href(target_path=plot_path, from_dir=summary_dir)

            caption_bits = [
                bit for bit in [
                    record.get("workflow", ""),
                    record.get("target_label", ""),
                    record.get("metric", ""),
                    record.get("plot_kind", ""),
                ] if bit
            ]
            caption = " | ".join(caption_bits)

            plots_html.append(
                f"""
                <figure class="plot-card">
                <a href="{html_escape(href)}" target="_blank">
                    <img src="{html_escape(href)}" alt="{html_escape(caption)}">
                </a>
                <figcaption>{html_escape(caption)}</figcaption>
                </figure>
                """
            )

        sections_html.append(
            f"""
            <section class="report-section">
              <h2>{html_escape(section_name)}</h2>
              <div class="plot-grid">
                {''.join(plots_html)}
              </div>
            </section>
            """
        )

    return "".join(sections_html)


def build_krona_section_html(krona_df: pd.DataFrame, summary_dir: Path) -> str:
    """Build the optional Krona section.

    Args:
        krona_df: Discovered Krona outputs.
        summary_dir: Summary directory.

    Returns:
        HTML fragment.
    """
    if krona_df.empty:
        return ""

    rows_html: list[str] = []
    for _, row in krona_df.iterrows():
        path = Path(str(row["path"]))
        href = relative_href(target_path=path, from_dir=summary_dir)
        link_html = (
            f'<a href="{html_escape(href)}" target="_blank">{html_escape(row["name"])}</a>'
        )
        rows_html.append(
            "<tr>"
            f"<td>{link_html}</td>"
            f"<td>{html_escape(row['type'])}</td>"
            f"<td>{html_escape(row['path_compact'])}</td>"
            "</tr>"
        )

    return f"""
    <section class="report-section">
      <h2>Krona outputs</h2>
      <p class="muted">
        Optional Krona outputs were discovered and are linked below where present.
      </p>
      <table class="data-table">
        <thead>
          <tr>
            <th>Name</th>
            <th>Type</th>
            <th>Path</th>
          </tr>
        </thead>
        <tbody>
          {''.join(rows_html)}
        </tbody>
      </table>
    </section>
    """


def build_html_report(
    title: str,
    run_manifest: pd.DataFrame,
    combined_wide: pd.DataFrame,
    combined_long: pd.DataFrame,
    summary_stats_df: pd.DataFrame,
    detection_overview: pd.DataFrame,
    top_signals_df: pd.DataFrame,
    issues_df: pd.DataFrame,
    plot_sections: dict[str, list[dict[str, str]]],
    krona_df: pd.DataFrame,
    summary_dir: Path,
) -> str:
    """Build the complete HTML report.

    Args:
        title: Report title.
        run_manifest: Normalised run manifest.
        combined_wide: Combined wide table.
        combined_long: Normalised long-format table.
        summary_stats_df: Summary statistics table.
        detection_overview: Detection overview table.
        top_signals_df: Top signal table.
        issues_df: Issues table.
        plot_sections: Grouped plot records.
        krona_df: Discovered Krona outputs.
        summary_dir: Report output directory.

    Returns:
        Complete HTML document.
    """
    workflow_counts = (
        run_manifest.groupby(["workflow", "is_shuffled_control_recomputed"])
        .size()
        .reset_index(name="n_runs")
        .rename(
            columns={
                "workflow": "Workflow",
                "is_shuffled_control_recomputed": "Shuffled control",
                "n_runs": "Runs",
            }
        )
    ) if not run_manifest.empty else pd.DataFrame()

    runs_display = run_manifest.copy()
    if not runs_display.empty:
        columns_to_keep = [
            column for column in [
                "run_name",
                "workflow",
                "target_label_display",
                "is_shuffled_control_recomputed",
                "run_path_compact",
            ] if column in runs_display.columns
        ]
        runs_display = runs_display[columns_to_keep].rename(
            columns={
                "run_name": "Run",
                "workflow": "Workflow",
                "target_label_display": "Target label",
                "is_shuffled_control_recomputed": "Shuffled control",
                "run_path_compact": "Run path",
            }
        )

    issues_display = issues_df.copy()
    if not issues_display.empty and "file_path" in issues_display.columns:
        issues_display["file_path"] = issues_display["file_path"].map(compact_path)

    combined_long_preview = combined_long.copy()
    if "run_path" in combined_long_preview.columns:
        combined_long_preview["run_path"] = combined_long_preview["run_path"].map(compact_path)

    combined_wide_preview = combined_wide.copy()
    for candidate in ["run_path", "run_dir", "run_directory", "source_dir"]:
        if candidate in combined_wide_preview.columns:
            combined_wide_preview[candidate] = combined_wide_preview[candidate].map(compact_path)

    style = """
    <style>
      body {
        font-family: Arial, Helvetica, sans-serif;
        margin: 0;
        padding: 0;
        color: #1a1a1a;
        background: #f7f9fc;
      }
      .container {
        max-width: 1280px;
        margin: 0 auto;
        padding: 32px;
        background: #ffffff;
      }
      h1, h2, h3 {
        color: #1f4e79;
      }
      h1 {
        margin-top: 0;
        font-size: 30px;
      }
      .lede {
        font-size: 16px;
        line-height: 1.5;
      }
      .muted {
        color: #666666;
      }
      .summary-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
        gap: 14px;
        margin: 24px 0 24px 0;
      }
      .summary-card {
        border: 1px solid #d6dfeb;
        border-radius: 10px;
        padding: 16px;
        background: #f8fbff;
      }
      .summary-label {
        font-size: 13px;
        color: #4a5d73;
        margin-bottom: 8px;
      }
      .summary-value {
        font-size: 28px;
        font-weight: bold;
        color: #1f4e79;
      }
      .report-section {
        margin-top: 34px;
      }
      .data-table {
        width: 100%;
        border-collapse: collapse;
        margin-top: 14px;
        font-size: 14px;
      }
      .data-table thead th {
        background: #1f4e79;
        color: white;
        text-align: left;
        padding: 10px 8px;
        position: sticky;
        top: 0;
      }
      .data-table td {
        border-bottom: 1px solid #e6edf5;
        padding: 8px;
        vertical-align: top;
      }
      .data-table tbody tr:nth-child(even) {
        background: #fbfdff;
      }
      .plot-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
        gap: 18px;
      }
      .plot-card {
        border: 1px solid #d6dfeb;
        border-radius: 10px;
        padding: 12px;
        margin: 0;
        background: #ffffff;
      }
      .plot-card img {
        width: 100%;
        height: auto;
        display: block;
        border-radius: 6px;
      }
      .plot-card figcaption {
        margin-top: 8px;
        font-size: 13px;
        color: #4a5d73;
      }
      .note-box {
        border-left: 4px solid #1f4e79;
        background: #f4f8fc;
        padding: 12px 16px;
        margin: 18px 0;
      }
      @media print {
        body {
          background: white;
        }
        .container {
          max-width: none;
          margin: 0;
          padding: 16px;
        }
        a {
          color: inherit;
          text-decoration: none;
        }
      }
    </style>
    """

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html_escape(title)}</title>
  {style}
</head>
<body>
  <div class="container">
    <h1>{html_escape(title)}</h1>
    <p class="lede">
      Stakeholder-facing summary of in silico ONT spike-in benchmarking runs.
      This report keeps the stronger table-based overview while adding
      experiment-level plots designed to show quickly whether signal increased
      above the zero-spike baseline.
    </p>

    <div class="note-box">
      <strong>How to use this report:</strong>
      Start with the executive summary and detection overview. Then use the
      experiment-level plots to judge whether observed signal tracks the number
      of spiked reads. Detailed tables and discovered outputs are provided
      further below for audit and interpretation.
    </div>

    <section class="report-section">
      <h2>Executive summary</h2>
      {build_executive_summary(run_manifest=run_manifest, detection_overview=detection_overview, issues_df=issues_df)}
    </section>

    <section class="report-section">
      <h2>Workflow overview</h2>
      {make_html_table(workflow_counts)}
    </section>

    <section class="report-section">
      <h2>Detection overview</h2>
      <p class="muted">
        For each workflow, target, and metric, this table shows the baseline at
        zero spiked reads, the strongest observed signal, the maximum rise above
        baseline, and whether any non-zero spike level exceeded the baseline.
      </p>
      {make_html_table(detection_overview, max_rows=100)}
    </section>

    <section class="report-section">
      <h2>Summary statistics</h2>
      {make_html_table(summary_stats_df, max_rows=100)}
    </section>

    <section class="report-section">
      <h2>Runs discovered</h2>
      <p class="muted">
        Technical run paths have been compacted for readability.
      </p>
      {make_html_table(runs_display, max_rows=100)}
    </section>

    {make_plot_sections_html(plot_sections=plot_sections, summary_dir=summary_dir)}

    {build_krona_section_html(krona_df=krona_df, summary_dir=summary_dir)}

    <section class="report-section">
      <h2>Top observed signals</h2>
      {make_html_table(top_signals_df, max_rows=30)}
    </section>

    <section class="report-section">
      <h2>Combined long table preview</h2>
      {make_html_table(combined_long_preview, max_rows=40)}
    </section>

    <section class="report-section">
      <h2>Combined wide table preview</h2>
      {make_html_table(combined_wide_preview, max_rows=25)}
    </section>

    <section class="report-section">
      <h2>Logged issues</h2>
      <p class="muted">
        These are files or run outputs that were missing, malformed, or still
        incomplete at the time the summary was generated.
      </p>
      {make_html_table(issues_display, max_rows=100)}
    </section>
  </div>
</body>
</html>
"""


def main() -> None:
    """Run report generation."""
    args = parse_args()

    summary_dir = Path(args.summary_dir).resolve()
    output_html = (
        Path(args.output_html).resolve()
        if args.output_html
        else summary_dir / "spikein_report_v5.html"
    )

    run_manifest_path = summary_dir / "run_manifest.tsv"
    combined_wide_path = summary_dir / "combined_wide.tsv"
    combined_long_path = summary_dir / "combined_long.tsv"
    summary_stats_path = summary_dir / "summary_statistics.tsv"
    issues_path = summary_dir / "missing_or_problematic_files.tsv"
    plot_manifest_path = summary_dir / "plot_manifest.tsv"

    run_manifest = coerce_numeric_columns(read_tsv(run_manifest_path))
    combined_long = coerce_numeric_columns(read_tsv(combined_long_path))
    combined_wide = (
        coerce_numeric_columns(read_tsv(combined_wide_path))
        if combined_wide_path.exists()
        else pd.DataFrame()
    )
    summary_stats_df = (
        coerce_numeric_columns(read_tsv(summary_stats_path))
        if summary_stats_path.exists()
        else pd.DataFrame()
    )
    issues_df = (
        coerce_numeric_columns(read_tsv(issues_path))
        if issues_path.exists()
        else pd.DataFrame(columns=["file_path", "issue"])
    )
    discovered_plot_manifest = (
        coerce_numeric_columns(read_tsv(plot_manifest_path))
        if plot_manifest_path.exists()
        else pd.DataFrame(columns=["plot_path", "workflow", "target_label", "metric", "plot_kind"])
    )

    manifest_norm = normalise_manifest(run_manifest=run_manifest)
    long_norm = normalise_long_table(combined_long=combined_long)

    if summary_stats_df.empty:
        summary_stats_df = summarise_key_metrics(combined_long=long_norm)

    detection_overview = build_detection_overview(combined_long=long_norm)
    top_signals_df = top_signal_rows(combined_long=long_norm, n_rows=30)

    report_plot_dir = summary_dir / "report_plots"
    generated_plot_manifest = create_experiment_plots(
        combined_long=long_norm,
        output_dir=report_plot_dir,
    )

    if discovered_plot_manifest.empty:
        plot_manifest = generated_plot_manifest.copy()
    elif generated_plot_manifest.empty:
        plot_manifest = discovered_plot_manifest.copy()
    else:
        plot_manifest = pd.concat(
            [generated_plot_manifest, discovered_plot_manifest],
            ignore_index=True,
            sort=False,
        ).drop_duplicates(subset=["plot_path"])

    # Force report-specific plots into a dedicated section.
    if not plot_manifest.empty and "plot_kind" in plot_manifest.columns:
        mask = plot_manifest["plot_kind"].astype(str).str.startswith("report_")
        if mask.any():
            report_plots = plot_manifest.loc[mask].copy()
            other_plots = plot_manifest.loc[~mask].copy()
            report_plots["section"] = "Experiment-level detection plots"
            plot_sections = select_plot_groups(
                plot_manifest=other_plots,
                max_plots_per_section=args.max_plots_per_section,
            )
            plot_sections["Experiment-level detection plots"] = (
                report_plots.head(args.max_plots_per_section).to_dict(orient="records")
            )
        else:
            plot_sections = select_plot_groups(
                plot_manifest=plot_manifest,
                max_plots_per_section=args.max_plots_per_section,
            )
    else:
        plot_sections = select_plot_groups(
            plot_manifest=plot_manifest,
            max_plots_per_section=args.max_plots_per_section,
        )

    krona_df = find_krona_outputs(summary_dir=summary_dir)

    html_text = build_html_report(
        title=args.title,
        run_manifest=manifest_norm,
        combined_wide=combined_wide,
        combined_long=long_norm,
        summary_stats_df=summary_stats_df,
        detection_overview=detection_overview,
        top_signals_df=top_signals_df,
        issues_df=issues_df,
        plot_sections=plot_sections,
        krona_df=krona_df,
        summary_dir=output_html.parent,
    )

    output_html.write_text(html_text, encoding="utf-8")
    print(f"[INFO] Wrote {output_html}")


if __name__ == "__main__":
    main()
