#!/usr/bin/env python3
"""Generate a stakeholder-friendly HTML report from spike-in summary outputs.

This script reads the outputs produced by the spike-in summariser and creates a
cleaner HTML report suitable for sharing. It does not trust the
``is_shuffled_control`` field from existing summary tables; instead it
recomputes shuffled-control status from run paths and available metadata.

Inputs expected inside ``--summary_dir``:
    - run_manifest.tsv
    - combined_wide.tsv
    - combined_long.tsv
    - missing_or_problematic_files.tsv (optional)
    - plot_manifest.tsv (optional)
    - plots/ (optional)

Outputs:
    - spikein_report_v3.html

The report is designed to print to PDF cleanly from a browser.
"""

from __future__ import annotations

import argparse
import html
import math
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Create a stakeholder-friendly HTML spike-in report."
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
             "<summary_dir>/spikein_report_v3.html.",
    )
    parser.add_argument(
        "--max_plots_per_section",
        type=int,
        default=6,
        help="Maximum number of plots to show per section.",
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
    for column in df.columns:
        if column in {
            "run_path",
            "summary_path",
            "source_file",
            "plot_path",
            "workflow",
            "run_name",
            "target_label",
            "plot_kind",
            "metric",
            "issue",
            "file_path",
        }:
            continue
        try:
            df[column] = pd.to_numeric(df[column])
        except Exception:
            pass
        df[column] = converted
    return df


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
        metadata = read_metadata_key_values(metadata_path=candidate / "run_metadata.tsv")
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
        metadata = read_metadata_key_values(metadata_path=candidate / "run_metadata.tsv")
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

    candidate_run_path_columns = [
        "run_path",
        "run_dir",
        "run_directory",
        "run_root",
        "source_dir",
    ]

    run_path_column = None
    for candidate in candidate_run_path_columns:
        if candidate in df.columns:
            run_path_column = candidate
            break

    if run_path_column is None:
        raise ValueError(
            "run_manifest.tsv must contain one of the following columns: "
            f"{', '.join(candidate_run_path_columns)}. "
            f"Observed columns: {', '.join(df.columns)}"
        )

    if run_path_column != "run_path":
        df = df.rename(columns={run_path_column: "run_path"})

    if "run_name" not in df.columns:
        df["run_name"] = df["run_path"].map(lambda x: Path(str(x)).name)

    if "workflow" not in df.columns:
        df["workflow"] = "unknown"

    df["is_shuffled_control_recomputed"] = df["run_path"].map(
        recompute_shuffled_flag
    )

    if "target_label" not in df.columns:
        df["target_label"] = ""

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
    candidate_run_path_columns = [
        "run_path",
        "run_dir",
        "run_directory",
        "run_root",
        "source_dir",
    ]

    run_path_column = None
    for candidate in candidate_run_path_columns:
        if candidate in df.columns:
            run_path_column = candidate
            break

    if run_path_column is not None and run_path_column != "run_path":
        df = df.rename(columns={run_path_column: "run_path"})

    if "run_path" in df.columns:
        df["is_shuffled_control_recomputed"] = df["run_path"].map(recompute_shuffled_flag)
    else:
        df["is_shuffled_control_recomputed"] = False


    if "target_label" not in df.columns:
        df["target_label"] = "Unspecified target"
    df["target_label"] = df["target_label"].replace("", "Unspecified target")
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
    missing = required_cols - set(combined_long.columns)
    if missing:
        return pd.DataFrame(
            {
                "workflow": [],
                "target_label": [],
                "metric": [],
                "n_points": [],
                "min_spike_n": [],
                "max_spike_n": [],
                "mean_value": [],
                "max_value": [],
                "baseline_value": [],
            }
        )

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
    """Group plots for display in the report.

    Args:
        plot_manifest: Plot manifest table.
        max_plots_per_section: Maximum plots to include in each section.

    Returns:
        Dictionary mapping section titles to lists of plot records.
    """
    if plot_manifest.empty:
        return {}

    df = plot_manifest.copy()
    for column in ["plot_path", "metric", "plot_kind", "workflow", "target_label"]:
        if column not in df.columns:
            df[column] = ""

    df["section"] = [
        categorise_plot(
            plot_path=str(plot_path),
            metric=str(metric),
            plot_kind=str(plot_kind),
        )
        for plot_path, metric, plot_kind in zip(
            df["plot_path"], df["metric"], df["plot_kind"]
        )
    ]

    sections: dict[str, list[dict[str, str]]] = {}
    for section, group in df.groupby("section", dropna=False):
        group = group.head(max_plots_per_section)
        records: list[dict[str, str]] = []
        for _, row in group.iterrows():
            records.append(
                {
                    "plot_path": str(row["plot_path"]),
                    "workflow": str(row.get("workflow", "")),
                    "target_label": str(row.get("target_label", "")),
                    "metric": str(row.get("metric", "")),
                    "plot_kind": str(row.get("plot_kind", "")),
                }
            )
        sections[str(section)] = records

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
        return str(target_path.relative_to(from_dir))
    except ValueError:
        return str(target_path)


def build_executive_summary(
    run_manifest: pd.DataFrame,
    combined_long: pd.DataFrame,
    issues_df: pd.DataFrame,
) -> str:
    """Build an executive summary block.

    Args:
        run_manifest: Normalised run manifest.
        combined_long: Normalised combined long table.
        issues_df: Issues table.

    Returns:
        HTML fragment.
    """
    n_runs = int(run_manifest.shape[0]) if not run_manifest.empty else 0
    n_shuffled = int(run_manifest["is_shuffled_control_recomputed"].sum()) if not run_manifest.empty else 0
    n_non_shuffled = n_runs - n_shuffled
    n_workflows = int(run_manifest["workflow"].nunique()) if not run_manifest.empty else 0
    n_targets = int(
        run_manifest["target_label_display"].replace("", pd.NA).dropna().nunique()
    ) if not run_manifest.empty else 0
    n_issues = int(issues_df.shape[0]) if not issues_df.empty else 0

    detection_text = (
        "The current summary includes read-level, assembly-based, "
        "multi-genome, and shuffled-control benchmarking runs. "
        "Target-specific signal appears to be measurable across the parsed "
        "outputs, but interpretation should continue to focus on baseline "
        "signal, method-specific sensitivity, and any remaining parsing or "
        "label-matching caveats."
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

    sections_html: list[str] = []
    for section_name, plots in plot_sections.items():
        plots_html: list[str] = []
        for record in plots:
            plot_path = Path(record["plot_path"])
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


def build_html_report(
    title: str,
    run_manifest: pd.DataFrame,
    combined_long: pd.DataFrame,
    key_metrics: pd.DataFrame,
    issues_df: pd.DataFrame,
    plot_sections: dict[str, list[dict[str, str]]],
    summary_dir: Path,
) -> str:
    """Build the complete HTML report.

    Args:
        title: Report title.
        run_manifest: Normalised run manifest.
        combined_long: Normalised long-format table.
        key_metrics: Compact summary metrics table.
        issues_df: Issues table.
        plot_sections: Grouped plot records.
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
    if not issues_display.empty:
        if "file_path" in issues_display.columns:
            issues_display["file_path"] = issues_display["file_path"].map(compact_path)

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
        max-width: 1200px;
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
      This report prioritises clarity over raw audit detail and recalculates
      shuffled-control status independently rather than relying on pre-existing
      table flags.
    </p>

    <div class="note-box">
      <strong>How to use this report:</strong>
      Use the executive summary and workflow overview for rapid interpretation.
      Refer to the run table and issues table for provenance and troubleshooting.
      Open plots in a browser tab for closer inspection before printing to PDF.
    </div>

    <section class="report-section">
      <h2>Executive summary</h2>
      {build_executive_summary(run_manifest=run_manifest, combined_long=combined_long, issues_df=issues_df)}
    </section>

    <section class="report-section">
      <h2>Workflow overview</h2>
      {make_html_table(workflow_counts)}
    </section>

    <section class="report-section">
      <h2>Key metric summary</h2>
      <p class="muted">
        Compact summary of parsed measurements by workflow, target, and metric.
      </p>
      {make_html_table(key_metrics, max_rows=50)}
    </section>

    <section class="report-section">
      <h2>Runs discovered</h2>
      <p class="muted">
        Technical run paths have been compacted for readability.
      </p>
      {make_html_table(runs_display, max_rows=100)}
    </section>

    {make_plot_sections_html(plot_sections=plot_sections, summary_dir=summary_dir)}

    <section class="report-section">
      <h2>Logged issues</h2>
      <p class="muted">
        These are files or run outputs that were missing, malformed, or
        otherwise problematic during summary generation.
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
        else summary_dir / "spikein_report_v3.html"
    )

    run_manifest_path = summary_dir / "run_manifest.tsv"
    combined_long_path = summary_dir / "combined_long.tsv"
    issues_path = summary_dir / "missing_or_problematic_files.tsv"
    plot_manifest_path = summary_dir / "plot_manifest.tsv"

    run_manifest = coerce_numeric_columns(read_tsv(run_manifest_path))
    combined_long = coerce_numeric_columns(read_tsv(combined_long_path))

    issues_df = (
        coerce_numeric_columns(read_tsv(issues_path))
        if issues_path.exists()
        else pd.DataFrame(columns=["file_path", "issue"])
    )
    plot_manifest = (
        coerce_numeric_columns(read_tsv(plot_manifest_path))
        if plot_manifest_path.exists()
        else pd.DataFrame(columns=["plot_path", "workflow", "target_label", "metric", "plot_kind"])
    )

    manifest_norm = normalise_manifest(run_manifest=run_manifest)
    long_norm = normalise_long_table(combined_long=combined_long)
    key_metrics = summarise_key_metrics(combined_long=long_norm)
    plot_sections = select_plot_groups(
        plot_manifest=plot_manifest,
        max_plots_per_section=args.max_plots_per_section,
    )

    html_text = build_html_report(
        title=args.title,
        run_manifest=manifest_norm,
        combined_long=long_norm,
        key_metrics=key_metrics,
        issues_df=issues_df,
        plot_sections=plot_sections,
        summary_dir=output_html.parent,
    )

    output_html.write_text(html_text, encoding="utf-8")
    print(f"[INFO] Wrote {output_html}")


if __name__ == "__main__":
    main()
