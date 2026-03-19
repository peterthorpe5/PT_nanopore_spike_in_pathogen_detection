#!/usr/bin/env python3
"""Generate an HTML report for ONT spike-in summary outputs.

This script combines the stronger layout from the newer stakeholder-focused
report with the richer tables and embedded plots from the earlier report. It is
designed to tolerate modest schema differences across summary files and to
degrade gracefully when optional inputs such as Krona outputs are absent.

Expected inputs inside ``--summary_dir``:
    - run_manifest.tsv
    - combined_wide.tsv (optional)
    - combined_long.tsv (optional)
    - summary_statistics.tsv (optional)
    - missing_or_problematic_files.tsv (optional)
    - plot_manifest.tsv (optional)
    - plots/ (optional)
    - krona_* directories or files (optional)

Output:
    - spikein_report_v4.html

The report is intended to open in a browser and print cleanly to PDF.
"""

from __future__ import annotations

import argparse
import html
import math
import os
from pathlib import Path
from typing import Any

import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Create a rich HTML report from spike-in summary outputs."
    )
    parser.add_argument(
        "--summary_dir",
        required=True,
        help="Directory containing summariser outputs.",
    )
    parser.add_argument(
        "--title",
        default="ONT spike-in summary report",
        help="Report title.",
    )
    parser.add_argument(
        "--output_html",
        default=None,
        help="Optional output HTML path. Defaults to "
             "<summary_dir>/spikein_report_v4.html.",
    )
    parser.add_argument(
        "--max_rows",
        type=int,
        default=100,
        help="Maximum number of rows to display in large tables.",
    )
    parser.add_argument(
        "--max_plots_per_section",
        type=int,
        default=12,
        help="Maximum number of plots to display per plot section.",
    )
    return parser.parse_args()


def read_tsv_if_exists(path: Path) -> pd.DataFrame:
    """Read a TSV file if present.

    Args:
        path: Path to the TSV file.

    Returns:
        DataFrame loaded from the file, or an empty DataFrame if absent.
    """
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def get_first_existing_column(
    dataframe: pd.DataFrame,
    candidates: list[str],
) -> str | None:
    """Return the first matching column name from a candidate list.

    Args:
        dataframe: Input DataFrame.
        candidates: Ordered candidate column names.

    Returns:
        First matching column name, or None if no match is found.
    """
    for candidate in candidates:
        if candidate in dataframe.columns:
            return candidate
    return None


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
    }

    for column in df.columns:
        if column in skip_columns:
            continue
        try:
            df[column] = pd.to_numeric(df[column])
        except Exception:
            pass

    return df


def read_metadata_key_values(metadata_path: Path) -> dict[str, str]:
    """Read a two-column metadata TSV into a dictionary.

    Args:
        metadata_path: Path to a metadata TSV.

    Returns:
        Metadata dictionary. Returns an empty dictionary if absent or malformed.
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
        True if the text strongly suggests a shuffled control.
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


def html_escape(value: Any) -> str:
    """Escape a value for HTML display.

    Args:
        value: Value to escape.

    Returns:
        HTML-safe string.
    """
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return html.escape(str(value))


def format_number(value: Any) -> str:
    """Format a number for HTML display.

    Args:
        value: Input value.

    Returns:
        Formatted number string or escaped text.
    """
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass

    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        if value.is_integer():
            return f"{int(value):,}"
        return f"{value:,.3f}"
    return html_escape(value)


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


def infer_target_label_from_run(run_path: str, fallback: str = "") -> str:
    """Infer a readable target label from run metadata or fallback text.

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

    run_path_column = get_first_existing_column(
        dataframe=df,
        candidates=["run_path", "run_dir", "run_directory", "run_root", "source_dir"],
    )
    workflow_column = get_first_existing_column(
        dataframe=df,
        candidates=["workflow", "workflow_type"],
    )
    target_label_column = get_first_existing_column(
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
    if "workflow" not in df.columns:
        df["workflow"] = "unknown"

    if target_label_column and target_label_column != "target_label":
        df = df.rename(columns={target_label_column: "target_label"})
    if "target_label" not in df.columns:
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
    df["run_path_compact"] = df["run_path"].map(compact_path)
    return df


def normalise_long_table(combined_long: pd.DataFrame) -> pd.DataFrame:
    """Normalise the combined long-format table.

    Args:
        combined_long: Raw long-format summary table.

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

    if run_path_column and run_path_column != "run_path":
        df = df.rename(columns={run_path_column: "run_path"})
    if workflow_column and workflow_column != "workflow":
        df = df.rename(columns={workflow_column: "workflow"})
    if metric_column and metric_column != "metric":
        df = df.rename(columns={metric_column: "metric"})

    if "workflow" not in df.columns:
        df["workflow"] = "unknown"

    if "run_path" in df.columns:
        df["is_shuffled_control_recomputed"] = df["run_path"].map(
            recompute_shuffled_flag
        )
    else:
        df["is_shuffled_control_recomputed"] = False

    if "target_label" not in df.columns:
        df["target_label"] = "Unspecified target"
    df["target_label"] = df["target_label"].replace("", "Unspecified target")
    return df


def summarise_key_metrics(
    combined_long: pd.DataFrame,
    summary_statistics: pd.DataFrame,
) -> pd.DataFrame:
    """Create a compact summary table of key metrics.

    This prefers the dedicated summary statistics table when present, but can
    fall back to recomputing from the long-format table.

    Args:
        combined_long: Normalised long-format summary table.
        summary_statistics: Summary statistics table, if present.

    Returns:
        Summary DataFrame.
    """
    if not summary_statistics.empty:
        df = summary_statistics.copy()
        workflow_column = get_first_existing_column(
            dataframe=df,
            candidates=["workflow", "workflow_type"],
        )
        metric_column = get_first_existing_column(
            dataframe=df,
            candidates=["metric", "metric_name"],
        )
        if workflow_column and workflow_column != "workflow":
            df = df.rename(columns={workflow_column: "workflow"})
        if metric_column and metric_column != "metric":
            df = df.rename(columns={metric_column: "metric"})
        return df

    if combined_long.empty:
        return pd.DataFrame()

    required_cols = {"workflow", "target_label", "metric", "spike_n", "value"}
    if not required_cols.issubset(set(combined_long.columns)):
        return pd.DataFrame()

    df = combined_long.copy()
    df["spike_n"] = pd.to_numeric(df["spike_n"], errors="coerce")
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["spike_n", "value"])

    rows: list[dict[str, Any]] = []
    grouped = df.groupby(["workflow", "target_label", "metric"], dropna=False)
    for (workflow, target_label, metric), group in grouped:
        baseline_values = group.loc[group["spike_n"] == 0, "value"]
        baseline_value = (
            baseline_values.mean() if not baseline_values.empty else math.nan
        )
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
                "baseline_value": (
                    float(baseline_value) if not math.isnan(baseline_value)
                    else math.nan
                ),
            }
        )

    summary_df = pd.DataFrame(rows)
    if not summary_df.empty:
        summary_df = summary_df.sort_values(
            by=["workflow", "target_label", "metric"]
        ).reset_index(drop=True)
    return summary_df


def top_signal_rows(combined_long: pd.DataFrame, n_rows: int = 25) -> pd.DataFrame:
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

    if "run_path" in columns_to_keep:
        df["run_path"] = df["run_path"].map(compact_path)

    return (
        df[columns_to_keep]
        .sort_values(by="value", ascending=False)
        .head(n_rows)
        .reset_index(drop=True)
    )


def make_html_table(
    dataframe: pd.DataFrame,
    max_rows: int | None = None,
    class_name: str = "data-table",
) -> str:
    """Render a DataFrame as an HTML table.

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

    if "qc" in text or "completion" in text:
        return "Quality control"
    if "baseline" in text:
        return "Baseline-adjusted plots"
    if "efficiency" in text:
        return "Detection-efficiency plots"
    if "assembly" in text or "contig" in text:
        return "Assembly-based plots"
    if "minimap" in text:
        return "Read-level minimap2 plots"
    if "kraken" in text:
        return "Read-level Kraken2 plots"
    return "Other plots"


def select_plot_groups(
    plot_manifest: pd.DataFrame,
    max_plots_per_section: int,
) -> dict[str, list[dict[str, str]]]:
    """Group plots for display in the report.

    Args:
        plot_manifest: Plot manifest table.
        max_plots_per_section: Maximum number of plots per section.

    Returns:
        Dictionary mapping section titles to plot records.
    """
    if plot_manifest.empty:
        return {}

    df = plot_manifest.copy()

    plot_path_column = get_first_existing_column(
        dataframe=df,
        candidates=["plot_path", "path"],
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
        candidates=["plot_kind"],
    )

    if plot_path_column is None:
        return {}

    if plot_path_column != "plot_path":
        df = df.rename(columns={plot_path_column: "plot_path"})
    if workflow_column and workflow_column != "workflow":
        df = df.rename(columns={workflow_column: "workflow"})
    if target_column and target_column != "target_label":
        df = df.rename(columns={target_column: "target_label"})
    if metric_column and metric_column != "metric":
        df = df.rename(columns={metric_column: "metric"})
    if plot_kind_column and plot_kind_column != "plot_kind":
        df = df.rename(columns={plot_kind_column: "plot_kind"})

    for column in ["workflow", "target_label", "metric", "plot_kind"]:
        if column not in df.columns:
            df[column] = ""

    df = df[df["plot_path"].astype(str).str.strip() != ""].copy()
    df = df.drop_duplicates(subset=["plot_path"])

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

    ordered_sections = [
        "Quality control",
        "Read-level Kraken2 plots",
        "Read-level minimap2 plots",
        "Assembly-based plots",
        "Baseline-adjusted plots",
        "Detection-efficiency plots",
        "Other plots",
    ]

    sections: dict[str, list[dict[str, str]]] = {}
    for section_name in ordered_sections:
        group = df[df["section"] == section_name].head(max_plots_per_section)
        if group.empty:
            continue
        sections[section_name] = group.to_dict(orient="records")

    return sections


def discover_krona_outputs(summary_dir: Path) -> list[dict[str, str]]:
    """Discover optional Krona outputs near the summary directory.

    Args:
        summary_dir: Summary directory.

    Returns:
        List of Krona output records.
    """
    roots_to_search = [summary_dir, summary_dir.parent]
    results: list[dict[str, str]] = []
    seen: set[str] = set()

    for root in roots_to_search:
        if not root.exists():
            continue

        for path in root.rglob("*"):
            if not path.exists():
                continue

            path_name = path.name.lower()
            parent_name = path.parent.name.lower()
            if "krona" not in path_name and "krona" not in parent_name:
                continue

            if path.suffix.lower() in {".html", ".htm", ".txt", ".tsv"} or path.is_dir():
                resolved = str(path.resolve())
                if resolved in seen:
                    continue
                seen.add(resolved)
                results.append(
                    {
                        "name": path.name,
                        "path": str(path),
                        "type": "directory" if path.is_dir() else path.suffix.lower().lstrip("."),
                        "path_compact": compact_path(str(path)),
                    }
                )

    results = sorted(results, key=lambda x: (x["type"], x["name"]))
    return results


def build_executive_summary(
    run_manifest: pd.DataFrame,
    issues_df: pd.DataFrame,
) -> str:
    """Build an executive summary block.

    Args:
        run_manifest: Normalised run manifest.
        issues_df: Issues table.

    Returns:
        HTML fragment.
    """
    n_runs = int(run_manifest.shape[0]) if not run_manifest.empty else 0
    n_shuffled = (
        int(run_manifest["is_shuffled_control_recomputed"].sum())
        if not run_manifest.empty else 0
    )
    n_non_shuffled = n_runs - n_shuffled
    n_workflows = (
        int(run_manifest["workflow"].nunique()) if not run_manifest.empty else 0
    )
    n_targets = (
        int(
            run_manifest["target_label_display"].replace("", pd.NA).dropna().nunique()
        )
        if not run_manifest.empty else 0
    )
    n_issues = int(issues_df.shape[0]) if not issues_df.empty else 0

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
    <p>
      This report combines stakeholder-facing overview content with the richer
      technical detail needed to inspect workflow behaviour. It recalculates
      shuffled-control status independently, embeds available plots directly in
      the page, and includes optional Krona outputs where present.
    </p>
    """


def make_plot_sections_html(
    plot_sections: dict[str, list[dict[str, str]]],
    summary_dir: Path,
) -> str:
    """Render grouped plots into HTML.

    Args:
        plot_sections: Plot groups to render.
        summary_dir: Directory containing the report.

    Returns:
        HTML fragment for plot sections.
    """
    if not plot_sections:
        return '<p class="muted">No plots were found for display.</p>'

    sections_html: list[str] = []
    for section_name, plots in plot_sections.items():
        plots_html: list[str] = []
        for record in plots:
            plot_path = Path(str(record.get("plot_path", "")))
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


def make_krona_section_html(
    krona_records: list[dict[str, str]],
    summary_dir: Path,
) -> str:
    """Render optional Krona outputs into HTML.

    Args:
        krona_records: Discovered Krona records.
        summary_dir: Directory containing the report.

    Returns:
        HTML fragment for the Krona section.
    """
    if not krona_records:
        return """
        <section class="report-section">
          <h2>Krona outputs</h2>
          <p class="muted">
            No Krona outputs were discovered. This is not treated as an error.
          </p>
        </section>
        """

    rows = []
    preview_blocks: list[str] = []
    for record in krona_records:
        path = Path(record["path"])
        href = relative_href(target_path=path, from_dir=summary_dir)
        rows.append(
            {
                "Name": record["name"],
                "Type": record["type"],
                "Path": record["path_compact"],
            }
        )
        if path.suffix.lower() in {".html", ".htm"}:
            preview_blocks.append(
                f"""
                <div class="krona-card">
                  <div class="krona-title">
                    <a href="{html_escape(href)}" target="_blank">
                      {html_escape(record["name"])}
                    </a>
                  </div>
                  <div class="krona-path">{html_escape(record["path_compact"])}</div>
                </div>
                """
            )

    krona_table = make_html_table(pd.DataFrame(rows))
    krona_previews = (
        '<div class="krona-grid">' + "".join(preview_blocks) + "</div>"
        if preview_blocks else ""
    )

    return f"""
    <section class="report-section">
      <h2>Krona outputs</h2>
      <p class="muted">
        Krona artefacts were discovered automatically and are listed here when
        present. Links open the original Krona output in a new tab.
      </p>
      {krona_table}
      {krona_previews}
    </section>
    """


def build_html_report(
    title: str,
    run_manifest: pd.DataFrame,
    combined_wide: pd.DataFrame,
    combined_long: pd.DataFrame,
    summary_statistics: pd.DataFrame,
    issues_df: pd.DataFrame,
    plot_sections: dict[str, list[dict[str, str]]],
    krona_records: list[dict[str, str]],
    summary_dir: Path,
    max_rows: int,
) -> str:
    """Build the complete HTML report.

    Args:
        title: Report title.
        run_manifest: Normalised run manifest.
        combined_wide: Combined wide table.
        combined_long: Normalised combined long table.
        summary_statistics: Summary statistics table.
        issues_df: Issues table.
        plot_sections: Grouped plot records.
        krona_records: Discovered Krona records.
        summary_dir: Report output directory.
        max_rows: Maximum number of rows to display in large tables.

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

    summary_display = summary_statistics.copy()
    top_signals_display = top_signal_rows(combined_long=combined_long, n_rows=25)

    combined_long_display = combined_long.copy()
    if "run_path" in combined_long_display.columns:
        combined_long_display["run_path"] = combined_long_display["run_path"].map(
            compact_path
        )

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
        max-width: 1360px;
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
      .krona-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
        gap: 14px;
        margin-top: 16px;
      }
      .krona-card {
        border: 1px solid #d6dfeb;
        border-radius: 10px;
        background: #ffffff;
        padding: 12px;
      }
      .krona-title {
        font-weight: bold;
        margin-bottom: 6px;
      }
      .krona-path {
        font-size: 13px;
        color: #4a5d73;
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
      This report combines the richer embedded-plot behaviour from the earlier
      report with the cleaner layout and corrected shuffled-control logic from
      the newer version.
    </p>

    <div class="note-box">
      <strong>How to use this report:</strong>
      Start with the executive summary and workflow overview, then review the
      embedded plots and summary statistics. More detailed tables are provided
      below for technical follow-up. Optional Krona outputs are listed when
      available and are not required for report generation.
    </div>

    <section class="report-section">
      <h2>Executive summary</h2>
      {build_executive_summary(run_manifest=run_manifest, issues_df=issues_df)}
    </section>

    <section class="report-section">
      <h2>Workflow overview</h2>
      {make_html_table(workflow_counts)}
    </section>

    <section class="report-section">
      <h2>Runs discovered</h2>
      <p class="muted">
        Run paths are compacted for readability. Shuffled-control status is
        recomputed in this report rather than taken from any pre-existing flag.
      </p>
      {make_html_table(runs_display, max_rows=max_rows)}
    </section>

    <section class="report-section">
      <h2>Summary statistics</h2>
      {make_html_table(summary_display, max_rows=max_rows)}
    </section>

    <section class="report-section">
      <h2>Top observed signals</h2>
      <p class="muted">
        Highest observed values from the long-format summary table. This can be
        useful for quickly identifying the strongest signal points across
        workflows, targets, and metrics.
      </p>
      {make_html_table(top_signals_display, max_rows=max_rows)}
    </section>

    {make_plot_sections_html(plot_sections=plot_sections, summary_dir=summary_dir)}

    {make_krona_section_html(krona_records=krona_records, summary_dir=summary_dir)}

    <section class="report-section">
      <h2>Combined long table preview</h2>
      {make_html_table(combined_long_display, max_rows=max_rows)}
    </section>

    <section class="report-section">
      <h2>Combined wide table preview</h2>
      {make_html_table(combined_wide, max_rows=max_rows)}
    </section>

    <section class="report-section">
      <h2>Logged issues</h2>
      <p class="muted">
        Missing or incomplete outputs may reflect runs that are still in
        progress as well as genuinely problematic outputs.
      </p>
      {make_html_table(issues_display, max_rows=max_rows)}
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
        else summary_dir / "spikein_report_v4.html"
    )

    run_manifest_path = summary_dir / "run_manifest.tsv"
    combined_wide_path = summary_dir / "combined_wide.tsv"
    combined_long_path = summary_dir / "combined_long.tsv"
    summary_stats_path = summary_dir / "summary_statistics.tsv"
    issues_path = summary_dir / "missing_or_problematic_files.tsv"
    plot_manifest_path = summary_dir / "plot_manifest.tsv"

    run_manifest = coerce_numeric_columns(read_tsv_if_exists(run_manifest_path))
    combined_wide = coerce_numeric_columns(read_tsv_if_exists(combined_wide_path))
    combined_long = coerce_numeric_columns(read_tsv_if_exists(combined_long_path))
    summary_statistics = coerce_numeric_columns(
        read_tsv_if_exists(summary_stats_path)
    )
    issues_df = coerce_numeric_columns(read_tsv_if_exists(issues_path))
    plot_manifest = coerce_numeric_columns(read_tsv_if_exists(plot_manifest_path))

    if run_manifest.empty:
        raise ValueError(
            f"No run_manifest.tsv found in summary directory: {summary_dir}"
        )

    manifest_norm = normalise_manifest(run_manifest=run_manifest)
    long_norm = normalise_long_table(combined_long=combined_long)
    summary_stats_norm = summarise_key_metrics(
        combined_long=long_norm,
        summary_statistics=summary_statistics,
    )
    plot_sections = select_plot_groups(
        plot_manifest=plot_manifest,
        max_plots_per_section=args.max_plots_per_section,
    )
    krona_records = discover_krona_outputs(summary_dir=summary_dir)

    html_text = build_html_report(
        title=args.title,
        run_manifest=manifest_norm,
        combined_wide=combined_wide,
        combined_long=long_norm,
        summary_statistics=summary_stats_norm,
        issues_df=issues_df,
        plot_sections=plot_sections,
        krona_records=krona_records,
        summary_dir=output_html.parent,
        max_rows=args.max_rows,
    )

    output_html.write_text(html_text, encoding="utf-8")
    print(f"[INFO] Wrote {output_html}")


if __name__ == "__main__":
    main()
