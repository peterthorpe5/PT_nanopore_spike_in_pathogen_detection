#!/usr/bin/env python3
"""Create HTML and PDF reports from spike-in summary outputs.

This script reads the outputs produced by summarise_spikein_runs_v2.py and
builds two human-readable reports:
- spikein_report.html
- spikein_report.pdf

The HTML report is intended for quick browsing and sharing. The PDF report is a
compact static summary suitable for emailing or sharing with collaborators.
"""

from __future__ import annotations

import argparse
import html
import math
from pathlib import Path
from typing import Iterable

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


PLOT_GROUP_ORDER = ["qc", "raw", "baseline_adjusted", "efficiency"]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for report generation."""
    parser = argparse.ArgumentParser(
        description="Build HTML and PDF reports from spike-in summary outputs."
    )
    parser.add_argument(
        "--summary_dir",
        required=True,
        help="Directory produced by summarise_spikein_runs_v2.py.",
    )
    parser.add_argument(
        "--title",
        default="ONT spike-in summary report",
        help="Title to use in the HTML and PDF reports.",
    )
    parser.add_argument(
        "--max_plots_per_group",
        type=int,
        default=12,
        help="Maximum number of plots to embed per plot group in the reports.",
    )
    return parser.parse_args()



def safe_read_tsv(path: Path) -> pd.DataFrame:
    """Read a TSV file, returning an empty table if it is unavailable."""
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame()



def build_overview_lines(
    run_manifest: pd.DataFrame,
    combined_wide: pd.DataFrame,
    combined_long: pd.DataFrame,
    issues_df: pd.DataFrame,
    plot_manifest: pd.DataFrame,
) -> list[str]:
    """Build short overview lines for the report header."""
    workflows = sorted(run_manifest["workflow_type"].dropna().astype(str).unique()) if not run_manifest.empty else []
    shuffled_n = int(run_manifest["is_shuffled_control"].fillna(False).sum()) if not run_manifest.empty and "is_shuffled_control" in run_manifest.columns else 0
    lines = [
        f"Runs discovered: {len(run_manifest)}",
        f"Parsed summary rows: {len(combined_wide)}",
        f"Long-form metric rows: {len(combined_long)}",
        f"Problem records: {len(issues_df)}",
        f"Plots generated: {len(plot_manifest)}",
        f"Shuffled-control runs: {shuffled_n}",
        f"Workflow types: {', '.join(workflows) if workflows else 'none detected'}",
    ]
    return lines



def build_key_findings(summary_stats: pd.DataFrame, issues_df: pd.DataFrame) -> list[str]:
    """Build a small set of plain-language findings for the report."""
    findings: list[str] = []
    if not summary_stats.empty:
        top_rows = summary_stats.sort_values(
            by=["max_metric_value", "n_rows"],
            ascending=[False, False],
            kind="stable",
        ).head(5)
        for _, row in top_rows.iterrows():
            findings.append(
                (
                    f"{row.get('workflow_type')} / {row.get('target_label')} / "
                    f"{row.get('metric_name')}: max={row.get('max_metric_value')}, "
                    f"median={row.get('median_metric_value')}, rows={row.get('n_rows')}"
                )
            )
    if not issues_df.empty:
        issue_counts = issues_df["category"].value_counts().head(5)
        findings.append(
            "Most common issue categories: "
            + ", ".join([f"{idx} ({val})" for idx, val in issue_counts.items()])
        )
    if not findings:
        findings.append("No summary statistics were available for interpretation.")
    return findings



def dataframe_to_html_table(dataframe: pd.DataFrame, max_rows: int = 30) -> str:
    """Render a DataFrame as a simple HTML table snippet."""
    if dataframe.empty:
        return "<p>No data available.</p>"
    return dataframe.head(max_rows).to_html(index=False, border=0, classes=["dataframe", "compact"])



def select_plots(plot_manifest: pd.DataFrame, max_plots_per_group: int) -> dict[str, list[Path]]:
    """Select plot images for embedding in the reports."""
    selected: dict[str, list[Path]] = {}
    if plot_manifest.empty:
        return selected

    for plot_group in PLOT_GROUP_ORDER:
        group_df = plot_manifest.loc[plot_manifest["plot_group"] == plot_group].copy()
        if group_df.empty:
            continue
        group_df = group_df.sort_values(
            by=["workflow_type", "metric_name", "target_label"],
            kind="stable",
        ).head(max_plots_per_group)
        selected[plot_group] = [Path(path) for path in group_df["png_path"].dropna().tolist()]

    remaining_groups = sorted(set(plot_manifest["plot_group"].dropna().astype(str)) - set(PLOT_GROUP_ORDER))
    for plot_group in remaining_groups:
        group_df = plot_manifest.loc[plot_manifest["plot_group"] == plot_group].copy()
        if group_df.empty:
            continue
        group_df = group_df.sort_values(
            by=["workflow_type", "metric_name", "target_label"],
            kind="stable",
        ).head(max_plots_per_group)
        selected[plot_group] = [Path(path) for path in group_df["png_path"].dropna().tolist()]

    return selected



def write_html_report(
    out_path: Path,
    title: str,
    overview_lines: Iterable[str],
    key_findings: Iterable[str],
    run_manifest: pd.DataFrame,
    summary_stats: pd.DataFrame,
    issues_df: pd.DataFrame,
    selected_plots: dict[str, list[Path]],
) -> None:
    """Write the HTML report to disk."""
    sections: list[str] = []
    sections.append(f"<h1>{html.escape(title)}</h1>")
    sections.append("<h2>Overview</h2>")
    sections.append("<ul>" + "".join([f"<li>{html.escape(line)}</li>" for line in overview_lines]) + "</ul>")

    sections.append("<h2>Key findings</h2>")
    sections.append("<ul>" + "".join([f"<li>{html.escape(line)}</li>" for line in key_findings]) + "</ul>")

    sections.append("<h2>Run manifest</h2>")
    sections.append(dataframe_to_html_table(run_manifest, max_rows=50))

    sections.append("<h2>Summary statistics</h2>")
    sections.append(dataframe_to_html_table(summary_stats, max_rows=50))

    sections.append("<h2>Issues</h2>")
    sections.append(dataframe_to_html_table(issues_df, max_rows=50))

    sections.append("<h2>Plots</h2>")
    for plot_group, paths in selected_plots.items():
        sections.append(f"<h3>{html.escape(plot_group.replace('_', ' ').title())}</h3>")
        if not paths:
            sections.append("<p>No plots selected for this group.</p>")
            continue
        sections.append('<div class="plot-grid">')
        for plot_path in paths:
            rel_path = plot_path.name if plot_path.is_absolute() else str(plot_path)
            sections.append(
                "<figure>"
                f"<img src=\"plots/{html.escape(Path(rel_path).name)}\" alt=\"{html.escape(Path(rel_path).name)}\">"
                f"<figcaption>{html.escape(Path(rel_path).name)}</figcaption>"
                "</figure>"
            )
        sections.append("</div>")

    html_text = f"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{html.escape(title)}</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 2rem; color: #222; }}
h1, h2, h3 {{ color: #1F4E79; }}
table.dataframe {{ border-collapse: collapse; width: 100%; margin-bottom: 1.5rem; }}
table.dataframe th {{ background: #1F4E79; color: white; padding: 0.5rem; text-align: left; }}
table.dataframe td {{ border-bottom: 1px solid #ddd; padding: 0.4rem; vertical-align: top; }}
table.compact {{ font-size: 0.9rem; }}
.plot-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 1rem; }}
figure {{ margin: 0; border: 1px solid #ddd; padding: 0.5rem; background: #fafafa; }}
img {{ max-width: 100%; height: auto; display: block; }}
figcaption {{ font-size: 0.85rem; margin-top: 0.5rem; color: #444; }}
</style>
</head>
<body>
{''.join(sections)}
</body>
</html>
"""
    out_path.write_text(html_text, encoding="utf-8")



def add_text_page(pdf: PdfPages, title: str, lines: Iterable[str]) -> None:
    """Add a text-only page to the PDF report."""
    fig = plt.figure(figsize=(8.27, 11.69))
    plt.axis("off")
    y_pos = 0.96
    fig.text(0.08, y_pos, title, fontsize=18, fontweight="bold")
    y_pos -= 0.05
    for line in lines:
        fig.text(0.08, y_pos, f"• {line}", fontsize=11)
        y_pos -= 0.03
        if y_pos < 0.08:
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            fig = plt.figure(figsize=(8.27, 11.69))
            plt.axis("off")
            y_pos = 0.96
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)



def add_table_page(pdf: PdfPages, title: str, dataframe: pd.DataFrame, max_rows: int = 20) -> None:
    """Add a table page to the PDF report using a matplotlib table."""
    fig, ax = plt.subplots(figsize=(11.69, 8.27))
    ax.axis("off")
    ax.set_title(title, fontsize=16, fontweight="bold", pad=12)
    if dataframe.empty:
        ax.text(0.05, 0.9, "No data available.", fontsize=12)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
        return

    display_df = dataframe.head(max_rows).copy()
    table = ax.table(
        cellText=display_df.astype(str).values,
        colLabels=list(display_df.columns),
        cellLoc="left",
        colLoc="left",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.3)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)



def add_image_page(pdf: PdfPages, image_path: Path, title: str) -> None:
    """Add a full-page plot image to the PDF report."""
    fig, ax = plt.subplots(figsize=(11.69, 8.27))
    ax.axis("off")
    ax.set_title(title, fontsize=14, fontweight="bold", pad=10)
    image = mpimg.imread(image_path)
    ax.imshow(image)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)



def write_pdf_report(
    out_path: Path,
    title: str,
    overview_lines: Iterable[str],
    key_findings: Iterable[str],
    run_manifest: pd.DataFrame,
    summary_stats: pd.DataFrame,
    issues_df: pd.DataFrame,
    selected_plots: dict[str, list[Path]],
) -> None:
    """Write the PDF report to disk."""
    with PdfPages(out_path) as pdf:
        add_text_page(pdf, title, overview_lines)
        add_text_page(pdf, "Key findings", key_findings)
        add_table_page(pdf, "Run manifest", run_manifest, max_rows=20)
        add_table_page(pdf, "Summary statistics", summary_stats, max_rows=20)
        add_table_page(pdf, "Issues", issues_df, max_rows=20)
        for plot_group, paths in selected_plots.items():
            for image_path in paths:
                if image_path.exists():
                    add_image_page(
                        pdf=pdf,
                        image_path=image_path,
                        title=f"{plot_group.replace('_', ' ').title()} - {image_path.name}",
                    )



def main() -> None:
    """Run the HTML and PDF report generation workflow."""
    args = parse_args()
    summary_dir = Path(args.summary_dir).expanduser().resolve()

    run_manifest = safe_read_tsv(summary_dir / "run_manifest.tsv")
    combined_wide = safe_read_tsv(summary_dir / "combined_wide.tsv")
    combined_long = safe_read_tsv(summary_dir / "combined_long.tsv")
    summary_stats = safe_read_tsv(summary_dir / "summary_statistics.tsv")
    issues_df = safe_read_tsv(summary_dir / "missing_or_problematic_files.tsv")
    plot_manifest = safe_read_tsv(summary_dir / "plot_manifest.tsv")

    overview_lines = build_overview_lines(
        run_manifest=run_manifest,
        combined_wide=combined_wide,
        combined_long=combined_long,
        issues_df=issues_df,
        plot_manifest=plot_manifest,
    )
    key_findings = build_key_findings(
        summary_stats=summary_stats,
        issues_df=issues_df,
    )
    selected_plots = select_plots(
        plot_manifest=plot_manifest,
        max_plots_per_group=args.max_plots_per_group,
    )

    write_html_report(
        out_path=summary_dir / "spikein_report.html",
        title=args.title,
        overview_lines=overview_lines,
        key_findings=key_findings,
        run_manifest=run_manifest,
        summary_stats=summary_stats,
        issues_df=issues_df,
        selected_plots=selected_plots,
    )
    write_pdf_report(
        out_path=summary_dir / "spikein_report.pdf",
        title=args.title,
        overview_lines=overview_lines,
        key_findings=key_findings,
        run_manifest=run_manifest,
        summary_stats=summary_stats,
        issues_df=issues_df,
        selected_plots=selected_plots,
    )

    print(f"[INFO] HTML report: {summary_dir / 'spikein_report.html'}")
    print(f"[INFO] PDF report: {summary_dir / 'spikein_report.pdf'}")


if __name__ == "__main__":
    main()
