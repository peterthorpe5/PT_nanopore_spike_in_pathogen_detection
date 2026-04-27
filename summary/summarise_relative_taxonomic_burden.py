#!/usr/bin/env python3
"""Summarise expected versus off-target taxonomic signal burden.

This script is intended to run after the spike-in summary workflow has produced
``reported_taxa_long.tsv``. It estimates, for each observation, how much of the
reported taxonomic signal belongs to expected target taxa and how much belongs
to non-expected off-target taxa. The resulting tables are designed to complement
confusion-matrix and clean-sensitivity summaries by showing signal strength, not
only whether an off-target taxon was present.

The output values should be interpreted as report-level taxonomic burden rather
than exact read partitioning. Kraken-style clade counts can include hierarchical
signal, so these ratios are best used as practical signal-to-noise summaries for
report interpretation.
"""

from __future__ import annotations

import argparse
import html
import math
from pathlib import Path
from typing import Iterable

import pandas as pd


ID_COLUMNS = [
    "workflow",
    "classifier",
    "run_name",
    "replicate",
    "spike_n_per_genome",
    "total_spike_n",
    "n_genomes",
    "is_shuffled_control",
]

PREFERRED_METRIC_PATTERNS = [
    "reported_clade_reads",
    "reported_direct_reads",
    "target_reads",
    "target_contigs",
    "alignments",
    "unique_reads",
    "found",
]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Summarise relative expected and off-target taxonomic signal from "
            "reported_taxa_long.tsv."
        )
    )
    parser.add_argument(
        "--reported_taxa_long_tsv",
        required=True,
        help="Path to reported_taxa_long.tsv from the spike-in summary workflow.",
    )
    parser.add_argument(
        "--combined_long_tsv",
        required=False,
        default=None,
        help=(
            "Optional combined_long.tsv. Used only as a fallback source of "
            "expected targets when reported_taxa_long.tsv lacks matched labels."
        ),
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for TSV and HTML summaries.",
    )
    parser.add_argument(
        "--min_signal",
        type=float,
        default=0.0,
        help="Minimum metric value retained before summarisation. Default: 0.",
    )
    parser.add_argument(
        "--restrict_to_plasmodium_like_offtargets",
        action="store_true",
        help=(
            "If set, off-target signal is restricted to taxa whose names look "
            "Plasmodium-like. By default, all non-expected taxa are included."
        ),
    )
    return parser.parse_args()


def normalise_columns(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Normalise common column names in reported taxa tables.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Raw reported taxa table.

    Returns
    -------
    pd.DataFrame
        Normalised table.
    """
    df = dataframe.copy()
    rename_map = {}
    if "workflow_type" in df.columns and "workflow" not in df.columns:
        rename_map["workflow_type"] = "workflow"
    if "metric" in df.columns and "metric_name" not in df.columns:
        rename_map["metric"] = "metric_name"
    if "value" in df.columns and "metric_value" not in df.columns:
        rename_map["value"] = "metric_value"
    if "spike_n" in df.columns and "spike_n_per_genome" not in df.columns:
        rename_map["spike_n"] = "spike_n_per_genome"
    if rename_map:
        df = df.rename(columns=rename_map)

    defaults = {
        "workflow": "unknown_workflow",
        "classifier": "unknown_classifier",
        "run_name": "unknown_run",
        "replicate": pd.NA,
        "spike_n_per_genome": pd.NA,
        "total_spike_n": pd.NA,
        "n_genomes": pd.NA,
        "is_shuffled_control": False,
        "taxon_name": "unknown_taxon",
        "metric_name": "unknown_metric",
        "metric_value": 0.0,
        "matched_target_labels": "",
    }
    for column, value in defaults.items():
        if column not in df.columns:
            df[column] = value

    if df["total_spike_n"].isna().all():
        df["total_spike_n"] = df["spike_n_per_genome"]

    text_columns = [
        "workflow",
        "classifier",
        "run_name",
        "taxon_name",
        "metric_name",
        "matched_target_labels",
    ]
    for column in text_columns:
        df[column] = df[column].fillna("").astype(str)

    numeric_columns = [
        "replicate",
        "spike_n_per_genome",
        "total_spike_n",
        "n_genomes",
        "metric_value",
        "percent",
    ]
    for column in numeric_columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

    df["is_shuffled_control"] = (
        df["is_shuffled_control"]
        .fillna(False)
        .astype(str)
        .str.lower()
        .isin(["true", "1", "yes"])
    )
    df["is_positive"] = (~df["is_shuffled_control"]) & (df["total_spike_n"] > 0)
    df["is_negative"] = df["is_shuffled_control"] | (df["total_spike_n"] == 0)
    return df


def load_expected_targets_from_combined(path: str | None) -> dict[tuple[str, str], set[str]]:
    """Load expected targets from combined_long.tsv as a fallback.

    Parameters
    ----------
    path : str | None
        Optional path to combined_long.tsv.

    Returns
    -------
    dict[tuple[str, str], set[str]]
        Mapping from ``(workflow, run_name)`` to expected target labels.
    """
    if not path:
        return {}
    input_path = Path(path)
    if not input_path.exists():
        return {}

    df = pd.read_csv(input_path, sep="\t")
    if "workflow_type" in df.columns and "workflow" not in df.columns:
        df = df.rename(columns={"workflow_type": "workflow"})
    if "target_label" not in df.columns:
        return {}
    if "run_name" not in df.columns:
        df["run_name"] = "unknown_run"
    if "workflow" not in df.columns:
        df["workflow"] = "unknown_workflow"

    spike_col = "total_spike_n" if "total_spike_n" in df.columns else None
    if spike_col is None:
        spike_col = "spike_n_per_genome" if "spike_n_per_genome" in df.columns else None
    if spike_col is not None:
        df[spike_col] = pd.to_numeric(df[spike_col], errors="coerce")
        df = df.loc[df[spike_col] > 0].copy()

    mapping: dict[tuple[str, str], set[str]] = {}
    for (workflow, run_name), group in df.groupby(["workflow", "run_name"], dropna=False):
        labels = set(group["target_label"].dropna().astype(str))
        labels = {x for x in labels if x and x != "single_target"}
        if labels:
            mapping[(str(workflow), str(run_name))] = labels
    return mapping


def split_label_string(value: str) -> set[str]:
    """Split matched target label strings into a set.

    Parameters
    ----------
    value : str
        Raw label string.

    Returns
    -------
    set[str]
        Labels found in the string.
    """
    if not value:
        return set()
    text = value.replace(",", ";")
    return {part.strip() for part in text.split(";") if part.strip()}


def looks_plasmodium_like(taxon_name: str) -> bool:
    """Return whether a taxon name appears to be Plasmodium-like.

    Parameters
    ----------
    taxon_name : str
        Taxon name.

    Returns
    -------
    bool
        True if the name appears Plasmodium-like.
    """
    lower = str(taxon_name).lower().strip()
    return "plasmodium" in lower or lower.startswith("p. ")


def choose_primary_metric(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Choose one primary signal metric per workflow/classifier/taxon context.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Normalised reported taxa table.

    Returns
    -------
    pd.DataFrame
        Subset containing one preferred metric family where possible.
    """
    df = dataframe.copy()
    metric_values = df["metric_name"].dropna().astype(str).unique().tolist()
    for pattern in PREFERRED_METRIC_PATTERNS:
        matched = [metric for metric in metric_values if metric.endswith(pattern)]
        if matched:
            return df.loc[df["metric_name"].isin(matched)].copy()
    return df


def classify_rows(
    dataframe: pd.DataFrame,
    expected_map: dict[tuple[str, str], set[str]],
    restrict_plasmodium_like: bool,
) -> pd.DataFrame:
    """Classify rows as expected or off-target signal.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Normalised reported taxa table.
    expected_map : dict[tuple[str, str], set[str]]
        Fallback expected target mapping.
    restrict_plasmodium_like : bool
        Whether to restrict off-targets to Plasmodium-like taxa.

    Returns
    -------
    pd.DataFrame
        Table with classification columns added.
    """
    df = dataframe.copy()
    expected_flags = []
    off_target_flags = []
    expected_label_strings = []

    for _, row in df.iterrows():
        matched_labels = split_label_string(row.get("matched_target_labels", ""))
        fallback_labels = expected_map.get(
            (str(row.get("workflow", "")), str(row.get("run_name", ""))),
            set(),
        )
        expected_labels = matched_labels or fallback_labels
        taxon_name = str(row.get("taxon_name", ""))

        is_expected = bool(matched_labels)
        if not is_expected and fallback_labels:
            is_expected = taxon_name in fallback_labels

        is_off_target = not is_expected
        if restrict_plasmodium_like and is_off_target:
            is_off_target = looks_plasmodium_like(taxon_name)

        expected_flags.append(is_expected)
        off_target_flags.append(is_off_target)
        expected_label_strings.append("; ".join(sorted(expected_labels)))

    df["expected_labels_for_observation"] = expected_label_strings
    df["is_expected_signal"] = expected_flags
    df["is_off_target_signal"] = off_target_flags
    return df


def safe_divide(numerator: float, denominator: float) -> float:
    """Safely divide two numbers.

    Parameters
    ----------
    numerator : float
        Numerator.
    denominator : float
        Denominator.

    Returns
    -------
    float
        Division result or NaN.
    """
    if denominator == 0 or pd.isna(denominator):
        return math.nan
    return numerator / denominator


def join_top_taxa(names: Iterable[str], values: Iterable[float], limit: int = 5) -> str:
    """Join the top taxa and values as a compact string.

    Parameters
    ----------
    names : Iterable[str]
        Taxon names.
    values : Iterable[float]
        Signal values.
    limit : int
        Maximum number of taxa to report.

    Returns
    -------
    str
        Compact summary string.
    """
    pairs = sorted(
        [(str(name), float(value)) for name, value in zip(names, values) if pd.notna(value)],
        key=lambda item: item[1],
        reverse=True,
    )
    return "; ".join([f"{name}={value:g}" for name, value in pairs[:limit]])


def summarise_by_observation(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Summarise expected and off-target signal by observation.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Classified reported taxa table.

    Returns
    -------
    pd.DataFrame
        One row per workflow/run/replicate/spike observation.
    """
    rows: list[dict[str, object]] = []
    group_cols = [col for col in ID_COLUMNS if col in dataframe.columns]

    for keys, group in dataframe.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        record = dict(zip(group_cols, keys))

        expected = group.loc[group["is_expected_signal"]].copy()
        off_target = group.loc[group["is_off_target_signal"]].copy()

        expected_signal = float(expected["metric_value"].sum())
        off_target_signal = float(off_target["metric_value"].sum())
        total_signal = expected_signal + off_target_signal

        dominant_off_target = ""
        dominant_off_target_signal = math.nan
        if not off_target.empty:
            off_by_taxon = (
                off_target.groupby("taxon_name", dropna=False)["metric_value"]
                .sum()
                .reset_index()
                .sort_values(by="metric_value", ascending=False)
            )
            dominant_off_target = str(off_by_taxon.iloc[0]["taxon_name"])
            dominant_off_target_signal = float(off_by_taxon.iloc[0]["metric_value"])

        record.update(
            {
                "metric_name": "; ".join(sorted(group["metric_name"].astype(str).unique())),
                "expected_labels": "; ".join(
                    sorted(
                        label
                        for value in group["expected_labels_for_observation"].astype(str).unique()
                        for label in split_label_string(value)
                    )
                ),
                "n_reported_taxa": int(group["taxon_name"].nunique()),
                "n_expected_taxa": int(expected["taxon_name"].nunique()),
                "n_off_target_taxa": int(off_target["taxon_name"].nunique()),
                "expected_signal": expected_signal,
                "off_target_signal": off_target_signal,
                "total_expected_plus_off_target_signal": total_signal,
                "off_target_to_expected_ratio": safe_divide(off_target_signal, expected_signal),
                "expected_signal_fraction": safe_divide(expected_signal, total_signal),
                "off_target_signal_fraction": safe_divide(off_target_signal, total_signal),
                "dominant_off_target_taxon": dominant_off_target,
                "dominant_off_target_signal": dominant_off_target_signal,
                "dominant_off_target_fraction_of_offtarget": safe_divide(
                    dominant_off_target_signal,
                    off_target_signal,
                ),
                "top_expected_taxa": join_top_taxa(
                    expected["taxon_name"],
                    expected["metric_value"],
                ),
                "top_off_target_taxa": join_top_taxa(
                    off_target["taxon_name"],
                    off_target["metric_value"],
                ),
            }
        )
        rows.append(record)

    return pd.DataFrame(rows)


def summarise_by_workflow(observation_df: pd.DataFrame) -> pd.DataFrame:
    """Summarise relative burden by workflow and classifier.

    Parameters
    ----------
    observation_df : pd.DataFrame
        Observation-level relative burden table.

    Returns
    -------
    pd.DataFrame
        Workflow-level summary table.
    """
    if observation_df.empty:
        return pd.DataFrame()

    df = observation_df.copy()
    df["is_positive"] = (~df["is_shuffled_control"].astype(bool)) & (df["total_spike_n"] > 0)
    pos_df = df.loc[df["is_positive"]].copy()
    if pos_df.empty:
        pos_df = df.copy()

    grouped = pos_df.groupby(["workflow", "classifier"], dropna=False)
    summary = grouped.agg(
        n_positive_observations=("workflow", "size"),
        mean_expected_signal=("expected_signal", "mean"),
        median_expected_signal=("expected_signal", "median"),
        mean_off_target_signal=("off_target_signal", "mean"),
        median_off_target_signal=("off_target_signal", "median"),
        mean_expected_signal_fraction=("expected_signal_fraction", "mean"),
        median_expected_signal_fraction=("expected_signal_fraction", "median"),
        mean_off_target_signal_fraction=("off_target_signal_fraction", "mean"),
        median_off_target_signal_fraction=("off_target_signal_fraction", "median"),
        mean_off_target_to_expected_ratio=("off_target_to_expected_ratio", "mean"),
        median_off_target_to_expected_ratio=("off_target_to_expected_ratio", "median"),
        mean_n_off_target_taxa=("n_off_target_taxa", "mean"),
        median_n_off_target_taxa=("n_off_target_taxa", "median"),
    ).reset_index()
    return summary.sort_values(
        by=["median_expected_signal_fraction", "median_off_target_signal_fraction"],
        ascending=[False, True],
    )


def summarise_by_spike(observation_df: pd.DataFrame) -> pd.DataFrame:
    """Summarise relative burden by workflow, classifier, and spike level.

    Parameters
    ----------
    observation_df : pd.DataFrame
        Observation-level relative burden table.

    Returns
    -------
    pd.DataFrame
        Spike-level summary table.
    """
    if observation_df.empty:
        return pd.DataFrame()

    df = observation_df.copy()
    df = df.loc[(~df["is_shuffled_control"].astype(bool)) & (df["total_spike_n"] > 0)].copy()
    if df.empty:
        return pd.DataFrame()

    grouped = df.groupby(["workflow", "classifier", "total_spike_n"], dropna=False)
    return grouped.agg(
        n_observations=("workflow", "size"),
        median_expected_signal=("expected_signal", "median"),
        median_off_target_signal=("off_target_signal", "median"),
        median_expected_signal_fraction=("expected_signal_fraction", "median"),
        median_off_target_signal_fraction=("off_target_signal_fraction", "median"),
        median_off_target_to_expected_ratio=("off_target_to_expected_ratio", "median"),
        median_n_off_target_taxa=("n_off_target_taxa", "median"),
    ).reset_index()


def write_html_report(
    observation_df: pd.DataFrame,
    workflow_df: pd.DataFrame,
    spike_df: pd.DataFrame,
    out_path: Path,
) -> None:
    """Write an HTML report for the relative burden summary.

    Parameters
    ----------
    observation_df : pd.DataFrame
        Observation-level table.
    workflow_df : pd.DataFrame
        Workflow-level table.
    spike_df : pd.DataFrame
        Spike-level table.
    out_path : Path
        Output HTML path.
    """
    style = """
    <style>
      body { font-family: Arial, Helvetica, sans-serif; margin: 28px; color: #1a1a1a; }
      h1, h2 { color: #1f4e79; }
      .note { max-width: 1000px; line-height: 1.45; }
      .table-wrap { overflow-x: auto; border: 1px solid #d9e2ef; margin: 18px 0; }
      table { border-collapse: collapse; font-size: 13px; width: 100%; }
      th { background: #1f4e79; color: white; padding: 7px; text-align: left; white-space: nowrap; }
      td { border-bottom: 1px solid #e6edf5; padding: 6px; vertical-align: top; }
      tr:nth-child(even) { background: #fbfdff; }
      .small { color: #555; font-size: 13px; }
    </style>
    """

    def to_html_table(df: pd.DataFrame, max_rows: int = 100) -> str:
        if df.empty:
            return "<p class='small'>No rows available.</p>"
        view = df.head(max_rows).copy()
        return "<div class='table-wrap'>" + view.to_html(index=False, escape=True) + "</div>"

    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Relative taxonomic burden summary</title>
{style}
</head>
<body>
<h1>Relative taxonomic burden summary</h1>
<p class="note">
This report summarises the balance between expected target signal and non-expected
reported taxa. It complements clean-sensitivity and off-target-count summaries by
asking whether off-target taxa are weak background signal or a substantial fraction
of the reported taxonomic signal. Ratios should be interpreted as report-level
signal-to-noise measures, not exact read partitions, because some taxonomic report
counts are hierarchical.
</p>
<h2>Workflow-level summary</h2>
{to_html_table(workflow_df, max_rows=200)}
<h2>Spike-level summary</h2>
{to_html_table(spike_df, max_rows=300)}
<h2>Observation-level detail</h2>
<p class="small">Showing the first 200 rows only. The full table is written as TSV.</p>
{to_html_table(observation_df, max_rows=200)}
</body>
</html>
"""
    out_path.write_text(html_text, encoding="utf-8")


def main() -> None:
    """Run the relative taxonomic burden workflow."""
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    reported_df = pd.read_csv(args.reported_taxa_long_tsv, sep="\t")
    reported_df = normalise_columns(reported_df)
    reported_df = reported_df.loc[reported_df["metric_value"] > args.min_signal].copy()
    reported_df = choose_primary_metric(reported_df)

    expected_map = load_expected_targets_from_combined(args.combined_long_tsv)
    classified_df = classify_rows(
        reported_df,
        expected_map=expected_map,
        restrict_plasmodium_like=args.restrict_to_plasmodium_like_offtargets,
    )

    observation_df = summarise_by_observation(classified_df)
    workflow_df = summarise_by_workflow(observation_df)
    spike_df = summarise_by_spike(observation_df)

    classified_df.to_csv(
        out_dir / "relative_taxonomic_burden_input_classified.tsv",
        sep="\t",
        index=False,
    )
    observation_df.to_csv(
        out_dir / "relative_taxonomic_burden_by_observation.tsv",
        sep="\t",
        index=False,
    )
    workflow_df.to_csv(
        out_dir / "relative_taxonomic_burden_by_workflow.tsv",
        sep="\t",
        index=False,
    )
    spike_df.to_csv(
        out_dir / "relative_taxonomic_burden_by_spike.tsv",
        sep="\t",
        index=False,
    )
    write_html_report(
        observation_df=observation_df,
        workflow_df=workflow_df,
        spike_df=spike_df,
        out_path=out_dir / "relative_taxonomic_burden_report.html",
    )

    manifest = pd.DataFrame(
        [
            {
                "output": "relative_taxonomic_burden_input_classified.tsv",
                "description": "Reported taxa table with expected/off-target classification columns added.",
            },
            {
                "output": "relative_taxonomic_burden_by_observation.tsv",
                "description": "One row per workflow/run/replicate/spike observation with expected and off-target signal ratios.",
            },
            {
                "output": "relative_taxonomic_burden_by_workflow.tsv",
                "description": "Workflow-level summary of expected signal fraction, off-target fraction, and signal ratios.",
            },
            {
                "output": "relative_taxonomic_burden_by_spike.tsv",
                "description": "Spike-level summary of relative signal burden across workflows.",
            },
            {
                "output": "relative_taxonomic_burden_report.html",
                "description": "HTML report for quick review.",
            },
        ]
    )
    manifest.to_csv(out_dir / "relative_taxonomic_burden_manifest.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
