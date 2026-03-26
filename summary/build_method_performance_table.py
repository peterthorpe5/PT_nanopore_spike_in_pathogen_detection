#!/usr/bin/env python3
"""
Build method-performance summary tables from the spike-in combined long table.

This script converts the long-format benchmarking table produced by the
spike-in summariser into compact method-level performance tables suitable
for inclusion in the HTML report.

It treats each workflow / metric / target-label combination as a separate
benchmarking method and calculates confusion-matrix-style performance
statistics using a configurable detection threshold.

Default interpretation
----------------------
- Positive observations:
  non-shuffled rows with spike_n > 0
- Negative observations:
  shuffled-control rows or rows with spike_n == 0

Default detection rule
----------------------
Detection is called when the observed metric value is greater than or equal
to the method-specific maximum observed negative value, with a floor set by
--min_detect_value. This is intentionally conservative.

Outputs
-------
All tables are written as tab-separated files.

1. method_performance.tsv
   One row per workflow / metric / target-label method.
2. method_performance_by_spike.tsv
   Detection rates by spike level for each method.
3. method_performance.html
   Standalone HTML page for review.
4. method_performance_fragment.html
   HTML fragment that can be inserted into the main report.

Example
-------
python build_method_performance_table.py \
    --combined_long_tsv spikein_summary_report/combined_long.tsv \
    --out_dir spikein_summary_report
"""

from __future__ import annotations

import argparse
import html
import math
from pathlib import Path
from typing import Optional

import pandas as pd


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Build method-level performance summary tables from a combined "
            "long-format spike-in benchmarking table."
        )
    )
    parser.add_argument(
        "--combined_long_tsv",
        required=True,
        help="Input combined_long.tsv file.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for performance tables.",
    )
    parser.add_argument(
        "--threshold_mode",
        choices=[
            "fixed",
            "baseline_mean",
            "baseline_max",
            "baseline_mean_plus_sd",
        ],
        default="baseline_max",
        help=("Detection threshold strategy. Default: baseline_max"),
    )
    parser.add_argument(
        "--min_detect_value",
        type=float,
        default=1.0,
        help="Minimum metric value that can count as a detection. Default: 1.0",
    )
    parser.add_argument(
        "--sd_multiplier",
        type=float,
        default=2.0,
        help=(
            "Standard deviation multiplier for baseline_mean_plus_sd mode. "
            "Default: 2.0"
        ),
    )
    return parser.parse_args()


def get_first_existing_column(
    dataframe: pd.DataFrame,
    candidates: list[str],
) -> Optional[str]:
    """
    Return the first candidate column present in a DataFrame.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Input table.
    candidates : list[str]
        Candidate column names.

    Returns
    -------
    Optional[str]
        First matching column name, otherwise None.
    """
    for candidate in candidates:
        if candidate in dataframe.columns:
            return candidate
    return None


def normalise_combined_long(combined_long: pd.DataFrame) -> pd.DataFrame:
    """
    Normalise the combined long-format table to a standard schema.

    Parameters
    ----------
    combined_long : pd.DataFrame
        Raw combined long-format table.

    Returns
    -------
    pd.DataFrame
        Normalised table.
    """
    dataframe = combined_long.copy()

    rename_map: dict[str, str] = {}

    workflow_column = get_first_existing_column(
        dataframe=dataframe,
        candidates=["workflow", "workflow_type"],
    )
    metric_column = get_first_existing_column(
        dataframe=dataframe,
        candidates=["metric", "metric_name"],
    )
    value_column = get_first_existing_column(
        dataframe=dataframe,
        candidates=["value", "metric_value"],
    )
    spike_column = get_first_existing_column(
        dataframe=dataframe,
        candidates=["spike_n", "spike_n_per_genome"],
    )
    shuffled_column = get_first_existing_column(
        dataframe=dataframe,
        candidates=[
            "is_shuffled_control_recomputed",
            "is_shuffled_control",
        ],
    )

    if workflow_column is None:
        raise ValueError("Could not find workflow column in combined long table.")
    if metric_column is None:
        raise ValueError("Could not find metric column in combined long table.")
    if value_column is None:
        raise ValueError("Could not find value column in combined long table.")
    if spike_column is None:
        raise ValueError("Could not find spike column in combined long table.")

    if workflow_column != "workflow":
        rename_map[workflow_column] = "workflow"
    if metric_column != "metric":
        rename_map[metric_column] = "metric"
    if value_column != "value":
        rename_map[value_column] = "value"
    if spike_column != "spike_n":
        rename_map[spike_column] = "spike_n"
    if shuffled_column and shuffled_column != "is_shuffled_control":
        rename_map[shuffled_column] = "is_shuffled_control"

    if rename_map:
        dataframe = dataframe.rename(columns=rename_map)

    if "target_label" not in dataframe.columns:
        dataframe["target_label"] = "single_target"
    if "is_shuffled_control" not in dataframe.columns:
        dataframe["is_shuffled_control"] = False
    if "replicate" not in dataframe.columns:
        dataframe["replicate"] = pd.NA

    dataframe["workflow"] = dataframe["workflow"].fillna("unknown").astype(str)
    dataframe["metric"] = dataframe["metric"].fillna("unknown_metric").astype(str)
    dataframe["target_label"] = (
        dataframe["target_label"]
        .fillna("single_target")
        .replace("", "single_target")
        .astype(str)
    )
    dataframe["value"] = pd.to_numeric(dataframe["value"], errors="coerce")
    dataframe["spike_n"] = pd.to_numeric(dataframe["spike_n"], errors="coerce")
    dataframe["is_shuffled_control"] = (
        dataframe["is_shuffled_control"]
        .fillna(False)
        .astype(str)
        .str.lower()
        .isin(["true", "1", "yes"])
    )

    dataframe = dataframe.dropna(subset=["value", "spike_n"]).copy()
    dataframe["is_positive"] = (
        (~dataframe["is_shuffled_control"]) & (dataframe["spike_n"] > 0)
    )
    dataframe["is_negative"] = (
        dataframe["is_shuffled_control"] | (dataframe["spike_n"] == 0)
    )
    return dataframe


def safe_divide(
    numerator: float,
    denominator: float,
) -> float:
    """
    Divide safely, returning NaN for zero denominators.

    Parameters
    ----------
    numerator : float
        Numerator.
    denominator : float
        Denominator.

    Returns
    -------
    float
        Division result, or NaN.
    """
    if denominator == 0:
        return math.nan
    return numerator / denominator


def choose_threshold(
    negative_values: pd.Series,
    threshold_mode: str,
    min_detect_value: float,
    sd_multiplier: float,
) -> float:
    """
    Choose the detection threshold for one method.

    Parameters
    ----------
    negative_values : pd.Series
        Negative-control metric values.
    threshold_mode : str
        Threshold strategy.
    min_detect_value : float
        Minimum allowed threshold.
    sd_multiplier : float
        Standard deviation multiplier.

    Returns
    -------
    float
        Chosen threshold.
    """
    negative_values = pd.to_numeric(negative_values, errors="coerce").dropna()

    if negative_values.empty:
        return float(min_detect_value)

    if threshold_mode == "fixed":
        return float(min_detect_value)
    if threshold_mode == "baseline_mean":
        return float(max(min_detect_value, negative_values.mean()))
    if threshold_mode == "baseline_max":
        return float(max(min_detect_value, negative_values.max()))
    if threshold_mode == "baseline_mean_plus_sd":
        threshold = negative_values.mean() + (
            sd_multiplier * negative_values.std(ddof=0)
        )
        return float(max(min_detect_value, threshold))

    return float(min_detect_value)


def compute_detection_calls(
    dataframe: pd.DataFrame,
    threshold_mode: str,
    min_detect_value: float,
    sd_multiplier: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute binary detection calls for each method and observation.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Normalised combined long table.
    threshold_mode : str
        Threshold strategy.
    min_detect_value : float
        Minimum allowed threshold.
    sd_multiplier : float
        Standard deviation multiplier.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Observation-level table with detection calls and method-level threshold
        summary table.
    """
    detection_rows: list[pd.DataFrame] = []
    threshold_rows: list[dict[str, object]] = []

    grouped = dataframe.groupby(
        ["workflow", "metric", "target_label"],
        dropna=False,
    )

    for (workflow, metric, target_label), group in grouped:
        negatives = group.loc[group["is_negative"], "value"]
        threshold_value = choose_threshold(
            negative_values=negatives,
            threshold_mode=threshold_mode,
            min_detect_value=min_detect_value,
            sd_multiplier=sd_multiplier,
        )

        group = group.copy()
        group["threshold_value"] = threshold_value
        group["detected"] = group["value"] >= threshold_value
        detection_rows.append(group)

        threshold_rows.append(
            {
                "workflow": workflow,
                "metric": metric,
                "target_label": target_label,
                "threshold_mode": threshold_mode,
                "threshold_value": threshold_value,
                "negative_mean": (
                    float(negatives.mean()) if not negatives.dropna().empty else math.nan
                ),
                "negative_max": (
                    float(negatives.max()) if not negatives.dropna().empty else math.nan
                ),
                "negative_sd": (
                    float(negatives.std(ddof=0))
                    if not negatives.dropna().empty
                    else math.nan
                ),
                "n_negative_for_threshold": int(negatives.dropna().shape[0]),
            }
        )

    detection_df = (
        pd.concat(detection_rows, ignore_index=True)
        if detection_rows
        else pd.DataFrame()
    )
    threshold_df = pd.DataFrame(threshold_rows)
    return detection_df, threshold_df


def first_spike_meeting_rate(
    dataframe: pd.DataFrame,
    min_rate: float,
) -> float:
    """
    Return the first spike level with detection rate at or above a threshold.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Positive-only detection table for one method.
    min_rate : float
        Minimum detection rate.

    Returns
    -------
    float
        First spike level meeting the requested rate, or NaN.
    """
    if dataframe.empty:
        return math.nan

    grouped = (
        dataframe.groupby("spike_n", dropna=False)["detected"]
        .mean()
        .reset_index()
        .sort_values(by="spike_n")
    )
    detected = grouped.loc[grouped["detected"] >= min_rate, "spike_n"]
    if detected.empty:
        return math.nan
    return float(detected.iloc[0])


def summarise_method_performance(
    detection_df: pd.DataFrame,
    threshold_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Summarise per-method confusion metrics and derived statistics.

    Parameters
    ----------
    detection_df : pd.DataFrame
        Observation-level detection table.
    threshold_df : pd.DataFrame
        Method-level threshold table.

    Returns
    -------
    pd.DataFrame
        Method-level summary table.
    """
    rows: list[dict[str, object]] = []

    grouped = detection_df.groupby(
        ["workflow", "metric", "target_label"],
        dropna=False,
    )

    for (workflow, metric, target_label), group in grouped:
        positives = group.loc[group["is_positive"]].copy()
        negatives = group.loc[group["is_negative"]].copy()

        tp = int((positives["detected"] == True).sum())
        fn = int((positives["detected"] == False).sum())
        fp = int((negatives["detected"] == True).sum())
        tn = int((negatives["detected"] == False).sum())

        sensitivity = safe_divide(numerator=tp, denominator=tp + fn)
        recall = sensitivity
        true_positive_rate = sensitivity
        specificity = safe_divide(numerator=tn, denominator=tn + fp)
        false_positive_rate = safe_divide(numerator=fp, denominator=fp + tn)
        false_negative_rate = safe_divide(numerator=fn, denominator=fn + tp)
        precision = safe_divide(numerator=tp, denominator=tp + fp)
        negative_predictive_value = safe_divide(
            numerator=tn,
            denominator=tn + fn,
        )
        accuracy = safe_divide(
            numerator=tp + tn,
            denominator=tp + tn + fp + fn,
        )
        balanced_accuracy = (
            (sensitivity + specificity) / 2
            if not math.isnan(sensitivity) and not math.isnan(specificity)
            else math.nan
        )
        f1_score = safe_divide(
            numerator=2 * tp,
            denominator=(2 * tp) + fp + fn,
        )
        youden_j = (
            sensitivity + specificity - 1
            if not math.isnan(sensitivity) and not math.isnan(specificity)
            else math.nan
        )

        lod_50 = first_spike_meeting_rate(dataframe=positives, min_rate=0.50)
        lod_95 = first_spike_meeting_rate(dataframe=positives, min_rate=0.95)

        threshold_row = threshold_df.loc[
            (
                (threshold_df["workflow"] == workflow)
                & (threshold_df["metric"] == metric)
                & (threshold_df["target_label"] == target_label)
            )
        ]
        threshold_record = (
            threshold_row.iloc[0].to_dict()
            if not threshold_row.empty
            else {}
        )

        rows.append(
            {
                "workflow": workflow,
                "metric": metric,
                "target_label": target_label,
                "n_positive": int(positives.shape[0]),
                "n_negative": int(negatives.shape[0]),
                "tp": tp,
                "fn": fn,
                "fp": fp,
                "tn": tn,
                "threshold_mode": threshold_record.get("threshold_mode", pd.NA),
                "threshold_value": threshold_record.get("threshold_value", math.nan),
                "negative_mean": threshold_record.get("negative_mean", math.nan),
                "negative_max": threshold_record.get("negative_max", math.nan),
                "negative_sd": threshold_record.get("negative_sd", math.nan),
                "sensitivity": sensitivity,
                "recall": recall,
                "true_positive_rate": true_positive_rate,
                "specificity": specificity,
                "false_positive_rate": false_positive_rate,
                "false_negative_rate": false_negative_rate,
                "precision": precision,
                "negative_predictive_value": negative_predictive_value,
                "accuracy": accuracy,
                "balanced_accuracy": balanced_accuracy,
                "f1_score": f1_score,
                "youden_j": youden_j,
                "lod50_spike_n": lod_50,
                "lod95_spike_n": lod_95,
            }
        )

    out_df = pd.DataFrame(rows)
    if not out_df.empty:
        out_df = out_df.sort_values(
            by=["workflow", "metric", "target_label"]
        ).reset_index(drop=True)
    return out_df


def summarise_detection_by_spike(
    detection_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Summarise detection rate by spike level for each method.

    Parameters
    ----------
    detection_df : pd.DataFrame
        Observation-level detection table.

    Returns
    -------
    pd.DataFrame
        Per-spike detection-rate summary table.
    """
    positive_df = detection_df.loc[detection_df["is_positive"]].copy()
    if positive_df.empty:
        return pd.DataFrame(
            columns=[
                "workflow",
                "metric",
                "target_label",
                "spike_n",
                "n_observations",
                "n_detected",
                "detection_rate",
            ]
        )

    grouped = (
        positive_df.groupby(
            ["workflow", "metric", "target_label", "spike_n"],
            dropna=False,
        )
        .agg(
            n_observations=("detected", "size"),
            n_detected=("detected", "sum"),
            detection_rate=("detected", "mean"),
        )
        .reset_index()
        .sort_values(by=["workflow", "metric", "target_label", "spike_n"])
        .reset_index(drop=True)
    )
    return grouped


def format_value(value: object) -> str:
    """
    Format values for HTML display.

    Parameters
    ----------
    value : object
        Input value.

    Returns
    -------
    str
        Formatted string.
    """
    if pd.isna(value):
        return ""
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        if value.is_integer():
            return str(int(value))
        return f"{value:.3f}"
    return str(value)


def dataframe_to_html_table(
    dataframe: pd.DataFrame,
    table_class: str = "data-table",
) -> str:
    """
    Convert a DataFrame to an HTML table.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Input table.
    table_class : str, optional
        HTML table class name, by default "data-table".

    Returns
    -------
    str
        HTML table string.
    """
    if dataframe.empty:
        return '<p class="muted">No rows available.</p>'

    header_cells = "".join(
        f"<th>{html.escape(str(column))}</th>"
        for column in dataframe.columns
    )
    body_rows = []
    for _, row in dataframe.iterrows():
        cells = "".join(
            f"<td>{html.escape(format_value(value))}</td>"
            for value in row.tolist()
        )
        body_rows.append(f"<tr>{cells}</tr>")

    return (
        f"<table class=\"{html.escape(table_class)}\">"
        f"<thead><tr>{header_cells}</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody>"
        f"</table>"
    )


def build_html_page(
    method_df: pd.DataFrame,
    by_spike_df: pd.DataFrame,
) -> str:
    """
    Build a standalone HTML page for the method-performance summaries.

    Parameters
    ----------
    method_df : pd.DataFrame
        Method-level performance table.
    by_spike_df : pd.DataFrame
        Per-spike detection-rate table.

    Returns
    -------
    str
        HTML page.
    """
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
      h1, h2 {
        color: #1f4e79;
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
      }
      .data-table td {
        border-bottom: 1px solid #e6edf5;
        padding: 8px;
        vertical-align: top;
      }
      .data-table tbody tr:nth-child(even) {
        background: #fbfdff;
      }
      .muted {
        color: #666666;
      }
    </style>
    """

    method_table = dataframe_to_html_table(dataframe=method_df)
    by_spike_table = dataframe_to_html_table(dataframe=by_spike_df)

    return f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <title>Method performance summary</title>
  {style}
</head>
<body>
  <div class=\"container\">
    <h1>Method performance summary</h1>
    <p class=\"muted\">
      Sensitivity, recall, and true positive rate are equivalent here and are
      all reported for convenience. The table also includes specificity,
      false-positive rate, false-negative rate, precision, F1 score, balanced
      accuracy, and simple limit-of-detection summaries.
    </p>

    <h2>Method-level performance</h2>
    {method_table}

    <h2>Detection rate by spike level</h2>
    {by_spike_table}
  </div>
</body>
</html>
"""


def build_html_fragment(
    method_df: pd.DataFrame,
) -> str:
    """
    Build an HTML fragment suitable for insertion into the main report.

    Parameters
    ----------
    method_df : pd.DataFrame
        Method-level performance table.

    Returns
    -------
    str
        HTML fragment.
    """
    table_html = dataframe_to_html_table(dataframe=method_df)
    return f"""
<section class=\"report-section\">
  <h2>Method performance summary</h2>
  <p class=\"muted\">
    Sensitivity, recall, and true positive rate are equivalent here and are
    all reported for convenience. Thresholding is method-specific and based on
    the selected negative-control rule.
  </p>
  {table_html}
</section>
"""


def main() -> None:
    """
    Run the method-performance summarisation workflow.
    """
    args = parse_args()

    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    combined_long = pd.read_csv(args.combined_long_tsv, sep="\t")
    combined_long = normalise_combined_long(combined_long=combined_long)

    detection_df, threshold_df = compute_detection_calls(
        dataframe=combined_long,
        threshold_mode=args.threshold_mode,
        min_detect_value=args.min_detect_value,
        sd_multiplier=args.sd_multiplier,
    )

    method_df = summarise_method_performance(
        detection_df=detection_df,
        threshold_df=threshold_df,
    )
    by_spike_df = summarise_detection_by_spike(detection_df=detection_df)

    method_tsv = out_dir / "method_performance.tsv"
    by_spike_tsv = out_dir / "method_performance_by_spike.tsv"
    page_html = out_dir / "method_performance.html"
    fragment_html = out_dir / "method_performance_fragment.html"

    method_df.to_csv(method_tsv, sep="\t", index=False)
    by_spike_df.to_csv(by_spike_tsv, sep="\t", index=False)
    page_html.write_text(
        build_html_page(method_df=method_df, by_spike_df=by_spike_df),
        encoding="utf-8",
    )
    fragment_html.write_text(
        build_html_fragment(method_df=method_df),
        encoding="utf-8",
    )

    print(f"[INFO] Wrote: {method_tsv}")
    print(f"[INFO] Wrote: {by_spike_tsv}")
    print(f"[INFO] Wrote: {page_html}")
    print(f"[INFO] Wrote: {fragment_html}")


if __name__ == "__main__":
    main()
