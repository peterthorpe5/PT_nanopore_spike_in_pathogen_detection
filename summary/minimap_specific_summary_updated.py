#!/usr/bin/env python3
"""Summarise minimap-only spike-in runs with minimap-specific logic.

This script is designed for the minimap-only benchmarking stream where the
reference-name format differs between the focused or specific masked database
and the full or shared database. It derives target-centred and taxon-centred
summaries directly from minimap BED outputs instead of relying on the generic
summariser that can miss true target signal when reference names do not match
spike FASTA accession headers exactly.

The script can scan one or more run roots recursively and will work with both:

- newer runs that contain ``minimap_sorted.filtered.bed``
- older minimap-fix runs that contain ``minimap_sorted.MASKED.bed``

Outputs are written as tab-separated files:

- ``minimap_combined_long.tsv``
- ``minimap_taxon_summary.tsv``
- ``minimap_reported_taxa_long.tsv``
- ``minimap_tracked_target_performance.tsv``
- ``minimap_real_world_observation_level.tsv``
- ``minimap_real_world_run_summary.tsv``
- ``minimap_real_world_reference_summary.tsv``
- ``minimap_summary_run_info.tsv``
- ``minimap_input_root_manifest.tsv``
- ``minimap_summary_counts.tsv``

The file ``minimap_summary_run_info.tsv`` is a true run manifest with one row per
mix directory, including the inferred reference scope. This makes it easy to
check whether ``specific_db`` and ``full_db`` runs were actually discovered.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Iterable

import pandas as pd


BED_FILENAMES = (
    "minimap_sorted.filtered.bed",
    "minimap_sorted.MASKED.bed",
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Summarise minimap-only spike-in benchmarking runs using a "
            "minimap-specific parser."
        )
    )
    parser.add_argument(
        "--input_roots",
        nargs="+",
        required=True,
        help="One or more directories to scan recursively for minimap mix runs.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Directory in which output TSV files will be written.",
    )
    parser.add_argument(
        "--single_target_label",
        default="Plasmodium vivax",
        help="Fallback single-target label for single-genome minimap runs.",
    )
    parser.add_argument(
        "--panel2_tsv",
        default="",
        help="Optional panel 2 TSV used when old runs do not have run_metadata.tsv.",
    )
    parser.add_argument(
        "--panel3_tsv",
        default="",
        help="Optional panel 3 TSV used when old runs do not have run_metadata.tsv.",
    )
    parser.add_argument(
        "--target_threshold_alignments",
        type=int,
        default=1,
        help="Detection threshold for alignment-count metrics.",
    )
    parser.add_argument(
        "--target_threshold_unique_reads",
        type=int,
        default=1,
        help="Detection threshold for unique-read-count metrics.",
    )
    parser.add_argument(
        "--real_world_threshold_alignments",
        type=int,
        default=1,
        help="Taxon-presence threshold for real-world summaries using alignments.",
    )
    parser.add_argument(
        "--real_world_threshold_unique_reads",
        type=int,
        default=1,
        help="Taxon-presence threshold for real-world summaries using unique reads.",
    )
    return parser.parse_args()


def read_simple_kv_tsv(path: Path) -> dict[str, str]:
    """Read a two-column key-value TSV.

    Parameters
    ----------
    path : Path
        Path to the TSV file.

    Returns
    -------
    dict[str, str]
        Mapping of key to value.
    """
    mapping: dict[str, str] = {}
    if not path.exists():
        return mapping
    with path.open(mode="r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        return mapping
    start_index = 1 if rows[0][:2] == ["parameter", "value"] else 0
    for row in rows[start_index:]:
        if len(row) < 2:
            continue
        mapping[str(row[0]).strip()] = str(row[1]).strip()
    return mapping


def load_panel_labels(path: Path) -> list[str]:
    """Load target labels from a pathogen panel TSV.

    Parameters
    ----------
    path : Path
        Path to the panel TSV.

    Returns
    -------
    list[str]
        Target labels from the second column.
    """
    if not path.exists():
        return []
    labels: list[str] = []
    dataframe = pd.read_csv(path, sep="\t")
    if dataframe.shape[1] < 2:
        return labels
    for value in dataframe.iloc[:, 1].dropna().tolist():
        labels.append(str(value).strip())
    return labels


def infer_reference_scope(path: Path) -> str:
    """Infer whether a run used the full or specific minimap reference.

    Parameters
    ----------
    path : Path
        Path to a run or mix directory.

    Returns
    -------
    str
        Inferred reference scope label.
    """
    path_text = str(path).lower()
    if any(token in path_text for token in ["maskedref", "focused", "specific"]):
        return "specific_db"
    if any(token in path_text for token in ["fullref", "shared"]):
        return "full_db"
    return "unknown_db"


def infer_workflow(run_name: str) -> str:
    """Infer workflow label from the run name.

    Parameters
    ----------
    run_name : str
        Run name.

    Returns
    -------
    str
        Workflow label.
    """
    name = run_name.lower()
    if "single" in name:
        return "single_read"
    if "multi" in name:
        return "multi_read"
    return "read_level"


def infer_expected_taxa(
    *,
    run_root: Path,
    run_name: str,
    metadata: dict[str, str],
    single_target_label: str,
    panel2_tsv: Path | None,
    panel3_tsv: Path | None,
) -> list[str]:
    """Infer the expected taxa for a run.

    Parameters
    ----------
    run_root : Path
        Run root directory.
    run_name : str
        Human-readable run name.
    metadata : dict[str, str]
        Parsed run metadata.
    single_target_label : str
        Fallback single-target label.
    panel2_tsv : Path | None
        Optional panel 2 TSV.
    panel3_tsv : Path | None
        Optional panel 3 TSV.

    Returns
    -------
    list[str]
        Expected target taxa for the run.
    """
    if "target_label" in metadata and metadata["target_label"]:
        return [metadata["target_label"]]
    if "pathogen_config_tsv" in metadata and metadata["pathogen_config_tsv"]:
        return load_panel_labels(Path(metadata["pathogen_config_tsv"]))

    lower_name = run_name.lower()
    if "panel2" in lower_name and panel2_tsv is not None and panel2_tsv.exists():
        return load_panel_labels(panel2_tsv)
    if "panel3" in lower_name and panel3_tsv is not None and panel3_tsv.exists():
        return load_panel_labels(panel3_tsv)
    if "single" in lower_name or "shuffle" in lower_name:
        return [single_target_label]

    if run_root.name == "main_run":
        parent_name = run_root.parent.name.lower()
        if "panel2" in parent_name and panel2_tsv is not None and panel2_tsv.exists():
            return load_panel_labels(panel2_tsv)
        if "panel3" in parent_name and panel3_tsv is not None and panel3_tsv.exists():
            return load_panel_labels(panel3_tsv)
        return [single_target_label]

    return []


def infer_is_shuffled(run_name: str, metadata: dict[str, str], run_root: Path) -> bool:
    """Infer whether a run is a shuffled control.

    Parameters
    ----------
    run_name : str
        Human-readable run name.
    metadata : dict[str, str]
        Parsed run metadata.
    run_root : Path
        Run root directory.

    Returns
    -------
    bool
        True if the run is a shuffled control, otherwise False.
    """
    if metadata.get("is_shuffled_control", "").strip().lower() == "true":
        return True
    path_text = str(run_root).lower()
    return "shuffle" in path_text or "shuffled" in path_text


def parse_mix_dir_name(mix_dir_name: str) -> tuple[int, int]:
    """Parse replicate and spike level from a mix directory name.

    Parameters
    ----------
    mix_dir_name : str
        Mix directory name.

    Returns
    -------
    tuple[int, int]
        Replicate number and spike level.
    """
    match = re.match(r"mix_rep(?P<rep>\d+)_n(?P<spike>\d+)$", mix_dir_name)
    if not match:
        raise ValueError(f"Could not parse replicate and spike from {mix_dir_name}")
    return int(match.group("rep")), int(match.group("spike"))


def canonicalise_taxon_name(raw_taxon: str) -> str:
    """Normalise a taxon name for downstream matching.

    Parameters
    ----------
    raw_taxon : str
        Raw taxon string.

    Returns
    -------
    str
        Canonical taxon string.
    """
    taxon = re.sub(r"\s+", " ", raw_taxon.replace("_", " ").strip())
    return taxon


def parse_taxon_from_reference(reference_name: str) -> str:
    """Parse a taxon label from a minimap reference name.

    Parameters
    ----------
    reference_name : str
        Reference name from the BED ``chrom`` field.

    Returns
    -------
    str
        Parsed taxon label.
    """
    parts = reference_name.split("_")
    if len(parts) >= 3 and parts[0].lower() == "plas":
        return canonicalise_taxon_name(" ".join(parts[2:]))
    if reference_name.startswith("Plasmodium_"):
        return canonicalise_taxon_name(reference_name)
    return canonicalise_taxon_name(reference_name)


def taxon_matches_expected(taxon_name: str, expected_taxon: str) -> bool:
    """Test whether a parsed taxon matches an expected target label.

    Parameters
    ----------
    taxon_name : str
        Parsed taxon name from the BED reference field.
    expected_taxon : str
        Expected target label.

    Returns
    -------
    bool
        True if the taxon should be considered a match for the expected label.
    """
    left = canonicalise_taxon_name(taxon_name).lower()
    right = canonicalise_taxon_name(expected_taxon).lower()
    return left == right or left.startswith(f"{right} ")


def choose_bed_path(mix_dir: Path) -> Path | None:
    """Choose the BED file present in a mix directory.

    Parameters
    ----------
    mix_dir : Path
        Mix directory.

    Returns
    -------
    Path | None
        BED path if present, otherwise None.
    """
    for filename in BED_FILENAMES:
        candidate = mix_dir / filename
        if candidate.exists() and candidate.stat().st_size > 0:
            return candidate
    return None


def parse_bed_counts(
    bed_path: Path,
) -> tuple[dict[str, int], dict[str, int], dict[str, set[str]], int, int]:
    """Parse taxon counts from a minimap BED file.

    Parameters
    ----------
    bed_path : Path
        Path to the BED file.

    Returns
    -------
    tuple[dict[str, int], dict[str, int], dict[str, set[str]], int, int]
        Reference-level alignment counts, taxon-level alignment counts,
        taxon-level unique read sets, total alignments, and total unique reads.
    """
    reference_counts: dict[str, int] = defaultdict(int)
    taxon_alignments: dict[str, int] = defaultdict(int)
    taxon_reads: dict[str, set[str]] = defaultdict(set)
    all_reads: set[str] = set()
    total_alignments = 0

    with bed_path.open(mode="r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            reference_name = fields[0]
            read_name = fields[3]
            taxon_name = parse_taxon_from_reference(reference_name)
            reference_counts[reference_name] += 1
            taxon_alignments[taxon_name] += 1
            taxon_reads[taxon_name].add(read_name)
            all_reads.add(read_name)
            total_alignments += 1

    return (
        dict(reference_counts),
        dict(taxon_alignments),
        dict(taxon_reads),
        total_alignments,
        len(all_reads),
    )


def first_spike_meeting_rate(dataframe: pd.DataFrame, min_rate: float) -> float:
    """Return the first spike level meeting a minimum detection rate.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Data frame containing ``spike_n`` and ``detected`` columns.
    min_rate : float
        Minimum detection rate.

    Returns
    -------
    float
        First spike level meeting the detection rate, or NaN if none do.
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
        Division result or NaN if the denominator is zero.
    """
    if denominator == 0:
        return math.nan
    return float(numerator / denominator)


def discover_mix_dirs(input_roots: Iterable[Path]) -> list[dict[str, object]]:
    """Discover mix directories containing minimap BED outputs.

    Parameters
    ----------
    input_roots : Iterable[Path]
        Input roots to scan.

    Returns
    -------
    list[dict[str, object]]
        Sorted unique mix-directory records containing the input root and BED.
    """
    records: dict[Path, dict[str, object]] = {}
    for input_root in input_roots:
        for bed_filename in BED_FILENAMES:
            for bed_path in input_root.rglob(bed_filename):
                mix_dir = bed_path.parent
                if mix_dir not in records:
                    records[mix_dir] = {
                        "mix_dir": mix_dir,
                        "input_root": input_root,
                        "bed_path": bed_path,
                    }
    return [records[path] for path in sorted(records)]


def build_long_tables(
    *,
    mix_records: list[dict[str, object]],
    single_target_label: str,
    panel2_tsv: Path | None,
    panel3_tsv: Path | None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build minimap target-centric, taxon-centric, and manifest tables.

    Parameters
    ----------
    mix_records : list[dict[str, object]]
        Mix-directory records.
    single_target_label : str
        Fallback single-target label.
    panel2_tsv : Path | None
        Optional panel 2 TSV.
    panel3_tsv : Path | None
        Optional panel 3 TSV.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        Combined long data frame, taxon summary data frame, and run manifest.
    """
    combined_rows: list[dict[str, object]] = []
    taxa_rows: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []

    for mix_record in mix_records:
        mix_dir = Path(mix_record["mix_dir"])
        input_root = Path(mix_record["input_root"])
        bed_path = Path(mix_record["bed_path"])

        run_root = mix_dir.parent
        raw_run_name = run_root.name
        run_name = run_root.parent.name if raw_run_name == "main_run" else raw_run_name
        metadata = read_simple_kv_tsv(run_root / "run_metadata.tsv")
        expected_taxa = infer_expected_taxa(
            run_root=run_root,
            run_name=run_name,
            metadata=metadata,
            single_target_label=single_target_label,
            panel2_tsv=panel2_tsv,
            panel3_tsv=panel3_tsv,
        )
        reference_scope = infer_reference_scope(run_root)
        workflow = infer_workflow(run_name)
        is_shuffled_control = infer_is_shuffled(
            run_name=run_name,
            metadata=metadata,
            run_root=run_root,
        )
        replicate, spike_n = parse_mix_dir_name(mix_dir.name)
        total_spiked_reads = spike_n * max(len(expected_taxa), 1)
        is_negative = bool(is_shuffled_control or spike_n == 0)

        (
            reference_counts,
            taxon_alignments,
            taxon_reads,
            total_alignments,
            total_unique_reads,
        ) = parse_bed_counts(bed_path=bed_path)

        combined_base = {
            "workflow": workflow,
            "workflow_type": workflow,
            "classifier": "minimap",
            "run_name": run_name,
            "run_path": str(run_root),
            "run_dir": str(run_root),
            "summary_path": str(bed_path),
            "mix_dir": mix_dir.name,
            "input_root": str(input_root),
            "reference_scope": reference_scope,
            "replicate": replicate,
            "spike_n": spike_n,
            "spike_n_per_genome": spike_n,
            "total_spiked_reads": total_spiked_reads,
            "total_spike_n": total_spiked_reads,
            "is_shuffled_control": int(is_shuffled_control),
            "is_negative": bool(is_negative),
            "n_expected_taxa": len(expected_taxa),
            "n_genomes": len(expected_taxa),
            "expected_taxa": "; ".join(expected_taxa),
            "bed_path": str(bed_path),
            "db_format": bed_path.name,
        }

        manifest_rows.append(
            {
                "input_root": str(input_root),
                "run_name": run_name,
                "run_dir": str(run_root),
                "mix_dir": mix_dir.name,
                "bed_path": str(bed_path),
                "reference_scope": reference_scope,
                "workflow": workflow,
                "classifier": "minimap",
                "replicate": replicate,
                "spike_n": spike_n,
                "spike_n_per_genome": spike_n,
                "total_spike_n": total_spiked_reads,
                "is_shuffled_control": int(is_shuffled_control),
                "is_negative": bool(is_negative),
                "n_expected_taxa": len(expected_taxa),
                "n_genomes": len(expected_taxa),
                "expected_taxa": "; ".join(expected_taxa),
                "db_format": bed_path.name,
                "total_alignments": total_alignments,
                "total_unique_reads": total_unique_reads,
                "n_taxa_reported": len(taxon_alignments),
            }
        )

        for metric_name, metric_value, target_label in [
            ("minimap_total_alignments", total_alignments, "all_taxa"),
            ("minimap_total_unique_reads", total_unique_reads, "all_taxa"),
        ]:
            combined_rows.append(
                {
                    **combined_base,
                    "metric": metric_name,
                    "metric_name": metric_name,
                    "target_label": target_label,
                    "value": metric_value,
                    "metric_value": metric_value,
                }
            )

        total_target_alignments = 0
        total_target_read_ids: set[str] = set()

        for expected_taxon in expected_taxa:
            expected_alignments = 0
            expected_read_ids: set[str] = set()
            matched_taxa: list[str] = []
            for taxon_name, alignments in taxon_alignments.items():
                if taxon_matches_expected(
                    taxon_name=taxon_name,
                    expected_taxon=expected_taxon,
                ):
                    expected_alignments += alignments
                    expected_read_ids.update(taxon_reads.get(taxon_name, set()))
                    matched_taxa.append(taxon_name)
            matched_taxa = sorted(set(matched_taxa))
            total_target_alignments += expected_alignments
            total_target_read_ids.update(expected_read_ids)

            for metric_name, metric_value in [
                ("minimap_target_alignments", expected_alignments),
                ("minimap_target_unique_reads", len(expected_read_ids)),
            ]:
                combined_rows.append(
                    {
                        **combined_base,
                        "metric": metric_name,
                        "metric_name": metric_name,
                        "target_label": expected_taxon,
                        "value": metric_value,
                        "metric_value": metric_value,
                        "matched_taxa": "; ".join(matched_taxa),
                    }
                )

        if expected_taxa:
            for metric_name, metric_value, target_label in [
                ("minimap_target_alignments", total_target_alignments, "total"),
                (
                    "minimap_target_unique_reads",
                    len(total_target_read_ids),
                    "total unique",
                ),
            ]:
                combined_rows.append(
                    {
                        **combined_base,
                        "metric": metric_name,
                        "metric_name": metric_name,
                        "target_label": target_label,
                        "value": metric_value,
                        "metric_value": metric_value,
                    }
                )

        for taxon_name, alignments in sorted(
            taxon_alignments.items(),
            key=lambda item: (-item[1], item[0]),
        ):
            unique_reads = len(taxon_reads.get(taxon_name, set()))
            matched_expected_labels = [
                expected_taxon
                for expected_taxon in expected_taxa
                if taxon_matches_expected(
                    taxon_name=taxon_name,
                    expected_taxon=expected_taxon,
                )
            ]
            reference_count = sum(
                1
                for reference_name in reference_counts
                if parse_taxon_from_reference(reference_name) == taxon_name
            )
            taxa_rows.append(
                {
                    **combined_base,
                    "taxon_name": taxon_name,
                    "reference_count": reference_count,
                    "alignments": alignments,
                    "unique_reads": unique_reads,
                    "matched_target_labels": "; ".join(matched_expected_labels),
                    "is_expected_taxon": int(bool(matched_expected_labels)),
                    "is_off_target_taxon": int(
                        bool(expected_taxa) and not bool(matched_expected_labels)
                    ),
                }
            )

    combined_df = pd.DataFrame(combined_rows)
    taxa_df = pd.DataFrame(taxa_rows)
    manifest_df = pd.DataFrame(manifest_rows)
    return combined_df, taxa_df, manifest_df


def compute_tracked_target_performance(
    combined_df: pd.DataFrame,
    threshold_map: dict[str, float],
) -> pd.DataFrame:
    """Compute tracked-target performance metrics.

    Parameters
    ----------
    combined_df : pd.DataFrame
        Target-centric long table.
    threshold_map : dict[str, float]
        Threshold values keyed by metric.

    Returns
    -------
    pd.DataFrame
        Tracked-target performance table.
    """
    subset = combined_df.loc[
        combined_df["metric"].isin(
            ["minimap_target_alignments", "minimap_target_unique_reads"]
        )
    ].copy()

    rows: list[dict[str, object]] = []
    grouped = subset.groupby(
        ["reference_scope", "workflow", "run_name", "metric", "target_label"],
        dropna=False,
    )

    for group_key, group in grouped:
        (
            reference_scope,
            workflow,
            run_name,
            metric,
            target_label,
        ) = group_key
        threshold_value = float(threshold_map.get(metric, 1))
        group = group.copy()
        group["detected"] = group["value"] >= threshold_value

        positives = group.loc[~group["is_negative"]].copy()
        negatives = group.loc[group["is_negative"]].copy()

        tp = int(positives["detected"].sum())
        fn = int((~positives["detected"]).sum())
        fp = int(negatives["detected"].sum())
        tn = int((~negatives["detected"]).sum())

        sensitivity = safe_divide(numerator=tp, denominator=(tp + fn))
        specificity = safe_divide(numerator=tn, denominator=(tn + fp))
        precision = safe_divide(numerator=tp, denominator=(tp + fp))
        negative_predictive_value = safe_divide(
            numerator=tn,
            denominator=(tn + fn),
        )
        accuracy = safe_divide(
            numerator=(tp + tn),
            denominator=(tp + tn + fp + fn),
        )
        if math.isnan(sensitivity) or math.isnan(specificity):
            balanced_accuracy = math.nan
        else:
            balanced_accuracy = float((sensitivity + specificity) / 2.0)
        false_positive_rate = safe_divide(numerator=fp, denominator=(fp + tn))
        false_negative_rate = safe_divide(numerator=fn, denominator=(fn + tp))
        f1_score = safe_divide(
            numerator=(2 * tp),
            denominator=((2 * tp) + fp + fn),
        )
        if math.isnan(sensitivity) or math.isnan(specificity):
            youden_j = math.nan
        else:
            youden_j = float(sensitivity + specificity - 1.0)

        lod50 = first_spike_meeting_rate(dataframe=positives, min_rate=0.50)
        lod95 = first_spike_meeting_rate(dataframe=positives, min_rate=0.95)
        lod100 = first_spike_meeting_rate(dataframe=positives, min_rate=1.00)

        rows.append(
            {
                "reference_scope": reference_scope,
                "workflow": workflow,
                "workflow_type": workflow,
                "run_name": run_name,
                "metric": metric,
                "metric_name": metric,
                "target_label": target_label,
                "threshold_value": threshold_value,
                "n_positive": int(positives.shape[0]),
                "n_negative": int(negatives.shape[0]),
                "tp": tp,
                "fn": fn,
                "fp": fp,
                "tn": tn,
                "sensitivity": sensitivity,
                "specificity": specificity,
                "precision": precision,
                "negative_predictive_value": negative_predictive_value,
                "accuracy": accuracy,
                "balanced_accuracy": balanced_accuracy,
                "false_positive_rate": false_positive_rate,
                "false_negative_rate": false_negative_rate,
                "f1_score": f1_score,
                "youden_j": youden_j,
                "lod50_spike_n": lod50,
                "lod95_spike_n": lod95,
                "lod100_spike_n": lod100,
            }
        )

    return pd.DataFrame(rows)


def build_observation_level_taxon_presence(
    taxa_df: pd.DataFrame,
    *,
    presence_metric: str,
    threshold_value: int,
) -> pd.DataFrame:
    """Collapse taxon rows to observation-level real-world flags.

    Parameters
    ----------
    taxa_df : pd.DataFrame
        Taxon-centric long table.
    presence_metric : str
        Either ``alignments`` or ``unique_reads``.
    threshold_value : int
        Threshold above which a taxon is considered present.

    Returns
    -------
    pd.DataFrame
        Observation-level taxonomic burden table.
    """
    rows: list[dict[str, object]] = []
    grouping_columns = [
        "reference_scope",
        "workflow",
        "workflow_type",
        "classifier",
        "run_name",
        "run_path",
        "run_dir",
        "mix_dir",
        "replicate",
        "spike_n",
        "spike_n_per_genome",
        "total_spiked_reads",
        "total_spike_n",
        "is_shuffled_control",
        "is_negative",
        "expected_taxa",
        "n_expected_taxa",
        "n_genomes",
    ]

    for group_key, group in taxa_df.groupby(grouping_columns, dropna=False):
        key_values = dict(
            zip(
                grouping_columns,
                group_key if isinstance(group_key, tuple) else (group_key,),
            )
        )
        positives = group.loc[~group["is_negative"]].copy()
        negatives = group.loc[group["is_negative"]].copy()

        present = group.loc[group[presence_metric] >= threshold_value].copy()
        expected_present = present.loc[present["is_expected_taxon"] == 1].copy()
        off_target_present = present.loc[present["is_off_target_taxon"] == 1].copy()

        expected_labels_present = sorted(
            {
                label.strip()
                for value in expected_present["matched_target_labels"].dropna().tolist()
                for label in str(value).split(";")
                if label.strip()
            }
        )

        n_reported_taxa = int(present.shape[0])
        n_off_target_taxa = int(off_target_present.shape[0])
        n_expected_taxa_present = len(expected_labels_present)
        any_expected = int(n_expected_taxa_present > 0)
        all_expected = int(
            n_expected_taxa_present == int(key_values["n_expected_taxa"])
            and int(key_values["n_expected_taxa"]) > 0
        )
        clean_positive = int(bool(all_expected) and n_off_target_taxa == 0)
        any_taxon_negative = int(bool(key_values["is_negative"]) and n_reported_taxa > 0)
        off_target_negative = int(
            bool(key_values["is_negative"]) and n_off_target_taxa > 0
        )
        clean_negative = int(bool(key_values["is_negative"]) and n_reported_taxa == 0)
        positive_off_target = int((not bool(key_values["is_negative"])) and n_off_target_taxa > 0)
        taxonomic_precision = safe_divide(
            numerator=n_expected_taxa_present,
            denominator=n_reported_taxa,
        )

        rows.append(
            {
                **key_values,
                "presence_metric": presence_metric,
                "threshold_value": threshold_value,
                "n_reported_taxa": n_reported_taxa,
                "n_expected_taxa_present": n_expected_taxa_present,
                "n_off_target_taxa": n_off_target_taxa,
                "expected_labels_present": "; ".join(expected_labels_present),
                "any_expected": any_expected,
                "all_expected": all_expected,
                "clean_positive": clean_positive,
                "positive_off_target": positive_off_target,
                "clean_negative": clean_negative,
                "any_taxon_negative": any_taxon_negative,
                "off_target_negative": off_target_negative,
                "taxonomic_precision": taxonomic_precision,
                "n_positive_observations": int(positives.shape[0]),
                "n_negative_observations": int(negatives.shape[0]),
            }
        )

    return pd.DataFrame(rows)


def summarise_real_world(
    observation_df: pd.DataFrame,
    grouping_columns: list[str],
) -> pd.DataFrame:
    """Summarise observation-level real-world metrics.

    Parameters
    ----------
    observation_df : pd.DataFrame
        Observation-level real-world table.
    grouping_columns : list[str]
        Grouping columns to use.

    Returns
    -------
    pd.DataFrame
        Real-world summary table.
    """
    rows: list[dict[str, object]] = []

    for group_key, group in observation_df.groupby(grouping_columns, dropna=False):
        key_values = dict(
            zip(
                grouping_columns,
                group_key if isinstance(group_key, tuple) else (group_key,),
            )
        )
        positives = group.loc[~group["is_negative"]].copy()
        negatives = group.loc[group["is_negative"]].copy()

        n_positive_observations = int(positives.shape[0])
        n_negative_observations = int(negatives.shape[0])
        n_observations = int(group.shape[0])

        n_positive_any_expected = int(positives["any_expected"].sum())
        n_positive_all_expected = int(positives["all_expected"].sum())
        n_positive_clean = int(positives["clean_positive"].sum())
        n_positive_with_off_target = int(positives["positive_off_target"].sum())

        n_negative_clean = int(negatives["clean_negative"].sum())
        n_negative_with_any_taxon = int(negatives["any_taxon_negative"].sum())
        n_negative_with_off_target = int(negatives["off_target_negative"].sum())

        any_expected_sensitivity = safe_divide(
            numerator=n_positive_any_expected,
            denominator=n_positive_observations,
        )
        all_expected_sensitivity = safe_divide(
            numerator=n_positive_all_expected,
            denominator=n_positive_observations,
        )
        clean_sensitivity = safe_divide(
            numerator=n_positive_clean,
            denominator=n_positive_observations,
        )
        positive_off_target_rate = safe_divide(
            numerator=n_positive_with_off_target,
            denominator=n_positive_observations,
        )
        taxonomic_specificity = safe_divide(
            numerator=n_negative_clean,
            denominator=n_negative_observations,
        )
        taxonomic_false_positive_rate = safe_divide(
            numerator=n_negative_with_any_taxon,
            denominator=n_negative_observations,
        )
        off_target_negative_rate = safe_divide(
            numerator=n_negative_with_off_target,
            denominator=n_negative_observations,
        )
        clean_precision = safe_divide(
            numerator=n_positive_clean,
            denominator=(n_positive_clean + n_negative_with_any_taxon),
        )
        if math.isnan(clean_precision) or math.isnan(clean_sensitivity):
            clean_f1 = math.nan
        else:
            denominator = clean_precision + clean_sensitivity
            clean_f1 = (
                math.nan
                if denominator == 0
                else float(2 * clean_precision * clean_sensitivity / denominator)
            )

        mean_reported_taxa_positive = (
            float(positives["n_reported_taxa"].mean())
            if not positives.empty
            else math.nan
        )
        mean_off_target_taxa_positive = (
            float(positives["n_off_target_taxa"].mean())
            if not positives.empty
            else math.nan
        )
        mean_taxonomic_precision_positive = (
            float(positives["taxonomic_precision"].mean())
            if not positives.empty
            else math.nan
        )
        median_off_target_taxa_positive = (
            float(median(positives["n_off_target_taxa"].tolist()))
            if not positives.empty
            else math.nan
        )

        row = {
            **key_values,
            "n_observations": n_observations,
            "n_positive_observations": n_positive_observations,
            "n_negative_observations": n_negative_observations,
            "n_positive_any_expected": n_positive_any_expected,
            "n_positive_all_expected": n_positive_all_expected,
            "n_positive_clean": n_positive_clean,
            "n_positive_with_off_target": n_positive_with_off_target,
            "n_negative_clean": n_negative_clean,
            "n_negative_with_any_taxon": n_negative_with_any_taxon,
            "n_negative_with_off_target": n_negative_with_off_target,
            "any_expected_sensitivity": any_expected_sensitivity,
            "all_expected_sensitivity": all_expected_sensitivity,
            "clean_sensitivity": clean_sensitivity,
            "positive_off_target_rate": positive_off_target_rate,
            "taxonomic_specificity": taxonomic_specificity,
            "taxonomic_false_positive_rate": taxonomic_false_positive_rate,
            "off_target_negative_rate": off_target_negative_rate,
            "clean_precision": clean_precision,
            "clean_f1": clean_f1,
            "mean_reported_taxa_positive": mean_reported_taxa_positive,
            "mean_off_target_taxa_positive": mean_off_target_taxa_positive,
            "mean_taxonomic_precision_positive": mean_taxonomic_precision_positive,
            "median_off_target_taxa_positive": median_off_target_taxa_positive,
        }
        rows.append(row)

    return pd.DataFrame(rows)


def build_reported_taxa_long(taxa_df: pd.DataFrame) -> pd.DataFrame:
    """Expand compact minimap taxon rows into long reported-taxa rows.

    Parameters
    ----------
    taxa_df : pd.DataFrame
        Compact taxon summary table.

    Returns
    -------
    pd.DataFrame
        Long-form reported-taxa table compatible with the broader pipeline.
    """
    rows: list[dict[str, object]] = []
    for _, row in taxa_df.iterrows():
        common = {
            "workflow": row.get("workflow"),
            "workflow_type": row.get("workflow_type"),
            "classifier": row.get("classifier", "minimap"),
            "run_name": row.get("run_name"),
            "run_dir": row.get("run_dir"),
            "summary_path": row.get("summary_path"),
            "is_shuffled_control": row.get("is_shuffled_control"),
            "replicate": row.get("replicate"),
            "spike_n_per_genome": row.get("spike_n_per_genome"),
            "total_spike_n": row.get("total_spike_n"),
            "n_genomes": row.get("n_genomes"),
            "assembly_n_contigs": pd.NA,
            "assembly_total_bases": pd.NA,
            "reported_taxa_source_tsv": row.get("bed_path"),
            "matched_from_column": pd.NA,
            "report_tsv": row.get("bed_path"),
            "filter_prefix": pd.NA,
            "matched_target_labels": row.get("matched_target_labels", ""),
            "taxid": pd.NA,
            "taxon_name": row.get("taxon_name"),
            "rank_code": pd.NA,
            "percent": pd.NA,
            "reference_scope": row.get("reference_scope"),
            "input_root": row.get("input_root"),
        }
        rows.append(
            {
                **common,
                "metric_name": "minimap_reported_clade_reads",
                "metric_value": row.get("alignments"),
            }
        )
        rows.append(
            {
                **common,
                "metric_name": "minimap_reported_unique_reads",
                "metric_value": row.get("unique_reads"),
            }
        )
    return pd.DataFrame(rows)


def build_input_root_manifest(
    input_roots: list[Path],
    manifest_df: pd.DataFrame,
) -> pd.DataFrame:
    """Build one row per input root for debugging discovery.

    Parameters
    ----------
    input_roots : list[Path]
        Input roots supplied by the user.
    manifest_df : pd.DataFrame
        Run manifest with one row per mix directory.

    Returns
    -------
    pd.DataFrame
        Input-root manifest.
    """
    rows: list[dict[str, object]] = []
    for input_root in input_roots:
        root_str = str(input_root)
        if manifest_df.empty:
            subset = pd.DataFrame()
        else:
            subset = manifest_df.loc[manifest_df["input_root"] == root_str].copy()
        rows.append(
            {
                "input_root": root_str,
                "exists": int(input_root.exists()),
                "n_mix_dirs_found": int(subset.shape[0]),
                "n_specific_db_mix_dirs": int(
                    subset.loc[subset["reference_scope"] == "specific_db"].shape[0]
                ),
                "n_full_db_mix_dirs": int(
                    subset.loc[subset["reference_scope"] == "full_db"].shape[0]
                ),
                "n_unknown_db_mix_dirs": int(
                    subset.loc[subset["reference_scope"] == "unknown_db"].shape[0]
                ),
                "workflows_present": "; ".join(
                    sorted(pd.unique(subset.get("workflow", pd.Series(dtype=str)).astype(str)))
                ) if not subset.empty else "",
                "reference_scopes_present": "; ".join(
                    sorted(
                        pd.unique(
                            subset.get("reference_scope", pd.Series(dtype=str)).astype(str)
                        )
                    )
                ) if not subset.empty else "",
            }
        )
    return pd.DataFrame(rows)


def write_tsv(dataframe: pd.DataFrame, path: Path) -> None:
    """Write a data frame as TSV.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Data frame to write.
    path : Path
        Output path.
    """
    dataframe.to_csv(path, sep="\t", index=False)


def main() -> None:
    """Run the minimap-specific summarisation workflow."""
    args = parse_args()

    input_roots = [Path(path).expanduser().resolve() for path in args.input_roots]
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    panel2_tsv = Path(args.panel2_tsv).expanduser().resolve() if args.panel2_tsv else None
    panel3_tsv = Path(args.panel3_tsv).expanduser().resolve() if args.panel3_tsv else None

    mix_records = discover_mix_dirs(input_roots=input_roots)
    if not mix_records:
        raise SystemExit("No minimap BED files were found under the input roots.")

    combined_df, taxa_df, manifest_df = build_long_tables(
        mix_records=mix_records,
        single_target_label=args.single_target_label,
        panel2_tsv=panel2_tsv,
        panel3_tsv=panel3_tsv,
    )
    reported_taxa_long_df = build_reported_taxa_long(taxa_df=taxa_df)

    threshold_map = {
        "minimap_target_alignments": float(args.target_threshold_alignments),
        "minimap_target_unique_reads": float(args.target_threshold_unique_reads),
    }
    tracked_df = compute_tracked_target_performance(
        combined_df=combined_df,
        threshold_map=threshold_map,
    )

    observation_rows: list[pd.DataFrame] = []
    real_world_threshold_map = {
        "alignments": int(args.real_world_threshold_alignments),
        "unique_reads": int(args.real_world_threshold_unique_reads),
    }
    for presence_metric, threshold_value in real_world_threshold_map.items():
        observation_rows.append(
            build_observation_level_taxon_presence(
                taxa_df=taxa_df,
                presence_metric=presence_metric,
                threshold_value=threshold_value,
            )
        )
    observation_df = pd.concat(observation_rows, ignore_index=True)

    real_world_run_df = summarise_real_world(
        observation_df=observation_df,
        grouping_columns=[
            "reference_scope",
            "workflow",
            "workflow_type",
            "classifier",
            "run_name",
            "presence_metric",
            "threshold_value",
        ],
    )
    real_world_reference_df = summarise_real_world(
        observation_df=observation_df,
        grouping_columns=[
            "reference_scope",
            "workflow",
            "workflow_type",
            "classifier",
            "presence_metric",
            "threshold_value",
        ],
    )

    input_root_manifest_df = build_input_root_manifest(
        input_roots=input_roots,
        manifest_df=manifest_df,
    )
    summary_counts_df = pd.DataFrame(
        {
            "metric": [
                "input_roots_supplied",
                "mix_dirs_found",
                "combined_long_rows",
                "taxon_summary_rows",
                "reported_taxa_long_rows",
                "tracked_target_rows",
                "real_world_observation_rows",
                "real_world_run_rows",
                "real_world_reference_rows",
            ],
            "value": [
                len(input_roots),
                int(manifest_df.shape[0]),
                int(combined_df.shape[0]),
                int(taxa_df.shape[0]),
                int(reported_taxa_long_df.shape[0]),
                int(tracked_df.shape[0]),
                int(observation_df.shape[0]),
                int(real_world_run_df.shape[0]),
                int(real_world_reference_df.shape[0]),
            ],
        }
    )

    write_tsv(dataframe=combined_df, path=out_dir / "minimap_combined_long.tsv")
    write_tsv(dataframe=taxa_df, path=out_dir / "minimap_taxon_summary.tsv")
    write_tsv(
        dataframe=reported_taxa_long_df,
        path=out_dir / "minimap_reported_taxa_long.tsv",
    )
    write_tsv(
        dataframe=tracked_df,
        path=out_dir / "minimap_tracked_target_performance.tsv",
    )
    write_tsv(
        dataframe=observation_df,
        path=out_dir / "minimap_real_world_observation_level.tsv",
    )
    write_tsv(
        dataframe=real_world_run_df,
        path=out_dir / "minimap_real_world_run_summary.tsv",
    )
    write_tsv(
        dataframe=real_world_reference_df,
        path=out_dir / "minimap_real_world_reference_summary.tsv",
    )
    write_tsv(dataframe=manifest_df, path=out_dir / "minimap_summary_run_info.tsv")
    write_tsv(
        dataframe=input_root_manifest_df,
        path=out_dir / "minimap_input_root_manifest.tsv",
    )
    write_tsv(
        dataframe=summary_counts_df,
        path=out_dir / "minimap_summary_counts.tsv",
    )


if __name__ == "__main__":
    main()
