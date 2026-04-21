#!/usr/bin/env python3
"""Summarise minimap2 BED hits against a reference FASTA with taxon metadata.

This script parses a filtered BED file produced from minimap2 alignments,
aggregates alignments by expected target and by reported taxon, and writes
summary tables that can be folded into the existing spike-in HTML reporting
workflow.

Targets are matched by exact reference identifiers taken from one or more
FASTA files used for simulation. Reported taxa are inferred from the main
reference FASTA header lines, including ``kraken:taxid|...`` metadata when
present.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


TAXID_RE = re.compile(r"kraken:taxid\|([^,\s]+)")
SPECIES_TOKEN_RE = re.compile(r"(?:^|_)([A-Z]\.[A-Za-z0-9._-]+)$")


@dataclass(frozen=True)
class ReferenceRecord:
    """Store parsed metadata for a single reference sequence."""

    ref_id: str
    header: str
    taxid: str
    taxon_name: str


@dataclass(frozen=True)
class TargetDefinition:
    """Store the label and exact reference identifiers for one target."""

    label: str
    fasta_path: Path
    reference_ids: frozenset[str]


@dataclass
class TaxonCounter:
    """Store aggregated hit counts for a reported taxon."""

    taxid: str
    taxon_name: str
    reference_ids: set[str]
    qnames: set[str]
    alignment_count: int = 0


@dataclass
class TargetCounter:
    """Store aggregated hit counts for an expected target."""

    label: str
    fasta_path: Path
    reference_ids: set[str]
    qnames: set[str]
    alignment_count: int = 0
    matched_taxa: set[str] | None = None

    def __post_init__(self) -> None:
        """Initialise mutable containers after dataclass creation."""
        if self.matched_taxa is None:
            self.matched_taxa = set()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Summarise minimap2 BED hits by exact target FASTA and by reported "
            "taxon derived from the main reference FASTA headers."
        )
    )
    parser.add_argument("--bed", required=True, help="Filtered BED file.")
    parser.add_argument(
        "--reference_fasta",
        required=True,
        help=(
            "Main reference FASTA used for minimap2. Headers may include "
            "kraken:taxid metadata."
        ),
    )
    parser.add_argument(
        "--target_fasta",
        action="append",
        default=[],
        help=(
            "Target FASTA used for simulation. Repeat together with "
            "--target_label for multiple targets."
        ),
    )
    parser.add_argument(
        "--target_label",
        action="append",
        default=[],
        help="Human-readable label for the corresponding --target_fasta.",
    )
    parser.add_argument(
        "--target_tsv",
        default=None,
        help=(
            "Optional TSV containing at least two columns: FASTA path and "
            "target label. Used mainly for multi-genome panels."
        ),
    )
    parser.add_argument(
        "--out_target_tsv",
        required=True,
        help="Output TSV with target-level counts.",
    )
    parser.add_argument(
        "--out_reported_taxa_tsv",
        required=True,
        help="Output TSV with taxon-level counts.",
    )
    parser.add_argument(
        "--out_meta_tsv",
        default=None,
        help="Optional output TSV with overall mapped-read totals.",
    )
    return parser.parse_args()


def read_fasta_headers(fasta_path: Path) -> dict[str, str]:
    """Read FASTA headers keyed by exact reference identifier.

    Parameters
    ----------
    fasta_path : Path
        FASTA file path.

    Returns
    -------
    dict[str, str]
        Mapping from exact reference identifier to full header line without the
        leading ``>``.
    """
    headers: dict[str, str] = {}
    with fasta_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            ref_id = header.split()[0]
            headers[ref_id] = header
    return headers


def derive_taxon_name(ref_id: str, header: str, taxid: str) -> str:
    """Derive a readable taxon label from a reference header.

    Parameters
    ----------
    ref_id : str
        Exact reference identifier.
    header : str
        Full FASTA header.
    taxid : str
        Parsed taxonomic identifier or fallback placeholder.

    Returns
    -------
    str
        Readable taxon label.
    """
    species_match = SPECIES_TOKEN_RE.search(ref_id)
    if species_match:
        token = species_match.group(1).replace("_", " ")
        return token

    header_tokens = header.split()
    if len(header_tokens) > 1:
        token = header_tokens[1].strip(",")
        if token:
            return token

    if taxid:
        return f"taxid_{taxid}"
    return ref_id


def parse_reference_metadata(reference_fasta: Path) -> dict[str, ReferenceRecord]:
    """Parse main reference FASTA metadata.

    Parameters
    ----------
    reference_fasta : Path
        Main minimap2 reference FASTA.

    Returns
    -------
    dict[str, ReferenceRecord]
        Mapping from exact reference identifier to parsed metadata.
    """
    reference_records: dict[str, ReferenceRecord] = {}
    headers = read_fasta_headers(reference_fasta)
    for ref_id, header in headers.items():
        taxid_match = TAXID_RE.search(header)
        taxid = taxid_match.group(1) if taxid_match else "unknown_taxid"
        taxon_name = derive_taxon_name(ref_id=ref_id, header=header, taxid=taxid)
        reference_records[ref_id] = ReferenceRecord(
            ref_id=ref_id,
            header=header,
            taxid=taxid,
            taxon_name=taxon_name,
        )
    return reference_records


def load_targets_from_tsv(target_tsv: Path) -> list[tuple[Path, str]]:
    """Load target FASTA definitions from a TSV file.

    Parameters
    ----------
    target_tsv : Path
        TSV file with at least two columns: FASTA path and label.

    Returns
    -------
    list[tuple[Path, str]]
        Ordered FASTA path and label pairs.
    """
    targets: list[tuple[Path, str]] = []
    with target_tsv.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row_index, row in enumerate(reader, start=1):
            if not row:
                continue
            if row_index == 1 and row[0].strip().lower() in {"fasta", "path", "fasta_path"}:
                continue
            if len(row) < 2:
                continue
            fasta_str = row[0].strip()
            label = row[1].strip()
            if not fasta_str or not label:
                continue
            targets.append((Path(fasta_str).expanduser(), label))
    return targets


def build_target_definitions(args: argparse.Namespace) -> list[TargetDefinition]:
    """Build target definitions from command-line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    list[TargetDefinition]
        Target definitions with exact reference identifiers.
    """
    target_pairs: list[tuple[Path, str]] = []
    if args.target_tsv:
        target_pairs.extend(load_targets_from_tsv(Path(args.target_tsv).expanduser()))

    if args.target_fasta or args.target_label:
        if len(args.target_fasta) != len(args.target_label):
            raise ValueError(
                "The number of --target_fasta and --target_label arguments must match."
            )
        target_pairs.extend(
            (Path(fasta_path).expanduser(), label)
            for fasta_path, label in zip(args.target_fasta, args.target_label, strict=True)
        )

    definitions: list[TargetDefinition] = []
    for fasta_path, label in target_pairs:
        headers = read_fasta_headers(fasta_path)
        definitions.append(
            TargetDefinition(
                label=label,
                fasta_path=fasta_path,
                reference_ids=frozenset(headers.keys()),
            )
        )
    return definitions


def iter_bed_rows(bed_path: Path) -> Iterable[tuple[str, str]]:
    """Yield reference identifiers and query names from a BED file.

    Parameters
    ----------
    bed_path : Path
        BED file path.

    Yields
    ------
    tuple[str, str]
        Reference identifier and query name.
    """
    with bed_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            yield fields[0], fields[3]


def summarise_hits(
    bed_path: Path,
    reference_records: dict[str, ReferenceRecord],
    target_definitions: list[TargetDefinition],
) -> tuple[list[TargetCounter], list[TaxonCounter], int, int]:
    """Summarise BED hits by target and by reported taxon.

    Parameters
    ----------
    bed_path : Path
        Filtered BED file.
    reference_records : dict[str, ReferenceRecord]
        Parsed reference metadata.
    target_definitions : list[TargetDefinition]
        Exact target reference definitions.

    Returns
    -------
    tuple[list[TargetCounter], list[TaxonCounter], int, int]
        Target counters, taxon counters, total alignments, and total unique
        mapped reads.
    """
    target_counters = {
        target.label: TargetCounter(
            label=target.label,
            fasta_path=target.fasta_path,
            reference_ids=set(target.reference_ids),
            qnames=set(),
        )
        for target in target_definitions
    }

    taxon_counters: dict[str, TaxonCounter] = {}
    all_qnames: set[str] = set()
    total_alignments = 0

    ref_to_targets: dict[str, set[str]] = defaultdict(set)
    for target in target_definitions:
        for ref_id in target.reference_ids:
            ref_to_targets[ref_id].add(target.label)

    for ref_id, qname in iter_bed_rows(bed_path):
        if ref_id not in reference_records:
            continue

        total_alignments += 1
        all_qnames.add(qname)
        record = reference_records[ref_id]

        if record.taxid not in taxon_counters:
            taxon_counters[record.taxid] = TaxonCounter(
                taxid=record.taxid,
                taxon_name=record.taxon_name,
                reference_ids=set(),
                qnames=set(),
            )
        taxon_counter = taxon_counters[record.taxid]
        taxon_counter.alignment_count += 1
        taxon_counter.reference_ids.add(ref_id)
        taxon_counter.qnames.add(qname)

        matched_labels = ref_to_targets.get(ref_id, set())
        for label in matched_labels:
            target_counter = target_counters[label]
            target_counter.alignment_count += 1
            target_counter.qnames.add(qname)
            target_counter.matched_taxa.add(record.taxon_name)

    target_counter_list = [target_counters[target.label] for target in target_definitions]
    taxon_counter_list = sorted(
        taxon_counters.values(),
        key=lambda item: (-item.alignment_count, item.taxon_name, item.taxid),
    )
    return target_counter_list, taxon_counter_list, total_alignments, len(all_qnames)


def write_target_summary(target_counters: list[TargetCounter], out_path: Path) -> None:
    """Write the target-level summary TSV.

    Parameters
    ----------
    target_counters : list[TargetCounter]
        Aggregated target counters.
    out_path : Path
        Output TSV path.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "target_label",
                "target_fasta",
                "target_reference_count",
                "target_alignments",
                "target_unique_reads",
                "matched_taxa",
            ]
        )
        for counter in target_counters:
            writer.writerow(
                [
                    counter.label,
                    str(counter.fasta_path),
                    len(counter.reference_ids),
                    counter.alignment_count,
                    len(counter.qnames),
                    "; ".join(sorted(counter.matched_taxa or set())),
                ]
            )


def write_reported_taxa_summary(
    taxon_counters: list[TaxonCounter],
    target_counters: list[TargetCounter],
    out_path: Path,
) -> None:
    """Write the reported-taxa summary TSV.

    Parameters
    ----------
    taxon_counters : list[TaxonCounter]
        Aggregated taxon counters.
    target_counters : list[TargetCounter]
        Aggregated target counters.
    out_path : Path
        Output TSV path.
    """
    target_taxon_lookup: dict[str, set[str]] = defaultdict(set)
    for counter in target_counters:
        for taxon_name in counter.matched_taxa or set():
            target_taxon_lookup[taxon_name].add(counter.label)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "taxon_name",
                "taxid",
                "reference_count",
                "alignments",
                "unique_reads",
                "matched_target_labels",
            ]
        )
        for counter in taxon_counters:
            writer.writerow(
                [
                    counter.taxon_name,
                    counter.taxid,
                    len(counter.reference_ids),
                    counter.alignment_count,
                    len(counter.qnames),
                    "; ".join(sorted(target_taxon_lookup.get(counter.taxon_name, set()))),
                ]
            )


def write_meta_summary(
    out_path: Path,
    total_alignments: int,
    total_unique_reads: int,
) -> None:
    """Write overall mapped totals.

    Parameters
    ----------
    out_path : Path
        Output TSV path.
    total_alignments : int
        Total alignment records in the BED file.
    total_unique_reads : int
        Total unique mapped read names in the BED file.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["total_alignments", total_alignments])
        writer.writerow(["total_unique_reads", total_unique_reads])


def main() -> None:
    """Run the BED-to-summary conversion workflow."""
    args = parse_args()

    bed_path = Path(args.bed).expanduser().resolve()
    reference_fasta = Path(args.reference_fasta).expanduser().resolve()
    out_target_tsv = Path(args.out_target_tsv).expanduser().resolve()
    out_reported_taxa_tsv = Path(args.out_reported_taxa_tsv).expanduser().resolve()
    out_meta_tsv = (
        Path(args.out_meta_tsv).expanduser().resolve()
        if args.out_meta_tsv
        else None
    )

    if not bed_path.exists():
        raise FileNotFoundError(f"BED file not found: {bed_path}")
    if not reference_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference_fasta}")

    target_definitions = build_target_definitions(args)
    reference_records = parse_reference_metadata(reference_fasta)
    target_counters, taxon_counters, total_alignments, total_unique_reads = summarise_hits(
        bed_path=bed_path,
        reference_records=reference_records,
        target_definitions=target_definitions,
    )

    write_target_summary(target_counters=target_counters, out_path=out_target_tsv)
    write_reported_taxa_summary(
        taxon_counters=taxon_counters,
        target_counters=target_counters,
        out_path=out_reported_taxa_tsv,
    )
    if out_meta_tsv is not None:
        write_meta_summary(
            out_path=out_meta_tsv,
            total_alignments=total_alignments,
            total_unique_reads=total_unique_reads,
        )


if __name__ == "__main__":
    main()
