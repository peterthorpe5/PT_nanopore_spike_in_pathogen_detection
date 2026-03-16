#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spikein_utils.py

Utilities for nanopore spike-in experiments:
- subsample FASTQ.GZ to N reads (reservoir sampling)
- combine NanoSim outputs into one FASTQ
- summarise Kraken2 report TSV for a target label
- count target hits from minimap-derived BED
- guess a reasonable target label from a FASTA filename/header

All outputs are tab-separated (TSV). No comma-separated outputs are produced.
"""

from __future__ import annotations

import argparse
import gzip
import os
import random
import re
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Tuple


def _open_text(path: str, mode: str):
    """Open a text file path, transparently handling gzip where relevant."""
    if path.endswith(".gz"):
        return gzip.open(path, mode=mode, encoding="utf-8", errors="replace")
    return open(path, mode=mode, encoding="utf-8", errors="replace")


def reservoir_sample_fastq(
    *,
    fastq_gz: str,
    n_reads: int,
    seed: int,
) -> List[Tuple[str, str, str, str]]:
    """
    Reservoir-sample N reads from a gzipped FASTQ without loading all reads.

    Parameters
    ----------
    fastq_gz
        Input FASTQ.GZ path.
    n_reads
        Number of reads to sample. If 0, returns an empty list.
    seed
        Random seed.

    Returns
    -------
    list
        List of FASTQ 4-line tuples (header, seq, plus, qual).
    """
    if n_reads < 0:
        raise ValueError("n_reads must be >= 0")

    rng = random.Random(seed)
    reservoir: List[Tuple[str, str, str, str]] = []

    if n_reads == 0:
        return reservoir

    with _open_text(fastq_gz, "rt") as fh:
        i = 0
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                break

            record = (header, seq, plus, qual)

            if i < n_reads:
                reservoir.append(record)
            else:
                j = rng.randint(0, i)
                if j < n_reads:
                    reservoir[j] = record
            i += 1

    return reservoir


def write_fastq_gz(
    *,
    records: Iterable[Tuple[str, str, str, str]],
    out_fastq_gz: str,
) -> None:
    """
    Write FASTQ records to a gzipped FASTQ.

    Parameters
    ----------
    records
        Iterable of FASTQ 4-line tuples.
    out_fastq_gz
        Output FASTQ.GZ path.
    """
    out_path = Path(out_fastq_gz)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(out_fastq_gz, mode="wt", encoding="utf-8") as out_fh:
        for header, seq, plus, qual in records:
            out_fh.write(header)
            out_fh.write(seq)
            out_fh.write(plus)
            out_fh.write(qual)


def sample_fastq_cli(args: argparse.Namespace) -> None:
    """CLI wrapper for reservoir FASTQ sampling."""
    records = reservoir_sample_fastq(
        fastq_gz=args.fastq_gz,
        n_reads=args.n_reads,
        seed=args.seed,
    )
    write_fastq_gz(records=records, out_fastq_gz=args.out_fastq_gz)


def combine_nanosim_fastq_cli(args: argparse.Namespace) -> None:
    """
    Combine NanoSim-produced FASTQ files into a single FASTQ.

    NanoSim may emit multiple FASTQ outputs depending on version/settings
    (e.g. aligned/unaligned categories). This function concatenates any FASTQ
    files with the given prefix.

    Parameters
    ----------
    sim_prefix
        Prefix path used with NanoSim --output.
    out_fastq
        Output FASTQ path (plain text).
    """
    sim_prefix = args.sim_prefix
    out_fastq = args.out_fastq

    candidates = []
    prefix_dir = str(Path(sim_prefix).parent)
    prefix_base = Path(sim_prefix).name

    for fn in os.listdir(prefix_dir or "."):
        if not fn.startswith(prefix_base):
            continue
        if fn.endswith(".fastq") or fn.endswith(".fq"):
            candidates.append(str(Path(prefix_dir) / fn))

    if not candidates:
        raise SystemExit(
            f"ERROR: no NanoSim fastq outputs found for prefix: {sim_prefix}"
        )

    candidates = sorted(candidates)

    with open(out_fastq, mode="wt", encoding="utf-8") as out_fh:
        for fp in candidates:
            with open(fp, mode="rt", encoding="utf-8", errors="replace") as in_fh:
                for line in in_fh:
                    out_fh.write(line)


def guess_label_from_fasta(*, fasta: str) -> str:
    """
    Guess a simple target label from a FASTA filename or first header.

    This is used to parse Kraken2 reports and minimap BED outputs with a
    pragmatic heuristic. You can override this later once we decide the most
    reliable label (e.g. taxon name vs accession).

    Parameters
    ----------
    fasta
        Path to FASTA.

    Returns
    -------
    str
        Guessed label string.
    """
    name = Path(fasta).name
    name = re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", name, flags=re.IGNORECASE)

    # Common pattern: GCF_..._Plas_inui_San_Antonio_1_V1_genomic
    m = re.search(r"(Plas[_\.][A-Za-z0-9_\.]+)", name)
    if m:
        return m.group(1).replace(".", "_")

    # Fallback to first header token
    opener = gzip.open if fasta.endswith(".gz") else open
    with opener(fasta, mode="rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                token = line[1:].strip().split()[0]
                return token

    return name


def guess_label_cli(args: argparse.Namespace) -> None:
    """CLI wrapper to print guessed label."""
    print(guess_label_from_fasta(fasta=args.fasta))


def summarise_kraken_report(
    *,
    kraken_report_tsv: str,
    target_label: str,
) -> Tuple[int, int]:
    """
    Summarise a Kraken2 report for a given target label.

    The Kraken2 report format is whitespace-delimited; we treat it as TSV-like
    by splitting on whitespace.

    We return:
      - classified reads
      - target reads (reads assigned directly to rows whose name contains
        target_label)

    Parameters
    ----------
    kraken_report_tsv
        Kraken2 --report output path.
    target_label
        Substring match against the taxon name column.

    Returns
    -------
    tuple
        (classified_reads, target_reads)
    """
    total_reads = None
    unclassified_reads = 0
    target_reads = 0
    target_label_lower = target_label.lower()

    with open(kraken_report_tsv, mode="rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split()
            if len(parts) < 6:
                continue

            try:
                clade_reads = int(parts[1])
                taxon_reads = int(parts[2])
            except ValueError:
                continue

            rank = parts[3]
            name = " ".join(parts[5:]).strip()

            if rank == "R" and name.lower() == "root":
                total_reads = clade_reads

            if rank == "U":
                unclassified_reads = clade_reads

            if target_label_lower in name.lower():
                target_reads += taxon_reads

    if total_reads is None:
        raise ValueError(
            f"Could not determine total reads from Kraken report: {kraken_report_tsv}"
        )

    classified_reads = total_reads - unclassified_reads
    return classified_reads, target_reads



def summarise_kraken_cli(args: argparse.Namespace) -> None:
    """CLI wrapper to write a 1-row TSV kraken summary."""
    classified_reads, target_reads = summarise_kraken_report(
        kraken_report_tsv=args.kraken_report_tsv,
        target_label=args.target_label,
    )

    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, mode="wt", encoding="utf-8") as out_fh:
        out_fh.write("kraken_report\tclassified_reads\ttarget_reads\ttarget_label\n")
        out_fh.write(
            f"{args.kraken_report_tsv}\t{classified_reads}\t{target_reads}\t{args.target_label}\n"
        )


def count_minimap_target_alignments(
    *,
    bed: str,
    target_label: str,
) -> int:
    """
    Count alignments in a BED file whose reference (chrom) contains target_label.

    Parameters
    ----------
    bed
        BED path output by bamToBed.
    target_label
        Substring used to match against the BED chrom column.

    Returns
    -------
    int
        Count of matching BED rows.
    """
    count = 0
    target_label_lower = target_label.lower()

    with open(bed, mode="rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom = line.split("\t", maxsplit=1)[0]
            if target_label_lower in chrom.lower():
                count += 1
    return count


def count_minimap_target_cli(args: argparse.Namespace) -> None:
    """CLI wrapper to count target alignments from a BED file."""
    print(
        count_minimap_target_alignments(
            bed=args.bed,
            target_label=args.target_label,
        )
    )


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""
    parser = argparse.ArgumentParser(prog="spikein_utils.py")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_sample = sub.add_parser("sample-fastq")
    p_sample.add_argument("--fastq_gz", required=True)
    p_sample.add_argument("--n_reads", required=True, type=int)
    p_sample.add_argument("--seed", required=True, type=int)
    p_sample.add_argument("--out_fastq_gz", required=True)
    p_sample.set_defaults(func=sample_fastq_cli)

    p_combine = sub.add_parser("combine-nanosim-fastq")
    p_combine.add_argument("--sim_prefix", required=True)
    p_combine.add_argument("--out_fastq", required=True)
    p_combine.set_defaults(func=combine_nanosim_fastq_cli)

    p_guess = sub.add_parser("guess-label")
    p_guess.add_argument("--fasta", required=True)
    p_guess.set_defaults(func=guess_label_cli)

    p_k = sub.add_parser("summarise-kraken")
    p_k.add_argument("--kraken_report_tsv", required=True)
    p_k.add_argument("--target_label", required=True)
    p_k.add_argument("--out_tsv", required=True)
    p_k.set_defaults(func=summarise_kraken_cli)

    p_m = sub.add_parser("count-minimap-target")
    p_m.add_argument("--bed", required=True)
    p_m.add_argument("--target_label", required=True)
    p_m.set_defaults(func=count_minimap_target_cli)

    return parser


def main() -> None:
    """Entry point."""
    parser = build_arg_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
