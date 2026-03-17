#!/usr/bin/env python3
"""Combine NanoSim FASTQ outputs into a single FASTQ file."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Combine NanoSim FASTQ outputs into one FASTQ file."
    )
    parser.add_argument(
        "--sim_prefix",
        required=True,
        help="NanoSim simulation output prefix.",
    )
    parser.add_argument(
        "--out_fastq",
        required=True,
        help="Output plain FASTQ file.",
    )
    return parser.parse_args()


def candidate_paths(sim_prefix: str) -> list[Path]:
    """Return candidate NanoSim FASTQ output files for a prefix."""
    prefix = Path(sim_prefix)
    parent = prefix.parent
    stem = prefix.name
    patterns = [
        f"{stem}*fastq",
        f"{stem}*fastq.gz",
        f"{stem}*fq",
        f"{stem}*fq.gz",
    ]
    hits: list[Path] = []
    for pattern in patterns:
        for path in sorted(parent.glob(pattern)):
            if path not in hits:
                hits.append(path)
    return hits


def is_fastq_path(path: Path) -> bool:
    """Return True if the path looks like a FASTQ-like file."""
    suffixes = "".join(path.suffixes).lower()
    return any(
        suffixes.endswith(token)
        for token in [".fastq", ".fastq.gz", ".fq", ".fq.gz"]
    )


def open_maybe_gzip(path: Path):
    """Open plain or gzipped text file for reading."""
    if "".join(path.suffixes).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def copy_fastq_records(src_paths: list[Path], out_fastq: str) -> int:
    """Copy FASTQ records from source files to a single output FASTQ file."""
    n_records = 0
    Path(out_fastq).parent.mkdir(parents=True, exist_ok=True)
    with open(out_fastq, "wt", encoding="utf-8") as out_handle:
        for src_path in src_paths:
            if not is_fastq_path(src_path):
                continue
            with open_maybe_gzip(src_path) as in_handle:
                while True:
                    header = in_handle.readline()
                    if not header:
                        break
                    seq = in_handle.readline()
                    plus = in_handle.readline()
                    qual = in_handle.readline()
                    if not qual:
                        raise ValueError(
                            f"Truncated FASTQ record in {src_path}."
                        )
                    out_handle.write(header)
                    out_handle.write(seq)
                    out_handle.write(plus)
                    out_handle.write(qual)
                    n_records += 1
    return n_records


def main() -> None:
    """Run the NanoSim FASTQ combination workflow."""
    args = parse_args()
    paths = candidate_paths(sim_prefix=args.sim_prefix)
    if not paths:
        raise SystemExit(f"No NanoSim output files found for prefix: {args.sim_prefix}")
    copied = copy_fastq_records(src_paths=paths, out_fastq=args.out_fastq)
    if copied == 0:
        raise SystemExit(
            f"No FASTQ records found for NanoSim prefix: {args.sim_prefix}"
        )


if __name__ == "__main__":
    main()
