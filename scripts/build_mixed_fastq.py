#!/usr/bin/env python3
"""Build a mixed gzipped FASTQ from one background and one or more spike files."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build a mixed gzipped FASTQ from gzipped FASTQ inputs."
    )
    parser.add_argument(
        "--background_fastq_gz",
        required=True,
        help="Background gzipped FASTQ file.",
    )
    parser.add_argument(
        "--spike_fastq_gz",
        nargs="+",
        required=True,
        help="One or more spike gzipped FASTQ files.",
    )
    parser.add_argument(
        "--out_fastq_gz",
        required=True,
        help="Output mixed gzipped FASTQ file.",
    )
    return parser.parse_args()


def validate_fastq_gz(path: str) -> None:
    """Validate gzip readability for a gzipped FASTQ file."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(f"Truncated FASTQ record found in {path}")
            if not header.startswith("@"):
                raise ValueError(f"Invalid FASTQ header found in {path}: {header.rstrip()}")


def copy_fastq(src_path: str, out_handle) -> None:
    """Copy FASTQ records from gzipped source to open text handle."""
    with gzip.open(src_path, "rt", encoding="utf-8", errors="replace") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(f"Truncated FASTQ record found in {src_path}")
            out_handle.write(header)
            out_handle.write(seq)
            out_handle.write(plus)
            out_handle.write(qual)


def main() -> None:
    """Run the mixed FASTQ creation workflow."""
    args = parse_args()
    validate_fastq_gz(args.background_fastq_gz)
    for path in args.spike_fastq_gz:
        validate_fastq_gz(path)

    out_path = Path(args.out_fastq_gz)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(out_path, "wt", encoding="utf-8") as out_handle:
        copy_fastq(src_path=args.background_fastq_gz, out_handle=out_handle)
        for spike_path in args.spike_fastq_gz:
            copy_fastq(src_path=spike_path, out_handle=out_handle)


if __name__ == "__main__":
    main()
