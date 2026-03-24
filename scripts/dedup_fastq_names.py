#!/usr/bin/env python3
"""Deduplicate read names in a gzipped FASTQ while preserving record order."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Deduplicate FASTQ read identifiers from a gzipped FASTQ."
    )
    parser.add_argument(
        "--input_fastq_gz",
        required=True,
        help="Input gzipped FASTQ file.",
    )
    parser.add_argument(
        "--output_fastq",
        required=True,
        help="Output plain FASTQ file.",
    )
    return parser.parse_args()


def main() -> None:
    """Run FASTQ header deduplication."""
    args = parse_args()
    seen: dict[str, int] = {}
    record_n = 0
    Path(args.output_fastq).parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(args.input_fastq_gz, "rt", encoding="utf-8", errors="replace") as infh, \
            open(args.output_fastq, "wt", encoding="utf-8") as outfh:
        while True:
            header = infh.readline()
            if not header:
                break
            seq = infh.readline()
            plus = infh.readline()
            qual = infh.readline()

            if not qual:
                raise ValueError("Truncated FASTQ record encountered.")

            record_n += 1
            header = header.rstrip("\n")
            seq = seq.rstrip("\n")
            plus = plus.rstrip("\n")
            qual = qual.rstrip("\n")

            if not header.startswith("@"):
                raise ValueError(
                    f"Invalid FASTQ header at record {record_n}: {header}"
                )

            core = header[1:].split()[0]
            extra = header[len(core) + 1:] if len(header) > len(core) + 1 else ""

            count = seen.get(core, 0) + 1
            seen[core] = count

            new_core = core if count == 1 else f"{core}_dup{count}"
            new_header = f"@{new_core}{extra}" if extra else f"@{new_core}"

            outfh.write(new_header + "\n")
            outfh.write(seq + "\n")
            outfh.write(plus + "\n")
            outfh.write(qual + "\n")


if __name__ == "__main__":
    main()
