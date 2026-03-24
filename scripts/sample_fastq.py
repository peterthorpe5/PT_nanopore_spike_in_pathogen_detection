#!/usr/bin/env python3
"""Sample FASTQ records from a gzipped FASTQ file using reservoir sampling."""

from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path
from typing import Iterator, Tuple


FastqRecord = Tuple[str, str, str, str]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Sample records from a gzipped FASTQ file."
    )
    parser.add_argument(
        "--fastq_gz",
        required=True,
        help="Input gzipped FASTQ file.",
    )
    parser.add_argument(
        "--n_reads",
        required=True,
        type=int,
        help="Number of reads to sample.",
    )
    parser.add_argument(
        "--seed",
        required=True,
        type=int,
        help="Random seed.",
    )
    parser.add_argument(
        "--out_fastq_gz",
        required=True,
        help="Output gzipped FASTQ file.",
    )
    return parser.parse_args()


def iter_fastq_records(path: str) -> Iterator[FastqRecord]:
    """Yield FASTQ records from a gzipped FASTQ file."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        record_n = 0
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(
                    f"Truncated FASTQ record encountered at record {record_n + 1}."
                )
            record_n += 1
            if not header.startswith("@"):
                raise ValueError(
                    f"Invalid FASTQ header at record {record_n}: {header.rstrip()}"
                )
            yield header, seq, plus, qual


def reservoir_sample(path: str, n_reads: int, seed: int) -> list[FastqRecord]:
    """Sample FASTQ records using reservoir sampling."""
    if n_reads < 0:
        raise ValueError("--n_reads must be non-negative.")
    rng = random.Random(seed)
    reservoir: list[FastqRecord] = []
    for idx, record in enumerate(iter_fastq_records(path)):
        if idx < n_reads:
            reservoir.append(record)
        else:
            j = rng.randint(0, idx)
            if j < n_reads:
                reservoir[j] = record
    return reservoir


def write_fastq_gz(records: list[FastqRecord], out_path: str) -> None:
    """Write FASTQ records to a gzipped FASTQ file."""
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_path, "wt", encoding="utf-8") as handle:
        for header, seq, plus, qual in records:
            handle.write(header)
            handle.write(seq)
            handle.write(plus)
            handle.write(qual)


def main() -> None:
    """Run the FASTQ sampling workflow."""
    args = parse_args()
    sampled_records = reservoir_sample(
        path=args.fastq_gz,
        n_reads=args.n_reads,
        seed=args.seed,
    )
    write_fastq_gz(records=sampled_records, out_path=args.out_fastq_gz)


if __name__ == "__main__":
    main()
