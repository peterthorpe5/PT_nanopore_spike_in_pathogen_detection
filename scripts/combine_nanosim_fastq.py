#!/usr/bin/env python3
"""
Combine NanoSim FASTQ outputs into a single FASTQ file.

This script is tolerant to different NanoSim output naming conventions,
including both:

- <prefix>_aligned_reads.fastq
- <prefix>_aligned_reads0.fastq
- <prefix>_aligned_reads1.fastq
- ...

It validates that each FASTQ record has four lines and raises a clear
error if a truncated record is encountered.
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path
from typing import Iterable, List


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Combine NanoSim FASTQ files for a given prefix."
    )
    parser.add_argument(
        "--sim_prefix",
        required=True,
        help="NanoSim simulation output prefix.",
    )
    parser.add_argument(
        "--out_fastq",
        required=True,
        help="Output combined FASTQ file.",
    )
    return parser.parse_args()


def discover_fastq_files(sim_prefix: str) -> List[str]:
    """
    Discover NanoSim FASTQ output files for a prefix.

    Parameters
    ----------
    sim_prefix : str
        NanoSim output prefix.

    Returns
    -------
    List[str]
        Sorted list of FASTQ files to combine.
    """
    patterns = [
        f"{sim_prefix}_aligned_reads.fastq",
        f"{sim_prefix}_aligned_reads*.fastq",
    ]

    found: List[str] = []
    for pattern in patterns:
        found.extend(glob.glob(pattern))

    unique_sorted = sorted(set(found))
    return [path for path in unique_sorted if os.path.isfile(path)]


def iter_fastq_records(fastq_path: str) -> Iterable[str]:
    """
    Yield FASTQ records from a file as 4-line strings.

    Parameters
    ----------
    fastq_path : str
        Input FASTQ path.

    Yields
    ------
    Iterable[str]
        FASTQ records as concatenated 4-line strings.

    Raises
    ------
    ValueError
        If a truncated FASTQ record is encountered.
    """
    with open(fastq_path, mode="r", encoding="utf-8") as handle:
        while True:
            header = handle.readline()
            if header == "":
                break

            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()

            if not seq or not plus or not qual:
                raise ValueError(
                    f"Truncated FASTQ record in {fastq_path}."
                )

            yield f"{header}{seq}{plus}{qual}"


def main() -> None:
    """Run FASTQ combination."""
    args = parse_args()

    sim_prefix = args.sim_prefix
    out_fastq = args.out_fastq

    fastq_files = discover_fastq_files(sim_prefix)

    if not fastq_files:
        raise FileNotFoundError(
            f"No NanoSim FASTQ files found for prefix: {sim_prefix}"
        )

    out_path = Path(out_fastq)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, mode="w", encoding="utf-8") as out_handle:
        for fastq_path in fastq_files:
            print(f"[INFO] Adding {fastq_path}")
            for record in iter_fastq_records(fastq_path):
                out_handle.write(record)

    print(f"[INFO] Wrote combined FASTQ: {out_fastq}")


if __name__ == "__main__":
    main()
