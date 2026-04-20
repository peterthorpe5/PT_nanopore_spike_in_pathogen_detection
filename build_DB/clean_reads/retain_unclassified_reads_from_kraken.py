#!/usr/bin/env python3
"""
Retain only Kraken-unclassified reads from a FASTQ file.

This script uses a Kraken classifications file to identify reads marked as
unclassified (status ``U``) and writes a filtered FASTQ containing only those
reads.

It is intended for stringent background clean-up workflows, where the user
wants to remove all reads that Kraken could classify against the chosen
database and keep only reads that remained unclassified.

Outputs
-------
All tabular outputs are written as tab-separated files.

1. A text file containing one retained read ID per line.
2. A tab-separated summary of retained and removed read counts.
3. A filtered FASTQ file containing only Kraken-unclassified reads.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


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
            "Retain only Kraken-unclassified reads from a FASTQ file."
        )
    )
    parser.add_argument(
        "--kraken_classifications",
        required=True,
        help="Path to Kraken classifications TSV.",
    )
    parser.add_argument(
        "--input_fastq",
        required=True,
        help="Path to the input FASTQ file.",
    )
    parser.add_argument(
        "--output_fastq",
        required=True,
        help="Path to the filtered output FASTQ file.",
    )
    parser.add_argument(
        "--retained_read_ids_out",
        required=True,
        help="Path to write retained unclassified read IDs.",
    )
    parser.add_argument(
        "--summary_tsv",
        required=True,
        help="Path to write a tab-separated summary file.",
    )
    return parser.parse_args()


def open_maybe_gzip(path: Path, mode: str):
    """
    Open a plain-text or gzip-compressed file.

    Parameters
    ----------
    path : Path
        File path to open.
    mode : str
        File mode, for example ``"rt"`` or ``"wt"``.

    Returns
    -------
    IO
        Open file handle.
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode)


def load_unclassified_read_ids(classifications_path: Path) -> set[str]:
    """
    Load read IDs classified as unclassified by Kraken.

    Parameters
    ----------
    classifications_path : Path
        Path to the Kraken classifications TSV.

    Returns
    -------
    set[str]
        Set of unclassified read IDs.
    """
    unclassified_ids: set[str] = set()

    with open_maybe_gzip(path=classifications_path, mode="rt") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < 2:
                continue

            status = fields[0]
            read_id = fields[1]

            if status == "U":
                unclassified_ids.add(read_id)

    return unclassified_ids


def filter_fastq_by_read_ids(
    input_fastq: Path,
    output_fastq: Path,
    retained_read_ids: set[str],
) -> tuple[int, int]:
    """
    Retain only FASTQ reads whose IDs are present in the supplied set.

    Parameters
    ----------
    input_fastq : Path
        Path to the input FASTQ file.
    output_fastq : Path
        Path to the filtered FASTQ file.
    retained_read_ids : set[str]
        Read IDs to retain.

    Returns
    -------
    tuple[int, int]
        Counts of reads kept and removed.
    """
    kept = 0
    removed = 0

    with open_maybe_gzip(path=input_fastq, mode="rt") as in_handle,             open_maybe_gzip(path=output_fastq, mode="wt") as out_handle:
        while True:
            header = in_handle.readline()
            if not header:
                break

            sequence = in_handle.readline()
            plus = in_handle.readline()
            quality = in_handle.readline()

            if not quality:
                raise ValueError("Input FASTQ appears truncated.")

            read_id = header.strip().split()[0].lstrip("@")

            if read_id in retained_read_ids:
                out_handle.write(header)
                out_handle.write(sequence)
                out_handle.write(plus)
                out_handle.write(quality)
                kept += 1
            else:
                removed += 1

    return kept, removed


def write_summary(
    summary_tsv: Path,
    kept: int,
    removed: int,
    total_unclassified_ids: int,
) -> None:
    """
    Write a tab-separated summary of filtering results.

    Parameters
    ----------
    summary_tsv : Path
        Output path for the summary TSV.
    kept : int
        Number of reads retained in the FASTQ.
    removed : int
        Number of reads removed from the FASTQ.
    total_unclassified_ids : int
        Number of unclassified read IDs seen in the Kraken file.
    """
    with summary_tsv.open("w") as handle:
        handle.write(
            "metric\tvalue\n"
            f"kraken_unclassified_read_ids\t{total_unclassified_ids}\n"
            f"fastq_reads_retained\t{kept}\n"
            f"fastq_reads_removed\t{removed}\n"
        )


def main() -> None:
    """
    Run the Kraken-unclassified FASTQ filtering workflow.
    """
    args = parse_args()

    retained_read_ids = load_unclassified_read_ids(
        classifications_path=Path(args.kraken_classifications)
    )

    Path(args.retained_read_ids_out).write_text(
        "\n".join(sorted(retained_read_ids)) + "\n"
    )

    kept, removed = filter_fastq_by_read_ids(
        input_fastq=Path(args.input_fastq),
        output_fastq=Path(args.output_fastq),
        retained_read_ids=retained_read_ids,
    )

    write_summary(
        summary_tsv=Path(args.summary_tsv),
        kept=kept,
        removed=removed,
        total_unclassified_ids=len(retained_read_ids),
    )

    print(f"Kraken-unclassified read IDs: {len(retained_read_ids)}")
    print(f"Reads retained in FASTQ: {kept}")
    print(f"Reads removed from FASTQ: {removed}")


if __name__ == "__main__":
    main()
