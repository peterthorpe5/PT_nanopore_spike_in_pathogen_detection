#!/usr/bin/env python3
"""Summarise basic assembly statistics from a FASTA file."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarise basic assembly statistics from a FASTA file."
    )
    parser.add_argument(
        "--assembly_fasta",
        required=True,
        help="Path to the assembly FASTA file.",
    )
    parser.add_argument(
        "--out_tsv",
        required=False,
        default=None,
        help="Optional output TSV path. If omitted, write to stdout.",
    )
    return parser.parse_args()


def read_fasta_lengths(assembly_fasta: str) -> list[int]:
    """Read contig lengths from a FASTA file.

    Args:
        assembly_fasta: Path to the assembly FASTA file.

    Returns:
        A list of contig lengths.
    """
    lengths: list[int] = []
    current_length = 0

    with open(assembly_fasta, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_length > 0:
                    lengths.append(current_length)
                current_length = 0
            else:
                current_length += len(line)

    if current_length > 0:
        lengths.append(current_length)

    return lengths


def calculate_n50(lengths: list[int]) -> int:
    """Calculate N50 from a list of contig lengths.

    Args:
        lengths: List of contig lengths.

    Returns:
        N50 value.
    """
    if not lengths:
        return 0

    sorted_lengths = sorted(lengths, reverse=True)
    half_total = sum(sorted_lengths) / 2
    running_total = 0

    for length in sorted_lengths:
        running_total += length
        if running_total >= half_total:
            return length

    return 0


def summarise_lengths(lengths: list[int]) -> dict[str, int]:
    """Generate summary statistics for contig lengths.

    Args:
        lengths: List of contig lengths.

    Returns:
        Dictionary of summary statistics.
    """
    if not lengths:
        return {
            "n_contigs": 0,
            "total_bases": 0,
            "largest_contig": 0,
            "n50": 0,
        }

    return {
        "n_contigs": len(lengths),
        "total_bases": sum(lengths),
        "largest_contig": max(lengths),
        "n50": calculate_n50(lengths),
    }


def format_tsv(stats: dict[str, int], assembly_fasta: str) -> str:
    """Format assembly statistics as TSV.

    Args:
        stats: Summary statistics dictionary.
        assembly_fasta: Path to the input FASTA.

    Returns:
        TSV-formatted string.
    """
    header = "assembly_fasta\tn_contigs\ttotal_bases\tlargest_contig\tn50"
    row = (
        f"{assembly_fasta}\t"
        f"{stats['n_contigs']}\t"
        f"{stats['total_bases']}\t"
        f"{stats['largest_contig']}\t"
        f"{stats['n50']}"
    )
    return f"{header}\n{row}\n"


def main() -> None:
    """Run the assembly statistics summariser."""
    args = parse_args()
    lengths = read_fasta_lengths(assembly_fasta=args.assembly_fasta)
    stats = summarise_lengths(lengths=lengths)
    output_text = format_tsv(
        stats=stats,
        assembly_fasta=args.assembly_fasta,
    )

    if args.out_tsv:
        out_path = Path(args.out_tsv)
        out_path.write_text(output_text, encoding="utf-8")
    else:
        print(output_text, end="")


if __name__ == "__main__":
    main()