#!/usr/bin/env python3
"""Calculate simple contig and base-count statistics for an assembly FASTA."""

from __future__ import annotations

import argparse


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate assembly contig and base statistics."
    )
    parser.add_argument(
        "--assembly_fasta",
        required=True,
        help="Input assembly FASTA file.",
    )
    return parser.parse_args()


def main() -> None:
    """Run assembly statistics calculation."""
    args = parse_args()
    n_contigs = 0
    n_bases = 0
    seq_len = 0

    with open(args.assembly_fasta, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if seq_len > 0:
                    n_contigs += 1
                    n_bases += seq_len
                seq_len = 0
            else:
                seq_len += len(line)

    if seq_len > 0:
        n_contigs += 1
        n_bases += seq_len

    print(f"{n_contigs}\t{n_bases}")


if __name__ == "__main__":
    main()
