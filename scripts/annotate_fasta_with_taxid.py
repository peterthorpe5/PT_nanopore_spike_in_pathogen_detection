#!/usr/bin/env python3
"""Rewrite FASTA headers with taxon annotations for MetaMaps database builds."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import TextIO


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Rewrite FASTA headers so each sequence header begins with "
            "kraken:taxid|<taxid>|, which can be used in downstream "
            "MetaMaps database construction workflows."
        )
    )
    parser.add_argument(
        "--input_fasta",
        required=True,
        help="Input FASTA file, plain text or gzipped.",
    )
    parser.add_argument(
        "--taxid",
        required=True,
        type=int,
        help="NCBI taxon ID to attach to all sequences in the FASTA.",
    )
    parser.add_argument(
        "--out_fasta",
        required=True,
        help="Output FASTA with rewritten headers.",
    )
    return parser.parse_args()


def open_maybe_gzip(path: str) -> TextIO:
    """Open a plain-text or gzipped FASTA file for reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def rewrite_headers(input_fasta: str, taxid: int, out_fasta: str) -> int:
    """Rewrite FASTA headers and return the number of sequences processed."""
    n_sequences = 0
    output_path = Path(out_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open_maybe_gzip(input_fasta) as in_handle, open(
        output_path, "wt", encoding="utf-8"
    ) as out_handle:
        for line in in_handle:
            if line.startswith(">"):
                n_sequences += 1
                header = line[1:].strip()
                seq_id = header.split()[0]
                out_handle.write(f">kraken:taxid|{taxid}|{seq_id}\n")
            else:
                out_handle.write(line)

    return n_sequences


def main() -> None:
    """Run FASTA header annotation."""
    args = parse_args()
    n_sequences = rewrite_headers(
        input_fasta=args.input_fasta,
        taxid=args.taxid,
        out_fasta=args.out_fasta,
    )
    if n_sequences == 0:
        raise SystemExit(f"No FASTA sequences found in {args.input_fasta}")
    print(f"[INFO] Rewrote {n_sequences} FASTA headers into {args.out_fasta}")


if __name__ == "__main__":
    main()
