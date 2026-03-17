#!/usr/bin/env python3
"""Create a mononucleotide-shuffled FASTA with short, safe contig names."""

from __future__ import annotations

import argparse
import csv
import gzip
import random
from collections import Counter
from pathlib import Path
from typing import Iterator, Tuple


FastaRecord = Tuple[str, str]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Shuffle FASTA contigs while preserving mononucleotide composition."
    )
    parser.add_argument(
        "--input_fasta",
        required=True,
        help="Input FASTA file, plain or gzipped.",
    )
    parser.add_argument(
        "--output_fasta",
        required=True,
        help="Output shuffled FASTA file.",
    )
    parser.add_argument(
        "--metadata_tsv",
        required=True,
        help="Output metadata TSV.",
    )
    parser.add_argument(
        "--seed",
        required=True,
        type=int,
        help="Random seed.",
    )
    return parser.parse_args()


def open_maybe_gzip(path: str, mode: str):
    """Open a plain or gzipped text file."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", errors="replace")
    return open(path, mode, encoding="utf-8", errors="replace")


def read_fasta(path: str) -> Iterator[FastaRecord]:
    """Yield FASTA records as header and sequence tuples."""
    header = None
    seq_parts: list[str] = []
    with open_maybe_gzip(path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"): 
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
        if header is not None:
            yield header, "".join(seq_parts)


def shuffle_seq(seq: str, rng: random.Random) -> str:
    """Shuffle A, C, G and T positions while preserving other characters."""
    seq_list = list(seq)
    positions = [idx for idx, base in enumerate(seq_list) if base in {"A", "C", "G", "T"}]
    bases = [seq_list[idx] for idx in positions]
    rng.shuffle(bases)
    for idx, base in zip(positions, bases):
        seq_list[idx] = base
    return "".join(seq_list)


def wrap_sequence(seq: str, width: int = 80) -> Iterator[str]:
    """Yield fixed-width sequence chunks."""
    for start in range(0, len(seq), width):
        yield seq[start:start + width]


def main() -> None:
    """Run the FASTA shuffling workflow."""
    args = parse_args()
    rng = random.Random(args.seed)
    records = list(read_fasta(args.input_fasta))
    if not records:
        raise SystemExit("No FASTA records found.")

    Path(args.output_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(args.metadata_tsv).parent.mkdir(parents=True, exist_ok=True)

    with open(args.output_fasta, "wt", encoding="utf-8") as out_handle, \
            open(args.metadata_tsv, "wt", encoding="utf-8", newline="") as meta_handle:
        writer = csv.writer(meta_handle, delimiter="\t")
        writer.writerow(
            [
                "record_index",
                "original_header",
                "shuffled_header",
                "length",
                "A",
                "C",
                "G",
                "T",
                "N_or_other",
            ]
        )

        for idx, (header, seq) in enumerate(records, start=1):
            shuffled_header = f"shuffle_contig_{idx}"
            shuffled_seq = shuffle_seq(seq=seq, rng=rng)
            out_handle.write(f">{shuffled_header}\n")
            for chunk in wrap_sequence(seq=shuffled_seq):
                out_handle.write(chunk + "\n")

            counts = Counter(seq)
            n_other = sum(
                count for base, count in counts.items() if base not in {"A", "C", "G", "T"}
            )
            writer.writerow(
                [
                    idx,
                    header,
                    shuffled_header,
                    len(seq),
                    counts.get("A", 0),
                    counts.get("C", 0),
                    counts.get("G", 0),
                    counts.get("T", 0),
                    n_other,
                ]
            )


if __name__ == "__main__":
    main()
