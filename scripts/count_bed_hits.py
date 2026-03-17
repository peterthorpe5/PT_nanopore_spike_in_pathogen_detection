#!/usr/bin/env python3
"""Count BED alignments that match one or more target substrings."""

from __future__ import annotations

import argparse
import csv
import sys


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Count BED hits matching target labels or aliases."
    )
    parser.add_argument(
        "--bed",
        required=True,
        help="Input BED file.",
    )
    parser.add_argument(
        "--target_label",
        action="append",
        default=[],
        help="Target label substring to search for in the BED reference name.",
    )
    parser.add_argument(
        "--aliases_tsv",
        default="",
        help=(
            "Optional TSV mapping file with columns target_label and alias. "
            "Aliases linked to any requested target label are also used."
        ),
    )
    parser.add_argument(
        "--print_matching_refs",
        action="store_true",
        help="Print a TSV of matching reference names and counts to stdout.",
    )
    return parser.parse_args()


def read_aliases(path: str, target_labels: list[str]) -> list[str]:
    """Read alias strings for requested target labels."""
    aliases: list[str] = []
    if not path:
        return aliases
    target_set = {label.casefold() for label in target_labels}
    with open(path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"target_label", "alias"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(
                "Alias TSV must contain target_label and alias columns."
            )
        for row in reader:
            if row["target_label"].casefold() in target_set:
                aliases.append(row["alias"])
    return aliases


def main() -> None:
    """Run BED hit counting."""
    args = parse_args()
    labels = list(args.target_label)
    labels.extend(read_aliases(path=args.aliases_tsv, target_labels=args.target_label))
    labels = [label.casefold() for label in labels if label]

    hit_count = 0
    ref_counts: dict[str, int] = {}

    with open(args.bed, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            ref_name = fields[0]
            ref_name_cf = ref_name.casefold()
            if labels and not any(label in ref_name_cf for label in labels):
                continue
            hit_count += 1
            ref_counts[ref_name] = ref_counts.get(ref_name, 0) + 1

    if args.print_matching_refs:
        writer = csv.writer(sys.stdout, delimiter="\t")
        writer.writerow(["reference_name", "n_hits"])
        for ref_name, n_hits in sorted(
            ref_counts.items(),
            key=lambda item: (-item[1], item[0]),
        ):
            writer.writerow([ref_name, n_hits])
    else:
        print(hit_count)


if __name__ == "__main__":
    main()
