#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_kraken_taxid_map.py

Build a TSV mapping of local genome FASTA files to NCBI taxonomic IDs using
NCBI assembly summary tables.

This script:
- Scans an input directory for FASTA/FA/FNA (optionally gzipped).
- Extracts assembly accessions (GCA_... or GCF_...) from filenames.
- Looks up the assembly accession in GenBank and RefSeq assembly summary tables.
- Writes a tab-separated mapping file suitable for driving kraken2-build
  --add-to-library with --taxid.

Outputs are TSV, not comma-separated.

Example:
python make_kraken_taxid_map.py \
  --genomes_dir genomes \
  --assembly_summary_genbank ncbi_metadata/assembly_summary_genbank.txt \
  --assembly_summary_refseq ncbi_metadata/assembly_summary_refseq.txt \
  --out_tsv genome_to_taxid.tsv
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes_dir", required=True)
    parser.add_argument("--assembly_summary_genbank", required=True)
    parser.add_argument("--assembly_summary_refseq", required=True)
    parser.add_argument("--out_tsv", required=True)
    parser.add_argument(
        "--extra_manual_tsv",
        required=False,
        default=None,
        help=(
            "Optional TSV with columns: genome_file<TAB>taxid<TAB>species. "
            "Used for files without GCA_/GCF_ in the name."
        ),
    )
    return parser.parse_args()


def iter_genome_files(*, genomes_dir: str) -> Iterable[Path]:
    """
    Yield genome files in a directory.

    Parameters
    ----------
    genomes_dir
        Directory containing genome FASTA files.

    Yields
    ------
    Path
        Path objects for FASTA-like files (including gzipped).
    """
    exts = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")
    for fp in sorted(Path(genomes_dir).iterdir()):
        if fp.is_file() and fp.name.lower().endswith(exts):
            yield fp


def extract_assembly_accession(*, filename: str) -> Optional[str]:
    """
    Extract an NCBI assembly accession from a filename.

    Parameters
    ----------
    filename
        Filename to parse.

    Returns
    -------
    str or None
        Assembly accession like GCA_000725905.1, or None if not found.
    """
    m = re.search(r"(GCA_\d+\.\d+|GCF_\d+\.\d+)", filename)
    return m.group(1) if m else None


def load_assembly_summary(*, path: str) -> Dict[str, Tuple[str, str]]:
    """
    Load an NCBI assembly summary file into a dict.

    Parameters
    ----------
    path
        Path to assembly_summary_*.txt

    Returns
    -------
    dict
        Mapping: assembly_accession -> (taxid, organism_name)
    """
    out: Dict[str, Tuple[str, str]] = {}
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            assembly_acc = parts[0]
            taxid = parts[5]
            organism = parts[7]
            out[assembly_acc] = (taxid, organism)
    return out


def load_manual_map(*, path: Optional[str]) -> Dict[str, Tuple[str, str]]:
    """
    Load an optional manual mapping TSV.

    Parameters
    ----------
    path
        TSV path with columns: genome_file, taxid, species

    Returns
    -------
    dict
        Mapping: genome_file -> (taxid, species)
    """
    out: Dict[str, Tuple[str, str]] = {}
    if path is None:
        return out
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header[:3] != ["genome_file", "taxid", "species"]:
            raise ValueError(
                "extra_manual_tsv must have header: genome_file<TAB>taxid<TAB>species"
            )
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            out[parts[0]] = (parts[1], parts[2])
    return out


def main() -> None:
    """Entry point."""
    args = parse_args()

    genbank = load_assembly_summary(path=args.assembly_summary_genbank)
    refseq = load_assembly_summary(path=args.assembly_summary_refseq)
    manual = load_manual_map(path=args.extra_manual_tsv)

    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "wt", encoding="utf-8") as out_fh:
        out_fh.write(
            "genome_file\tassembly_accession\ttaxid\tspecies\tsource\n"
        )

        for fp in iter_genome_files(genomes_dir=args.genomes_dir):
            assembly = extract_assembly_accession(filename=fp.name)

            if assembly is None:
                if fp.name in manual:
                    taxid, species = manual[fp.name]
                    out_fh.write(f"{fp}\tNA\t{taxid}\t{species}\tmanual\n")
                else:
                    out_fh.write(f"{fp}\tNA\tNA\tNA\tmissing_assembly\n")
                continue

            if assembly in refseq:
                taxid, species = refseq[assembly]
                out_fh.write(f"{fp}\t{assembly}\t{taxid}\t{species}\trefseq\n")
            elif assembly in genbank:
                taxid, species = genbank[assembly]
                out_fh.write(f"{fp}\t{assembly}\t{taxid}\t{species}\tgenbank\n")
            else:
                out_fh.write(f"{fp}\t{assembly}\tNA\tNA\tnot_found\n")


if __name__ == "__main__":
    main()
