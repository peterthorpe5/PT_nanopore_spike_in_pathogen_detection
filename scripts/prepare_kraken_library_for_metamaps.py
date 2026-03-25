#!/usr/bin/env python3
"""
Prepare Kraken library FASTA files for MetaMaps database building.

This script converts one or more Kraken `library.fna` files into a combined
FASTA with synthetic taxon IDs suitable for MetaMaps-style database building.

It creates an extended taxonomy directory by copying an existing NCBI taxonomy
directory and then appending synthetic child taxon IDs beneath the original
parent taxon IDs found in the Kraken FASTA headers.

This version uses a two-pass workflow and performs explicit validation to
ensure that every synthetic taxon ID written into the combined FASTA is also
present in both `nodes.dmp` and `names.dmp`.

Expected Kraken FASTA header format
-----------------------------------
>kraken:taxid|5855|NC_009906.1 Plasmodium vivax chromosome 1 ...

Outputs
-------
- combined_input.fa
- taxonomy/
- manifest.tsv
- excluded_records.tsv
"""

from __future__ import annotations

import argparse
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple


HEADER_RE = re.compile(r"^>(?:kraken:taxid\|(?P<taxid>\d+)\|)(?P<rest>.+)$")


@dataclass
class SyntheticRecord:
    """Container for a converted Kraken FASTA record."""

    source_fasta: str
    old_taxid: int
    new_taxid: int
    accession: str
    old_header: str
    new_header: str
    sequence: str
    synthetic_name: str


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
            "Convert Kraken library FASTA files into a MetaMaps-ready FASTA "
            "and extended taxonomy."
        )
    )
    parser.add_argument(
        "--input_fastas",
        nargs="+",
        required=True,
        help="One or more Kraken library FASTA files.",
    )
    parser.add_argument(
        "--taxonomy_dir",
        required=True,
        help="Input taxonomy directory containing nodes.dmp and names.dmp.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for combined FASTA, extended taxonomy, and logs.",
    )
    parser.add_argument(
        "--taxid_offset",
        type=int,
        default=1,
        help=(
            "Offset to add to the maximum existing taxon ID when allocating "
            "synthetic taxon IDs. Default: 1"
        ),
    )
    return parser.parse_args()


def require_file(file_path: Path) -> None:
    """
    Ensure a file exists.

    Parameters
    ----------
    file_path : Path
        File path to validate.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    if not file_path.is_file():
        raise FileNotFoundError(f"Required file not found: {file_path}")


def require_dir(dir_path: Path) -> None:
    """
    Ensure a directory exists.

    Parameters
    ----------
    dir_path : Path
        Directory path to validate.

    Raises
    ------
    FileNotFoundError
        If the directory does not exist.
    """
    if not dir_path.is_dir():
        raise FileNotFoundError(f"Required directory not found: {dir_path}")


def iter_fasta_records(fasta_path: Path) -> Iterator[Tuple[str, str]]:
    """
    Iterate through FASTA records.

    Parameters
    ----------
    fasta_path : Path
        Input FASTA file.

    Yields
    ------
    Iterator[Tuple[str, str]]
        Tuples of (header, sequence_string_without_wrapping).
    """
    header: Optional[str] = None
    seq_lines: List[str] = []

    with fasta_path.open(mode="r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        if header is not None:
            yield header, "".join(seq_lines)


def wrap_sequence(sequence: str, width: int = 80) -> str:
    """
    Wrap a sequence string to a fixed width.

    Parameters
    ----------
    sequence : str
        Sequence string.
    width : int, optional
        Wrap width, by default 80.

    Returns
    -------
    str
        Wrapped sequence with trailing newline.
    """
    return "\n".join(
        sequence[i:i + width] for i in range(0, len(sequence), width)
    ) + "\n"


def read_valid_taxids(nodes_path: Path) -> Tuple[set[int], int]:
    """
    Read valid taxon IDs and maximum existing taxon ID from nodes.dmp.

    Parameters
    ----------
    nodes_path : Path
        Path to taxonomy/nodes.dmp.

    Returns
    -------
    Tuple[set[int], int]
        Set of valid taxon IDs and maximum observed taxon ID.
    """
    valid_taxids: set[int] = set()
    max_taxid = 0

    with nodes_path.open(mode="r", encoding="utf-8") as handle:
        for line in handle:
            parts = [field.strip() for field in line.split("|")]
            if not parts or not parts[0]:
                continue
            taxid = int(parts[0])
            valid_taxids.add(taxid)
            if taxid > max_taxid:
                max_taxid = taxid

    return valid_taxids, max_taxid


def copy_taxonomy_dir(src_dir: Path, dst_dir: Path) -> None:
    """
    Copy taxonomy directory contents into a destination directory.

    Parameters
    ----------
    src_dir : Path
        Source taxonomy directory.
    dst_dir : Path
        Destination taxonomy directory.
    """
    if dst_dir.exists():
        shutil.rmtree(dst_dir)
    shutil.copytree(src=src_dir, dst=dst_dir)


def extract_accession(rest: str) -> str:
    """
    Extract a leading accession-like token from a Kraken header remainder.

    Parameters
    ----------
    rest : str
        Header text after the Kraken taxid block.

    Returns
    -------
    str
        Leading accession-like token, or a fallback label.
    """
    token = rest.split(maxsplit=1)[0].strip()
    return token if token else "unknown_record"


def build_synthetic_name(accession: str, parent_taxid: int, rest: str) -> str:
    """
    Build a scientific name for a synthetic child taxon.

    Parameters
    ----------
    accession : str
        Sequence accession token.
    parent_taxid : int
        Parent NCBI taxon ID.
    rest : str
        Original Kraken header remainder.

    Returns
    -------
    str
        Synthetic scientific name.
    """
    return (
        f"{accession} synthetic child of taxid {parent_taxid}: "
        f"{rest.strip()}"
    )


def prepare_outputs(out_dir: Path) -> Dict[str, Path]:
    """
    Prepare output directory structure.

    Parameters
    ----------
    out_dir : Path
        Base output directory.

    Returns
    -------
    Dict[str, Path]
        Dictionary of key output paths.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    return {
        "taxonomy_out": out_dir / "taxonomy",
        "combined_fasta": out_dir / "combined_input.fa",
        "manifest_tsv": out_dir / "manifest.tsv",
        "excluded_tsv": out_dir / "excluded_records.tsv",
    }


def append_node_line(new_taxid: int, parent_taxid: int) -> str:
    """
    Create a nodes.dmp line for a synthetic child taxon.

    Parameters
    ----------
    new_taxid : int
        Synthetic taxon ID.
    parent_taxid : int
        Parent taxon ID.

    Returns
    -------
    str
        Formatted nodes.dmp line.
    """
    return (
        f"{new_taxid}\t|\t{parent_taxid}\t|\tno rank\t|\t\t|\t0\t|\t0\t|"
        f"\t11\t|\t1\t|\t1\t|\t0\t|\t0\t|\t0\t|\t\t|\n"
    )


def append_name_line(new_taxid: int, scientific_name: str) -> str:
    """
    Create a names.dmp line for a synthetic scientific name.

    Parameters
    ----------
    new_taxid : int
        Synthetic taxon ID.
    scientific_name : str
        Scientific name.

    Returns
    -------
    str
        Formatted names.dmp line.
    """
    return (
        f"{new_taxid}\t|\t{scientific_name}\t|\t\t|\tscientific name\t|\n"
    )


def write_combined_fasta(
    combined_fasta: Path,
    records: List[SyntheticRecord],
) -> None:
    """
    Write the combined FASTA file.

    Parameters
    ----------
    combined_fasta : Path
        Output FASTA path.
    records : List[SyntheticRecord]
        Converted records.
    """
    with combined_fasta.open(mode="w", encoding="utf-8") as handle:
        for record in records:
            handle.write(record.new_header + "\n")
            handle.write(wrap_sequence(sequence=record.sequence))


def write_manifest(
    manifest_tsv: Path,
    records: List[SyntheticRecord],
) -> None:
    """
    Write the manifest TSV.

    Parameters
    ----------
    manifest_tsv : Path
        Output manifest path.
    records : List[SyntheticRecord]
        Converted records.
    """
    with manifest_tsv.open(mode="w", encoding="utf-8") as handle:
        handle.write(
            "source_fasta\told_taxid\tnew_taxid\taccession\told_header\tnew_header\n"
        )
        for record in records:
            handle.write(
                f"{record.source_fasta}\t{record.old_taxid}\t{record.new_taxid}\t"
                f"{record.accession}\t{record.old_header}\t{record.new_header}\n"
            )


def write_excluded(
    excluded_tsv: Path,
    excluded_rows: List[Tuple[str, str, str, str]],
) -> None:
    """
    Write the excluded-record log TSV.

    Parameters
    ----------
    excluded_tsv : Path
        Output excluded-record log path.
    excluded_rows : List[Tuple[str, str, str, str]]
        Rows in the form (source_fasta, reason, old_taxid, old_header).
    """
    with excluded_tsv.open(mode="w", encoding="utf-8") as handle:
        handle.write("source_fasta\treason\told_taxid\told_header\n")
        for row in excluded_rows:
            handle.write("\t".join(row) + "\n")


def extend_taxonomy(
    taxonomy_out: Path,
    records: List[SyntheticRecord],
) -> None:
    """
    Extend copied taxonomy files with synthetic nodes and names.

    Parameters
    ----------
    taxonomy_out : Path
        Output taxonomy directory.
    records : List[SyntheticRecord]
        Converted records.
    """
    nodes_out = taxonomy_out / "nodes.dmp"
    names_out = taxonomy_out / "names.dmp"

    with nodes_out.open(mode="a", encoding="utf-8") as nodes_handle:
        for record in records:
            nodes_handle.write(
                append_node_line(
                    new_taxid=record.new_taxid,
                    parent_taxid=record.old_taxid,
                )
            )

    with names_out.open(mode="a", encoding="utf-8") as names_handle:
        for record in records:
            names_handle.write(
                append_name_line(
                    new_taxid=record.new_taxid,
                    scientific_name=record.synthetic_name,
                )
            )


def read_taxids_from_nodes(nodes_path: Path) -> set[int]:
    """
    Read all taxon IDs present in nodes.dmp.

    Parameters
    ----------
    nodes_path : Path
        Path to nodes.dmp.

    Returns
    -------
    set[int]
        Taxon IDs present in nodes.dmp.
    """
    taxids: set[int] = set()
    with nodes_path.open(mode="r", encoding="utf-8") as handle:
        for line in handle:
            parts = [field.strip() for field in line.split("|")]
            if parts and parts[0]:
                taxids.add(int(parts[0]))
    return taxids


def read_taxids_from_names(names_path: Path) -> set[int]:
    """
    Read all taxon IDs present in names.dmp.

    Parameters
    ----------
    names_path : Path
        Path to names.dmp.

    Returns
    -------
    set[int]
        Taxon IDs present in names.dmp.
    """
    taxids: set[int] = set()
    with names_path.open(mode="r", encoding="utf-8") as handle:
        for line in handle:
            parts = [field.strip() for field in line.split("|")]
            if parts and parts[0]:
                taxids.add(int(parts[0]))
    return taxids


def read_taxids_from_combined_fasta(fasta_path: Path) -> set[int]:
    """
    Read all taxon IDs present in the combined FASTA headers.

    Parameters
    ----------
    fasta_path : Path
        Combined FASTA path.

    Returns
    -------
    set[int]
        Taxon IDs found in FASTA headers.
    """
    taxids: set[int] = set()
    with fasta_path.open(mode="r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            match = HEADER_RE.match(line.rstrip("\n"))
            if match:
                taxids.add(int(match.group("taxid")))
    return taxids


def validate_outputs(
    taxonomy_out: Path,
    combined_fasta: Path,
    records: List[SyntheticRecord],
) -> None:
    """
    Validate consistency across FASTA, nodes.dmp, names.dmp, and manifest data.

    Parameters
    ----------
    taxonomy_out : Path
        Output taxonomy directory.
    combined_fasta : Path
        Combined FASTA path.
    records : List[SyntheticRecord]
        Converted records.

    Raises
    ------
    RuntimeError
        If any synthetic taxon IDs are missing from outputs.
    """
    nodes_taxids = read_taxids_from_nodes(nodes_path=taxonomy_out / "nodes.dmp")
    names_taxids = read_taxids_from_names(names_path=taxonomy_out / "names.dmp")
    fasta_taxids = read_taxids_from_combined_fasta(fasta_path=combined_fasta)
    manifest_taxids = {record.new_taxid for record in records}

    missing_in_nodes = sorted(fasta_taxids - nodes_taxids)
    missing_in_names = sorted(fasta_taxids - names_taxids)
    missing_manifest_in_nodes = sorted(manifest_taxids - nodes_taxids)
    missing_manifest_in_names = sorted(manifest_taxids - names_taxids)

    if (
        missing_in_nodes
        or missing_in_names
        or missing_manifest_in_nodes
        or missing_manifest_in_names
    ):
        error_parts: List[str] = []
        if missing_in_nodes:
            error_parts.append(
                f"Missing FASTA taxids in nodes.dmp: {missing_in_nodes[:10]}"
            )
        if missing_in_names:
            error_parts.append(
                f"Missing FASTA taxids in names.dmp: {missing_in_names[:10]}"
            )
        if missing_manifest_in_nodes:
            error_parts.append(
                "Missing manifest taxids in nodes.dmp: "
                f"{missing_manifest_in_nodes[:10]}"
            )
        if missing_manifest_in_names:
            error_parts.append(
                "Missing manifest taxids in names.dmp: "
                f"{missing_manifest_in_names[:10]}"
            )
        raise RuntimeError("Validation failed. " + " | ".join(error_parts))


def main() -> None:
    """
    Run Kraken-to-MetaMaps preparation.
    """
    args = parse_args()

    input_fastas = [Path(path) for path in args.input_fastas]
    taxonomy_dir = Path(args.taxonomy_dir)
    out_dir = Path(args.out_dir)

    for fasta_path in input_fastas:
        require_file(file_path=fasta_path)

    require_dir(dir_path=taxonomy_dir)
    require_file(file_path=taxonomy_dir / "nodes.dmp")
    require_file(file_path=taxonomy_dir / "names.dmp")

    outputs = prepare_outputs(out_dir=out_dir)
    copy_taxonomy_dir(
        src_dir=taxonomy_dir,
        dst_dir=outputs["taxonomy_out"],
    )

    valid_taxids, max_taxid = read_valid_taxids(
        nodes_path=taxonomy_dir / "nodes.dmp"
    )
    next_taxid = max_taxid + args.taxid_offset

    records: List[SyntheticRecord] = []
    excluded_rows: List[Tuple[str, str, str, str]] = []

    total_records = 0

    for fasta_path in input_fastas:
        print(f"[INFO] Processing {fasta_path}")
        for header, sequence in iter_fasta_records(fasta_path=fasta_path):
            total_records += 1

            match = HEADER_RE.match(header)
            if not match:
                excluded_rows.append(
                    (str(fasta_path), "invalid_header_format", "", header)
                )
                continue

            old_taxid = int(match.group("taxid"))
            rest = match.group("rest").strip()

            if old_taxid not in valid_taxids:
                excluded_rows.append(
                    (
                        str(fasta_path),
                        "missing_taxid_in_taxonomy",
                        str(old_taxid),
                        header,
                    )
                )
                continue

            accession = extract_accession(rest=rest)
            new_taxid = next_taxid
            next_taxid += 1

            synthetic_name = build_synthetic_name(
                accession=accession,
                parent_taxid=old_taxid,
                rest=rest,
            )
            new_header = f">kraken:taxid|{new_taxid}|{accession}"
            #new_header = f">kraken:taxid|{new_taxid}|{rest}"

            records.append(
                SyntheticRecord(
                    source_fasta=str(fasta_path),
                    old_taxid=old_taxid,
                    new_taxid=new_taxid,
                    accession=accession,
                    old_header=header,
                    new_header=new_header,
                    sequence=sequence,
                    synthetic_name=synthetic_name,
                )
            )

    write_combined_fasta(
        combined_fasta=outputs["combined_fasta"],
        records=records,
    )
    write_manifest(
        manifest_tsv=outputs["manifest_tsv"],
        records=records,
    )
    write_excluded(
        excluded_tsv=outputs["excluded_tsv"],
        excluded_rows=excluded_rows,
    )
    extend_taxonomy(
        taxonomy_out=outputs["taxonomy_out"],
        records=records,
    )
    validate_outputs(
        taxonomy_out=outputs["taxonomy_out"],
        combined_fasta=outputs["combined_fasta"],
        records=records,
    )

    print(f"[INFO] Total input records: {total_records}")
    print(f"[INFO] Kept records: {len(records)}")
    print(f"[INFO] Excluded records: {len(excluded_rows)}")
    print(f"[INFO] Combined FASTA: {outputs['combined_fasta']}")
    print(f"[INFO] Extended taxonomy: {outputs['taxonomy_out']}")
    print(f"[INFO] Manifest: {outputs['manifest_tsv']}")
    print(f"[INFO] Excluded log: {outputs['excluded_tsv']}")
    print("[INFO] Validation passed")


if __name__ == "__main__":
    main()