#!/usr/bin/env python3
"""
Convert Kraken-style library FASTA files into Metabuli-ready inputs.

This script merges one or more Kraken-style FASTA files into a single FASTA
suitable for Metabuli custom database building, and writes a matching
NCBI-style accession2taxid table.

Expected Kraken-style headers
-----------------------------
Examples:
>kraken:taxid|5855|NC_009906.1 Plasmodium vivax chromosome 1
>kraken:taxid|12345|CP001849.1 some description
>kraken:taxid|67890|NZ_ABC12345.1

Metabuli-ready outputs
----------------------
1. merged_metabuli_input.fa
   FASTA headers rewritten as:
   >ACCESSION.VERSION

2. accession2taxid.tsv
   Four-column NCBI-style mapping:
   accession<TAB>accession.version<TAB>taxid<TAB>gi

3. fasta_list.txt
   A one-line file containing the absolute path to the merged FASTA.

4. manifest.tsv
   Detailed record-level mapping from original headers to rewritten headers.

5. skipped_records.tsv
   Records skipped because the header could not be parsed or an accession
   could not be extracted.

Notes
-----
- Duplicate accession.version identifiers are only written once.
- The first occurrence is kept and later duplicates are logged.
- Output files are tab-separated.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from typing import Iterator, Optional


HEADER_RE = re.compile(
    r"^>kraken:taxid\|(?P<taxid>\d+)\|(?P<rest>.+)$"
)
ACCESSION_RE = re.compile(
    r"^(?P<accession>[A-Za-z0-9_]+)(?:\.(?P<version>\d+))?$"
)


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
            "Merge Kraken-style FASTA libraries and generate Metabuli-ready "
            "FASTA and accession2taxid files."
        )
    )
    parser.add_argument(
        "--input_fastas",
        nargs="+",
        required=True,
        help="One or more Kraken-style FASTA or FASTA.GZ input files.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for Metabuli-ready files.",
    )
    parser.add_argument(
        "--merged_fasta_name",
        default="merged_metabuli_input.fa",
        help="Output FASTA filename. Default: merged_metabuli_input.fa",
    )
    parser.add_argument(
        "--accession2taxid_name",
        default="accession2taxid.tsv",
        help="Output accession2taxid filename. Default: accession2taxid.tsv",
    )
    parser.add_argument(
        "--fasta_list_name",
        default="fasta_list.txt",
        help="Output FASTA list filename. Default: fasta_list.txt",
    )
    parser.add_argument(
        "--manifest_name",
        default="manifest.tsv",
        help="Output manifest filename. Default: manifest.tsv",
    )
    parser.add_argument(
        "--skipped_name",
        default="skipped_records.tsv",
        help="Output skipped-records filename. Default: skipped_records.tsv",
    )
    return parser.parse_args()


def open_maybe_gzip(file_path: Path):
    """
    Open a plain text or gzip-compressed file.

    Parameters
    ----------
    file_path : Path
        Input file path.

    Returns
    -------
    IO[str]
        Text-mode file handle.
    """
    if str(file_path).endswith(".gz"):
        return gzip.open(file_path, mode="rt", encoding="utf-8", errors="replace")
    return open(file_path, mode="rt", encoding="utf-8", errors="replace")


def iter_fasta_records(fasta_path: Path) -> Iterator[tuple[str, str]]:
    """
    Iterate over FASTA records.

    Parameters
    ----------
    fasta_path : Path
        Input FASTA path.

    Yields
    ------
    Iterator[tuple[str, str]]
        Header and sequence string.
    """
    header: Optional[str] = None
    seq_lines: list[str] = []

    with open_maybe_gzip(fasta_path) as handle:
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
    Wrap a sequence string to a fixed line width.

    Parameters
    ----------
    sequence : str
        Sequence string.
    width : int, optional
        Output line width, by default 80.

    Returns
    -------
    str
        Wrapped sequence string.
    """
    return "\n".join(
        sequence[i:i + width] for i in range(0, len(sequence), width)
    )


def extract_accession_token(rest: str) -> Optional[str]:
    """
    Extract the accession token from the remainder of a Kraken header.

    Parameters
    ----------
    rest : str
        Header text after the Kraken taxid block.

    Returns
    -------
    Optional[str]
        Accession token if found, otherwise None.
    """
    token = rest.split(maxsplit=1)[0].strip()
    if not token:
        return None
    return token


def split_accession(accession_token: str) -> tuple[str, str]:
    """
    Split an accession token into accession and accession.version values.

    Parameters
    ----------
    accession_token : str
        Accession token, for example NC_009906.1.

    Returns
    -------
    tuple[str, str]
        accession, accession.version
    """
    match = ACCESSION_RE.match(accession_token)
    if match is None:
        return accession_token, accession_token

    accession = match.group("accession")
    version = match.group("version")

    if version is None:
        accession_version = accession
    else:
        accession_version = f"{accession}.{version}"

    return accession, accession_version


def main() -> None:
    """
    Run the Kraken-to-Metabuli conversion workflow.
    """
    args = parse_args()

    input_fastas = [Path(path).expanduser().resolve() for path in args.input_fastas]
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    for fasta_path in input_fastas:
        if not fasta_path.is_file():
            raise FileNotFoundError(f"Required input FASTA not found: {fasta_path}")

    merged_fasta = out_dir / args.merged_fasta_name
    accession2taxid_tsv = out_dir / args.accession2taxid_name
    fasta_list_txt = out_dir / args.fasta_list_name
    manifest_tsv = out_dir / args.manifest_name
    skipped_tsv = out_dir / args.skipped_name

    seen_accessions: set[str] = set()

    n_total = 0
    n_written = 0
    n_skipped = 0
    n_duplicate = 0

    with open(merged_fasta, mode="wt", encoding="utf-8") as fasta_handle, \
        open(accession2taxid_tsv, mode="wt", encoding="utf-8") as map_handle, \
        open(manifest_tsv, mode="wt", encoding="utf-8") as manifest_handle, \
        open(skipped_tsv, mode="wt", encoding="utf-8") as skipped_handle:

        manifest_handle.write(
            "source_fasta\told_header\tnew_header\ttaxid\taccession\t"
            "accession_version\n"
        )
        skipped_handle.write(
            "source_fasta\treason\told_header\n"
        )

        for fasta_path in input_fastas:
            for header, sequence in iter_fasta_records(fasta_path=fasta_path):
                n_total += 1

                header_match = HEADER_RE.match(header)
                if header_match is None:
                    skipped_handle.write(
                        f"{fasta_path}\tinvalid_header_format\t{header}\n"
                    )
                    n_skipped += 1
                    continue

                taxid = header_match.group("taxid")
                rest = header_match.group("rest").strip()

                accession_token = extract_accession_token(rest=rest)
                if accession_token is None:
                    skipped_handle.write(
                        f"{fasta_path}\tmissing_accession_token\t{header}\n"
                    )
                    n_skipped += 1
                    continue

                accession, accession_version = split_accession(
                    accession_token=accession_token
                )
                if not accession_version:
                    skipped_handle.write(
                        f"{fasta_path}\tinvalid_accession\t{header}\n"
                    )
                    n_skipped += 1
                    continue

                if accession_version in seen_accessions:
                    skipped_handle.write(
                        f"{fasta_path}\tduplicate_accession_version\t{header}\n"
                    )
                    n_duplicate += 1
                    continue

                seen_accessions.add(accession_version)

                new_header = f">{accession_version}"
                fasta_handle.write(new_header + "\n")
                fasta_handle.write(wrap_sequence(sequence=sequence) + "\n")

                map_handle.write(
                    f"{accession}\t{accession_version}\t{taxid}\t0\n"
                )

                manifest_handle.write(
                    f"{fasta_path}\t{header}\t{new_header}\t{taxid}\t"
                    f"{accession}\t{accession_version}\n"
                )

                n_written += 1

    with open(fasta_list_txt, mode="wt", encoding="utf-8") as handle:
        handle.write(str(merged_fasta) + "\n")

    print(f"[INFO] Input FASTA files: {len(input_fastas)}")
    print(f"[INFO] Total records seen: {n_total}")
    print(f"[INFO] Records written: {n_written}")
    print(f"[INFO] Duplicates skipped: {n_duplicate}")
    print(f"[INFO] Other skipped records: {n_skipped}")
    print(f"[INFO] Merged FASTA: {merged_fasta}")
    print(f"[INFO] accession2taxid: {accession2taxid_tsv}")
    print(f"[INFO] FASTA list: {fasta_list_txt}")
    print(f"[INFO] Manifest: {manifest_tsv}")
    print(f"[INFO] Skipped log: {skipped_tsv}")


if __name__ == "__main__":
    main()
