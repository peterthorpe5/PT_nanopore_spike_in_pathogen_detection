#!/usr/bin/env python3
"""Build a validated shared minimap2 reference FASTA.

This script merges the three Kraken-derived source FASTA files used in the
current ONT spike-in benchmarking project into one normalised FASTA suitable
for minimap2.

The script is designed to be conservative and streaming-based so it can cope
with very large input files. It:

- reads each input FASTA record by record
- skips malformed preamble text before the first header
- removes internal whitespace from sequence lines
- converts sequences to upper case
- keeps only records with valid nucleotide IUPAC characters
- renames duplicate reference names so minimap2 receives unique targets
- writes a tab-separated discard log and summary table
- writes to a temporary output and renames only on success
- can optionally build a minimap2 index after writing the FASTA

All tabular outputs are written as tab-separated files.
"""

from __future__ import annotations

import argparse
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Optional, TextIO


VALID_SEQ_RE = re.compile(r"^[ACGTURYSWKMBDHVN\-]+$")
WHITESPACE_RE = re.compile(r"\s+")


@dataclass
class FastaRecord:
    """Represent one parsed FASTA record.

    Attributes
    ----------
    source_path : Path
        Source FASTA file path.
    header : str
        Full FASTA header, including the leading ``>``.
    sequence : str
        Sequence text after whitespace removal and upper-casing.
    """

    source_path: Path
    header: str
    sequence: str


@dataclass
class SourceStats:
    """Track per-source FASTA statistics.

    Attributes
    ----------
    source_path : Path
        Source FASTA path.
    records_seen : int
        Number of headered records encountered.
    records_kept : int
        Number of records written to output.
    records_discarded : int
        Number of records discarded.
    preamble_lines_discarded : int
        Number of non-header lines seen before the first header.
    duplicate_names_renamed : int
        Number of duplicate reference names renamed.
    """

    source_path: Path
    records_seen: int = 0
    records_kept: int = 0
    records_discarded: int = 0
    preamble_lines_discarded: int = 0
    duplicate_names_renamed: int = 0


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Merge the Kraken source FASTA files into a validated minimap2 "
            "reference FASTA."
        )
    )
    parser.add_argument(
        "--custom_fasta",
        default=(
            "/home/pthorpe001/data/project_back_up_2024/"
            "kraken_bact_virus_plasmo_fungal/library/library.fna"
        ),
        help="Path to the top-level custom/library FASTA.",
    )
    parser.add_argument(
        "--viral_fasta",
        default=(
            "/home/pthorpe001/data/project_back_up_2024/"
            "kraken_bact_virus_plasmo_fungal/library/viral/library.fna"
        ),
        help="Path to the viral library FASTA.",
    )
    parser.add_argument(
        "--bacterial_fasta",
        default=(
            "/home/pthorpe001/data/project_back_up_2024/"
            "kraken_bact_virus_plasmo_fungal/library/bacteria/library.fna"
        ),
        help="Path to the bacterial library FASTA.",
    )
    parser.add_argument(
        "--output_fasta",
        default=(
            "/home/pthorpe001/data/databases/minimap2/"
            "shared_bact_viral_plasmo_refs.fa"
        ),
        help="Output FASTA path.",
    )
    parser.add_argument(
        "--discarded_tsv",
        default=(
            "/home/pthorpe001/data/databases/minimap2/"
            "shared_bact_viral_plasmo_refs.discarded.tsv"
        ),
        help="TSV path for discarded records and reasons.",
    )
    parser.add_argument(
        "--summary_tsv",
        default=(
            "/home/pthorpe001/data/databases/minimap2/"
            "shared_bact_viral_plasmo_refs.summary.tsv"
        ),
        help="TSV path for the overall summary.",
    )
    parser.add_argument(
        "--output_mmi",
        default=(
            "/home/pthorpe001/data/databases/minimap2/"
            "shared_bact_viral_plasmo_refs.mmi"
        ),
        help="Optional minimap2 index path.",
    )
    parser.add_argument(
        "--build_index",
        action="store_true",
        help="Build a minimap2 index after writing the FASTA.",
    )
    parser.add_argument(
        "--minimap2_exe",
        default="minimap2",
        help="Minimap2 executable name or absolute path.",
    )
    return parser.parse_args()


def ensure_parent_dir(*, path: Path) -> None:
    """Create the parent directory for a path if needed.

    Parameters
    ----------
    path : Path
        Target path whose parent directory should exist.
    """
    path.parent.mkdir(parents=True, exist_ok=True)


def split_header_name(*, header: str) -> tuple[str, str]:
    """Split a FASTA header into reference name and remainder.

    Parameters
    ----------
    header : str
        Header line including the leading ``>``.

    Returns
    -------
    tuple[str, str]
        Reference name and the remainder of the header text.
    """
    header_body = header[1:].strip()
    if not header_body:
        return "unnamed_record", ""

    parts = header_body.split(maxsplit=1)
    name = parts[0]
    remainder = parts[1] if len(parts) > 1 else ""
    return name, remainder


def make_unique_header(
    *,
    header: str,
    used_names: dict[str, int],
) -> tuple[str, bool]:
    """Make a FASTA header name unique if necessary.

    Parameters
    ----------
    header : str
        Input FASTA header including the leading ``>``.
    used_names : dict[str, int]
        Mapping of already-used reference names to counts.

    Returns
    -------
    tuple[str, bool]
        Updated header and a flag indicating whether renaming occurred.
    """
    name, remainder = split_header_name(header=header)

    if name not in used_names:
        used_names[name] = 1
        return header, False

    used_names[name] += 1
    new_name = f"{name}__dup{used_names[name]}"
    if remainder:
        return f">{new_name} {remainder}", True
    return f">{new_name}", True


def write_discard_record(
    *,
    handle: TextIO,
    source_path: Path,
    record_index: int,
    header: str,
    reason: str,
) -> None:
    """Write one discarded-entry row.

    Parameters
    ----------
    handle : TextIO
        Open TSV handle.
    source_path : Path
        Source FASTA path.
    record_index : int
        Running record index within the source.
    header : str
        FASTA header or placeholder text.
    reason : str
        Discard reason.
    """
    handle.write(
        f"{source_path}\t{record_index}\t{header}\t{reason}\n"
    )


def iter_fasta_records(
    *,
    source_path: Path,
    discard_handle: TextIO,
    stats: SourceStats,
) -> Iterator[FastaRecord]:
    """Yield valid FASTA-like records from a source file.

    Parameters
    ----------
    source_path : Path
        Input FASTA path.
    discard_handle : TextIO
        TSV handle for discard logging.
    stats : SourceStats
        Mutable statistics tracker.

    Yields
    ------
    FastaRecord
        Parsed FASTA record with cleaned sequence.
    """
    current_header: Optional[str] = None
    current_sequence_chunks: list[str] = []

    def flush_current_record() -> Optional[FastaRecord]:
        """Validate and package the current record.

        Returns
        -------
        Optional[FastaRecord]
            Parsed record if valid, otherwise ``None``.
        """
        if current_header is None:
            return None

        if not current_sequence_chunks:
            stats.records_discarded += 1
            write_discard_record(
                handle=discard_handle,
                source_path=source_path,
                record_index=stats.records_seen,
                header=current_header,
                reason="empty_sequence",
            )
            return None

        sequence = "".join(current_sequence_chunks).upper()
        if not VALID_SEQ_RE.match(sequence):
            stats.records_discarded += 1
            write_discard_record(
                handle=discard_handle,
                source_path=source_path,
                record_index=stats.records_seen,
                header=current_header,
                reason="invalid_sequence_characters",
            )
            return None

        return FastaRecord(
            source_path=source_path,
            header=current_header,
            sequence=sequence,
        )

    with source_path.open(mode="r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            stripped = line.strip()

            if not stripped:
                continue

            if stripped.startswith(">"):  # FASTA header
                if current_header is not None:
                    record = flush_current_record()
                    if record is not None:
                        yield record

                stats.records_seen += 1
                current_header = stripped
                current_sequence_chunks = []
                continue

            if current_header is None:
                stats.preamble_lines_discarded += 1
                write_discard_record(
                    handle=discard_handle,
                    source_path=source_path,
                    record_index=0,
                    header="PREAMBLE",
                    reason=(
                        "non_header_before_first_record:"
                        f"{stripped[:120]}"
                    ),
                )
                continue

            cleaned = WHITESPACE_RE.sub("", stripped).upper()
            current_sequence_chunks.append(cleaned)

    if current_header is not None:
        record = flush_current_record()
        if record is not None:
            yield record


def write_wrapped_sequence(*, handle: TextIO, sequence: str, width: int = 80) -> None:
    """Write a sequence in wrapped FASTA format.

    Parameters
    ----------
    handle : TextIO
        Open FASTA handle.
    sequence : str
        Sequence text.
    width : int, optional
        Line width for wrapping, by default 80.
    """
    for start in range(0, len(sequence), width):
        handle.write(f"{sequence[start:start + width]}\n")


def merge_fastas(
    *,
    source_paths: Iterable[Path],
    output_fasta: Path,
    discarded_tsv: Path,
    summary_tsv: Path,
) -> dict[str, int]:
    """Merge FASTA sources into one validated output FASTA.

    Parameters
    ----------
    source_paths : Iterable[Path]
        Source FASTA paths in the desired merge order.
    output_fasta : Path
        Output merged FASTA path.
    discarded_tsv : Path
        TSV path for discard logging.
    summary_tsv : Path
        TSV path for summary metrics.

    Returns
    -------
    dict[str, int]
        Overall summary metrics.
    """
    ensure_parent_dir(path=output_fasta)
    ensure_parent_dir(path=discarded_tsv)
    ensure_parent_dir(path=summary_tsv)

    tmp_output = output_fasta.with_suffix(output_fasta.suffix + ".tmp")
    used_names: dict[str, int] = {}
    per_source_stats: list[SourceStats] = []

    total_records_seen = 0
    total_records_kept = 0
    total_records_discarded = 0
    total_preamble_lines_discarded = 0
    total_duplicate_names_renamed = 0

    with (
        tmp_output.open(mode="w", encoding="utf-8") as out_handle,
        discarded_tsv.open(mode="w", encoding="utf-8") as discard_handle,
    ):
        discard_handle.write("source_path\trecord_index\theader\treason\n")

        for source_path in source_paths:
            if not source_path.exists():
                raise FileNotFoundError(f"Required source FASTA not found: {source_path}")

            stats = SourceStats(source_path=source_path)
            per_source_stats.append(stats)

            for record in iter_fasta_records(
                source_path=source_path,
                discard_handle=discard_handle,
                stats=stats,
            ):
                updated_header, renamed = make_unique_header(
                    header=record.header,
                    used_names=used_names,
                )
                if renamed:
                    stats.duplicate_names_renamed += 1

                out_handle.write(f"{updated_header}\n")
                write_wrapped_sequence(handle=out_handle, sequence=record.sequence)
                stats.records_kept += 1

            total_records_seen += stats.records_seen
            total_records_kept += stats.records_kept
            total_records_discarded += stats.records_discarded
            total_preamble_lines_discarded += stats.preamble_lines_discarded
            total_duplicate_names_renamed += stats.duplicate_names_renamed

    tmp_output.replace(output_fasta)

    with summary_tsv.open(mode="w", encoding="utf-8") as summary_handle:
        summary_handle.write("scope\tmetric\tvalue\n")
        summary_handle.write(f"overall\toutput_fasta\t{output_fasta}\n")
        summary_handle.write(f"overall\trecords_seen\t{total_records_seen}\n")
        summary_handle.write(f"overall\trecords_kept\t{total_records_kept}\n")
        summary_handle.write(
            f"overall\trecords_discarded\t{total_records_discarded}\n"
        )
        summary_handle.write(
            "overall\tpreamble_lines_discarded"
            f"\t{total_preamble_lines_discarded}\n"
        )
        summary_handle.write(
            "overall\tduplicate_names_renamed"
            f"\t{total_duplicate_names_renamed}\n"
        )

        for stats in per_source_stats:
            summary_handle.write(
                f"{stats.source_path}\trecords_seen\t{stats.records_seen}\n"
            )
            summary_handle.write(
                f"{stats.source_path}\trecords_kept\t{stats.records_kept}\n"
            )
            summary_handle.write(
                f"{stats.source_path}\trecords_discarded\t{stats.records_discarded}\n"
            )
            summary_handle.write(
                "{0}\tpreamble_lines_discarded\t{1}\n".format(
                    stats.source_path,
                    stats.preamble_lines_discarded,
                )
            )
            summary_handle.write(
                "{0}\tduplicate_names_renamed\t{1}\n".format(
                    stats.source_path,
                    stats.duplicate_names_renamed,
                )
            )

    return {
        "records_seen": total_records_seen,
        "records_kept": total_records_kept,
        "records_discarded": total_records_discarded,
        "preamble_lines_discarded": total_preamble_lines_discarded,
        "duplicate_names_renamed": total_duplicate_names_renamed,
    }


def build_minimap_index(
    *,
    minimap2_exe: str,
    output_fasta: Path,
    output_mmi: Path,
) -> None:
    """Build a minimap2 index for the merged FASTA.

    Parameters
    ----------
    minimap2_exe : str
        Minimap2 executable name or path.
    output_fasta : Path
        Merged FASTA path.
    output_mmi : Path
        Output minimap2 index path.
    """
    ensure_parent_dir(path=output_mmi)
    subprocess.run(
        [
            minimap2_exe,
            "-d",
            str(output_mmi),
            str(output_fasta),
        ],
        check=True,
    )


def main() -> None:
    """Run FASTA merging and optional minimap2 indexing."""
    args = parse_args()

    source_paths = [
        Path(args.custom_fasta).expanduser().resolve(),
        Path(args.viral_fasta).expanduser().resolve(),
        Path(args.bacterial_fasta).expanduser().resolve(),
    ]
    output_fasta = Path(args.output_fasta).expanduser().resolve()
    discarded_tsv = Path(args.discarded_tsv).expanduser().resolve()
    summary_tsv = Path(args.summary_tsv).expanduser().resolve()
    output_mmi = Path(args.output_mmi).expanduser().resolve()

    summary = merge_fastas(
        source_paths=source_paths,
        output_fasta=output_fasta,
        discarded_tsv=discarded_tsv,
        summary_tsv=summary_tsv,
    )

    if summary["records_kept"] == 0:
        raise RuntimeError(
            "No valid FASTA records were kept. Output FASTA should not be used."
        )

    if args.build_index:
        build_minimap_index(
            minimap2_exe=args.minimap2_exe,
            output_fasta=output_fasta,
            output_mmi=output_mmi,
        )

    print(f"[INFO] Output FASTA: {output_fasta}")
    print(f"[INFO] Discard log: {discarded_tsv}")
    print(f"[INFO] Summary TSV: {summary_tsv}")
    print(f"[INFO] Records kept: {summary['records_kept']}")
    print(f"[INFO] Records discarded: {summary['records_discarded']}")
    print(
        "[INFO] Preamble lines discarded: "
        f"{summary['preamble_lines_discarded']}"
    )
    print(
        "[INFO] Duplicate names renamed: "
        f"{summary['duplicate_names_renamed']}"
    )
    if args.build_index:
        print(f"[INFO] Minimap2 index: {output_mmi}")


if __name__ == "__main__":
    main()
