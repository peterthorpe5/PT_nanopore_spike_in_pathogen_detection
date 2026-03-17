#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N SENS_shuffle_ctrl

set -euo pipefail

# run_spikein_shuffled_control.sh
#
# Negative control for the single-genome spike-in workflow.
#
# Strategy:
#   1. Shuffle each pathogen contig independently (mononucleotide shuffle)
#   2. Simulate ONT-like reads from the shuffled reference using the same
#      NanoSim model as the real-pathogen workflow
#   3. Run the standard one-sample spike-in pipeline using the shuffled FASTA
#
# Notes:
#   - Headers in the shuffled FASTA are deliberately shortened to reduce the
#     risk of very long downstream read names.
#   - The script performs gzip integrity checks and writes provenance metadata.
#   - The main pipeline is expected to be run_spikein_one_sample.sh.

###############################################################################
# User-configurable inputs
###############################################################################

REAL_FASTQ="${REAL_FASTQ:-/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC1023_AmM008WB.fastq.gz}"

MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna.gz}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna}"

MONKEY_DEPLETION_FASTA="${MONKEY_DEPLETION_FASTA:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/M.nemestrina_tonkeana_nigra.fasta}"
PLASMO_MASKED_GZ="${PLASMO_MASKED_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"

DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/host_plus_plasmo_depletion.fasta}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/MRC1023_AmM008WB.depleted.fastq.gz}"

PATHOGEN_FASTA="${PATHOGEN_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/GCA_014843685.1_ASM1484368v1_genomic.fna}"
TARGET_LABEL="${TARGET_LABEL:-Plasmodium vivax}"

MINIMAP_DB_GZ="${MINIMAP_DB_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"
KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"

SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
MAIN_SCRIPT="${MAIN_SCRIPT:-${SCRIPT_DIR}/run_spikein_one_sample.sh}"

OUT_DIR="${OUT_DIR:-spikein_shuffled_control_out}"
THREADS="${THREADS:-12}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-3}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"

SHUFFLE_SEED="${SHUFFLE_SEED:-123}"
SHUFFLE_MODE="${SHUFFLE_MODE:-mononucleotide}"

###############################################################################
# Helper functions
###############################################################################

log_info() {
    printf '[INFO] %s\n' "$*"
}

log_warn() {
    printf '[WARN] %s\n' "$*" >&2
}

log_error() {
    printf '[ERROR] %s\n' "$*" >&2
}

require_file() {
    local file_path="$1"
    if [[ ! -f "${file_path}" ]]; then
        log_error "Required file not found: ${file_path}"
        exit 1
    fi
}

require_dir() {
    local dir_path="$1"
    if [[ ! -d "${dir_path}" ]]; then
        log_error "Required directory not found: ${dir_path}"
        exit 1
    fi
}

require_exe() {
    local exe_name="$1"
    if ! command -v "${exe_name}" >/dev/null 2>&1; then
        log_error "Required executable not found on PATH: ${exe_name}"
        exit 1
    fi
}

gzip_test_if_gz() {
    local file_path="$1"
    if [[ "${file_path}" == *.gz ]]; then
        log_info "Testing gzip integrity: ${file_path}"
        gzip -t "${file_path}"
    fi
}

###############################################################################
# Checks
###############################################################################

mkdir -p "${OUT_DIR}"

require_file "${REAL_FASTQ}"
require_file "${PATHOGEN_FASTA}"
require_file "${MONKEY_DEPLETION_FASTA}"
require_file "${PLASMO_MASKED_GZ}"
require_file "${MINIMAP_DB_GZ}"
require_file "${MAIN_SCRIPT}"
require_dir "${KRAKEN_DB_DIR}"

if [[ ! -f "${MONKEY_SMALL_GZ}" && ! -f "${MONKEY_SMALL_FASTA}" ]]; then
    log_error "Neither MONKEY_SMALL_GZ nor MONKEY_SMALL_FASTA exists."
    exit 1
fi

require_exe python3
require_exe gzip
require_exe bash

gzip_test_if_gz "${REAL_FASTQ}"
gzip_test_if_gz "${PATHOGEN_FASTA}"
gzip_test_if_gz "${PLASMO_MASKED_GZ}"
gzip_test_if_gz "${MINIMAP_DB_GZ}"
if [[ -f "${MONKEY_SMALL_GZ}" ]]; then
    gzip_test_if_gz "${MONKEY_SMALL_GZ}"
fi

###############################################################################
# Output paths
###############################################################################

SHUFFLED_FASTA="${OUT_DIR}/pathogen.shuffled.fasta"
SHUFFLE_META_TSV="${OUT_DIR}/shuffle_metadata.tsv"
RUN_META_TSV="${OUT_DIR}/run_metadata.tsv"
RUN_LOG="${OUT_DIR}/shuffle_control_run.log"

###############################################################################
# Create shuffled FASTA
###############################################################################

log_info "Creating shuffled pathogen FASTA"

python3 - "${PATHOGEN_FASTA}" "${SHUFFLED_FASTA}" "${SHUFFLE_META_TSV}" "${SHUFFLE_SEED}" "${SHUFFLE_MODE}" <<'PY'
"""Create a mononucleotide-shuffled FASTA with short, safe headers."""

import gzip
import random
import sys
from collections import Counter
from pathlib import Path


def open_maybe_gzip(path_str, mode):
    """Open plain-text or gzipped file depending on suffix."""
    path = Path(path_str)
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", errors="replace")
    return open(path, mode, encoding="utf-8", errors="replace")


def read_fasta(path_str):
    """Yield FASTA records as (header, sequence) tuples."""
    header = None
    seq_parts = []

    with open_maybe_gzip(path_str, "rt") as handle:
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


def wrap_sequence(seq, width=80):
    """Yield fixed-width chunks from a sequence string."""
    for start in range(0, len(seq), width):
        yield seq[start:start + width]


def shuffle_mono(seq, rng):
    """Shuffle A/C/G/T positions while leaving other characters in place."""
    seq_list = list(seq)
    acgt_positions = [
        idx for idx, base in enumerate(seq_list)
        if base in {"A", "C", "G", "T"}
    ]
    acgt_bases = [seq_list[idx] for idx in acgt_positions]
    rng.shuffle(acgt_bases)

    for idx, base in zip(acgt_positions, acgt_bases):
        seq_list[idx] = base

    return "".join(seq_list)


input_fasta = sys.argv[1]
output_fasta = sys.argv[2]
metadata_tsv = sys.argv[3]
seed = int(sys.argv[4])
shuffle_mode = sys.argv[5]

if shuffle_mode != "mononucleotide":
    raise SystemExit(f"Unsupported SHUFFLE_MODE: {shuffle_mode}")

rng = random.Random(seed)
records = list(read_fasta(input_fasta))

if not records:
    raise SystemExit("No FASTA records found in PATHOGEN_FASTA.")

with open(output_fasta, "wt", encoding="utf-8") as out_handle, \
        open(metadata_tsv, "wt", encoding="utf-8") as meta_handle:
    meta_handle.write(
        "record_index\toriginal_header\tshuffled_header\tlength\t"
        "A\tC\tG\tT\tN_or_other\n"
    )

    for idx, (header, seq) in enumerate(records, start=1):
        shuffled_seq = shuffle_mono(seq=seq, rng=rng)
        shuffled_header = f"shuffle_contig_{idx}"

        out_handle.write(f">{shuffled_header}\n")
        for chunk in wrap_sequence(seq=shuffled_seq, width=80):
            out_handle.write(chunk + "\n")

        counts = Counter(seq)
        n_other = sum(
            value for key, value in counts.items()
            if key not in {"A", "C", "G", "T"}
        )

        meta_handle.write(
            f"{idx}\t{header}\t{shuffled_header}\t{len(seq)}\t"
            f"{counts.get('A', 0)}\t{counts.get('C', 0)}\t"
            f"{counts.get('G', 0)}\t{counts.get('T', 0)}\t{n_other}\n"
        )
PY

require_file "${SHUFFLED_FASTA}"
require_file "${SHUFFLE_META_TSV}"

###############################################################################
# Write provenance
###############################################################################

log_info "Writing run metadata"

{
    printf 'parameter\tvalue\n'
    printf 'real_fastq\t%s\n' "${REAL_FASTQ}"
    printf 'original_pathogen_fasta\t%s\n' "${PATHOGEN_FASTA}"
    printf 'shuffled_fasta\t%s\n' "${SHUFFLED_FASTA}"
    printf 'target_label\t%s\n' "${TARGET_LABEL}"
    printf 'monkey_small_gz\t%s\n' "${MONKEY_SMALL_GZ}"
    printf 'monkey_small_fasta\t%s\n' "${MONKEY_SMALL_FASTA}"
    printf 'monkey_depletion_fasta\t%s\n' "${MONKEY_DEPLETION_FASTA}"
    printf 'plasmo_masked_gz\t%s\n' "${PLASMO_MASKED_GZ}"
    printf 'depletion_ref_fasta\t%s\n' "${DEPLETION_REF_FASTA}"
    printf 'depleted_fastq\t%s\n' "${DEPLETED_FASTQ}"
    printf 'minimap_db_gz\t%s\n' "${MINIMAP_DB_GZ}"
    printf 'kraken_db_dir\t%s\n' "${KRAKEN_DB_DIR}"
    printf 'threads\t%s\n' "${THREADS}"
    printf 'train_reads_n\t%s\n' "${TRAIN_READS_N}"
    printf 'sim_pool_n\t%s\n' "${SIM_POOL_N}"
    printf 'spike_levels\t%s\n' "${SPIKE_LEVELS}"
    printf 'replicates\t%s\n' "${REPLICATES}"
    printf 'do_host_depletion\t%s\n' "${DO_HOST_DEPLETION}"
    printf 'shuffle_seed\t%s\n' "${SHUFFLE_SEED}"
    printf 'shuffle_mode\t%s\n' "${SHUFFLE_MODE}"
    printf 'main_script\t%s\n' "${MAIN_SCRIPT}"
} > "${RUN_META_TSV}"

###############################################################################
# Run the main one-sample pipeline with shuffled FASTA
###############################################################################

log_info "Launching main one-sample pipeline using shuffled pathogen FASTA"
log_info "Main script: ${MAIN_SCRIPT}"
log_info "Main output directory: ${OUT_DIR}/main_pipeline_run"
log_info "Run log: ${RUN_LOG}"

REAL_FASTQ="${REAL_FASTQ}" \
MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ}" \
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA}" \
MONKEY_DEPLETION_FASTA="${MONKEY_DEPLETION_FASTA}" \
PLASMO_MASKED_GZ="${PLASMO_MASKED_GZ}" \
DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA}" \
DEPLETED_FASTQ="${DEPLETED_FASTQ}" \
PATHOGEN_FASTA="${SHUFFLED_FASTA}" \
TARGET_LABEL="${TARGET_LABEL}" \
MINIMAP_DB_GZ="${MINIMAP_DB_GZ}" \
KRAKEN_DB_DIR="${KRAKEN_DB_DIR}" \
SCRIPT_DIR="${SCRIPT_DIR}" \
OUT_DIR="${OUT_DIR}/main_pipeline_run" \
THREADS="${THREADS}" \
TRAIN_READS_N="${TRAIN_READS_N}" \
SIM_POOL_N="${SIM_POOL_N}" \
SPIKE_LEVELS="${SPIKE_LEVELS}" \
REPLICATES="${REPLICATES}" \
DO_HOST_DEPLETION="${DO_HOST_DEPLETION}" \
bash "${MAIN_SCRIPT}" 2>&1 | tee "${RUN_LOG}"

log_info "Shuffled control completed successfully"