#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 64G
#$ -adds l_hard h_vmem 64G
#$ -N melon_mixed

set -euo pipefail

log_info() {
  printf '[INFO] %s\n' "$*"
}

log_error() {
  printf '[ERROR] %s\n' "$*" >&2
}

require_file() {
  [[ -f "$1" ]] || {
    log_error "Required file not found: $1"
    exit 1
  }
}

require_dir() {
  [[ -d "$1" ]] || {
    log_error "Required directory not found: $1"
    exit 1
  }
}

require_exe() {
  command -v "$1" >/dev/null 2>&1 || {
    log_error "Required executable not found on PATH: $1"
    exit 1
  }
}

THREADS="${THREADS:-12}"
SAMPLE_NAME="${SAMPLE_NAME:-mixed_genome_test}"

REPO_DIR="${REPO_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection}"
REAL_FASTQ="${REAL_FASTQ:-/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC1023_AmM008WB.fastq.gz}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/MRC1023_AmM008WB.depleted.fastq.gz}"

READS="${READS:-${DEPLETED_FASTQ}}"
MELON_DB="${MELON_DB:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/databases/melon/database}"
KRAKEN_DB="${KRAKEN_DB:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"

OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/melon_mixed_$(date +%Y%m%d_%H%M%S)}"

USE_PREFILTER="${USE_PREFILTER:-false}"
MIN_READ_LENGTH="${MIN_READ_LENGTH:-1000}"
MIN_READ_QUALITY="${MIN_READ_QUALITY:-10}"

for x in melon nanoq gzip find; do
  require_exe "${x}"
done

require_file "${READS}"
require_dir "${KRAKEN_DB}"

if [[ ! -d "${MELON_DB}" ]]; then
  log_error "Melon database directory not found: ${MELON_DB}"
  log_error "Please set MELON_DB to the Melon database folder."
  exit 1
fi

mkdir -p "${OUT_DIR}"
mkdir -p "${OUT_DIR}/qc"
mkdir -p "${OUT_DIR}/logs"

log_info "Sample name: ${SAMPLE_NAME}"
log_info "Input reads: ${READS}"
log_info "Melon database: ${MELON_DB}"
log_info "Output directory: ${OUT_DIR}"
log_info "Threads: ${THREADS}"

QC_READS="${OUT_DIR}/qc/${SAMPLE_NAME}.q${MIN_READ_QUALITY}.min${MIN_READ_LENGTH}.fastq.gz"

log_info "Running NanoQ filtering"
gzip -dc "${READS}" | \
  nanoq \
    --quality "${MIN_READ_QUALITY}" \
    --length "${MIN_READ_LENGTH}" | \
  gzip > "${QC_READS}"

require_file "${QC_READS}"

log_info "Running Melon"
if [[ "${USE_PREFILTER}" == "true" ]]; then
  melon \
    "${QC_READS}" \
    -d "${MELON_DB}" \
    -o "${OUT_DIR}" \
    -k "${KRAKEN_DB}" \
    > "${OUT_DIR}/logs/melon.stdout.log" \
    2> "${OUT_DIR}/logs/melon.stderr.log"
else
  melon \
    "${QC_READS}" \
    -d "${MELON_DB}" \
    -o "${OUT_DIR}" \
    > "${OUT_DIR}/logs/melon.stdout.log" \
    2> "${OUT_DIR}/logs/melon.stderr.log"
fi

log_info "Melon run complete"
log_info "Listing output files"
find "${OUT_DIR}" -maxdepth 2 -type f | sort

log_info "Done: ${OUT_DIR}"
