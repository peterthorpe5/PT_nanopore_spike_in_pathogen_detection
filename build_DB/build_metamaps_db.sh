#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 4
#$ -jc long
#$ -mods l_hard mfree 64G
#$ -adds l_hard h_vmem 64G
#$ -N build_metamaps_db

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

THREADS="${THREADS:-4}"
KRAKEN_DB="${KRAKEN_DB:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"
OUT_BASE="${OUT_BASE:-/home/pthorpe001/data/databases}"
OUT_DIR="${OUT_DIR:-${OUT_BASE}/metamaps/custom_$(date +%Y%m%d_%H%M%S)}"
DB_NAME="${DB_NAME:-custom_metamaps_db}"
NCBI_TAXONOMY_DIR="${NCBI_TAXONOMY_DIR:-${KRAKEN_DB}/taxonomy}"

require_exe perl
require_exe buildDB.pl
require_exe find
require_exe sort
require_dir "${KRAKEN_DB}"
require_dir "${KRAKEN_DB}/library"
require_dir "${NCBI_TAXONOMY_DIR}"

mkdir -p "${OUT_DIR}"

COMBINED_FASTA="${OUT_DIR}/combined_input.fa"
MANIFEST_TSV="${OUT_DIR}/database_manifest.tsv"
DB_DIR="${OUT_DIR}/${DB_NAME}"

log_info "Kraken DB: ${KRAKEN_DB}"
log_info "MetaMaps output dir: ${OUT_DIR}"
log_info "Taxonomy dir: ${NCBI_TAXONOMY_DIR}"

mapfile -t LIB_FASTAS < <(
    find "${KRAKEN_DB}/library" -type f -name 'library.fna' | sort
)

if [[ "${#LIB_FASTAS[@]}" -eq 0 ]]; then
    log_error "No library.fna files found under ${KRAKEN_DB}/library"
    exit 1
fi

: > "${COMBINED_FASTA}"
printf 'source_fasta\n' > "${MANIFEST_TSV}"

for fasta in "${LIB_FASTAS[@]}"; do
    require_file "${fasta}"
    log_info "Adding ${fasta}"
    cat "${fasta}" >> "${COMBINED_FASTA}"
    printf '%s\n' "${fasta}" >> "${MANIFEST_TSV}"
done

[[ -s "${COMBINED_FASTA}" ]] || {
    log_error "Combined FASTA is empty: ${COMBINED_FASTA}"
    exit 1
}

perl "$(command -v buildDB.pl)" \
    --DB "${DB_DIR}" \
    --FASTAs "${COMBINED_FASTA}" \
    --taxonomy "${NCBI_TAXONOMY_DIR}"

log_info "Done"
log_info "Database dir: ${DB_DIR}"
log_info "Database FASTA: ${DB_DIR}/DB.fa"
log_info "Manifest: ${MANIFEST_TSV}"
log_info "Combined input FASTA: ${COMBINED_FASTA}"
