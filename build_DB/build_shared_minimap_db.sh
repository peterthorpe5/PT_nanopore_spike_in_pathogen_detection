#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 8
#$ -jc long
#$ -mods l_hard mfree 64G
#$ -adds l_hard h_vmem 64G
#$ -N build_shared_minimap_db

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

THREADS="${THREADS:-8}"
KRAKEN_DB="${KRAKEN_DB:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"
OUT_BASE="${OUT_BASE:-/home/pthorpe001/data/databases}"
SHARED_DIR="${SHARED_DIR:-${OUT_BASE}/shared_refs}"
MINIMAP_DIR="${MINIMAP_DIR:-${OUT_BASE}/minimap2}"

require_exe minimap2
require_exe find
require_exe sort
require_exe cp
require_dir "${KRAKEN_DB}"
require_dir "${KRAKEN_DB}/library"

mkdir -p "${SHARED_DIR}"
mkdir -p "${MINIMAP_DIR}"

COMBINED_FASTA="${SHARED_DIR}/shared_bact_viral_plasmo_refs.fa"
MANIFEST_TSV="${SHARED_DIR}/shared_bact_viral_plasmo_refs.manifest.tsv"
MINIMAP_FASTA="${MINIMAP_DIR}/shared_bact_viral_plasmo_refs.fa"
MINIMAP_INDEX="${MINIMAP_DIR}/shared_bact_viral_plasmo_refs.mmi"

log_info "Kraken DB: ${KRAKEN_DB}"
log_info "Shared FASTA: ${COMBINED_FASTA}"
log_info "Minimap FASTA: ${MINIMAP_FASTA}"
log_info "Minimap index: ${MINIMAP_INDEX}"

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

cp "${COMBINED_FASTA}" "${MINIMAP_FASTA}"

log_info "Building minimap2 index"
minimap2 \
    -x map-ont \
    -d "${MINIMAP_INDEX}" \
    "${MINIMAP_FASTA}"

log_info "Done"
log_info "Combined FASTA: ${COMBINED_FASTA}"
log_info "Manifest: ${MANIFEST_TSV}"
log_info "Minimap FASTA: ${MINIMAP_FASTA}"
log_info "Minimap index: ${MINIMAP_INDEX}"
