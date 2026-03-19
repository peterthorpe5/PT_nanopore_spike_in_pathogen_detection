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

log_info() { printf '[INFO] %s\n' "$*"; }
log_error() { printf '[ERROR] %s\n' "$*" >&2; }
require_file() { [[ -f "$1" ]] || { log_error "Required file not found: $1"; exit 1; }; }
require_dir() { [[ -d "$1" ]] || { log_error "Required directory not found: $1"; exit 1; }; }
require_exe() { command -v "$1" >/dev/null 2>&1 || { log_error "Required executable not found on PATH: $1"; exit 1; }; }

if [[ -n "${REPO_DIR:-}" ]]; then
    REPO_DIR="${REPO_DIR}"
elif [[ -n "${SGE_O_WORKDIR:-}" ]]; then
    REPO_DIR="${SGE_O_WORKDIR}"
else
    REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection"
fi

PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR:-${REPO_DIR}/scripts}"
PATHOGEN_CONFIG_TSV="${PATHOGEN_CONFIG_TSV:-${REPO_DIR}/configs/pathogen_panel_3_with_taxid.tsv}"
NCBI_TAXONOMY_DIR="${NCBI_TAXONOMY_DIR:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal/taxonomy}"
ANNOTATE_FASTA_PY="${ANNOTATE_FASTA_PY:-${PY_SCRIPTS_DIR}/annotate_fasta_with_taxid.py}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/metamaps_db/custom_$(date +%Y%m%d_%H%M%S)}"
DB_NAME="${DB_NAME:-custom_metamaps_db}"

require_exe python3
require_exe perl
require_exe buildDB.pl
require_file "${PATHOGEN_CONFIG_TSV}"
require_dir "${NCBI_TAXONOMY_DIR}"
require_file "${ANNOTATE_FASTA_PY}"

mkdir -p "${OUT_DIR}/annotated_fastas"
COMBINED_FASTA="${OUT_DIR}/combined_input.fa"
MANIFEST_TSV="${OUT_DIR}/database_manifest.tsv"
DB_DIR="${OUT_DIR}/${DB_NAME}"

: > "${COMBINED_FASTA}"
printf 'label\ttaxid\tinput_fasta\tannotated_fasta\n' > "${MANIFEST_TSV}"

while IFS=$'\t' read -r fasta label taxid extra; do
    [[ -n "${fasta}" ]] || continue
    [[ "${fasta}" == "fasta" ]] && continue
    require_file "${fasta}"
    [[ -n "${label}" ]] || { log_error "Missing label for ${fasta}"; exit 1; }
    [[ -n "${taxid}" ]] || { log_error "Missing taxid for ${fasta}"; exit 1; }

    safe_label="$(printf '%s' "${label}" | tr ' ' '_' | tr '/' '_')"
    annotated_fasta="${OUT_DIR}/annotated_fastas/${safe_label}.taxid_${taxid}.fa"
    python3 "${ANNOTATE_FASTA_PY}" \
        --input_fasta "${fasta}" \
        --taxid "${taxid}" \
        --out_fasta "${annotated_fasta}"

    cat "${annotated_fasta}" >> "${COMBINED_FASTA}"
    printf '%s\t%s\t%s\t%s\n' "${label}" "${taxid}" "${fasta}" "${annotated_fasta}" >> "${MANIFEST_TSV}"
done < "${PATHOGEN_CONFIG_TSV}"

[[ -s "${COMBINED_FASTA}" ]] || { log_error "Combined FASTA is empty: ${COMBINED_FASTA}"; exit 1; }

perl "$(command -v buildDB.pl)" \
    --DB "${DB_DIR}" \
    --FASTAs "${COMBINED_FASTA}" \
    --taxonomy "${NCBI_TAXONOMY_DIR}"

log_info "Done: ${DB_DIR}"
log_info "Database FASTA: ${DB_DIR}/DB.fa"
log_info "Manifest: ${MANIFEST_TSV}"
