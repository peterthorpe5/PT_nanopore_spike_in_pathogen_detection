#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N spike_single_mm_full

set -euo pipefail

log_error() {
    printf '[ERROR] %s\n' "$*" >&2
}

log_info() {
    printf '[INFO] %s\n' "$*"
}

require_file() {
    [[ -f "$1" ]] || {
        log_error "Required file not found: $1"
        exit 1
    }
}

if [[ -n "${REPO_DIR:-}" ]]; then
    REPO_DIR="${REPO_DIR}"
elif [[ -n "${SGE_O_WORKDIR:-}" ]]; then
    REPO_DIR="${SGE_O_WORKDIR}"
else
    REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection"
fi

SHELLS_DIR="${SHELLS_DIR:-${REPO_DIR}/shells}"
CONFIG_DIR="${CONFIG_DIR:-${REPO_DIR}/configs}"
source "${CONFIG_DIR}/pipeline_paths.sh"

FULL_MINIMAP_DB_FASTA="${FULL_MINIMAP_DB_FASTA:-${FULL_MINIMAP_DB_FASTA_DEFAULT:-${MINIMAP_DB_FASTA_DEFAULT:-}}}"
[[ -n "${FULL_MINIMAP_DB_FASTA}" ]] || {
    log_error "FULL_MINIMAP_DB_FASTA is empty and no full-reference default was set in pipeline_paths.sh"
    exit 1
}
require_file "${FULL_MINIMAP_DB_FASTA}"

export MINIMAP_DB_FASTA="${FULL_MINIMAP_DB_FASTA}"
export OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_single_minimap_fullref_$(date +%Y%m%d_%H%M%S)}"

log_info "Single minimap full-reference run using: ${MINIMAP_DB_FASTA}"

bash "${SHELLS_DIR}/run_spikein_single_minimap_only.sh"
