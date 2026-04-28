#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N rerun_single_flye_medaka

set -euo pipefail

log_info() { printf '[INFO] %s\n' "$*"; }
log_error() { printf '[ERROR] %s\n' "$*" >&2; }
require_file() { [[ -f "$1" ]] || { log_error "Required file not found: $1"; exit 1; }; }
require_dir() { [[ -d "$1" ]] || { log_error "Required directory not found: $1"; exit 1; }; }

PROJECT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity"
REPO_DIR_DEFAULT="${PROJECT_DIR}/PT_nanopore_spike_in_pathogen_detection"
REPO_DIR="${REPO_DIR:-${REPO_DIR_DEFAULT}}"
SHELLS_DIR="${SHELLS_DIR:-${REPO_DIR}/shells}"
RUN_SCRIPT="${RUN_SCRIPT:-${SHELLS_DIR}/run_spikein_single_flye_medaka.sh}"

require_dir "${PROJECT_DIR}"
require_dir "${REPO_DIR}"
require_file "${RUN_SCRIPT}"

cd "${PROJECT_DIR}"

if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate kraken16s || log_info "conda activate kraken16s failed; continuing with current environment"
fi

export REPO_DIR
export SHELLS_DIR
export CONFIG_DIR="${CONFIG_DIR:-${REPO_DIR}/configs}"
export PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR:-${REPO_DIR}/scripts}"
export MEDAKA_MODEL="${MEDAKA_MODEL:-auto}"
export OUT_DIR="${OUT_DIR:-${PROJECT_DIR}/runs/spikein_single_flye_medaka_rerun_$(date +%Y%m%d_%H%M%S)}"

log_info "Project dir: ${PROJECT_DIR}"
log_info "Repo dir: ${REPO_DIR}"
log_info "Using run script: ${RUN_SCRIPT}"
log_info "Medaka model: ${MEDAKA_MODEL}"
log_info "Output dir: ${OUT_DIR}"
log_info "Submitting single-genome Flye plus Medaka assembly rerun"

bash "${RUN_SCRIPT}"

log_info "Single-genome Flye plus Medaka rerun finished"
log_info "Output dir: ${OUT_DIR}"
