#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N spike_shuffle_read

set -euo pipefail
log_info(){ printf '[INFO] %s\n' "$*"; }
log_error(){ printf '[ERROR] %s\n' "$*" >&2; }
require_file(){ [[ -f "$1" ]] || { log_error "Required file not found: $1"; exit 1; }; }
if [[ -n "${REPO_DIR:-}" ]]; then REPO_DIR="${REPO_DIR}"; elif [[ -n "${SGE_O_WORKDIR:-}" ]]; then REPO_DIR="${SGE_O_WORKDIR}"; else REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection"; fi
PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR:-${REPO_DIR}/scripts}"
SHUFFLE_FASTA_PY="${PY_SCRIPTS_DIR}/shuffle_fasta.py"
PATHOGEN_FASTA="${PATHOGEN_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/GCA_014843685.1_ASM1484368v1_genomic.fna}"
SHUFFLE_SEED="${SHUFFLE_SEED:-123}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_shuffle_read_$(date +%Y%m%d_%H%M%S)}"
require_file "${SHUFFLE_FASTA_PY}"; require_file "${PATHOGEN_FASTA}"; mkdir -p "${OUT_DIR}"
SHUFFLED_FASTA="${OUT_DIR}/pathogen.shuffled.fasta"; SHUFFLE_META_TSV="${OUT_DIR}/shuffle_metadata.tsv"
python3 "${SHUFFLE_FASTA_PY}" --input_fasta "${PATHOGEN_FASTA}" --output_fasta "${SHUFFLED_FASTA}" --metadata_tsv "${SHUFFLE_META_TSV}" --seed "${SHUFFLE_SEED}"
REPO_DIR="${REPO_DIR}" PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR}" PATHOGEN_FASTA="${SHUFFLED_FASTA}" OUT_DIR="${OUT_DIR}/main_run" bash "${REPO_DIR}/run_spikein_single_readlevel.sh"
