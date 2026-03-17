#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 24
#$ -jc long
#$ -mods l_hard mfree 450G
#$ -adds l_hard h_vmem 450G
#$ -N spike_shuffle_flye

set -euo pipefail

default_repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-${default_repo_dir}}"
PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR:-${REPO_DIR}/scripts}"
PATHOGEN_FASTA="${PATHOGEN_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/GCA_014843685.1_ASM1484368v1_genomic.fna}"
SHUFFLE_SEED="${SHUFFLE_SEED:-123}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_shuffle_flye_$(date +%Y%m%d_%H%M%S)}"

mkdir -p "${OUT_DIR}"
python3 "${PY_SCRIPTS_DIR}/shuffle_fasta.py" --input_fasta "${PATHOGEN_FASTA}" --output_fasta "${OUT_DIR}/pathogen.shuffled.fasta" --metadata_tsv "${OUT_DIR}/shuffle_metadata.tsv" --seed "${SHUFFLE_SEED}"
REPO_DIR="${REPO_DIR}" PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR}" PATHOGEN_FASTA="${OUT_DIR}/pathogen.shuffled.fasta" OUT_DIR="${OUT_DIR}/main_run" bash "${REPO_DIR}/run_spikein_single_flye.sh"
