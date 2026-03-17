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

if [[ -n "${REPO_DIR:-}" ]]; then
    REPO_DIR="${REPO_DIR}"
elif [[ -n "${SGE_O_WORKDIR:-}" ]]; then
    REPO_DIR="${SGE_O_WORKDIR}"
else
    REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection"
fi

PY_DIR="${REPO_DIR}/scripts"
CONFIG_DIR="${REPO_DIR}/configs"

PATHOGEN_CONFIG_TSV="${CONFIG_DIR}/pathogen_panel_3.tsv"

SHUFFLE_FASTA_PY="${PY_DIR}/shuffle_fasta.py"
SAMPLE_FASTQ_PY="${PY_DIR}/sample_fastq.py"
BUILD_MIXED_FASTQ_PY="${PY_DIR}/build_mixed_fastq.py"
COMBINE_NANOSIM_FASTQ_PY="${PY_DIR}/combine_nanosim_fastq.py"
SUMMARISE_KRAKEN_PY="${PY_DIR}/summarise_kraken_report.py"
COUNT_BED_HITS_PY="${PY_DIR}/count_bed_hits.py"
DEDUP_FASTQ_NAMES_PY="${PY_DIR}/dedup_fastq_names.py"
ASSEMBLY_STATS_PY="${PY_DIR}/assembly_stats.py"

PATHOGEN_FASTA="${PATHOGEN_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/GCA_014843685.1_ASM1484368v1_genomic.fna}"
SHUFFLE_SEED="${SHUFFLE_SEED:-123}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_shuffle_flye_$(date +%Y%m%d_%H%M%S)}"

mkdir -p "${OUT_DIR}"
python3 "${PY_SCRIPTS_DIR}/shuffle_fasta.py" --input_fasta "${PATHOGEN_FASTA}" --output_fasta "${OUT_DIR}/pathogen.shuffled.fasta" --metadata_tsv "${OUT_DIR}/shuffle_metadata.tsv" --seed "${SHUFFLE_SEED}"
REPO_DIR="${REPO_DIR}" PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR}" PATHOGEN_FASTA="${OUT_DIR}/pathogen.shuffled.fasta" OUT_DIR="${OUT_DIR}/main_run" bash "${REPO_DIR}/run_spikein_single_flye.sh"
