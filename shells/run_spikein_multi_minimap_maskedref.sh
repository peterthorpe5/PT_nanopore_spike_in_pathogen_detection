#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N spike_multi_mm_masked

set -euo pipefail

log_error() {
    printf '[ERROR] %s\n' "$*" >&2
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

MASKED_MINIMAP_DB_FASTA="${MASKED_MINIMAP_DB_FASTA:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta}"
require_file "${MASKED_MINIMAP_DB_FASTA}"

export MINIMAP_DB_FASTA="${MASKED_MINIMAP_DB_FASTA}"
export OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_minimap_maskedref_$(date +%Y%m%d_%H%M%S)}"

bash "${SHELLS_DIR}/run_spikein_multi_minimap_only.sh"
