#!/usr/bin/env bash

set -euo pipefail

cd /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity

conda activate kraken16s

export REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection"
export SHELLS_DIR="${REPO_DIR}/shells"
export CONFIG_DIR="${REPO_DIR}/configs"

log_info() {
    printf '[INFO] %s\n' "$*"
}

submit_job() {
    local script_path="$1"
    if [[ ! -f "${script_path}" ]]; then
        printf '[ERROR] Missing script: %s\n' "${script_path}" >&2
        exit 1
    fi
    log_info "Submitting ${script_path}"
    qsub -V "${script_path}"
}

log_info "Repository: ${REPO_DIR}"
log_info "Shells dir: ${SHELLS_DIR}"
log_info "Configs dir: ${CONFIG_DIR}"


# Single-genome workflows
submit_job "${SHELLS_DIR}/run_spikein_single_readlevel.sh"
submit_job "${SHELLS_DIR}/run_spikein_single_flye.sh"
submit_job "${SHELLS_DIR}/run_spikein_single_flye_medaka.sh"


# Shuffled control workflows
unset PATHOGEN_CONFIG_TSV

export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_shuffle_read_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_shuffled_readlevel.sh"

export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_shuffle_flye_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_shuffled_flye.sh"

export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_shuffle_flye_medaka_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_shuffled_flye_medaka_clean.sh"


# Multi-genome workflows: 2-genome panel
export PATHOGEN_CONFIG_TSV="${CONFIG_DIR}/pathogen_panel_2.tsv"
log_info "Using config: ${PATHOGEN_CONFIG_TSV}"

export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_read_panel2_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_multi_readlevel.sh"

export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_flye_panel2_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_multi_flye_clean.sh"

export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_flye_medaka_panel2_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_multi_flye_medaka_clean.sh"


# Multi-genome workflows: 3-genome panel
export PATHOGEN_CONFIG_TSV="${CONFIG_DIR}/pathogen_panel_3.tsv"
log_info "Using config: ${PATHOGEN_CONFIG_TSV}"

export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_read_panel3_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_multi_readlevel.sh"

export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_flye_panel3_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_multi_flye_clean.sh"

export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_flye_medaka_panel3_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_multi_flye_medaka_clean.sh"



####################
# metabuli workflows

conda activate kraken16s


export METABULI_DB_DIR="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/databases/metabuli/custom_metabuli_db"

export PATHOGEN_CONFIG_TSV="${CONFIG_DIR}/pathogen_panel_1.tsv"
export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_single_metabuli_panel1_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_readlevel_metabuli.sh"

export PATHOGEN_CONFIG_TSV="${CONFIG_DIR}/pathogen_panel_2.tsv"
export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_metabuli_panel2_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_readlevel_metabuli.sh"

export PATHOGEN_CONFIG_TSV="${CONFIG_DIR}/pathogen_panel_3.tsv"
export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_metabuli_panel3_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_readlevel_metabuli.sh"

##### REPEAT THOSE THAT FAILED


qsub "${SHELLS_DIR}/rerun_minimap_shared_db_q15_len500.sh"
qsub "${SHELLS_DIR}/rerun_minimap_focused_plasdb_q15_len500.sh"

qsub "${SHELLS_DIR}/rerun_single_genome_flye_assembly.sh"
qsub "${SHELLS_DIR}/rerun_single_genome_flye_medaka_assembly.sh"
qsub "${SHELLS_DIR}/rerun_shuffled_single_genome_flye_assembly.sh"
qsub "${SHELLS_DIR}/rerun_shuffled_single_genome_flye_medaka_assembly.sh"



################################
#   THESE DONT WROK, DONT RUN

####################
# Metamaps workflows

conda activate metamaps

export PATHOGEN_CONFIG_TSV="${CONFIG_DIR}/pathogen_panel_1.tsv"
export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_single_metamaps_panel1_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_single_readlevel_metamaps.sh"


export PATHOGEN_CONFIG_TSV="${CONFIG_DIR}/pathogen_panel_2.tsv"
export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_metamaps_panel2_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_multi_readlevel_metamaps.sh"


export PATHOGEN_CONFIG_TSV="${CONFIG_DIR}/pathogen_panel_3.tsv"
export OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_metamaps_panel3_$(date +%Y%m%d_%H%M%S)"
submit_job "${SHELLS_DIR}/run_spikein_multi_readlevel_metamaps.sh"


log_info "All qsub commands submitted"