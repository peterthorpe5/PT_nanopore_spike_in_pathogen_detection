#!/usr/bin/env bash
set -euo pipefail

log_info() { printf '[INFO] %s\n' "$*"; }

PROJECT_DIR="${PROJECT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity}"
REPO_DIR="${REPO_DIR:-${PROJECT_DIR}/PT_nanopore_spike_in_pathogen_detection}"
SHELLS_DIR="${SHELLS_DIR:-${REPO_DIR}/shells}"
CONFIG_DIR="${CONFIG_DIR:-${REPO_DIR}/configs}"

# Kraken2 does not have one universal recommended confidence threshold.
# This run uses a common real-world screening threshold of 0.10 to reduce
# database-related false positives while retaining sensitivity.
export KRAKEN_CONFIDENCE="${KRAKEN_CONFIDENCE:-0.10}"

# Metabuli's Nanopore long-read recommended min-score is 0.008.
export MIN_SCORE="${MIN_SCORE:-0.008}"
export METABULI_DB_DIR="${METABULI_DB_DIR:-/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/databases/metabuli/custom_metabuli_db}"

export SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 3000 4000 5000 7000 10000}"
export REPLICATES="${REPLICATES:-12}"

log_info "Submitting recommended-threshold benchmark jobs"
log_info "Kraken2 confidence: ${KRAKEN_CONFIDENCE}"
log_info "Metabuli min-score: ${MIN_SCORE}"
log_info "Spike levels: ${SPIKE_LEVELS}"

# Kraken2 confidence-filtered read-level runs.
unset OUT_DIR
export OUT_DIR="${PROJECT_DIR}/runs/spikein_single_read_kraken_confidence_${KRAKEN_CONFIDENCE}_$(date +%Y%m%d_%H%M%S)"
qsub -V "${SHELLS_DIR}/run_spikein_single_readlevel_kraken_confidence.sh"

unset OUT_DIR
export PATHOGEN_CONFIG_TSV="${CONFIG_DIR}/pathogen_panel_2.tsv"
export OUT_DIR="${PROJECT_DIR}/runs/spikein_multi_read_panel2_kraken_confidence_${KRAKEN_CONFIDENCE}_$(date +%Y%m%d_%H%M%S)"
qsub -V "${SHELLS_DIR}/run_spikein_multi_readlevel_kraken_confidence.sh"

unset OUT_DIR
export PATHOGEN_CONFIG_TSV="${CONFIG_DIR}/pathogen_panel_3.tsv"
export OUT_DIR="${PROJECT_DIR}/runs/spikein_multi_read_panel3_kraken_confidence_${KRAKEN_CONFIDENCE}_$(date +%Y%m%d_%H%M%S)"
qsub -V "${SHELLS_DIR}/run_spikein_multi_readlevel_kraken_confidence.sh"

# Metabuli recommended Nanopore min-score runs.
unset OUT_DIR
export PATHOGEN_CONFIG_TSV="${CONFIG_DIR}/pathogen_panel_1.tsv"
export OUT_DIR="${PROJECT_DIR}/runs/spikein_single_metabuli_panel1_min_score_${MIN_SCORE}_$(date +%Y%m%d_%H%M%S)"
qsub -V "${SHELLS_DIR}/run_spikein_readlevel_metabuli.sh"

unset OUT_DIR
export PATHOGEN_CONFIG_TSV="${CONFIG_DIR}/pathogen_panel_2.tsv"
export OUT_DIR="${PROJECT_DIR}/runs/spikein_multi_metabuli_panel2_min_score_${MIN_SCORE}_$(date +%Y%m%d_%H%M%S)"
qsub -V "${SHELLS_DIR}/run_spikein_readlevel_metabuli.sh"

unset OUT_DIR
export PATHOGEN_CONFIG_TSV="${CONFIG_DIR}/pathogen_panel_3.tsv"
export OUT_DIR="${PROJECT_DIR}/runs/spikein_multi_metabuli_panel3_min_score_${MIN_SCORE}_$(date +%Y%m%d_%H%M%S)"
qsub -V "${SHELLS_DIR}/run_spikein_readlevel_metabuli.sh"

log_info "Submitted recommended-threshold jobs"
