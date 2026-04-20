#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 16
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N minimap_rerun_focus

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

source "${REPO_DIR}/configs/pipeline_paths.sh"

FOCUSED_DB_FASTA="${FOCUSED_DB_FASTA:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta}"
require_file "${FOCUSED_DB_FASTA}"
require_dir "${REPO_DIR}"
require_dir "${REPO_DIR}/shells"
require_dir "${REPO_DIR}/summary"
require_exe python3
require_exe minimap2
require_exe samtools
require_exe kraken2
require_exe bamToBed
require_exe read_analysis.py
require_exe simulator.py

STAMP="$(date +%Y%m%d_%H%M%S)"
THREADS="${THREADS:-16}"
MIN_MAPQ="${MIN_MAPQ:-15}"
MIN_ALIGN_LEN="${MIN_ALIGN_LEN:-500}"

PROJECT_BASE="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity"
NEW_RUNS_DIR="${PROJECT_BASE}/runs_minimapfix_focused_q15_len500_${STAMP}"
PATCHED_SHELLS_DIR="${NEW_RUNS_DIR}/patched_shells"
SUMMARY_STAGE_DIR="${PROJECT_BASE}/summary_input_minimapfix_focused_q15_len500_${STAMP}"
SUMMARY_OUT_DIR="${PROJECT_BASE}/spikein_summary_minimapfix_focused_q15_len500_${STAMP}"

mkdir -p "${NEW_RUNS_DIR}" "${PATCHED_SHELLS_DIR}" "${SUMMARY_STAGE_DIR}" "${SUMMARY_OUT_DIR}"

log_info "Using focused minimap database"
log_info "MINIMAP_DB_FASTA=${FOCUSED_DB_FASTA}"
log_info "MIN_MAPQ=${MIN_MAPQ}"
log_info "MIN_ALIGN_LEN=${MIN_ALIGN_LEN}"
log_info "New runs will be written under ${NEW_RUNS_DIR}"

for shell_name in \
    run_spikein_single_readlevel.sh \
    run_spikein_multi_readlevel.sh \
    run_spikein_shuffled_readlevel.sh; do
    cp "${REPO_DIR}/shells/${shell_name}" "${PATCHED_SHELLS_DIR}/${shell_name}"
    require_file "${PATCHED_SHELLS_DIR}/${shell_name}"
done

python3 - <<'PY'
from pathlib import Path
import os

patched_dir = Path(os.environ['PATCHED_SHELLS_DIR'])
replacements = {
    'samtools view -Sh -q "${MIN_MAPQ}" -e "rlen >= ${MIN_ALIGN_LEN}" -': 'samtools view -@ "${THREADS}" -Sh -q "${MIN_MAPQ}" -e "length(seq) >= ${MIN_ALIGN_LEN}" -',
    'samtools sort -O bam -@ "${THREADS}" -': 'samtools sort -@ "${THREADS}" -O bam -',
}
for script_name in ["run_spikein_single_readlevel.sh", "run_spikein_multi_readlevel.sh"]:
    path = patched_dir / script_name
    text = path.read_text()
    original = text
    for old, new in replacements.items():
        text = text.replace(old, new)
    if text == original:
        raise SystemExit(f"No minimap pipeline replacement made in {path}")
    path.write_text(text)
PY

export REPO_DIR
export PATCHED_SHELLS_DIR
export SHELLS_DIR="${PATCHED_SHELLS_DIR}"
export THREADS
export MIN_MAPQ
export MIN_ALIGN_LEN
export MINIMAP_DB_FASTA="${FOCUSED_DB_FASTA}"

log_info "Running single-genome read-level benchmark with focused minimap database"
export OUT_DIR="${NEW_RUNS_DIR}/spikein_single_read_minimapfix_focused_${STAMP}"
bash "${PATCHED_SHELLS_DIR}/run_spikein_single_readlevel.sh"

log_info "Running shuffled negative-control read-level benchmark with focused minimap database"
export OUT_DIR="${NEW_RUNS_DIR}/spikein_shuffle_read_minimapfix_focused_${STAMP}"
bash "${PATCHED_SHELLS_DIR}/run_spikein_shuffled_readlevel.sh"

log_info "Running two-genome read-level benchmark with focused minimap database"
export PATHOGEN_CONFIG_TSV="${REPO_DIR}/configs/pathogen_panel_2.tsv"
export OUT_DIR="${NEW_RUNS_DIR}/spikein_multi_read_panel2_minimapfix_focused_${STAMP}"
bash "${PATCHED_SHELLS_DIR}/run_spikein_multi_readlevel.sh"

log_info "Running three-genome read-level benchmark with focused minimap database"
export PATHOGEN_CONFIG_TSV="${REPO_DIR}/configs/pathogen_panel_3.tsv"
export OUT_DIR="${NEW_RUNS_DIR}/spikein_multi_read_panel3_minimapfix_focused_${STAMP}"
bash "${PATCHED_SHELLS_DIR}/run_spikein_multi_readlevel.sh"

log_info "Staging run directories for merged summary"
find "${PROJECT_BASE}/runs" -mindepth 1 -maxdepth 1 -type d \
    ! -name 'spikein_single_read_*' \
    ! -name 'spikein_shuffle_read_*' \
    ! -name 'spikein_multi_read_panel2_*' \
    ! -name 'spikein_multi_read_panel3_*' \
    -print0 | while IFS= read -r -d '' existing_run; do
    ln -s "${existing_run}" "${SUMMARY_STAGE_DIR}/$(basename "${existing_run}")"
done

find "${NEW_RUNS_DIR}" -mindepth 1 -maxdepth 1 -type d ! -name patched_shells -print0 | while IFS= read -r -d '' new_run; do
    ln -s "${new_run}" "${SUMMARY_STAGE_DIR}/$(basename "${new_run}")"
done

log_info "Running merged summary build"
python3 "${REPO_DIR}/summary/summarise_spikein_runs_v4_reported_taxa_full.py" \
    --input_dirs "${SUMMARY_STAGE_DIR}" \
    --out_dir "${SUMMARY_OUT_DIR}" \
    --verbose

python3 "${REPO_DIR}/summary/build_method_performance_table_with_reported_taxa.py" \
    --combined_long_tsv "${SUMMARY_OUT_DIR}/combined_long.tsv" \
    --reported_taxa_long_tsv "${SUMMARY_OUT_DIR}/reported_taxa_long.tsv" \
    --out_dir "${SUMMARY_OUT_DIR}" \
    --threshold_mode fixed \
    --min_detect_value 1

python3 "${REPO_DIR}/summary/make_spikein_report_v5.py" \
    --summary_dir "${SUMMARY_OUT_DIR}" \
    --title "ONT spike-in summary report: minimap rerun focused DB Q15 len500"

python3 "${REPO_DIR}/summary/make_spikein_replicate_report_v2.py" \
    --summary_dir "${SUMMARY_OUT_DIR}" \
    --out_dir "${SUMMARY_OUT_DIR}/replicate_resolved_report" \
    --title "ONT spike-in replicate-resolved report: minimap rerun focused DB Q15 len500" \
    --threshold_mode fixed \
    --min_detect_value 1

python3 "${REPO_DIR}/summary/make_spikein_threshold_calibration_report_v3.py" \
    --summary_dir "${SUMMARY_OUT_DIR}" \
    --out_dir "${SUMMARY_OUT_DIR}/threshold_calibration_report_v3" \
    --title "ONT spike-in threshold calibration report: minimap rerun focused DB Q15 len500" \
    --min_detect_value 1 \
    --target_fpr 0.05

log_info "Finished focused-database minimap rerun"
log_info "Run directory: ${NEW_RUNS_DIR}"
log_info "Summary directory: ${SUMMARY_OUT_DIR}"
