#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 2
#$ -jc long
#$ -N summary_old_new

set -euo pipefail

PROJECT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity"
REPO_DIR="${PROJECT_DIR}/PT_nanopore_spike_in_pathogen_detection"
SUMMARY_DIR="${REPO_DIR}/summary"
OLD_SUMMARY_DIR="${SUMMARY_DIR}/old"
RUNS_DIR="${PROJECT_DIR}/runs"
OLD_OUT_DIR="${PROJECT_DIR}/spikein_summary_report_old"
NEW_OUT_DIR="${PROJECT_DIR}/spikein_summary_report_new"

OLD_SUMMARISER="${OLD_SUMMARY_DIR}/summarise_spikein_runs_v4_reported_taxa_full.py"
NEW_SUMMARISER="${SUMMARY_DIR}/summarise_spikein_runs_v4_reported_taxa_full_minimap_patched.py"
METHOD_TABLE_SCRIPT="${OLD_SUMMARY_DIR}/build_method_performance_table_with_reported_taxa.py"
MAIN_REPORT_SCRIPT="${OLD_SUMMARY_DIR}/make_spikein_report_v5.py"
REPLICATE_REPORT_SCRIPT="${OLD_SUMMARY_DIR}/make_spikein_replicate_report_v2.py"
THRESHOLD_REPORT_SCRIPT="${OLD_SUMMARY_DIR}/make_spikein_threshold_calibration_report_v3.py"
REAL_WORLD_REPORT_SCRIPT="${OLD_SUMMARY_DIR}/build_combined_real_world_report.py"

log_info() {
    printf '[INFO] %s\n' "$*"
}

log_error() {
    printf '[ERROR] %s\n' "$*" >&2
}

require_file() {
    [[ -f "$1" ]] || {
        log_error "Required file not found: $1"
        exit 1
    }
}

cd "${PROJECT_DIR}"

mkdir -p "${RUNS_DIR}"

shopt -s nullglob
for run_dir in runs_*; do
    [[ -d "${run_dir}" ]] || continue
    log_info "Moving ${run_dir} into ${RUNS_DIR}"
    mv "${run_dir}" "${RUNS_DIR}/"
done
shopt -u nullglob

require_file "${OLD_SUMMARISER}"
require_file "${NEW_SUMMARISER}"
require_file "${METHOD_TABLE_SCRIPT}"
require_file "${MAIN_REPORT_SCRIPT}"
require_file "${REPLICATE_REPORT_SCRIPT}"
require_file "${THRESHOLD_REPORT_SCRIPT}"
require_file "${REAL_WORLD_REPORT_SCRIPT}"

log_info "Removing previous output directories"
rm -rf "${OLD_OUT_DIR}" "${NEW_OUT_DIR}"

log_info "Running OLD summary workflow"
python "${OLD_SUMMARISER}" \
    --input_dirs "${RUNS_DIR}" \
    --out_dir "${OLD_OUT_DIR}" \
    --verbose

python "${METHOD_TABLE_SCRIPT}" \
    --combined_long_tsv "${OLD_OUT_DIR}/combined_long.tsv" \
    --out_dir "${OLD_OUT_DIR}" \
    --reported_taxa_long_tsv "${OLD_OUT_DIR}/reported_taxa_long.tsv" \
    --threshold_mode fixed \
    --min_detect_value 1

python "${MAIN_REPORT_SCRIPT}" \
    --summary_dir "${OLD_OUT_DIR}" \
    --title "ONT spike-in summary report (old workflow)"

python "${REPLICATE_REPORT_SCRIPT}" \
    --summary_dir "${OLD_OUT_DIR}" \
    --out_dir "${OLD_OUT_DIR}/replicate_resolved_report" \
    --title "ONT spike-in replicate-resolved report (old workflow)" \
    --threshold_mode fixed \
    --min_detect_value 1

python "${THRESHOLD_REPORT_SCRIPT}" \
    --summary_dir "${OLD_OUT_DIR}" \
    --out_dir "${OLD_OUT_DIR}/threshold_calibration_report_v3" \
    --title "ONT spike-in threshold calibration report (old workflow)" \
    --min_detect_value 1 \
    --target_fpr 0.05

python "${REAL_WORLD_REPORT_SCRIPT}" \
    --method_performance_xlsx "${OLD_OUT_DIR}/method_performance.xlsx" \
    --replicate_report_xlsx "${OLD_OUT_DIR}/replicate_resolved_report/replicate_resolved_report.xlsx" \
    --threshold_report_xlsx "${OLD_OUT_DIR}/threshold_calibration_report_v3/threshold_calibration_report.xlsx" \
    --out_dir "${OLD_OUT_DIR}/combined_real_world_report" \
    --report_title "Combined ONT spike-in benchmark report with real-world taxonomic burden (old workflow)"


##########################################################################
log_info "Running NEW summary workflow"


##########################################################################
log_info "Running NEW summary workflow"

MINIMAP_INPUT_ROOTS=()
while IFS= read -r d; do
    MINIMAP_INPUT_ROOTS+=("${d}")
done < <(
    find "${RUNS_DIR}" -maxdepth 1 -mindepth 1 -type d \
        \( -name "*minimap*maskedref*" -o -name "*minimap*focused*" -o -name "*minimap*shared*" -o -name "*minimap*fullref*" \) \
        ! -name "*broken*" \
        | sort
)

python "${SUMMARY_DIR}/minimap_specific_summary_updated.py" \
    --input_roots "${MINIMAP_INPUT_ROOTS[@]}" \
    --out_dir "${NEW_OUT_DIR}/minimap_specific_summary" \
    --single_target_label "Plasmodium vivax" \
    --panel2_tsv "${REPO_DIR}/configs/pathogen_panel_2.tsv" \
    --panel3_tsv "${REPO_DIR}/configs/pathogen_panel_3.tsv" \
    --target_threshold_alignments 1 \
    --target_threshold_unique_reads 1 \
    --real_world_threshold_alignments 1 \
    --real_world_threshold_unique_reads 1

python "${NEW_SUMMARISER}" \
    --input_dirs "${RUNS_DIR}" \
    --out_dir "${NEW_OUT_DIR}" \
    --verbose

python "${METHOD_TABLE_SCRIPT}" \
    --combined_long_tsv "${NEW_OUT_DIR}/combined_long.tsv" \
    --out_dir "${NEW_OUT_DIR}" \
    --reported_taxa_long_tsv "${NEW_OUT_DIR}/reported_taxa_long.tsv" \
    --threshold_mode fixed \
    --min_detect_value 1

python "${MAIN_REPORT_SCRIPT}" \
    --summary_dir "${NEW_OUT_DIR}" \
    --title "ONT spike-in summary report (new workflow)"

python "${REPLICATE_REPORT_SCRIPT}" \
    --summary_dir "${NEW_OUT_DIR}" \
    --out_dir "${NEW_OUT_DIR}/replicate_resolved_report" \
    --title "ONT spike-in replicate-resolved report (new workflow)" \
    --threshold_mode fixed \
    --min_detect_value 1

python "${THRESHOLD_REPORT_SCRIPT}" \
    --summary_dir "${NEW_OUT_DIR}" \
    --out_dir "${NEW_OUT_DIR}/threshold_calibration_report_v3" \
    --title "ONT spike-in threshold calibration report (new workflow)" \
    --min_detect_value 1 \
    --target_fpr 0.05

python "${REAL_WORLD_REPORT_SCRIPT}" \
    --method_performance_xlsx "${NEW_OUT_DIR}/method_performance.xlsx" \
    --replicate_report_xlsx "${NEW_OUT_DIR}/replicate_resolved_report/replicate_resolved_report.xlsx" \
    --threshold_report_xlsx "${NEW_OUT_DIR}/threshold_calibration_report_v3/threshold_calibration_report.xlsx" \
    --out_dir "${NEW_OUT_DIR}/combined_real_world_report" \
    --report_title "Combined ONT spike-in benchmark report with real-world taxonomic burden (new workflow)"

log_info "Done"
log_info "Old outputs: ${OLD_OUT_DIR}"
log_info "New outputs: ${NEW_OUT_DIR}"


PLOT_PUBLICATION_SCRIPT="${SUMMARY_DIR}/plot_publication_summary_figures.py"


python "${PLOT_PUBLICATION_SCRIPT}" \
    --method_performance_xlsx "${NEW_OUT_DIR}/method_performance.xlsx" \
    --real_world_xlsx "${NEW_OUT_DIR}/combined_real_world_report/combined_real_world_report.xlsx" \
    --replicate_report_xlsx "${NEW_OUT_DIR}/replicate_resolved_report/replicate_resolved_report.xlsx" \
    --threshold_report_xlsx "${NEW_OUT_DIR}/threshold_calibration_report_v3/threshold_calibration_report.xlsx" \
    --minimap_tracked_tsv "${NEW_OUT_DIR}/minimap_specific_summary/minimap_tracked_target_performance.tsv" \
    --minimap_real_world_tsv "${NEW_OUT_DIR}/minimap_specific_summary/minimap_real_world_reference_summary.tsv" \
    --out_dir "${NEW_OUT_DIR}/publication_summary_plots"
