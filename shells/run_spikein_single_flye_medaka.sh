#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N spike_single_flye

set -euo pipefail

log_info() { printf '[INFO] %s\n' "$*"; }
log_error() { printf '[ERROR] %s\n' "$*" >&2; }
require_file() { [[ -f "$1" ]] || { log_error "Required file not found: $1"; exit 1; }; }
require_dir() { [[ -d "$1" ]] || { log_error "Required directory not found: $1"; exit 1; }; }
require_exe() { command -v "$1" >/dev/null 2>&1 || { log_error "Required executable not found on PATH: $1"; exit 1; }; }

run_medaka_consensus() {
  local reads_fastq="$1"
  local draft_fasta="$2"
  local medaka_out_dir="$3"
  local threads="$4"
  local model="$5"

  mkdir -p "${medaka_out_dir}"

  if [[ "${model}" == "auto" ]]; then
    log_info "Running Medaka with automatic model selection"
    medaka_consensus \
      -i "${reads_fastq}" \
      -d "${draft_fasta}" \
      -o "${medaka_out_dir}" \
      -t "${threads}"
  else
    log_info "Running Medaka with explicit model: ${model}"
    medaka_consensus \
      -i "${reads_fastq}" \
      -d "${draft_fasta}" \
      -o "${medaka_out_dir}" \
      -t "${threads}" \
      -m "${model}"
  fi
}

if [[ -n "${REPO_DIR:-}" ]]; then REPO_DIR="${REPO_DIR}"; elif [[ -n "${SGE_O_WORKDIR:-}" ]]; then REPO_DIR="${SGE_O_WORKDIR}"; else REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection"; fi

source "${REPO_DIR}/configs/pipeline_paths.sh"

PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR:-${REPO_DIR}/scripts}"
SAMPLE_FASTQ_PY="${PY_SCRIPTS_DIR}/sample_fastq.py"
BUILD_MIXED_FASTQ_PY="${PY_SCRIPTS_DIR}/build_mixed_fastq.py"
COMBINE_NANOSIM_FASTQ_PY="${PY_SCRIPTS_DIR}/combine_nanosim_fastq.py"
SUMMARISE_KRAKEN_PY="${PY_SCRIPTS_DIR}/summarise_kraken_report_with_reported_taxa.py"
DEDUP_FASTQ_NAMES_PY="${PY_SCRIPTS_DIR}/dedup_fastq_names.py"
ASSEMBLY_STATS_PY="${PY_SCRIPTS_DIR}/assembly_stats.py"

REAL_FASTQ="${REAL_FASTQ:-${REAL_FASTQ_DEFAULT}}"
MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-${MONKEY_SMALL_GZ_DEFAULT}}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-${MONKEY_SMALL_FASTA_DEFAULT}}"
DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-${DEPLETION_REF_FASTA_DEFAULT}}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-${DEPLETED_FASTQ_DEFAULT}}"
PATHOGEN_FASTA="${PATHOGEN_FASTA:-${DEFAULT_PATHOGEN_FASTA}}"
TARGET_LABEL="${TARGET_LABEL:-${DEFAULT_TARGET_LABEL}}"
KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-${KRAKEN_DB_DIR_DEFAULT}}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_single_flye_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-12}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-12}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"
MEDAKA_MODEL="${MEDAKA_MODEL:-auto}"

for x in python3 read_analysis.py simulator.py flye kraken2 gzip medaka_consensus; do require_exe "$x"; done
for f in "${SAMPLE_FASTQ_PY}" "${BUILD_MIXED_FASTQ_PY}" "${COMBINE_NANOSIM_FASTQ_PY}" "${SUMMARISE_KRAKEN_PY}" "${DEDUP_FASTQ_NAMES_PY}" "${ASSEMBLY_STATS_PY}" "${REAL_FASTQ}" "${PATHOGEN_FASTA}"; do require_file "$f"; done
require_dir "${KRAKEN_DB_DIR}"
mkdir -p "${OUT_DIR}"


########################################################################
# Job-local working directory
########################################################################
# The permanent output path requested by the caller is kept as FINAL_OUT_DIR.
# Heavy temporary files are written under $TMPDIR, then selected result files
# are copied back when the script exits. Set COPY_WORKING_FILES=true to copy
# everything, including FASTQ, BAM, FASTA, Flye, and Medaka intermediates.
# Set USE_TMPDIR=false to bypass this behaviour for local debugging.
FINAL_OUT_DIR="${OUT_DIR}"
COPY_WORKING_FILES="${COPY_WORKING_FILES:-false}"
USE_TMPDIR="${USE_TMPDIR:-true}"

copy_results_back() {
    local exit_code="$?"

    if [[ "${USE_TMPDIR}" != "true" ]]; then
        exit "${exit_code}"
    fi

    if [[ -z "${WORK_OUT_DIR:-}" || -z "${FINAL_OUT_DIR:-}" ]]; then
        exit "${exit_code}"
    fi

    mkdir -p "${FINAL_OUT_DIR}"

    if [[ -d "${WORK_OUT_DIR}" ]]; then
        log_info "Copying results from ${WORK_OUT_DIR} to ${FINAL_OUT_DIR}"

        if [[ "${COPY_WORKING_FILES}" == "true" ]]; then
            rsync -a "${WORK_OUT_DIR}/" "${FINAL_OUT_DIR}/" || true
        else
            rsync -a \
                --include='*/' \
                --include='*.tsv' \
                --include='*.txt' \
                --include='*.log' \
                --include='*.err' \
                --include='*.out' \
                --include='*.json' \
                --include='*.html' \
                --include='*.xlsx' \
                --include='*.pdf' \
                --include='*.svg' \
                --include='*.png' \
                --include='*.report' \
                --include='*.report.tsv' \
                --include='*.summary.tsv' \
                --include='*.reported_taxa.tsv' \
                --include='*.reported_taxa_long.tsv' \
                --include='*.meta.tsv' \
                --include='*.target_summary.tsv' \
                --include='*.classifications.tsv' \
                --exclude='*.fastq' \
                --exclude='*.fq' \
                --exclude='*.fastq.gz' \
                --exclude='*.fq.gz' \
                --exclude='*.bam' \
                --exclude='*.bai' \
                --exclude='*.sam' \
                --exclude='*.fasta' \
                --exclude='*.fa' \
                --exclude='*.fna' \
                --exclude='flye_out/***' \
                --exclude='medaka_out/***' \
                --exclude='medaka_round1_out/***' \
                --exclude='medaka_round2_out/***' \
                --exclude='simulated_pathogen*/***' \
                --exclude='sim_pool*' \
                --exclude='nanosim_training*' \
                --exclude='*' \
                "${WORK_OUT_DIR}/" "${FINAL_OUT_DIR}/" || true
        fi

        # Rewrite result paths so downstream summary scripts see permanent paths,
        # not job-local $TMPDIR paths that disappear after the job exits.
        if command -v python3 >/dev/null 2>&1; then
            python3 - "${FINAL_OUT_DIR}" "${WORK_OUT_DIR}" "${FINAL_OUT_DIR}" <<'PYTHON_REWRITE_PATHS'
from pathlib import Path
import sys

final_dir = Path(sys.argv[1])
old = sys.argv[2]
new = sys.argv[3]

for path in final_dir.rglob('*'):
    if not path.is_file():
        continue
    if path.suffix not in {'.tsv', '.txt', '.log', '.json', '.html'}:
        continue
    try:
        text = path.read_text(encoding='utf-8')
    except UnicodeDecodeError:
        continue
    updated = text.replace(old, new)
    if updated != text:
        path.write_text(updated, encoding='utf-8')
PYTHON_REWRITE_PATHS
        fi
    fi

    exit "${exit_code}"
}

if [[ "${USE_TMPDIR}" == "true" ]]; then
    if [[ -z "${TMPDIR:-}" || ! -d "${TMPDIR}" ]]; then
        log_error "TMPDIR is not set or does not exist. Set USE_TMPDIR=false only for local debugging."
        exit 1
    fi

    RUN_STAMP="$(basename "${FINAL_OUT_DIR}")"
    WORK_OUT_DIR="${TMPDIR}/${USER:-user}_${JOB_ID:-manual}_${RUN_STAMP}"
    mkdir -p "${WORK_OUT_DIR}"
    mkdir -p "${FINAL_OUT_DIR}"
    OUT_DIR="${WORK_OUT_DIR}"
    trap copy_results_back EXIT

    log_info "Permanent output directory: ${FINAL_OUT_DIR}"
    log_info "Job-local working directory: ${WORK_OUT_DIR}"
fi
########################################################################

TRAIN_FASTQ="${OUT_DIR}/train_reads.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"
SIM_PREFIX="${OUT_DIR}/simulated_pathogen"
SIM_POOL_FASTQ="${OUT_DIR}/sim_pool.fastq"
SIM_POOL_FASTQ_GZ="${OUT_DIR}/sim_pool.fastq.gz"
SUMMARY_TSV="${OUT_DIR}/spikein_flye_summary.tsv"
WORK_FASTQ="${DEPLETED_FASTQ}"

python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${REAL_FASTQ}" --n_reads "${TRAIN_READS_N}" --seed 1 --out_fastq_gz "${TRAIN_FASTQ}"
if [[ ! -s "${MONKEY_SMALL_FASTA}" ]]; then gzip -cd "${MONKEY_SMALL_GZ}" > "${MONKEY_SMALL_FASTA}"; fi
read_analysis.py genome --read "${TRAIN_FASTQ}" --ref_g "${MONKEY_SMALL_FASTA}" --output "${NS_MODEL_PREFIX}" --num_threads "${THREADS}" --fastq
simulator.py genome --ref_g "${PATHOGEN_FASTA}" --model_prefix "${NS_MODEL_PREFIX}" --output "${SIM_PREFIX}" --number "${SIM_POOL_N}" -dna_type linear --seed 42 --num_threads "${THREADS}" --fastq
python3 "${COMBINE_NANOSIM_FASTQ_PY}" --sim_prefix "${SIM_PREFIX}" --out_fastq "${SIM_POOL_FASTQ}"
gzip -c "${SIM_POOL_FASTQ}" > "${SIM_POOL_FASTQ_GZ}"

if [[ -s "${WORK_FASTQ}" ]]; then :; elif [[ "${DO_HOST_DEPLETION}" == "true" ]]; then minimap2 -t "${THREADS}" -a -x map-ont "${DEPLETION_REF_FASTA}" "${REAL_FASTQ}" | samtools view -h -T "${DEPLETION_REF_FASTA}" -b -f 4 -@ "${THREADS}" - | samtools fastq -@ "${THREADS}" - | gzip -c > "${WORK_FASTQ}"; else cp "${REAL_FASTQ}" "${WORK_FASTQ}"; fi



printf 'replicate\tspike_n\tmixed_fastq_gz\tflye_out_dir\tflye_assembly_fasta\tpolished_assembly_fasta\tcontig_count\ttotal_bases\tkraken_report\tkraken_classified_contigs\tkraken_target_contigs\n' > "${SUMMARY_TSV}"

for rep in $(seq 1 "${REPLICATES}"); do
  for spike_n in ${SPIKE_LEVELS}; do
    MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"
    mkdir -p "${MIX_DIR}"

    SPIKE_FASTQ_GZ="${MIX_DIR}/spike.fastq.gz"
    MIX_FASTQ_GZ="${MIX_DIR}/mixed.fastq.gz"
    DEDUP_FASTQ="${MIX_DIR}/mixed.dedup.fastq"
    FLYE_OUT_DIR="${MIX_DIR}/flye_out"
    ASSEMBLY_FASTA="${FLYE_OUT_DIR}/assembly.fasta"
    MEDAKA_OUT_DIR="${MIX_DIR}/medaka_out"
    POLISHED_FASTA="${MEDAKA_OUT_DIR}/consensus.fasta"
    ASM_STATS_TSV="${MIX_DIR}/assembly.stats.tsv"
    KRAKEN_REPORT="${MIX_DIR}/assembly.kraken.report.tsv"
    KRAKEN_OUT="${MIX_DIR}/assembly.kraken.classifications.tsv"
    KRAKEN_SUM_TSV="${MIX_DIR}/assembly.kraken.summary.tsv"

    SPIKE_SEED=$(( rep * 100000 + spike_n ))
    python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${SIM_POOL_FASTQ_GZ}" --n_reads "${spike_n}" --seed "${SPIKE_SEED}" --out_fastq_gz "${SPIKE_FASTQ_GZ}"
    python3 "${BUILD_MIXED_FASTQ_PY}" --background_fastq_gz "${WORK_FASTQ}" --spike_fastq_gz "${SPIKE_FASTQ_GZ}" --out_fastq_gz "${MIX_FASTQ_GZ}"
    python3 "${DEDUP_FASTQ_NAMES_PY}" --input_fastq_gz "${MIX_FASTQ_GZ}" --output_fastq "${DEDUP_FASTQ}"

    flye --meta --nano-hq "${DEDUP_FASTQ}" --out-dir "${FLYE_OUT_DIR}" --threads "${THREADS}"
    run_medaka_consensus "${DEDUP_FASTQ}" "${ASSEMBLY_FASTA}" "${MEDAKA_OUT_DIR}" "${THREADS}" "${MEDAKA_MODEL}"

    require_file "${POLISHED_FASTA}"
    python3 "${ASSEMBLY_STATS_PY}" --assembly_fasta "${POLISHED_FASTA}" --out_tsv "${ASM_STATS_TSV}"
    kraken2 --db "${KRAKEN_DB_DIR}" --threads "${THREADS}" --report "${KRAKEN_REPORT}" --output "${KRAKEN_OUT}" "${POLISHED_FASTA}" >/dev/null
    python3 "${SUMMARISE_KRAKEN_PY}" --kraken_report_tsv "${KRAKEN_REPORT}" --target_label "${TARGET_LABEL}" --out_tsv "${KRAKEN_SUM_TSV}"

    CONTIG_COUNT="$(awk 'NR==2{print $2}' "${ASM_STATS_TSV}")"
    TOTAL_BASES="$(awk 'NR==2{print $3}' "${ASM_STATS_TSV}")"
    K_CLASSIFIED="$(awk 'NR==2{print $2}' "${KRAKEN_SUM_TSV}")"
    K_TARGET="$(awk 'NR==2{print $3}' "${KRAKEN_SUM_TSV}")"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${rep}" "${spike_n}" "${MIX_FASTQ_GZ}" "${FLYE_OUT_DIR}" "${ASSEMBLY_FASTA}" "${POLISHED_FASTA}" \
      "${CONTIG_COUNT}" "${TOTAL_BASES}" "${KRAKEN_REPORT}" "${K_CLASSIFIED}" "${K_TARGET}" >> "${SUMMARY_TSV}"
  done
done

log_info "Done: ${SUMMARY_TSV}"
