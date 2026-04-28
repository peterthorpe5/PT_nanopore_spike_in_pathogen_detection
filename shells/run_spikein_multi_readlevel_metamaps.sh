#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 240G
#$ -adds l_hard h_vmem 240G
#$ -N MMspike_multi_mm

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

PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR:-${REPO_DIR}/scripts}"
CONFIG_DIR="${CONFIG_DIR:-${REPO_DIR}/configs}"
PATHOGEN_CONFIG_TSV="${PATHOGEN_CONFIG_TSV:-${DEFAULT_PATHOGEN_PANEL_3_WITH_TAXID}}"
SAMPLE_FASTQ_PY="${SAMPLE_FASTQ_PY:-${PY_SCRIPTS_DIR}/sample_fastq.py}"
BUILD_MIXED_FASTQ_PY="${BUILD_MIXED_FASTQ_PY:-${PY_SCRIPTS_DIR}/build_mixed_fastq.py}"
COMBINE_NANOSIM_FASTQ_PY="${COMBINE_NANOSIM_FASTQ_PY:-${PY_SCRIPTS_DIR}/combine_nanosim_fastq.py}"
SUMMARISE_METAMAPS_WIMP_PY="${SUMMARISE_METAMAPS_WIMP_PY:-${PY_SCRIPTS_DIR}/summarise_metamaps_wimp.py}"

REAL_FASTQ="${REAL_FASTQ:-${REAL_FASTQ_DEFAULT}}"
MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-${MONKEY_SMALL_GZ_DEFAULT}}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-${MONKEY_SMALL_FASTA_DEFAULT}}"
DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-${DEPLETION_REF_FASTA_DEFAULT}}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-${DEPLETED_FASTQ_DEFAULT}}"
METAMAPS_DB_DIR="${METAMAPS_DB_DIR:-${METAMAPS_DB_DIR_DEFAULT}}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_metamaps_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-12}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-12}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"
MAX_MEMORY_GB="${MAX_MEMORY_GB:-140}"

for x in python3 minimap2 samtools read_analysis.py simulator.py metamaps; do require_exe "$x"; done
for f in "${SAMPLE_FASTQ_PY}" "${BUILD_MIXED_FASTQ_PY}" "${COMBINE_NANOSIM_FASTQ_PY}" "${SUMMARISE_METAMAPS_WIMP_PY}" "${REAL_FASTQ}" "${PATHOGEN_CONFIG_TSV}"; do require_file "$f"; done
require_dir "${METAMAPS_DB_DIR}"
require_file "${METAMAPS_DB_DIR}/DB.fa"
require_dir "${METAMAPS_DB_DIR}/taxonomy"
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

mapfile -t PANEL_LINES < <(awk 'BEGIN{FS="\t"} NR>1 && NF>=2 {print $1"\t"$2}' "${PATHOGEN_CONFIG_TSV}")
[[ ${#PANEL_LINES[@]} -gt 0 ]] || { log_error "No pathogen entries found in ${PATHOGEN_CONFIG_TSV}"; exit 1; }

TRAIN_FASTQ="${OUT_DIR}/train_reads.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"
SUMMARY_TSV="${OUT_DIR}/spikein_multi_metamaps_summary.tsv"
WORK_FASTQ="${DEPLETED_FASTQ}"

python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${REAL_FASTQ}" --n_reads "${TRAIN_READS_N}" --seed 1 --out_fastq_gz "${TRAIN_FASTQ}"
if [[ ! -s "${MONKEY_SMALL_FASTA}" ]]; then gzip -cd "${MONKEY_SMALL_GZ}" > "${MONKEY_SMALL_FASTA}"; fi
read_analysis.py genome --read "${TRAIN_FASTQ}" --ref_g "${MONKEY_SMALL_FASTA}" --output "${NS_MODEL_PREFIX}" --num_threads "${THREADS}" --fastq

GENOME_INDEX=0
for line in "${PANEL_LINES[@]}"; do
    GENOME_INDEX=$((GENOME_INDEX + 1))
    pathogen_fasta="$(printf '%s' "${line}" | cut -f1)"
    target_label="$(printf '%s' "${line}" | cut -f2)"
    require_file "${pathogen_fasta}"
    sim_prefix="${OUT_DIR}/simulated_pathogen_${GENOME_INDEX}"
    sim_pool_fastq="${OUT_DIR}/sim_pool_${GENOME_INDEX}.fastq"
    sim_pool_fastq_gz="${OUT_DIR}/sim_pool_${GENOME_INDEX}.fastq.gz"
    simulator.py genome --ref_g "${pathogen_fasta}" --model_prefix "${NS_MODEL_PREFIX}" --output "${sim_prefix}" --number "${SIM_POOL_N}" -dna_type linear --seed "$((40 + GENOME_INDEX))" --num_threads "${THREADS}" --fastq
    python3 "${COMBINE_NANOSIM_FASTQ_PY}" --sim_prefix "${sim_prefix}" --out_fastq "${sim_pool_fastq}"
    gzip -c "${sim_pool_fastq}" > "${sim_pool_fastq_gz}"
done

if [[ -s "${WORK_FASTQ}" ]]; then :; elif [[ "${DO_HOST_DEPLETION}" == "true" ]]; then minimap2 -t "${THREADS}" -a -x map-ont "${DEPLETION_REF_FASTA}" "${REAL_FASTQ}" | samtools view -h -T "${DEPLETION_REF_FASTA}" -b -f 4 -@ "${THREADS}" - | samtools fastq -@ "${THREADS}" - | gzip -c > "${WORK_FASTQ}"; else cp "${REAL_FASTQ}" "${WORK_FASTQ}"; fi

header='replicate\tspike_n\ttotal_spiked_reads\tmixed_fastq_gz\tmetamaps_prefix\tmetamaps_wimp\tmetamaps_total_reads_definedGenomes'
for line in "${PANEL_LINES[@]}"; do target_label="$(printf '%s' "${line}" | cut -f2)"; safe_label="$(printf '%s' "${target_label}" | tr ' ' '_' | tr '/' '_')"; header+="\tmetamaps_${safe_label}_reads"; done
printf '%b\n' "${header}" > "${SUMMARY_TSV}"

unset MALLOC_ARENA_MAX || true

for rep in $(seq 1 "${REPLICATES}"); do
  for spike_n in ${SPIKE_LEVELS}; do
    MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"
    mkdir -p "${MIX_DIR}"

    log_info "Processing ${MIX_DIR}"

    tmp_spikes=()
    total_spiked=0
    GENOME_INDEX=0

    for line in "${PANEL_LINES[@]}"; do
      GENOME_INDEX=$((GENOME_INDEX + 1))
      sim_pool_fastq_gz="${OUT_DIR}/sim_pool_${GENOME_INDEX}.fastq.gz"
      spike_fastq_gz="${MIX_DIR}/spike_${GENOME_INDEX}.fastq.gz"
      spike_seed=$(( rep * 100000 + GENOME_INDEX * 1000 + spike_n ))

      python3 "${SAMPLE_FASTQ_PY}" \
        --fastq_gz "${sim_pool_fastq_gz}" \
        --n_reads "${spike_n}" \
        --seed "${spike_seed}" \
        --out_fastq_gz "${spike_fastq_gz}"

      tmp_spikes+=("${spike_fastq_gz}")
      total_spiked=$(( total_spiked + spike_n ))
    done

    combined_spike="${MIX_DIR}/spike_combined.fastq.gz"

    mix_args=(--background_fastq_gz "${tmp_spikes[0]}")
    for s in "${tmp_spikes[@]:1}"; do
      mix_args+=(--spike_fastq_gz "${s}")
    done
    mix_args+=(--out_fastq_gz "${combined_spike}")

    python3 "${BUILD_MIXED_FASTQ_PY}" "${mix_args[@]}"

    mix_fastq_gz="${MIX_DIR}/mixed.fastq.gz"
    python3 "${BUILD_MIXED_FASTQ_PY}" \
      --background_fastq_gz "${WORK_FASTQ}" \
      --spike_fastq_gz "${combined_spike}" \
      --out_fastq_gz "${mix_fastq_gz}"

    metamaps_prefix="${MIX_DIR}/metamaps"
    wimp_tsv="${metamaps_prefix}.EM.WIMP"
    wimp_total="NA"

    target_values=()

    if metamaps mapDirectly \
        -t "${THREADS}" \
        --all \
        --maxmemory "${MAX_MEMORY_GB}" \
        -r "${METAMAPS_DB_DIR}/DB.fa" \
        -q "${mix_fastq_gz}" \
        -o "${metamaps_prefix}"; then

        unset MALLOC_ARENA_MAX

        if metamaps classify \
            -t 1 \
            --mappings "${metamaps_prefix}" \
            --DB "${METAMAPS_DB_DIR}"; then

            if [[ -s "${wimp_tsv}" ]]; then
                for line in "${PANEL_LINES[@]}"; do
                    target_label="$(printf '%s' "${line}" | cut -f2)"
                    safe_label="$(printf '%s' "${target_label}" | tr ' ' '_' | tr '/' '_')"
                    wsum="${MIX_DIR}/metamaps.${safe_label}.summary.tsv"

                    python3 "${SUMMARISE_METAMAPS_WIMP_PY}" \
                        --wimp_tsv "${wimp_tsv}" \
                        --analysis_level definedGenomes \
                        --target_label "${target_label}" \
                        --out_tsv "${wsum}" || true

                    if [[ "${wimp_total}" == "NA" && -s "${wsum}" ]]; then
                        wimp_total="$(awk 'NR==2{print $3}' "${wsum}")"
                    fi

                    if [[ -s "${wsum}" ]]; then
                        target_values+=("$(awk 'NR==2{print $5}' "${wsum}")")
                    else
                        target_values+=("NA")
                    fi
                done
            else
                log_error "Missing or empty WIMP output for ${MIX_DIR}"
                for line in "${PANEL_LINES[@]}"; do
                    target_values+=("NA")
                done
            fi
        else
            log_error "MetaMaps classify failed for ${MIX_DIR}"
            for line in "${PANEL_LINES[@]}"; do
                target_values+=("NA")
            done
        fi
    else
        log_error "MetaMaps mapDirectly failed for ${MIX_DIR}"
        for line in "${PANEL_LINES[@]}"; do
            target_values+=("NA")
        done
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s' \
        "${rep}" \
        "${spike_n}" \
        "${total_spiked}" \
        "${mix_fastq_gz}" \
        "${metamaps_prefix}" \
        "${wimp_tsv}" \
        "${wimp_total}" \
        >> "${SUMMARY_TSV}"

    for value in "${target_values[@]}"; do
        printf '\t%s' "${value}" >> "${SUMMARY_TSV}"
    done
    printf '\n' >> "${SUMMARY_TSV}"
  done
done

log_info "Done: ${SUMMARY_TSV}"
