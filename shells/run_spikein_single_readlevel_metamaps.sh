#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 200G
#$ -adds l_hard h_vmem 200G
#$ -N MMspike_single_mm

set -euo pipefail

log_info() { printf '[INFO] %s\n' "$*"; }
log_warn() { printf '[WARN] %s\n' "$*" >&2; }
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
SAMPLE_FASTQ_PY="${SAMPLE_FASTQ_PY:-${PY_SCRIPTS_DIR}/sample_fastq.py}"
BUILD_MIXED_FASTQ_PY="${BUILD_MIXED_FASTQ_PY:-${PY_SCRIPTS_DIR}/build_mixed_fastq.py}"
COMBINE_NANOSIM_FASTQ_PY="${COMBINE_NANOSIM_FASTQ_PY:-${PY_SCRIPTS_DIR}/combine_nanosim_fastq.py}"
SUMMARISE_METAMAPS_WIMP_PY="${SUMMARISE_METAMAPS_WIMP_PY:-${PY_SCRIPTS_DIR}/summarise_metamaps_wimp.py}"

REAL_FASTQ="${REAL_FASTQ:-${REAL_FASTQ_DEFAULT}}"
MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-${MONKEY_SMALL_GZ_DEFAULT}}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-${MONKEY_SMALL_FASTA_DEFAULT}}"
DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-${DEPLETION_REF_FASTA_DEFAULT}}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-${DEPLETED_FASTQ_DEFAULT}}"
PATHOGEN_FASTA="${PATHOGEN_FASTA:-${DEFAULT_PATHOGEN_FASTA}}"
TARGET_LABEL="${TARGET_LABEL:-${DEFAULT_TARGET_LABEL}}"
METAMAPS_DB_DIR="${METAMAPS_DB_DIR:-${METAMAPS_DB_DIR_DEFAULT}}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_single_metamaps_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-12}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-3}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"
MAX_MEMORY_GB="${MAX_MEMORY_GB:-120}"

require_exe python3
require_exe minimap2
require_exe samtools
require_exe read_analysis.py
require_exe simulator.py
require_exe metamaps
require_file "${SAMPLE_FASTQ_PY}"
require_file "${BUILD_MIXED_FASTQ_PY}"
require_file "${COMBINE_NANOSIM_FASTQ_PY}"
require_file "${SUMMARISE_METAMAPS_WIMP_PY}"
require_file "${REAL_FASTQ}"
require_file "${PATHOGEN_FASTA}"
require_dir "${METAMAPS_DB_DIR}"
require_file "${METAMAPS_DB_DIR}/DB.fa"
require_dir "${METAMAPS_DB_DIR}/taxonomy"
mkdir -p "${OUT_DIR}"

TRAIN_FASTQ="${OUT_DIR}/train_reads.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"
SIM_PREFIX="${OUT_DIR}/simulated_pathogen"
SIM_POOL_FASTQ="${OUT_DIR}/sim_pool.fastq"
SIM_POOL_FASTQ_GZ="${OUT_DIR}/sim_pool.fastq.gz"
SUMMARY_TSV="${OUT_DIR}/spikein_metamaps_summary.tsv"
RUN_META_TSV="${OUT_DIR}/run_metadata.tsv"
WORK_FASTQ="${DEPLETED_FASTQ}"

python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${REAL_FASTQ}" --n_reads "${TRAIN_READS_N}" --seed 1 --out_fastq_gz "${TRAIN_FASTQ}"

if [[ ! -s "${MONKEY_SMALL_FASTA}" ]]; then
    gzip -cd "${MONKEY_SMALL_GZ}" > "${MONKEY_SMALL_FASTA}"
fi

read_analysis.py genome --read "${TRAIN_FASTQ}" --ref_g "${MONKEY_SMALL_FASTA}" --output "${NS_MODEL_PREFIX}" --num_threads "${THREADS}" --fastq
simulator.py genome --ref_g "${PATHOGEN_FASTA}" --model_prefix "${NS_MODEL_PREFIX}" --output "${SIM_PREFIX}" --number "${SIM_POOL_N}" -dna_type linear --seed 42 --num_threads "${THREADS}" --fastq
python3 "${COMBINE_NANOSIM_FASTQ_PY}" --sim_prefix "${SIM_PREFIX}" --out_fastq "${SIM_POOL_FASTQ}"
gzip -c "${SIM_POOL_FASTQ}" > "${SIM_POOL_FASTQ_GZ}"

if [[ -s "${WORK_FASTQ}" ]]; then
    log_info "Reusing depleted FASTQ: ${WORK_FASTQ}"
elif [[ "${DO_HOST_DEPLETION}" == "true" ]]; then
    minimap2 -t "${THREADS}" -a -x map-ont "${DEPLETION_REF_FASTA}" "${REAL_FASTQ}" | samtools view -h -T "${DEPLETION_REF_FASTA}" -b -f 4 -@ "${THREADS}" - | samtools fastq -@ "${THREADS}" - | gzip -c > "${WORK_FASTQ}"
else
    cp "${REAL_FASTQ}" "${WORK_FASTQ}"
fi

{
    printf 'parameter\tvalue\n'
    printf 'repo_dir\t%s\n' "${REPO_DIR}"
    printf 'pathogen_fasta\t%s\n' "${PATHOGEN_FASTA}"
    printf 'target_label\t%s\n' "${TARGET_LABEL}"
    printf 'metamaps_db_dir\t%s\n' "${METAMAPS_DB_DIR}"
    printf 'threads\t%s\n' "${THREADS}"
    printf 'replicates\t%s\n' "${REPLICATES}"
    printf 'spike_levels\t%s\n' "${SPIKE_LEVELS}"
    printf 'max_memory_gb\t%s\n' "${MAX_MEMORY_GB}"
} > "${RUN_META_TSV}"

printf 'replicate\tspike_n\tmixed_fastq_gz\tmetamaps_prefix\tmetamaps_wimp\tmetamaps_total_reads_definedGenomes\tmetamaps_target_reads\n' > "${SUMMARY_TSV}"

unset MALLOC_ARENA_MAX || true

for rep in $(seq 1 "${REPLICATES}"); do
    for spike_n in ${SPIKE_LEVELS}; do
        MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"
        mkdir -p "${MIX_DIR}"
        SPIKE_FASTQ_GZ="${MIX_DIR}/spike.fastq.gz"
        MIX_FASTQ_GZ="${MIX_DIR}/mixed.fastq.gz"
        METAMAPS_PREFIX="${MIX_DIR}/metamaps"
        WIMP_TSV="${METAMAPS_PREFIX}.EM.WIMP"
        WIMP_SUMMARY_TSV="${MIX_DIR}/metamaps.wimp.summary.tsv"
        SPIKE_SEED=$(( rep * 100000 + spike_n ))

        python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${SIM_POOL_FASTQ_GZ}" --n_reads "${spike_n}" --seed "${SPIKE_SEED}" --out_fastq_gz "${SPIKE_FASTQ_GZ}"
        python3 "${BUILD_MIXED_FASTQ_PY}" --background_fastq_gz "${WORK_FASTQ}" --spike_fastq_gz "${SPIKE_FASTQ_GZ}" --out_fastq_gz "${MIX_FASTQ_GZ}"

        metamaps mapDirectly \
            -t "${THREADS}" \
            --all \
            --maxmemory "${MAX_MEMORY_GB}" \
            -r "${METAMAPS_DB_DIR}/DB.fa" \
            -q "${MIX_FASTQ_GZ}" \
            -o "${METAMAPS_PREFIX}"

        mm_total="NA"
        mm_target="NA"

        unset MALLOC_ARENA_MAX

        if metamaps classify \
            -t 1 \
            --mappings "${METAMAPS_PREFIX}" \
            --DB "${METAMAPS_DB_DIR}"; then

            if [[ -s "${WIMP_TSV}" ]]; then
                python3 "${SUMMARISE_METAMAPS_WIMP_PY}" \
                    --wimp_tsv "${WIMP_TSV}" \
                    --analysis_level definedGenomes \
                    --target_label "${TARGET_LABEL}" \
                    --out_tsv "${WIMP_SUMMARY_TSV}" || true

                if [[ -s "${WIMP_SUMMARY_TSV}" ]]; then
                    mm_total="$(awk 'NR==2{print $3}' "${WIMP_SUMMARY_TSV}")"
                    mm_target="$(awk 'NR==2{print $5}' "${WIMP_SUMMARY_TSV}")"
                else
                    log_warn "Empty WIMP summary for ${MIX_DIR}"
                fi
            else
                log_warn "Missing or empty WIMP output for ${MIX_DIR}"
            fi
        else
            log_warn "MetaMaps classify failed for ${MIX_DIR}"
        fi

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${rep}" "${spike_n}" "${MIX_FASTQ_GZ}" "${METAMAPS_PREFIX}" \
            "${WIMP_TSV}" "${mm_total}" "${mm_target}" \
            >> "${SUMMARY_TSV}"



    done
done

log_info "Done: ${SUMMARY_TSV}"
