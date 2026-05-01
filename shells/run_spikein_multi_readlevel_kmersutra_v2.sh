#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 240G
#$ -adds l_hard h_vmem 240G
#$ -N KSspike_multi

set -euo pipefail

log_info() { printf '[INFO] %s\n' "$*"; }
log_warn() { printf '[WARN] %s\n' "$*" >&2; }
log_error() { printf '[ERROR] %s\n' "$*" >&2; }

require_file() {
    [[ -f "$1" ]] || {
        log_error "Required file not found: $1"
        exit 1
    }
}

require_dir() {
    [[ -d "$1" ]] || {
        log_error "Required directory not found: $1"
        exit 1
    }
}

require_exe() {
    command -v "$1" >/dev/null 2>&1 || {
        log_error "Required executable not found on PATH: $1"
        exit 1
    }
}

get_tsv_value() {
    local tsv_path="$1"
    local key_column="$2"
    local key_value="$3"
    local value_column="$4"

    awk -F '\t' \
        -v key_column_name="${key_column}" \
        -v key_value="${key_value}" \
        -v value_column_name="${value_column}" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == key_column_name) key_index = i
                if ($i == value_column_name) value_index = i
            }
            next
        }
        key_index && value_index && $key_index == key_value {
            print $value_index
            exit
        }
    ' "${tsv_path}"
}

run_kmersutra_screen() {
    local input_fastq="$1"
    local sample_id="$2"
    local out_dir="$3"
    local stdout_log="$4"
    local stderr_log="$5"
    local command_tsv="$6"

    local -a ks_command=()

    if command -v "${KMERSUTRA_SCREEN_CMD}" >/dev/null 2>&1; then
        ks_command=("${KMERSUTRA_SCREEN_CMD}")
    elif [[ -n "${KMERSUTRA_DIR:-}" && -f "${KMERSUTRA_DIR}/scripts/screen_reads_for_clade_kmers.py" ]]; then
        log_warn "${KMERSUTRA_SCREEN_CMD} was not found on PATH. Falling back to the repository script."
        export PYTHONPATH="${KMERSUTRA_DIR}:${PYTHONPATH:-}"
        ks_command=(python3 "${KMERSUTRA_DIR}/scripts/screen_reads_for_clade_kmers.py")
    else
        log_error "Could not find ${KMERSUTRA_SCREEN_CMD} or a fallback KmerSutra screening script."
        log_error "Install KmerSutra with pip install -e . or set KMERSUTRA_DIR correctly."
        exit 1
    fi

    ks_command+=(
        --input "${input_fastq}"
        --input_format fastq
        --panel "${KMERSUTRA_PANEL}"
        --sample_id "${sample_id}"
        --out_dir "${out_dir}"
        --threads "${KMERSUTRA_THREADS}"
        --chunk_size "${KMERSUTRA_CHUNK_SIZE}"
        --max_mismatches "${KMERSUTRA_MAX_MISMATCHES}"
        --fuzzy_min_k "${KMERSUTRA_FUZZY_MIN_K}"
    )

    if [[ "${KMERSUTRA_VERBOSE}" == "true" ]]; then
        ks_command+=(--verbose)
    fi

    {
        printf 'field\tvalue\n'
        printf 'sample_id\t%s\n' "${sample_id}"
        printf 'input_fastq\t%s\n' "${input_fastq}"
        printf 'out_dir\t%s\n' "${out_dir}"
        printf 'command\t'
        printf '%q ' "${ks_command[@]}"
        printf '\n'
    } > "${command_tsv}"

    log_info "Running KmerSutra for ${sample_id}"
    log_info "KmerSutra stderr log: ${stderr_log}"

    if ! "${ks_command[@]}" > "${stdout_log}" 2> "${stderr_log}"; then
        log_error "KmerSutra screening failed for ${sample_id}"
        log_error "Command details: ${command_tsv}"
        log_error "Stdout log: ${stdout_log}"
        log_error "Stderr log: ${stderr_log}"
        exit 1
    fi
}

resolve_repo_dir() {
    local candidate
    local -a candidates=()

    if [[ -n "${REPO_DIR:-}" ]]; then
        candidates+=("${REPO_DIR}")
    fi

    if [[ -n "${SGE_O_WORKDIR:-}" ]]; then
        candidates+=("${SGE_O_WORKDIR}/PT_nanopore_spike_in_pathogen_detection")
        candidates+=("${SGE_O_WORKDIR}")
    fi

    candidates+=("${PWD}/PT_nanopore_spike_in_pathogen_detection")
    candidates+=("${PWD}")
    candidates+=("/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection")
    candidates+=("/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection")

    for candidate in "${candidates[@]}"; do
        if [[ -f "${candidate}/configs/pipeline_paths.sh" ]]; then
            printf '%s\n' "${candidate}"
            return 0
        fi
    done

    log_error "Could not resolve REPO_DIR. Looked for configs/pipeline_paths.sh in:"
    for candidate in "${candidates[@]}"; do
        log_error "  ${candidate}"
    done
    return 1
}

REPO_DIR="$(resolve_repo_dir)"
log_info "Resolved repository directory: ${REPO_DIR}"

source "${REPO_DIR}/configs/pipeline_paths.sh"

PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR:-${REPO_DIR}/scripts}"
CONFIG_DIR="${CONFIG_DIR:-${REPO_DIR}/configs}"
PATHOGEN_CONFIG_TSV="${PATHOGEN_CONFIG_TSV:-${DEFAULT_PATHOGEN_PANEL_3}}"

SAMPLE_FASTQ_PY="${PY_SCRIPTS_DIR}/sample_fastq.py"
BUILD_MIXED_FASTQ_PY="${PY_SCRIPTS_DIR}/build_mixed_fastq.py"
COMBINE_NANOSIM_FASTQ_PY="${PY_SCRIPTS_DIR}/combine_nanosim_fastq.py"

REAL_FASTQ="${REAL_FASTQ:-${REAL_FASTQ_DEFAULT}}"
MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-${MONKEY_SMALL_GZ_DEFAULT}}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-${MONKEY_SMALL_FASTA_DEFAULT}}"
DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-${DEPLETION_REF_FASTA_DEFAULT}}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-${DEPLETED_FASTQ_DEFAULT}}"

OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_kmersutra_panel3_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-12}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-12}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"

KMERSUTRA_DIR="${KMERSUTRA_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/kmersutra}"
KMERSUTRA_PANEL="${KMERSUTRA_PANEL:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/kmersutra_data/kmersutra_plasmodium_panel3_build/species_kmer_panel.tsv.gz}"
KMERSUTRA_SCREEN_CMD="${KMERSUTRA_SCREEN_CMD:-kmersutra-screen}"
KMERSUTRA_THREADS="${KMERSUTRA_THREADS:-${THREADS}}"
KMERSUTRA_CHUNK_SIZE="${KMERSUTRA_CHUNK_SIZE:-1000}"
KMERSUTRA_MAX_MISMATCHES="${KMERSUTRA_MAX_MISMATCHES:-0}"
KMERSUTRA_FUZZY_MIN_K="${KMERSUTRA_FUZZY_MIN_K:-101}"
KMERSUTRA_VERBOSE="${KMERSUTRA_VERBOSE:-true}"

for executable in python3 read_analysis.py simulator.py minimap2 samtools; do
    require_exe "${executable}"
done

for required_file in \
    "${SAMPLE_FASTQ_PY}" \
    "${BUILD_MIXED_FASTQ_PY}" \
    "${COMBINE_NANOSIM_FASTQ_PY}" \
    "${REAL_FASTQ}" \
    "${PATHOGEN_CONFIG_TSV}" \
    "${KMERSUTRA_PANEL}"
do
    require_file "${required_file}"
done

if ! command -v "${KMERSUTRA_SCREEN_CMD}" >/dev/null 2>&1; then
    if [[ -n "${KMERSUTRA_DIR:-}" && -f "${KMERSUTRA_DIR}/scripts/screen_reads_for_clade_kmers.py" ]]; then
        log_warn "${KMERSUTRA_SCREEN_CMD} is not on PATH. The script will use the fallback repository script."
        export PYTHONPATH="${KMERSUTRA_DIR}:${PYTHONPATH:-}"
    else
        log_error "${KMERSUTRA_SCREEN_CMD} is not on PATH and no fallback script was found."
        log_error "Expected fallback: ${KMERSUTRA_DIR}/scripts/screen_reads_for_clade_kmers.py"
        exit 1
    fi
fi

########################################################################
# Job-local working directory
########################################################################

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
                --include='*.tsv.gz' \
                --include='*.txt' \
                --include='*.log' \
                --include='*.json' \
                --include='*.html' \
                --include='*.xlsx' \
                --include='*.pdf' \
                --include='*.svg' \
                --include='*.png' \
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
                --exclude='simulated_pathogen*/***' \
                --exclude='sim_pool*' \
                --exclude='nanosim_training*' \
                --exclude='*' \
                "${WORK_OUT_DIR}/" "${FINAL_OUT_DIR}/" || true
        fi

        python3 - "${FINAL_OUT_DIR}" "${WORK_OUT_DIR}" "${FINAL_OUT_DIR}" <<'PYTHON_REWRITE_PATHS'
from pathlib import Path
import sys

final_dir = Path(sys.argv[1])
old = sys.argv[2]
new = sys.argv[3]

for path in final_dir.rglob("*"):
    if not path.is_file():
        continue
    if path.suffix not in {".tsv", ".txt", ".log", ".json", ".html"}:
        continue
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        continue
    updated = text.replace(old, new)
    if updated != text:
        path.write_text(updated, encoding="utf-8")
PYTHON_REWRITE_PATHS
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

mapfile -t PANEL_LINES < <(
    awk 'BEGIN{FS="\t"} NR>1 && NF>=2 {print $1 "\t" $2}' "${PATHOGEN_CONFIG_TSV}"
)

[[ ${#PANEL_LINES[@]} -gt 0 ]] || {
    log_error "No pathogen entries found in ${PATHOGEN_CONFIG_TSV}"
    exit 1
}

TRAIN_FASTQ="${OUT_DIR}/train_reads.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"
SUMMARY_TSV="${OUT_DIR}/spikein_multi_kmersutra_summary.tsv"
RUN_META_TSV="${OUT_DIR}/run_metadata.tsv"
WORK_FASTQ="${DEPLETED_FASTQ}"

{
    printf 'parameter\tvalue\n'
    printf 'repo_dir\t%s\n' "${REPO_DIR}"
    printf 'pathogen_config_tsv\t%s\n' "${PATHOGEN_CONFIG_TSV}"
    printf 'kmersutra_dir\t%s\n' "${KMERSUTRA_DIR}"
    printf 'kmersutra_panel\t%s\n' "${KMERSUTRA_PANEL}"
    printf 'kmersutra_screen_cmd\t%s\n' "${KMERSUTRA_SCREEN_CMD}"
    printf 'kmersutra_threads\t%s\n' "${KMERSUTRA_THREADS}"
    printf 'kmersutra_chunk_size\t%s\n' "${KMERSUTRA_CHUNK_SIZE}"
    printf 'kmersutra_max_mismatches\t%s\n' "${KMERSUTRA_MAX_MISMATCHES}"
    printf 'kmersutra_fuzzy_min_k\t%s\n' "${KMERSUTRA_FUZZY_MIN_K}"
    printf 'threads\t%s\n' "${THREADS}"
    printf 'replicates\t%s\n' "${REPLICATES}"
    printf 'spike_levels\t%s\n' "${SPIKE_LEVELS}"
} > "${RUN_META_TSV}"

log_info "KmerSutra screen command: ${KMERSUTRA_SCREEN_CMD}"
log_info "KmerSutra threads: ${KMERSUTRA_THREADS}"
log_info "KmerSutra chunk size: ${KMERSUTRA_CHUNK_SIZE}"
log_info "KmerSutra max mismatches: ${KMERSUTRA_MAX_MISMATCHES}"
log_info "KmerSutra fuzzy minimum k: ${KMERSUTRA_FUZZY_MIN_K}"

python3 "${SAMPLE_FASTQ_PY}" \
    --fastq_gz "${REAL_FASTQ}" \
    --n_reads "${TRAIN_READS_N}" \
    --seed 1 \
    --out_fastq_gz "${TRAIN_FASTQ}"

if [[ ! -s "${MONKEY_SMALL_FASTA}" ]]; then
    gzip -cd "${MONKEY_SMALL_GZ}" > "${MONKEY_SMALL_FASTA}"
fi

read_analysis.py genome \
    --read "${TRAIN_FASTQ}" \
    --ref_g "${MONKEY_SMALL_FASTA}" \
    --output "${NS_MODEL_PREFIX}" \
    --num_threads "${THREADS}" \
    --fastq

GENOME_INDEX=0
for line in "${PANEL_LINES[@]}"; do
    GENOME_INDEX=$((GENOME_INDEX + 1))
    pathogen_fasta="$(printf '%s' "${line}" | cut -f1)"
    require_file "${pathogen_fasta}"

    sim_prefix="${OUT_DIR}/simulated_pathogen_${GENOME_INDEX}"
    sim_pool_fastq="${OUT_DIR}/sim_pool_${GENOME_INDEX}.fastq"
    sim_pool_fastq_gz="${OUT_DIR}/sim_pool_${GENOME_INDEX}.fastq.gz"

    simulator.py genome \
        --ref_g "${pathogen_fasta}" \
        --model_prefix "${NS_MODEL_PREFIX}" \
        --output "${sim_prefix}" \
        --number "${SIM_POOL_N}" \
        -dna_type linear \
        --seed "$((40 + GENOME_INDEX))" \
        --num_threads "${THREADS}" \
        --fastq

    python3 "${COMBINE_NANOSIM_FASTQ_PY}" \
        --sim_prefix "${sim_prefix}" \
        --out_fastq "${sim_pool_fastq}"

    gzip -c "${sim_pool_fastq}" > "${sim_pool_fastq_gz}"
done

if [[ -s "${WORK_FASTQ}" ]]; then
    log_info "Reusing depleted FASTQ: ${WORK_FASTQ}"
elif [[ "${DO_HOST_DEPLETION}" == "true" ]]; then
    minimap2 -t "${THREADS}" -a -x map-ont "${DEPLETION_REF_FASTA}" "${REAL_FASTQ}" \
        | samtools view -h -T "${DEPLETION_REF_FASTA}" -b -f 4 -@ "${THREADS}" - \
        | samtools fastq -@ "${THREADS}" - \
        | gzip -c > "${WORK_FASTQ}"
else
    cp "${REAL_FASTQ}" "${WORK_FASTQ}"
fi

header='replicate\tspike_n\ttotal_spiked_reads\tmixed_fastq_gz\tkmersutra_out_dir\tkmersutra_calls_tsv\tkmersutra_evidence_tsv\tkmersutra_read_hits_tsv_gz\tkmersutra_ml_features_tsv\tkmersutra_stdout_log\tkmersutra_stderr_log\tkmersutra_command_tsv\tkmersutra_runtime_seconds'
for line in "${PANEL_LINES[@]}"; do
    target_label="$(printf '%s' "${line}" | cut -f2)"
    safe_label="$(printf '%s' "${target_label}" | tr ' ' '_' | tr '/' '_' | tr -cd '[:alnum:]_')"
    header+="\tkmersutra_${safe_label}_call"
    header+="\tkmersutra_${safe_label}_unique_kmers"
    header+="\tkmersutra_${safe_label}_positive_reads"
    header+="\tkmersutra_${safe_label}_confidence"
done
printf '%b\n' "${header}" > "${SUMMARY_TSV}"

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

        mix_fastq_gz="${MIX_DIR}/mixed.fastq.gz"

        python3 "${BUILD_MIXED_FASTQ_PY}" \
            --background_fastq_gz "${WORK_FASTQ}" \
            --spike_fastq_gz "${tmp_spikes[@]}" \
            --out_fastq_gz "${mix_fastq_gz}"

        ks_out_dir="${MIX_DIR}/kmersutra"
        mkdir -p "${ks_out_dir}"

        ks_stdout_log="${ks_out_dir}/kmersutra_screen.stdout.log"
        ks_stderr_log="${ks_out_dir}/kmersutra_screen.stderr.log"
        ks_command_tsv="${ks_out_dir}/kmersutra_screen.command.tsv"

        ks_start_seconds="${SECONDS}"
        run_kmersutra_screen \
            "${mix_fastq_gz}" \
            "rep${rep}_n${spike_n}" \
            "${ks_out_dir}" \
            "${ks_stdout_log}" \
            "${ks_stderr_log}" \
            "${ks_command_tsv}"
        ks_runtime_seconds=$((SECONDS - ks_start_seconds))

        calls_tsv="${ks_out_dir}/species_detection_calls.tsv"
        evidence_tsv="${ks_out_dir}/sample_species_kmer_evidence.tsv"
        read_hits_tsv_gz="${ks_out_dir}/read_level_species_kmer_hits.tsv.gz"
        ml_features_tsv="${ks_out_dir}/sequence_ml_features.tsv"

        require_file "${calls_tsv}"
        require_file "${evidence_tsv}"
        require_file "${read_hits_tsv_gz}"

        if [[ ! -s "${ml_features_tsv}" ]]; then
            ml_features_tsv="NA"
        fi

        row="${rep}\t${spike_n}\t${total_spiked}\t${mix_fastq_gz}\t${ks_out_dir}\t${calls_tsv}\t${evidence_tsv}\t${read_hits_tsv_gz}\t${ml_features_tsv}\t${ks_stdout_log}\t${ks_stderr_log}\t${ks_command_tsv}\t${ks_runtime_seconds}"

        for line in "${PANEL_LINES[@]}"; do
            target_label="$(printf '%s' "${line}" | cut -f2)"

            call_value="$(get_tsv_value "${calls_tsv}" species_name "${target_label}" call)"
            unique_kmers="$(get_tsv_value "${evidence_tsv}" species_name "${target_label}" n_unique_diagnostic_kmers)"
            positive_reads="$(get_tsv_value "${evidence_tsv}" species_name "${target_label}" n_positive_reads)"
            confidence="$(get_tsv_value "${evidence_tsv}" species_name "${target_label}" confidence_score)"

            row+="\t${call_value:-NA}\t${unique_kmers:-0}\t${positive_reads:-0}\t${confidence:-0}"
        done

        printf '%b\n' "${row}" >> "${SUMMARY_TSV}"
        log_info "Finished rep=${rep} spike_n=${spike_n} KmerSutra runtime=${ks_runtime_seconds}s"
    done
done

log_info "Done: ${SUMMARY_TSV}"
