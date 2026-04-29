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
export KMERSUTRA_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/kmersutra"
export PYTHONPATH="${KMERSUTRA_DIR}:${PYTHONPATH:-}"


log_info() { printf '[INFO] %s\n' "$*"; }
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
KMERSUTRA_MAX_MISMATCHES="${KMERSUTRA_MAX_MISMATCHES:-0}"
KMERSUTRA_FUZZY_MIN_K="${KMERSUTRA_FUZZY_MIN_K:-101}"

for executable in python3 read_analysis.py simulator.py minimap2 samtools; do
    require_exe "${executable}"
done

for required_file in \
    "${SAMPLE_FASTQ_PY}" \
    "${BUILD_MIXED_FASTQ_PY}" \
    "${COMBINE_NANOSIM_FASTQ_PY}" \
    "${REAL_FASTQ}" \
    "${PATHOGEN_CONFIG_TSV}" \
    "${KMERSUTRA_PANEL}" \
    "${KMERSUTRA_DIR}/scripts/screen_reads_for_clade_kmers.py"
do
    require_file "${required_file}"
done

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
    printf 'kmersutra_max_mismatches\t%s\n' "${KMERSUTRA_MAX_MISMATCHES}"
    printf 'kmersutra_fuzzy_min_k\t%s\n' "${KMERSUTRA_FUZZY_MIN_K}"
    printf 'threads\t%s\n' "${THREADS}"
    printf 'replicates\t%s\n' "${REPLICATES}"
    printf 'spike_levels\t%s\n' "${SPIKE_LEVELS}"
} > "${RUN_META_TSV}"

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

header='replicate\tspike_n\ttotal_spiked_reads\tmixed_fastq_gz\tkmersutra_out_dir\tkmersutra_calls_tsv\tkmersutra_evidence_tsv\tkmersutra_read_hits_tsv_gz'
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

        python3 "${KMERSUTRA_DIR}/scripts/screen_reads_for_clade_kmers.py" \
            --input "${mix_fastq_gz}" \
            --input_format fastq \
            --panel "${KMERSUTRA_PANEL}" \
            --sample_id "rep${rep}_n${spike_n}" \
            --out_dir "${ks_out_dir}" \
            --max_mismatches "${KMERSUTRA_MAX_MISMATCHES}" \
            --fuzzy_min_k "${KMERSUTRA_FUZZY_MIN_K}" \
            --verbose

        calls_tsv="${ks_out_dir}/species_detection_calls.tsv"
        evidence_tsv="${ks_out_dir}/sample_species_kmer_evidence.tsv"
        read_hits_tsv_gz="${ks_out_dir}/read_level_species_kmer_hits.tsv.gz"

        require_file "${calls_tsv}"
        require_file "${evidence_tsv}"
        require_file "${read_hits_tsv_gz}"

        row="${rep}\t${spike_n}\t${total_spiked}\t${mix_fastq_gz}\t${ks_out_dir}\t${calls_tsv}\t${evidence_tsv}\t${read_hits_tsv_gz}"

        for line in "${PANEL_LINES[@]}"; do
            target_label="$(printf '%s' "${line}" | cut -f2)"

            call_value="$(
                awk -F '\t' -v species="${target_label}" '
                    NR == 1 {
                        for (i = 1; i <= NF; i++) {
                            if ($i == "species_name") species_col = i
                            if ($i == "call") call_col = i
                        }
                        next
                    }
                    $species_col == species {print $call_col}
                ' "${calls_tsv}" | head -n 1
            )"

            unique_kmers="$(
                awk -F '\t' -v species="${target_label}" '
                    NR == 1 {
                        for (i = 1; i <= NF; i++) {
                            if ($i == "species_name") species_col = i
                            if ($i == "n_unique_diagnostic_kmers") value_col = i
                        }
                        next
                    }
                    $species_col == species {print $value_col}
                ' "${evidence_tsv}" | head -n 1
            )"

            positive_reads="$(
                awk -F '\t' -v species="${target_label}" '
                    NR == 1 {
                        for (i = 1; i <= NF; i++) {
                            if ($i == "species_name") species_col = i
                            if ($i == "n_positive_reads") value_col = i
                        }
                        next
                    }
                    $species_col == species {print $value_col}
                ' "${evidence_tsv}" | head -n 1
            )"

            confidence="$(
                awk -F '\t' -v species="${target_label}" '
                    NR == 1 {
                        for (i = 1; i <= NF; i++) {
                            if ($i == "species_name") species_col = i
                            if ($i == "confidence_score") value_col = i
                        }
                        next
                    }
                    $species_col == species {print $value_col}
                ' "${evidence_tsv}" | head -n 1
            )"

            row+="\t${call_value:-NA}\t${unique_kmers:-0}\t${positive_reads:-0}\t${confidence:-0}"
        done

        printf '%b\n' "${row}" >> "${SUMMARY_TSV}"
    done
done

log_info "Done: ${SUMMARY_TSV}"
