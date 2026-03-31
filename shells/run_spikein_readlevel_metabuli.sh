#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 240G
#$ -adds l_hard h_vmem 240G
#$ -N MBspike_read

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
CONFIG_DIR="${CONFIG_DIR:-${REPO_DIR}/configs}"

PATHOGEN_CONFIG_TSV="${PATHOGEN_CONFIG_TSV:-${DEFAULT_PATHOGEN_PANEL_3_WITH_TAXID}}"
SAMPLE_FASTQ_PY="${SAMPLE_FASTQ_PY:-${PY_SCRIPTS_DIR}/sample_fastq.py}"
BUILD_MIXED_FASTQ_PY="${BUILD_MIXED_FASTQ_PY:-${PY_SCRIPTS_DIR}/build_mixed_fastq.py}"
COMBINE_NANOSIM_FASTQ_PY="${COMBINE_NANOSIM_FASTQ_PY:-${PY_SCRIPTS_DIR}/combine_nanosim_fastq.py}"
SUMMARISE_METABULI_REPORT_PY="${SUMMARISE_METABULI_REPORT_PY:-${PY_SCRIPTS_DIR}/summarise_metabuli_report_with_reported_taxa.py}"

REAL_FASTQ="${REAL_FASTQ:-${REAL_FASTQ_DEFAULT}}"
MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-${MONKEY_SMALL_GZ_DEFAULT}}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-${MONKEY_SMALL_FASTA_DEFAULT}}"
DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-${DEPLETION_REF_FASTA_DEFAULT}}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-${DEPLETED_FASTQ_DEFAULT}}"

METABULI_DB_DIR="${METABULI_DB_DIR:-}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_metabuli_$(date +%Y%m%d_%H%M%S)}"

THREADS="${THREADS:-12}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-3}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"

MAX_RAM_GB="${MAX_RAM_GB:-120}"
# best for nanopore
MIN_SCORE="${MIN_SCORE:-0.008}"

for x in python3 minimap2 samtools read_analysis.py simulator.py metabuli; do
    require_exe "$x"
done

for f in \
    "${SAMPLE_FASTQ_PY}" \
    "${BUILD_MIXED_FASTQ_PY}" \
    "${COMBINE_NANOSIM_FASTQ_PY}" \
    "${SUMMARISE_METABULI_REPORT_PY}" \
    "${REAL_FASTQ}" \
    "${PATHOGEN_CONFIG_TSV}"
do
    require_file "$f"
done

require_dir "${METABULI_DB_DIR}"
mkdir -p "${OUT_DIR}"

mapfile -t PANEL_LINES < <(awk 'BEGIN{FS="\t"} NR>1 && NF>=2 {print $1 "\t" $2 "\t" $3}' "${PATHOGEN_CONFIG_TSV}")
[[ ${#PANEL_LINES[@]} -gt 0 ]] || { log_error "No pathogen entries found in ${PATHOGEN_CONFIG_TSV}"; exit 1; }

TRAIN_FASTQ="${OUT_DIR}/train_reads.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"
SUMMARY_TSV="${OUT_DIR}/spikein_metabuli_summary.tsv"
WORK_FASTQ="${DEPLETED_FASTQ}"

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
    target_label="$(printf '%s' "${line}" | cut -f2)"
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
    :
elif [[ "${DO_HOST_DEPLETION}" == "true" ]]; then
    minimap2 -t "${THREADS}" -a -x map-ont "${DEPLETION_REF_FASTA}" "${REAL_FASTQ}" \
        | samtools view -h -T "${DEPLETION_REF_FASTA}" -b -f 4 -@ "${THREADS}" - \
        | samtools fastq -@ "${THREADS}" - \
        | gzip -c > "${WORK_FASTQ}"
else
    cp "${REAL_FASTQ}" "${WORK_FASTQ}"
fi

header='replicate\tspike_n\ttotal_spiked_reads\tn_genomes\tmixed_fastq_gz\tmetabuli_report_tsv\tmetabuli_classifications_tsv\tmetabuli_total_reads_in_report'
for line in "${PANEL_LINES[@]}"; do
    target_label="$(printf '%s' "${line}" | cut -f2)"
    safe_label="$(printf '%s' "${target_label}" | tr ' ' '_' | tr '/' '_' | tr -cd '[:alnum:]_')"
    header+="\tmetabuli_${safe_label}_clade_reads"
    header+="\tmetabuli_${safe_label}_direct_reads"
    header+="\tmetabuli_${safe_label}_found"
done
printf '%b\n' "${header}" > "${SUMMARY_TSV}"

for ((rep=1; rep<=REPLICATES; rep++)); do
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

        [[ ${#tmp_spikes[@]} -gt 0 ]] || { log_warn "No spike FASTQs for ${MIX_DIR}"; continue; }

        mix_fastq_gz="${MIX_DIR}/mixed.fastq.gz"
        python3 "${BUILD_MIXED_FASTQ_PY}" \
            --background_fastq_gz "${WORK_FASTQ}" \
            --spike_fastq_gz "${tmp_spikes[@]}" \
            --out_fastq_gz "${mix_fastq_gz}"

        job_id="metabuli"
        report_tsv="${MIX_DIR}/${job_id}_report.tsv"
        classifications_tsv="${MIX_DIR}/${job_id}_classifications.tsv"
        helper_tsv="${MIX_DIR}/${job_id}.summary.tsv"

        summary_values=""
        n_genomes="${#PANEL_LINES[@]}"

        if metabuli classify \
            --seq-mode 3 \
            --threads "${THREADS}" \
            --max-ram "${MAX_RAM_GB}" \
            --min-score "${MIN_SCORE}" \
            "${mix_fastq_gz}" \
            "${METABULI_DB_DIR}" \
            "${MIX_DIR}" \
            "${job_id}"; then

            if [[ -s "${report_tsv}" ]]; then
                helper_args=(
                    python3 "${SUMMARISE_METABULI_REPORT_PY}"
                    --report_tsv "${report_tsv}"
                    --out_tsv "${helper_tsv}"
                )

                for line in "${PANEL_LINES[@]}"; do
                    target_label="$(printf '%s' "${line}" | cut -f2)"
                    target_taxid="$(printf '%s' "${line}" | cut -f3)"
                    helper_args+=(--target-label "${target_label}")
                    if [[ -n "${target_taxid}" ]]; then
                        helper_args+=(--target-taxid "${target_taxid}")
                    fi
                done

                "${helper_args[@]}"

                summary_values="$(awk 'NR==2{print}' "${helper_tsv}")"
            else
                log_warn "Missing report file for ${MIX_DIR}"
            fi
        else
            log_warn "Metabuli classify failed for ${MIX_DIR}"
        fi

        if [[ -z "${summary_values}" ]]; then
            summary_values="0"
            for line in "${PANEL_LINES[@]}"; do
                summary_values+="\t0\t0\t0"
            done
        fi

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%b\n' \
            "${rep}" \
            "${spike_n}" \
            "${total_spiked}" \
            "${n_genomes}" \
            "${mix_fastq_gz}" \
            "${report_tsv}" \
            "${classifications_tsv}" \
            "${summary_values}" \
            >> "${SUMMARY_TSV}"
    done
done

log_info "Done: ${SUMMARY_TSV}"
