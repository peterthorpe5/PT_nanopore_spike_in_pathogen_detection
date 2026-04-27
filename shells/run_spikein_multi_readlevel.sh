#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N spike_multi_read

set -euo pipefail

log_info() { printf '[INFO] %s\n' "$*"; }
log_warn() { printf '[WARN] %s\n' "$*" >&2; }
log_error() { printf '[ERROR] %s\n' "$*" >&2; }
require_file() { [[ -f "$1" ]] || { log_error "Required file not found: $1"; exit 1; }; }
require_dir() { [[ -d "$1" ]] || { log_error "Required directory not found: $1"; exit 1; }; }
require_exe() { command -v "$1" >/dev/null 2>&1 || { log_error "Required executable not found on PATH: $1"; exit 1; }; }

resolve_repo_dir() {
    if [[ -n "${REPO_DIR:-}" ]]; then
        printf '%s\n' "${REPO_DIR}"
    elif [[ -n "${SGE_O_WORKDIR:-}" && -f "${SGE_O_WORKDIR}/configs/pipeline_paths.sh" ]]; then
        printf '%s\n' "${SGE_O_WORKDIR}"
    else
        printf '%s\n' "/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection"
    fi
}

run_minimap_to_bed() {
    local ref_fasta="$1"
    local query_fastq="$2"
    local bam_out="$3"
    local bed_out="$4"
    local split_prefix="$5"

    rm -f "${bam_out}" "${bed_out}"
    rm -f "${split_prefix}"*

    minimap2 \
        --split-prefix "${split_prefix}" \
        -t "${THREADS}" \
        -ax map-ont \
        "${ref_fasta}" \
        "${query_fastq}" \
        | samtools view -h -q "${MIN_MAPQ}" -e "rlen >= ${MIN_ALIGN_LEN}" - \
        | samtools sort -O bam -@ "${THREADS}" - \
        | tee "${bam_out}" \
        | bamToBed > "${bed_out}"

    samtools quickcheck -v "${bam_out}"
}

REPO_DIR="$(resolve_repo_dir)"
source "${REPO_DIR}/configs/pipeline_paths.sh"

PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR:-${REPO_DIR}/scripts}"
CONFIG_DIR="${CONFIG_DIR:-${REPO_DIR}/configs}"
PATHOGEN_CONFIG_TSV="${PATHOGEN_CONFIG_TSV:-${DEFAULT_PATHOGEN_PANEL_3}}"
SAMPLE_FASTQ_PY="${PY_SCRIPTS_DIR}/sample_fastq.py"
BUILD_MIXED_FASTQ_PY="${PY_SCRIPTS_DIR}/build_mixed_fastq.py"
COMBINE_NANOSIM_FASTQ_PY="${PY_SCRIPTS_DIR}/combine_nanosim_fastq.py"
SUMMARISE_KRAKEN_PY="${PY_SCRIPTS_DIR}/summarise_kraken_report_with_reported_taxa.py"
COUNT_BED_HITS_PY="${PY_SCRIPTS_DIR}/count_bed_hits.py"

REAL_FASTQ="${REAL_FASTQ:-${REAL_FASTQ_DEFAULT}}"
MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-${MONKEY_SMALL_GZ_DEFAULT}}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-${MONKEY_SMALL_FASTA_DEFAULT}}"
DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-${DEPLETION_REF_FASTA_DEFAULT}}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-${DEPLETED_FASTQ_DEFAULT}}"
MINIMAP_DB_FASTA="${MINIMAP_DB_FASTA:-${MINIMAP_DB_FASTA_DEFAULT}}"
KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-${KRAKEN_DB_DIR_DEFAULT}}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_read_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-12}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-${SPIKE_LEVELS_DEFAULT:-0 1 5 10 25 50 100 250 500 1000 2500 3000 4000 5000 7000 10000}}"
REPLICATES="${REPLICATES:-12}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"
MIN_MAPQ="${MIN_MAPQ:-15}"
MIN_ALIGN_LEN="${MIN_ALIGN_LEN:-500}"

for x in python3 minimap2 samtools kraken2 bamToBed read_analysis.py simulator.py; do
    require_exe "${x}"
done

for f in \
    "${SAMPLE_FASTQ_PY}" \
    "${BUILD_MIXED_FASTQ_PY}" \
    "${COMBINE_NANOSIM_FASTQ_PY}" \
    "${SUMMARISE_KRAKEN_PY}" \
    "${COUNT_BED_HITS_PY}" \
    "${REAL_FASTQ}" \
    "${PATHOGEN_CONFIG_TSV}" \
    "${MINIMAP_DB_FASTA}"; do
    require_file "${f}"
done
require_dir "${KRAKEN_DB_DIR}"
mkdir -p "${OUT_DIR}"

mapfile -t PANEL_LINES < <(awk 'BEGIN{FS="\t"} NR>1 && NF>=2 {print $1"\t"$2}' "${PATHOGEN_CONFIG_TSV}")
[[ ${#PANEL_LINES[@]} -gt 0 ]] || { log_error "No pathogen entries found in ${PATHOGEN_CONFIG_TSV}"; exit 1; }

TRAIN_FASTQ="${OUT_DIR}/train_reads.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"
SUMMARY_TSV="${OUT_DIR}/spikein_multi_summary.tsv"
RUN_META_TSV="${OUT_DIR}/run_metadata.tsv"
WORK_FASTQ="${DEPLETED_FASTQ}"

log_info "OUT_DIR=${OUT_DIR}"
log_info "PATHOGEN_CONFIG_TSV=${PATHOGEN_CONFIG_TSV}"
log_info "MINIMAP_DB_FASTA=${MINIMAP_DB_FASTA}"
log_info "SPIKE_LEVELS=${SPIKE_LEVELS}"

python3 "${SAMPLE_FASTQ_PY}" \
    --fastq_gz "${REAL_FASTQ}" \
    --n_reads "${TRAIN_READS_N}" \
    --seed 1 \
    --out_fastq_gz "${TRAIN_FASTQ}"

if [[ ! -s "${MONKEY_SMALL_FASTA}" ]]; then
    require_file "${MONKEY_SMALL_GZ}"
    gzip -cd "${MONKEY_SMALL_GZ}" > "${MONKEY_SMALL_FASTA}"
fi
require_file "${MONKEY_SMALL_FASTA}"

rm -f "${NS_MODEL_PREFIX}"_primary.bam "${NS_MODEL_PREFIX}"_primary.bam.bai
read_analysis.py genome \
    --read "${TRAIN_FASTQ}" \
    --ref_g "${MONKEY_SMALL_FASTA}" \
    --output "${NS_MODEL_PREFIX}" \
    --num_threads "${THREADS}" \
    --fastq

if [[ ! -s "${NS_MODEL_PREFIX}_primary.bam" ]]; then
    log_error "NanoSim training failed to create ${NS_MODEL_PREFIX}_primary.bam"
    log_error "Check the read_analysis.py output above. Do not reuse this OUT_DIR for other jobs."
    exit 1
fi
samtools quickcheck -v "${NS_MODEL_PREFIX}_primary.bam"

GENOME_INDEX=0
for line in "${PANEL_LINES[@]}"; do
    GENOME_INDEX=$((GENOME_INDEX + 1))
    pathogen_fasta="$(printf '%s' "${line}" | cut -f1)"
    target_label="$(printf '%s' "${line}" | cut -f2)"
    require_file "${pathogen_fasta}"
    sim_prefix="${OUT_DIR}/simulated_pathogen_${GENOME_INDEX}"
    sim_pool_fastq="${OUT_DIR}/sim_pool_${GENOME_INDEX}.fastq"
    sim_pool_fastq_gz="${OUT_DIR}/sim_pool_${GENOME_INDEX}.fastq.gz"
    log_info "Simulating ${target_label} from ${pathogen_fasta}"
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
    require_file "${DEPLETION_REF_FASTA}"
    minimap2 \
        --split-prefix "${OUT_DIR}/depletion_split_index" \
        -t "${THREADS}" \
        -a \
        -x map-ont \
        "${DEPLETION_REF_FASTA}" \
        "${REAL_FASTQ}" \
        | samtools view -h -T "${DEPLETION_REF_FASTA}" -b -f 4 -@ "${THREADS}" - \
        | samtools fastq -@ "${THREADS}" - \
        | gzip -c > "${WORK_FASTQ}"
else
    cp "${REAL_FASTQ}" "${WORK_FASTQ}"
fi

{
    printf 'parameter\tvalue\n'
    printf 'repo_dir\t%s\n' "${REPO_DIR}"
    printf 'pathogen_config_tsv\t%s\n' "${PATHOGEN_CONFIG_TSV}"
    printf 'threads\t%s\n' "${THREADS}"
    printf 'replicates\t%s\n' "${REPLICATES}"
    printf 'spike_levels\t%s\n' "${SPIKE_LEVELS}"
    printf 'minimap_db_fasta\t%s\n' "${MINIMAP_DB_FASTA}"
} > "${RUN_META_TSV}"

header='replicate\tspike_n\ttotal_spiked_reads\tmixed_fastq_gz\tkraken_report'
for line in "${PANEL_LINES[@]}"; do
    target_label="$(printf '%s' "${line}" | cut -f2)"
    safe_label="$(printf '%s' "${target_label}" | tr ' ' '_' | tr '/' '_')"
    header+="\tkraken_${safe_label}_reads\tminimap_${safe_label}_alignments"
done
printf '%b\n' "${header}" > "${SUMMARY_TSV}"

for rep in $(seq 1 "${REPLICATES}"); do
    for spike_n in ${SPIKE_LEVELS}; do
        MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"
        mkdir -p "${MIX_DIR}"
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
        python3 "${BUILD_MIXED_FASTQ_PY}" \
            --background_fastq_gz "${tmp_spikes[0]}" \
            $(for s in "${tmp_spikes[@]:1}"; do printf -- ' --spike_fastq_gz %q' "$s"; done) \
            --out_fastq_gz "${combined_spike}"

        mix_fastq_gz="${MIX_DIR}/mixed.fastq.gz"
        python3 "${BUILD_MIXED_FASTQ_PY}" \
            --background_fastq_gz "${WORK_FASTQ}" \
            --spike_fastq_gz "${combined_spike}" \
            --out_fastq_gz "${mix_fastq_gz}"

        kraken_report="${MIX_DIR}/kraken.report.tsv"
        kraken_out="${MIX_DIR}/kraken.classifications.tsv"
        kraken2 \
            --db "${KRAKEN_DB_DIR}" \
            --threads "${THREADS}" \
            --report "${kraken_report}" \
            --output "${kraken_out}" \
            "${mix_fastq_gz}" >/dev/null

        bam_out="${MIX_DIR}/minimap_sorted.bam"
        bed_out="${MIX_DIR}/minimap_sorted.MASKED.bed"
        run_minimap_to_bed \
            "${MINIMAP_DB_FASTA}" \
            "${mix_fastq_gz}" \
            "${bam_out}" \
            "${bed_out}" \
            "${MIX_DIR}/minimap_split_index"

        row="${rep}\t${spike_n}\t${total_spiked}\t${mix_fastq_gz}\t${kraken_report}"
        for line in "${PANEL_LINES[@]}"; do
            target_label="$(printf '%s' "${line}" | cut -f2)"
            safe_label="$(printf '%s' "${target_label}" | tr ' ' '_' | tr '/' '_')"
            ksum="${MIX_DIR}/kraken.${safe_label}.summary.tsv"
            python3 "${SUMMARISE_KRAKEN_PY}" \
                --kraken_report_tsv "${kraken_report}" \
                --target_label "${target_label}" \
                --out_tsv "${ksum}"
            ktarget="$(awk 'NR==2{print $3}' "${ksum}")"
            maltarget="$(python3 "${COUNT_BED_HITS_PY}" --bed "${bed_out}" --target_label "${target_label}")"
            row+="\t${ktarget}\t${maltarget}"
        done
        printf '%b\n' "${row}" >> "${SUMMARY_TSV}"
    done
done

log_info "Done: ${SUMMARY_TSV}"
