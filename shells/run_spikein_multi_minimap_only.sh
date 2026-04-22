#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N spike_multi_mm_only

set -euo pipefail

log_info() { printf '[INFO] %s\n' "$*"; }
log_warn() { printf '[WARN] %s\n' "$*" >&2; }
log_error() { printf '[ERROR] %s\n' "$*" >&2; }
require_file() { [[ -f "$1" ]] || { log_error "Required file not found: $1"; exit 1; }; }
require_exe() { command -v "$1" >/dev/null 2>&1 || { log_error "Required executable not found on PATH: $1"; exit 1; }; }

check_first_fasta_header() {
    local fasta="$1"
    python3 - "$fasta" <<'PY'
import sys
from pathlib import Path

path = Path(sys.argv[1])
with path.open("rt", errors="replace") as handle:
    for line_number, line in enumerate(handle, start=1):
        if not line.strip():
            continue
        if line.startswith(">"):
            print(line[1:].strip())
            raise SystemExit(0)
        print(
            f"First non-empty line does not begin with '>': line {line_number}",
            file=sys.stderr,
        )
        raise SystemExit(1)
print("No non-empty lines were found in the FASTA file.", file=sys.stderr)
raise SystemExit(1)
PY
}

validate_reference_fasta() {
    local fasta="$1"
    local header_preview

    require_file "${fasta}"
    [[ -s "${fasta}" ]] || {
        log_error "Reference FASTA exists but is empty: ${fasta}"
        exit 1
    }

    header_preview="$(check_first_fasta_header "${fasta}")" || {
        log_error "Reference FASTA does not appear to be valid FASTA: ${fasta}"
        exit 1
    }

    log_info "Using minimap reference FASTA: ${fasta}"
    log_info "Reference first header: ${header_preview}"
}

check_minimap_log_for_zero_targets() {
    local log_path="$1"

    if grep -q 'loaded/built the index for 0 target sequence(s)' "${log_path}"; then
        log_error "Minimap2 reported zero target sequences in the reference."
        log_error "Reference FASTA used: ${MINIMAP_DB_FASTA}"
        log_error "See minimap stderr log: ${log_path}"
        exit 1
    fi
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
SUMMARISE_MINIMAP_PY="${SUMMARISE_MINIMAP_PY:-${PY_SCRIPTS_DIR}/summarise_minimap_bed_hits_with_reported_taxa.py}"

REAL_FASTQ="${REAL_FASTQ:-${REAL_FASTQ_DEFAULT}}"
MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-${MONKEY_SMALL_GZ_DEFAULT}}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-${MONKEY_SMALL_FASTA_DEFAULT}}"
DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-${DEPLETION_REF_FASTA_DEFAULT}}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-${DEPLETED_FASTQ_DEFAULT}}"
MINIMAP_DB_FASTA="${MINIMAP_DB_FASTA:-${MINIMAP_DB_FASTA_DEFAULT}}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_minimap_only_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-${THREADS_DEFAULT}}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-${SPIKE_LEVELS_DEFAULT}}"
REPLICATES="${REPLICATES:-${REPLICATES_DEFAULT}}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"
MIN_MAPQ="${MIN_MAPQ:-15}"
MIN_READ_LEN="${MIN_READ_LEN:-500}"
IS_SHUFFLED_CONTROL="${IS_SHUFFLED_CONTROL:-false}"

require_exe python3
require_exe minimap2
require_exe samtools
require_exe bamToBed
require_exe read_analysis.py
require_exe simulator.py

require_file "${SAMPLE_FASTQ_PY}"
require_file "${BUILD_MIXED_FASTQ_PY}"
require_file "${COMBINE_NANOSIM_FASTQ_PY}"
require_file "${SUMMARISE_MINIMAP_PY}"
require_file "${REAL_FASTQ}"
require_file "${PATHOGEN_CONFIG_TSV}"
validate_reference_fasta "${MINIMAP_DB_FASTA}"
mkdir -p "${OUT_DIR}"

mapfile -t PANEL_LINES < <(awk 'BEGIN{FS="\t"} NR>1 && NF>=2 {print $1"\t"$2}' "${PATHOGEN_CONFIG_TSV}")
[[ ${#PANEL_LINES[@]} -gt 0 ]] || { log_error "No pathogen entries found in ${PATHOGEN_CONFIG_TSV}"; exit 1; }
N_GENOMES="${#PANEL_LINES[@]}"

TRAIN_FASTQ="${OUT_DIR}/train_reads.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"
SUMMARY_TSV="${OUT_DIR}/spikein_multi_minimap_only_summary.tsv"
RUN_META_TSV="${OUT_DIR}/run_metadata.tsv"
WORK_FASTQ="${DEPLETED_FASTQ}"

python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${REAL_FASTQ}" --n_reads "${TRAIN_READS_N}" --seed 1 --out_fastq_gz "${TRAIN_FASTQ}"

if [[ ! -s "${MONKEY_SMALL_FASTA}" ]]; then
    gzip -cd "${MONKEY_SMALL_GZ}" > "${MONKEY_SMALL_FASTA}"
fi

read_analysis.py genome --read "${TRAIN_FASTQ}" --ref_g "${MONKEY_SMALL_FASTA}" --output "${NS_MODEL_PREFIX}" --num_threads "${THREADS}" --fastq

GENOME_INDEX=0
for line in "${PANEL_LINES[@]}"; do
    GENOME_INDEX=$((GENOME_INDEX + 1))
    pathogen_fasta="$(printf '%s' "${line}" | cut -f1)"
    require_file "${pathogen_fasta}"
    sim_prefix="${OUT_DIR}/simulated_pathogen_${GENOME_INDEX}"
    sim_pool_fastq="${OUT_DIR}/sim_pool_${GENOME_INDEX}.fastq"
    sim_pool_fastq_gz="${OUT_DIR}/sim_pool_${GENOME_INDEX}.fastq.gz"
    simulator.py genome --ref_g "${pathogen_fasta}" --model_prefix "${NS_MODEL_PREFIX}" --output "${sim_prefix}" --number "${SIM_POOL_N}" -dna_type linear --seed "$((40 + GENOME_INDEX))" --num_threads "${THREADS}" --fastq
    python3 "${COMBINE_NANOSIM_FASTQ_PY}" --sim_prefix "${sim_prefix}" --out_fastq "${sim_pool_fastq}"
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

{
    printf 'parameter\tvalue\n'
    printf 'repo_dir\t%s\n' "${REPO_DIR}"
    printf 'pathogen_config_tsv\t%s\n' "${PATHOGEN_CONFIG_TSV}"
    printf 'minimap_db_fasta\t%s\n' "${MINIMAP_DB_FASTA}"
    printf 'threads\t%s\n' "${THREADS}"
    printf 'replicates\t%s\n' "${REPLICATES}"
    printf 'spike_levels\t%s\n' "${SPIKE_LEVELS}"
    printf 'min_mapq\t%s\n' "${MIN_MAPQ}"
    printf 'min_read_len\t%s\n' "${MIN_READ_LEN}"
    printf 'n_genomes\t%s\n' "${N_GENOMES}"
    printf 'is_shuffled_control\t%s\n' "${IS_SHUFFLED_CONTROL}"
} > "${RUN_META_TSV}"

header='replicate\tspike_n\ttotal_spiked_reads\tn_genomes\tis_shuffled_control\tmixed_fastq_gz\tminimap_bam\tminimap_bed\tminimap_target_summary_tsv\tminimap_reported_taxa_tsv\tminimap_meta_tsv\tminimap_total_alignments\tminimap_total_unique_reads'
for line in "${PANEL_LINES[@]}"; do
    target_label="$(printf '%s' "${line}" | cut -f2)"
    safe_label="$(printf '%s' "${target_label}" | tr ' ' '_' | tr '/' '_' | tr -cd '[:alnum:]_')"
    header+="\tminimap_${safe_label}_alignments\tminimap_${safe_label}_unique_reads"
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
            python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${sim_pool_fastq_gz}" --n_reads "${spike_n}" --seed "${spike_seed}" --out_fastq_gz "${spike_fastq_gz}"
            tmp_spikes+=("${spike_fastq_gz}")
            total_spiked=$(( total_spiked + spike_n ))
        done

        combined_spike="${MIX_DIR}/spike_combined.fastq.gz"
        mix_args=(--background_fastq_gz "${tmp_spikes[0]}")
        for spike_fastq_gz in "${tmp_spikes[@]:1}"; do
            mix_args+=(--spike_fastq_gz "${spike_fastq_gz}")
        done
        mix_args+=(--out_fastq_gz "${combined_spike}")
        python3 "${BUILD_MIXED_FASTQ_PY}" "${mix_args[@]}"

        MIX_FASTQ_GZ="${MIX_DIR}/mixed.fastq.gz"
        python3 "${BUILD_MIXED_FASTQ_PY}" --background_fastq_gz "${WORK_FASTQ}" --spike_fastq_gz "${combined_spike}" --out_fastq_gz "${MIX_FASTQ_GZ}"

        BAM_OUT="${MIX_DIR}/minimap_sorted.bam"
        BED_OUT="${MIX_DIR}/minimap_sorted.filtered.bed"
        TARGET_SUMMARY_TSV="${MIX_DIR}/minimap.target_summary.tsv"
        REPORTED_TAXA_TSV="${MIX_DIR}/minimap.reported_taxa.tsv"
        META_TSV="${MIX_DIR}/minimap.meta.tsv"
        MINIMAP_LOG="${MIX_DIR}/minimap.stderr.log"

        minimap2 -t "${THREADS}" -ax map-ont "${MINIMAP_DB_FASTA}" "${MIX_FASTQ_GZ}" 2> "${MINIMAP_LOG}" \
            | samtools view -h -q "${MIN_MAPQ}" -e "length(seq) >= ${MIN_READ_LEN}" - \
            | samtools sort -O bam -@ "${THREADS}" -o "${BAM_OUT}" -

        check_minimap_log_for_zero_targets "${MINIMAP_LOG}"

        bamToBed -i "${BAM_OUT}" > "${BED_OUT}"

        python3 "${SUMMARISE_MINIMAP_PY}" \
            --bed "${BED_OUT}" \
            --reference_fasta "${MINIMAP_DB_FASTA}" \
            --target_tsv "${PATHOGEN_CONFIG_TSV}" \
            --out_target_tsv "${TARGET_SUMMARY_TSV}" \
            --out_reported_taxa_tsv "${REPORTED_TAXA_TSV}" \
            --out_meta_tsv "${META_TSV}"

        require_file "${TARGET_SUMMARY_TSV}"
        require_file "${REPORTED_TAXA_TSV}"
        require_file "${META_TSV}"

        MINIMAP_TOTAL_ALIGNMENTS="$(awk -F '\t' '$1=="total_alignments"{print $2}' "${META_TSV}")"
        MINIMAP_TOTAL_UNIQUE_READS="$(awk -F '\t' '$1=="total_unique_reads"{print $2}' "${META_TSV}")"

        row="${rep}\t${spike_n}\t${total_spiked}\t${N_GENOMES}\t${IS_SHUFFLED_CONTROL}\t${MIX_FASTQ_GZ}\t${BAM_OUT}\t${BED_OUT}\t${TARGET_SUMMARY_TSV}\t${REPORTED_TAXA_TSV}\t${META_TSV}\t${MINIMAP_TOTAL_ALIGNMENTS}\t${MINIMAP_TOTAL_UNIQUE_READS}"
        for line in "${PANEL_LINES[@]}"; do
            target_label="$(printf '%s' "${line}" | cut -f2)"
            target_alignments="$(awk -F '\t' -v target_label="${target_label}" '$1==target_label{print $4}' "${TARGET_SUMMARY_TSV}")"
            target_unique_reads="$(awk -F '\t' -v target_label="${target_label}" '$1==target_label{print $5}' "${TARGET_SUMMARY_TSV}")"
            row+="\t${target_alignments}\t${target_unique_reads}"
        done
        printf '%b\n' "${row}" >> "${SUMMARY_TSV}"
    done
done

log_info "Done: ${SUMMARY_TSV}"
