#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N spike_single_read

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

SAMPLE_FASTQ_PY="${PY_SCRIPTS_DIR}/sample_fastq.py"
BUILD_MIXED_FASTQ_PY="${PY_SCRIPTS_DIR}/build_mixed_fastq.py"
COMBINE_NANOSIM_FASTQ_PY="${PY_SCRIPTS_DIR}/combine_nanosim_fastq.py"
SUMMARISE_KRAKEN_PY="${PY_SCRIPTS_DIR}/summarise_kraken_report_with_reported_taxa.py"
COUNT_BED_HITS_PY="${PY_SCRIPTS_DIR}/count_bed_hits.py"

REAL_FASTQ="${REAL_FASTQ:-${REAL_FASTQ_DEFAULT}}"
MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-${MONKEY_SMALL_GZ_DEFAULT}}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-${MONKEY_SMALL_FASTA_DEFAULT}}"
MONKEY_DEPLETION_FASTA="${MONKEY_DEPLETION_FASTA:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/M.nemestrina_tonkeana_nigra.fasta}"
PLASMO_MASKED_GZ="${PLASMO_MASKED_GZ:-${MINIMAP_DB_FASTA_DEFAULT}}"
DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-${DEPLETION_REF_FASTA_DEFAULT}}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-${DEPLETED_FASTQ_DEFAULT}}"
PATHOGEN_FASTA="${PATHOGEN_FASTA:-${DEFAULT_PATHOGEN_FASTA}}"
TARGET_LABEL="${TARGET_LABEL:-${DEFAULT_TARGET_LABEL}}"
MINIMAP_DB_FASTA="${MINIMAP_DB_FASTA:-${MINIMAP_DB_FASTA_DEFAULT}}"
KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-${KRAKEN_DB_DIR_DEFAULT}}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_single_read_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-12}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-12}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"
MIN_MAPQ="${MIN_MAPQ:-15}"
MIN_ALIGN_LEN="${MIN_ALIGN_LEN:-500}"

require_exe python3
require_exe minimap2
require_exe samtools
require_exe kraken2
require_exe bamToBed
require_exe read_analysis.py
require_exe simulator.py

require_file "${SAMPLE_FASTQ_PY}"
require_file "${BUILD_MIXED_FASTQ_PY}"
require_file "${COMBINE_NANOSIM_FASTQ_PY}"
require_file "${SUMMARISE_KRAKEN_PY}"
require_file "${COUNT_BED_HITS_PY}"
require_file "${REAL_FASTQ}"
require_file "${PATHOGEN_FASTA}"
require_file "${MONKEY_DEPLETION_FASTA}"
require_file "${PLASMO_MASKED_GZ}"
require_file "${MINIMAP_DB_FASTA}"
require_dir "${KRAKEN_DB_DIR}"
mkdir -p "${OUT_DIR}"

TRAIN_FASTQ="${OUT_DIR}/train_reads.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"
SIM_PREFIX="${OUT_DIR}/simulated_pathogen"
SIM_POOL_FASTQ="${OUT_DIR}/sim_pool.fastq"
SIM_POOL_FASTQ_GZ="${OUT_DIR}/sim_pool.fastq.gz"
# MINIMAP_DB_FASTA="${OUT_DIR}/plas_outgrps_genomes_Hard_MASKED.fasta"
SUMMARY_TSV="${OUT_DIR}/spikein_summary.tsv"
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
    printf 'threads\t%s\n' "${THREADS}"
    printf 'replicates\t%s\n' "${REPLICATES}"
    printf 'spike_levels\t%s\n' "${SPIKE_LEVELS}"
} > "${RUN_META_TSV}"

printf 'replicate\tspike_n\tmixed_fastq_gz\tkraken_report\tkraken_classified_reads\tkraken_target_reads\tminimap_bam\tminimap_target_alignments\n' > "${SUMMARY_TSV}"

for rep in $(seq 1 "${REPLICATES}"); do
    for spike_n in ${SPIKE_LEVELS}; do
        MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"
        mkdir -p "${MIX_DIR}"
        SPIKE_FASTQ_GZ="${MIX_DIR}/spike.fastq.gz"
        MIX_FASTQ_GZ="${MIX_DIR}/mixed.fastq.gz"
        KRAKEN_REPORT="${MIX_DIR}/kraken.report.tsv"
        KRAKEN_OUT="${MIX_DIR}/kraken.classifications.tsv"
        KRAKEN_SUM_TSV="${MIX_DIR}/kraken.summary.tsv"
        BAM_OUT="${MIX_DIR}/minimap_sorted.bam"
        BED_OUT="${MIX_DIR}/minimap_sorted.MASKED.bed"
        SPIKE_SEED=$(( rep * 100000 + spike_n ))

        python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${SIM_POOL_FASTQ_GZ}" --n_reads "${spike_n}" --seed "${SPIKE_SEED}" --out_fastq_gz "${SPIKE_FASTQ_GZ}"
        python3 "${BUILD_MIXED_FASTQ_PY}" --background_fastq_gz "${WORK_FASTQ}" --spike_fastq_gz "${SPIKE_FASTQ_GZ}" --out_fastq_gz "${MIX_FASTQ_GZ}"

        kraken2 --db "${KRAKEN_DB_DIR}" --threads "${THREADS}" --report "${KRAKEN_REPORT}" --output "${KRAKEN_OUT}" "${MIX_FASTQ_GZ}" >/dev/null
        python3 "${SUMMARISE_KRAKEN_PY}" --kraken_report_tsv "${KRAKEN_REPORT}" --target_label "${TARGET_LABEL}" --out_tsv "${KRAKEN_SUM_TSV}"
        K_CLASSIFIED="$(awk 'NR==2{print $2}' "${KRAKEN_SUM_TSV}")"
        K_TARGET="$(awk 'NR==2{print $3}' "${KRAKEN_SUM_TSV}")"

        minimap2 -t "${THREADS}" -ax map-ont "${MINIMAP_DB_FASTA}" "${MIX_FASTQ_GZ}" | samtools view -Sh -q "${MIN_MAPQ}" -e "rlen >= ${MIN_ALIGN_LEN}" - | samtools sort -O bam -@ "${THREADS}" - | tee "${BAM_OUT}" | bamToBed > "${BED_OUT}"
        MINIMAP_TARGET_ALN="$(python3 "${COUNT_BED_HITS_PY}" --bed "${BED_OUT}" --target_label "${TARGET_LABEL}")"

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${rep}" "${spike_n}" "${MIX_FASTQ_GZ}" "${KRAKEN_REPORT}" "${K_CLASSIFIED}" "${K_TARGET}" "${BAM_OUT}" "${MINIMAP_TARGET_ALN}" >> "${SUMMARY_TSV}"
    done
done

log_info "Done: ${SUMMARY_TSV}"
