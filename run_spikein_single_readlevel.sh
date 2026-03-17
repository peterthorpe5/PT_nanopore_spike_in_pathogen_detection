#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N spike1_read

set -euo pipefail

log_info() { printf '[INFO] %s\n' "$*"; }
log_error() { printf '[ERROR] %s\n' "$*" >&2; }
require_file() { [[ -f "$1" ]] || { log_error "Required file not found: $1"; exit 1; }; }
require_dir() { [[ -d "$1" ]] || { log_error "Required directory not found: $1"; exit 1; }; }
require_exe() { command -v "$1" >/dev/null 2>&1 || { log_error "Required executable not found: $1"; exit 1; }; }

default_repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-${default_repo_dir}}"
PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR:-${REPO_DIR}/scripts}"

REAL_FASTQ="${REAL_FASTQ:-/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC1023_AmM008WB.fastq.gz}"
MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna.gz}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna}"
MONKEY_DEPLETION_FASTA="${MONKEY_DEPLETION_FASTA:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/M.nemestrina_tonkeana_nigra.fasta}"
PLASMO_MASKED_GZ="${PLASMO_MASKED_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"
DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/host_plus_plasmo_depletion.fasta}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/MRC1023_AmM008WB.depleted.fastq.gz}"
PATHOGEN_FASTA="${PATHOGEN_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/GCA_014843685.1_ASM1484368v1_genomic.fna}"
TARGET_LABEL="${TARGET_LABEL:-Plasmodium vivax}"
MINIMAP_DB_GZ="${MINIMAP_DB_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"
KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"
ALIASES_TSV="${ALIASES_TSV:-}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_single_read_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-12}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-3}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
MIN_ALIGN_LEN="${MIN_ALIGN_LEN:-500}"
MIN_MAPQ="${MIN_MAPQ:-15}"

SAMPLE_FASTQ_PY="${PY_SCRIPTS_DIR}/sample_fastq.py"
COMBINE_NANOSIM_FASTQ_PY="${PY_SCRIPTS_DIR}/combine_nanosim_fastq.py"
BUILD_MIXED_FASTQ_PY="${PY_SCRIPTS_DIR}/build_mixed_fastq.py"
SUMMARISE_KRAKEN_PY="${PY_SCRIPTS_DIR}/summarise_kraken_report.py"
COUNT_BED_HITS_PY="${PY_SCRIPTS_DIR}/count_bed_hits.py"

mkdir -p "${OUT_DIR}" "$(dirname "${MONKEY_SMALL_FASTA}")" "$(dirname "${DEPLETION_REF_FASTA}")" "$(dirname "${DEPLETED_FASTQ}")"
for f in "${REAL_FASTQ}" "${PATHOGEN_FASTA}" "${MONKEY_DEPLETION_FASTA}" "${PLASMO_MASKED_GZ}" "${MINIMAP_DB_GZ}"; do require_file "$f"; done
[[ -f "${MONKEY_SMALL_GZ}" || -f "${MONKEY_SMALL_FASTA}" ]] || { log_error "Neither MONKEY_SMALL_GZ nor MONKEY_SMALL_FASTA exists."; exit 1; }
require_dir "${KRAKEN_DB_DIR}"
for f in "${SAMPLE_FASTQ_PY}" "${COMBINE_NANOSIM_FASTQ_PY}" "${BUILD_MIXED_FASTQ_PY}" "${SUMMARISE_KRAKEN_PY}" "${COUNT_BED_HITS_PY}"; do require_file "$f"; done
for x in python3 minimap2 samtools gzip kraken2 bamToBed read_analysis.py simulator.py; do require_exe "$x"; done

if [[ ! -s "${MONKEY_SMALL_FASTA}" ]]; then gzip -cd "${MONKEY_SMALL_GZ}" > "${MONKEY_SMALL_FASTA}"; fi
if [[ ! -s "${DEPLETION_REF_FASTA}" ]]; then cat "${MONKEY_DEPLETION_FASTA}" > "${DEPLETION_REF_FASTA}"; gzip -cd "${PLASMO_MASKED_GZ}" >> "${DEPLETION_REF_FASTA}"; fi

TRAIN_FASTQ="${OUT_DIR}/train_${TRAIN_READS_N}.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"
SIM_PREFIX="${OUT_DIR}/simulated_pathogen"
SIM_POOL_FASTQ="${OUT_DIR}/sim_pool.fastq"
SIM_POOL_FASTQ_GZ="${OUT_DIR}/sim_pool.fastq.gz"
MINIMAP_DB_FASTA="${OUT_DIR}/plas_outgrps_genomes_Hard_MASKED.fasta"
SUMMARY_TSV="${OUT_DIR}/spikein_summary.tsv"
RUN_META_TSV="${OUT_DIR}/run_metadata.tsv"

python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${REAL_FASTQ}" --n_reads "${TRAIN_READS_N}" --seed 1 --out_fastq_gz "${TRAIN_FASTQ}"
read_analysis.py genome --read "${TRAIN_FASTQ}" --ref_g "${MONKEY_SMALL_FASTA}" --output "${NS_MODEL_PREFIX}" --num_threads "${THREADS}" --fastq
simulator.py genome --ref_g "${PATHOGEN_FASTA}" --model_prefix "${NS_MODEL_PREFIX}" --output "${SIM_PREFIX}" --number "${SIM_POOL_N}" --dna_type linear --seed 42 --num_threads "${THREADS}" --fastq
python3 "${COMBINE_NANOSIM_FASTQ_PY}" --sim_prefix "${SIM_PREFIX}" --out_fastq "${SIM_POOL_FASTQ}"
gzip -c "${SIM_POOL_FASTQ}" > "${SIM_POOL_FASTQ_GZ}"

WORK_FASTQ="${DEPLETED_FASTQ}"
if [[ -s "${WORK_FASTQ}" ]]; then
  log_info "Reusing existing depleted background: ${WORK_FASTQ}"
elif [[ "${DO_HOST_DEPLETION}" == "true" ]]; then
  minimap2 -t "${THREADS}" -a -x map-ont "${DEPLETION_REF_FASTA}" "${REAL_FASTQ}" | samtools view -h -T "${DEPLETION_REF_FASTA}" -b -f 4 -@ "${THREADS}" - | samtools fastq -@ "${THREADS}" - | gzip -c > "${WORK_FASTQ}"
else
  cp "${REAL_FASTQ}" "${WORK_FASTQ}"
fi

if [[ ! -s "${MINIMAP_DB_FASTA}" ]]; then gzip -cd "${MINIMAP_DB_GZ}" > "${MINIMAP_DB_FASTA}"; fi

{
  printf 'parameter\tvalue\n'
  printf 'pipeline\tsingle_readlevel\n'
  printf 'pathogen_fasta\t%s\n' "${PATHOGEN_FASTA}"
  printf 'target_label\t%s\n' "${TARGET_LABEL}"
  printf 'aliases_tsv\t%s\n' "${ALIASES_TSV}"
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

    BED_ARGS=(--bed "${BED_OUT}" --target_label "${TARGET_LABEL}")
    if [[ -n "${ALIASES_TSV}" ]]; then BED_ARGS+=(--aliases_tsv "${ALIASES_TSV}"); fi
    MINIMAP_TARGET_ALN="$(python3 "${COUNT_BED_HITS_PY}" "${BED_ARGS[@]}")"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${rep}" "${spike_n}" "${MIX_FASTQ_GZ}" "${KRAKEN_REPORT}" "${K_CLASSIFIED}" "${K_TARGET}" "${BAM_OUT}" "${MINIMAP_TARGET_ALN}" >> "${SUMMARY_TSV}"
  done
done

log_info "Done: ${SUMMARY_TSV}"
