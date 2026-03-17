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
MINIMAP_DB_GZ="${MINIMAP_DB_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"
KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"
ALIASES_TSV="${ALIASES_TSV:-}"
PATHOGEN_CONFIG_TSV="${PATHOGEN_CONFIG_TSV:-${REPO_DIR}/configs/pathogen_panel_3.tsv}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_read_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-12}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000}"
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
for f in "${REAL_FASTQ}" "${MONKEY_DEPLETION_FASTA}" "${PLASMO_MASKED_GZ}" "${MINIMAP_DB_GZ}" "${PATHOGEN_CONFIG_TSV}"; do require_file "$f"; done
[[ -f "${MONKEY_SMALL_GZ}" || -f "${MONKEY_SMALL_FASTA}" ]] || { log_error "Neither MONKEY_SMALL_GZ nor MONKEY_SMALL_FASTA exists."; exit 1; }
require_dir "${KRAKEN_DB_DIR}"
for f in "${SAMPLE_FASTQ_PY}" "${COMBINE_NANOSIM_FASTQ_PY}" "${BUILD_MIXED_FASTQ_PY}" "${SUMMARISE_KRAKEN_PY}" "${COUNT_BED_HITS_PY}"; do require_file "$f"; done
for x in python3 minimap2 samtools gzip kraken2 bamToBed read_analysis.py simulator.py; do require_exe "$x"; done

mapfile -t PANEL_LINES < <(awk 'BEGIN{FS="\t"} NR>1 && NF>=2 {print $1"\t"$2}' "${PATHOGEN_CONFIG_TSV}")
[[ "${#PANEL_LINES[@]}" -ge 2 ]] || { log_error "PATHOGEN_CONFIG_TSV must contain at least two genomes."; exit 1; }

declare -a PATHOGEN_FASTA_ARR=()
declare -a TARGET_LABEL_ARR=()
for line in "${PANEL_LINES[@]}"; do
  fasta="${line%%$'\t'*}"
  label="${line#*$'\t'}"
  require_file "${fasta}"
  PATHOGEN_FASTA_ARR+=("${fasta}")
  TARGET_LABEL_ARR+=("${label}")
done
N_GENOMES="${#PATHOGEN_FASTA_ARR[@]}"

if [[ ! -s "${MONKEY_SMALL_FASTA}" ]]; then gzip -cd "${MONKEY_SMALL_GZ}" > "${MONKEY_SMALL_FASTA}"; fi
if [[ ! -s "${DEPLETION_REF_FASTA}" ]]; then cat "${MONKEY_DEPLETION_FASTA}" > "${DEPLETION_REF_FASTA}"; gzip -cd "${PLASMO_MASKED_GZ}" >> "${DEPLETION_REF_FASTA}"; fi

TRAIN_FASTQ="${OUT_DIR}/train_${TRAIN_READS_N}.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"
MINIMAP_DB_FASTA="${OUT_DIR}/plas_outgrps_genomes_Hard_MASKED.fasta"
SUMMARY_TSV="${OUT_DIR}/spikein_multi_summary.tsv"
RUN_META_TSV="${OUT_DIR}/run_metadata.tsv"

python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${REAL_FASTQ}" --n_reads "${TRAIN_READS_N}" --seed 1 --out_fastq_gz "${TRAIN_FASTQ}"
read_analysis.py genome --read "${TRAIN_FASTQ}" --ref_g "${MONKEY_SMALL_FASTA}" --output "${NS_MODEL_PREFIX}" --num_threads "${THREADS}" --fastq

SIM_POOL_FASTQ_GZ_ARR=()
for i in "${!PATHOGEN_FASTA_ARR[@]}"; do
  genome_idx=$(( i + 1 ))
  sim_prefix="${OUT_DIR}/simulated_pathogen_${genome_idx}"
  sim_pool_fastq="${OUT_DIR}/sim_pool_${genome_idx}.fastq"
  sim_pool_fastq_gz="${OUT_DIR}/sim_pool_${genome_idx}.fastq.gz"
  simulator.py genome --ref_g "${PATHOGEN_FASTA_ARR[$i]}" --model_prefix "${NS_MODEL_PREFIX}" --output "${sim_prefix}" --number "${SIM_POOL_N}" --dna_type linear --seed $((42 + i)) --num_threads "${THREADS}" --fastq
  python3 "${COMBINE_NANOSIM_FASTQ_PY}" --sim_prefix "${sim_prefix}" --out_fastq "${sim_pool_fastq}"
  gzip -c "${sim_pool_fastq}" > "${sim_pool_fastq_gz}"
  SIM_POOL_FASTQ_GZ_ARR+=("${sim_pool_fastq_gz}")
done

WORK_FASTQ="${DEPLETED_FASTQ}"
if [[ -s "${WORK_FASTQ}" ]]; then log_info "Reusing existing depleted background: ${WORK_FASTQ}";
elif [[ "${DO_HOST_DEPLETION}" == "true" ]]; then
  minimap2 -t "${THREADS}" -a -x map-ont "${DEPLETION_REF_FASTA}" "${REAL_FASTQ}" | samtools view -h -T "${DEPLETION_REF_FASTA}" -b -f 4 -@ "${THREADS}" - | samtools fastq -@ "${THREADS}" - | gzip -c > "${WORK_FASTQ}"
else
  cp "${REAL_FASTQ}" "${WORK_FASTQ}"
fi

if [[ ! -s "${MINIMAP_DB_FASTA}" ]]; then gzip -cd "${MINIMAP_DB_GZ}" > "${MINIMAP_DB_FASTA}"; fi

{
  printf 'parameter\tvalue\n'
  printf 'pipeline\tmulti_readlevel\n'
  printf 'pathogen_config_tsv\t%s\n' "${PATHOGEN_CONFIG_TSV}"
  printf 'n_genomes\t%s\n' "${N_GENOMES}"
} > "${RUN_META_TSV}"

header=(replicate spike_n_per_genome n_genomes mixed_fastq_gz kraken_report kraken_classified_reads)
for label in "${TARGET_LABEL_ARR[@]}"; do header+=("kraken_target_reads_${label// /_}"); done
header+=(minimap_bam)
for label in "${TARGET_LABEL_ARR[@]}"; do header+=("minimap_target_alignments_${label// /_}"); done
(IFS=$'\t'; printf '%s\n' "${header[*]}") > "${SUMMARY_TSV}"

for rep in $(seq 1 "${REPLICATES}"); do
  for spike_n in ${SPIKE_LEVELS}; do
    MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"
    mkdir -p "${MIX_DIR}"
    spike_args=()
    for i in "${!SIM_POOL_FASTQ_GZ_ARR[@]}"; do
      spike_fastq_gz="${MIX_DIR}/spike_genome_$((i + 1)).fastq.gz"
      seed=$(( rep * 100000 + (i + 1) * 1000 + spike_n ))
      python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${SIM_POOL_FASTQ_GZ_ARR[$i]}" --n_reads "${spike_n}" --seed "${seed}" --out_fastq_gz "${spike_fastq_gz}"
      spike_args+=(--spike_fastq_gz "${spike_fastq_gz}")
    done
    MIX_FASTQ_GZ="${MIX_DIR}/mixed.fastq.gz"
    KRAKEN_REPORT="${MIX_DIR}/kraken.report.tsv"
    KRAKEN_OUT="${MIX_DIR}/kraken.classifications.tsv"
    KRAKEN_SUM_TSV="${MIX_DIR}/kraken.summary.tsv"
    BAM_OUT="${MIX_DIR}/minimap_sorted.bam"
    BED_OUT="${MIX_DIR}/minimap_sorted.MASKED.bed"

    python3 "${BUILD_MIXED_FASTQ_PY}" --background_fastq_gz "${WORK_FASTQ}" "${spike_args[@]}" --out_fastq_gz "${MIX_FASTQ_GZ}"
    KRAKEN_LABEL_ARGS=()
    for label in "${TARGET_LABEL_ARR[@]}"; do KRAKEN_LABEL_ARGS+=(--target_label "${label}"); done
    kraken2 --db "${KRAKEN_DB_DIR}" --threads "${THREADS}" --report "${KRAKEN_REPORT}" --output "${KRAKEN_OUT}" "${MIX_FASTQ_GZ}" >/dev/null
    python3 "${SUMMARISE_KRAKEN_PY}" --kraken_report_tsv "${KRAKEN_REPORT}" "${KRAKEN_LABEL_ARGS[@]}" --out_tsv "${KRAKEN_SUM_TSV}"
    minimap2 -t "${THREADS}" -ax map-ont "${MINIMAP_DB_FASTA}" "${MIX_FASTQ_GZ}" | samtools view -Sh -q "${MIN_MAPQ}" -e "rlen >= ${MIN_ALIGN_LEN}" - | samtools sort -O bam -@ "${THREADS}" - | tee "${BAM_OUT}" | bamToBed > "${BED_OUT}"

    row=("${rep}" "${spike_n}" "${N_GENOMES}" "${MIX_FASTQ_GZ}" "${KRAKEN_REPORT}" "$(awk 'NR==2{print $2}' "${KRAKEN_SUM_TSV}")")
    field_idx=3
    for _ in "${TARGET_LABEL_ARR[@]}"; do row+=("$(awk -v idx="${field_idx}" 'NR==2{print $idx}' "${KRAKEN_SUM_TSV}")"); field_idx=$(( field_idx + 1 )); done
    row+=("${BAM_OUT}")
    for label in "${TARGET_LABEL_ARR[@]}"; do
      BED_ARGS=(--bed "${BED_OUT}" --target_label "${label}")
      if [[ -n "${ALIASES_TSV}" ]]; then BED_ARGS+=(--aliases_tsv "${ALIASES_TSV}"); fi
      row+=("$(python3 "${COUNT_BED_HITS_PY}" "${BED_ARGS[@]}")")
    done
    (IFS=$'\t'; printf '%s\n' "${row[*]}") >> "${SUMMARY_TSV}"
  done
done

log_info "Done: ${SUMMARY_TSV}"
