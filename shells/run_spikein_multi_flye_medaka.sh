#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N spike_multi_flye

set -euo pipefail

log_info(){ printf '[INFO] %s\n' "$*"; }
log_error(){ printf '[ERROR] %s\n' "$*" >&2; }
require_file(){ [[ -f "$1" ]] || { log_error "Required file not found: $1"; exit 1; }; }
require_dir(){ [[ -d "$1" ]] || { log_error "Required directory not found: $1"; exit 1; }; }
require_exe(){ command -v "$1" >/dev/null 2>&1 || { log_error "Required executable not found on PATH: $1"; exit 1; }; }

run_medaka_consensus(){
  local reads_fastq="$1"
  local draft_fasta="$2"
  local medaka_out_dir="$3"
  local threads="$4"
  local model="$5"

  mkdir -p "${medaka_out_dir}"

  if [[ "${model}" == "auto" ]]; then
    log_info "Running Medaka with automatic model selection"
    medaka_consensus \
      -i "${reads_fastq}" \
      -d "${draft_fasta}" \
      -o "${medaka_out_dir}" \
      -t "${threads}"
  else
    log_info "Running Medaka with explicit model: ${model}"
    medaka_consensus \
      -i "${reads_fastq}" \
      -d "${draft_fasta}" \
      -o "${medaka_out_dir}" \
      -t "${threads}" \
      -m "${model}"
  fi
}

if [[ -n "${REPO_DIR:-}" ]]; then REPO_DIR="${REPO_DIR}"; elif [[ -n "${SGE_O_WORKDIR:-}" ]]; then REPO_DIR="${SGE_O_WORKDIR}"; else REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection"; fi
source "${REPO_DIR}/configs/pipeline_paths.sh"
PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR:-${REPO_DIR}/scripts}"
CONFIG_DIR="${CONFIG_DIR:-${REPO_DIR}/configs}"
PATHOGEN_CONFIG_TSV="${PATHOGEN_CONFIG_TSV:-${DEFAULT_PATHOGEN_PANEL_3}}"
SAMPLE_FASTQ_PY="${PY_SCRIPTS_DIR}/sample_fastq.py"
BUILD_MIXED_FASTQ_PY="${PY_SCRIPTS_DIR}/build_mixed_fastq.py"
COMBINE_NANOSIM_FASTQ_PY="${PY_SCRIPTS_DIR}/combine_nanosim_fastq.py"
SUMMARISE_KRAKEN_PY="${PY_SCRIPTS_DIR}/summarise_kraken_report_with_reported_taxa.py"
DEDUP_FASTQ_NAMES_PY="${PY_SCRIPTS_DIR}/dedup_fastq_names.py"
ASSEMBLY_STATS_PY="${PY_SCRIPTS_DIR}/assembly_stats.py"
REAL_FASTQ="${REAL_FASTQ:-${REAL_FASTQ_DEFAULT}}"
MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-${MONKEY_SMALL_GZ_DEFAULT}}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-${MONKEY_SMALL_FASTA_DEFAULT}}"
DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-${DEPLETION_REF_FASTA_DEFAULT}}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-${DEPLETED_FASTQ_DEFAULT}}"
KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-${KRAKEN_DB_DIR_DEFAULT}}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_flye_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-12}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-12}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"
MEDAKA_MODEL="${MEDAKA_MODEL:-auto}"

for x in python3 read_analysis.py simulator.py flye kraken2 minimap2 samtools medaka_consensus; do require_exe "$x"; done
for f in "${SAMPLE_FASTQ_PY}" "${BUILD_MIXED_FASTQ_PY}" "${COMBINE_NANOSIM_FASTQ_PY}" "${SUMMARISE_KRAKEN_PY}" "${DEDUP_FASTQ_NAMES_PY}" "${ASSEMBLY_STATS_PY}" "${PATHOGEN_CONFIG_TSV}" "${REAL_FASTQ}"; do require_file "$f"; done
require_dir "${KRAKEN_DB_DIR}"
mkdir -p "${OUT_DIR}"

mapfile -t PANEL_LINES < <(awk 'BEGIN{FS="\t"} NR>1 && NF>=2 {print $1"\t"$2}' "${PATHOGEN_CONFIG_TSV}")
[[ ${#PANEL_LINES[@]} -gt 0 ]] || { log_error "No pathogen entries found in ${PATHOGEN_CONFIG_TSV}"; exit 1; }

TRAIN_FASTQ="${OUT_DIR}/train_reads.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"
SUMMARY_TSV="${OUT_DIR}/spikein_multi_flye_summary.tsv"
WORK_FASTQ="${DEPLETED_FASTQ}"

python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${REAL_FASTQ}" --n_reads "${TRAIN_READS_N}" --seed 1 --out_fastq_gz "${TRAIN_FASTQ}"
if [[ ! -s "${MONKEY_SMALL_FASTA}" ]]; then gzip -cd "${MONKEY_SMALL_GZ}" > "${MONKEY_SMALL_FASTA}"; fi
read_analysis.py genome --read "${TRAIN_FASTQ}" --ref_g "${MONKEY_SMALL_FASTA}" --output "${NS_MODEL_PREFIX}" --num_threads "${THREADS}" --fastq

GENOME_INDEX=0
for line in "${PANEL_LINES[@]}"; do
  GENOME_INDEX=$((GENOME_INDEX+1))
  pathogen_fasta="$(printf '%s' "${line}" | cut -f1)"
  require_file "${pathogen_fasta}"
  sim_prefix="${OUT_DIR}/simulated_pathogen_${GENOME_INDEX}"
  sim_pool_fastq="${OUT_DIR}/sim_pool_${GENOME_INDEX}.fastq"
  sim_pool_fastq_gz="${OUT_DIR}/sim_pool_${GENOME_INDEX}.fastq.gz"
  simulator.py genome --ref_g "${pathogen_fasta}" --model_prefix "${NS_MODEL_PREFIX}" --output "${sim_prefix}" --number "${SIM_POOL_N}" -dna_type linear --seed "$((40+GENOME_INDEX))" --num_threads "${THREADS}" --fastq
  python3 "${COMBINE_NANOSIM_FASTQ_PY}" --sim_prefix "${sim_prefix}" --out_fastq "${sim_pool_fastq}"
  gzip -c "${sim_pool_fastq}" > "${sim_pool_fastq_gz}"
done

if [[ -s "${WORK_FASTQ}" ]]; then :; elif [[ "${DO_HOST_DEPLETION}" == "true" ]]; then minimap2 -t "${THREADS}" -a -x map-ont "${DEPLETION_REF_FASTA}" "${REAL_FASTQ}" | samtools view -h -T "${DEPLETION_REF_FASTA}" -b -f 4 -@ "${THREADS}" - | samtools fastq -@ "${THREADS}" - | gzip -c > "${WORK_FASTQ}"; else cp "${REAL_FASTQ}" "${WORK_FASTQ}"; fi



header='replicate	spike_n	total_spiked_reads	mixed_fastq_gz	flye_out_dir	flye_assembly_fasta	polished_assembly_fasta	contig_count	total_bases	kraken_report'
for line in "${PANEL_LINES[@]}"; do target_label="$(printf '%s' "${line}" | cut -f2)"; safe_label="$(printf '%s' "${target_label}" | tr ' ' '_' | tr '/' '_')"; header+="	kraken_${safe_label}_contigs"; done
printf '%b\n' "${header}" > "${SUMMARY_TSV}"

for rep in $(seq 1 "${REPLICATES}"); do
 for spike_n in ${SPIKE_LEVELS}; do
  MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"
  mkdir -p "${MIX_DIR}"
  tmp_spikes=()
  total_spiked=0
  GENOME_INDEX=0

  for line in "${PANEL_LINES[@]}"; do
    GENOME_INDEX=$((GENOME_INDEX+1))
    sim_pool_fastq_gz="${OUT_DIR}/sim_pool_${GENOME_INDEX}.fastq.gz"
    spike_fastq_gz="${MIX_DIR}/spike_${GENOME_INDEX}.fastq.gz"
    spike_seed=$(( rep * 100000 + GENOME_INDEX * 1000 + spike_n ))
    python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${sim_pool_fastq_gz}" --n_reads "${spike_n}" --seed "${spike_seed}" --out_fastq_gz "${spike_fastq_gz}"
    tmp_spikes+=("${spike_fastq_gz}")
    total_spiked=$(( total_spiked + spike_n ))
  done


  mix_fastq_gz="${MIX_DIR}/mixed.fastq.gz"
  python3 "${BUILD_MIXED_FASTQ_PY}" \
    --background_fastq_gz "${WORK_FASTQ}" \
    --spike_fastq_gz "${tmp_spikes[@]}" \
    --out_fastq_gz "${mix_fastq_gz}"

  dedup_fastq="${MIX_DIR}/mixed.dedup.fastq"
  python3 "${DEDUP_FASTQ_NAMES_PY}" --input_fastq_gz "${mix_fastq_gz}" --output_fastq "${dedup_fastq}"

  flye_out_dir="${MIX_DIR}/flye_out"
  assembly_fasta="${flye_out_dir}/assembly.fasta"
  medaka_out_dir="${MIX_DIR}/medaka_out"
  polished_fasta="${medaka_out_dir}/consensus.fasta"
  asm_stats_tsv="${MIX_DIR}/assembly.stats.tsv"

  flye --meta --nano-hq "${dedup_fastq}" --out-dir "${flye_out_dir}" --threads "${THREADS}"
  run_medaka_consensus "${dedup_fastq}" "${assembly_fasta}" "${medaka_out_dir}" "${THREADS}" "${MEDAKA_MODEL}"
  require_file "${polished_fasta}"

  python3 "${ASSEMBLY_STATS_PY}" --assembly_fasta "${polished_fasta}" --out_tsv "${asm_stats_tsv}"

  kraken_report="${MIX_DIR}/assembly.kraken.report.tsv"
  kraken_out="${MIX_DIR}/assembly.kraken.classifications.tsv"
  kraken2 --db "${KRAKEN_DB_DIR}" --threads "${THREADS}" --report "${kraken_report}" --output "${kraken_out}" "${polished_fasta}" >/dev/null

  contig_count="$(awk 'NR==2{print $2}' "${asm_stats_tsv}")"
  total_bases="$(awk 'NR==2{print $3}' "${asm_stats_tsv}")"
  row="${rep}	${spike_n}	${total_spiked}	${mix_fastq_gz}	${flye_out_dir}	${assembly_fasta}	${polished_fasta}	${contig_count}	${total_bases}	${kraken_report}"

  for line in "${PANEL_LINES[@]}"; do
    target_label="$(printf '%s' "${line}" | cut -f2)"
    safe_label="$(printf '%s' "${target_label}" | tr ' ' '_' | tr '/' '_')"
    ksum="${MIX_DIR}/assembly.kraken.${safe_label}.summary.tsv"
    python3 "${SUMMARISE_KRAKEN_PY}" --kraken_report_tsv "${kraken_report}" --target_label "${target_label}" --out_tsv "${ksum}"
    ktarget="$(awk 'NR==2{print $3}' "${ksum}")"
    row+="	${ktarget}"
  done

  printf '%b\n' "${row}" >> "${SUMMARY_TSV}"
 done
done

log_info "Done: ${SUMMARY_TSV}"
