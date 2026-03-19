#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 240G
#$ -adds l_hard h_vmem 240G
#$ -N MMspike_multi_mm

set -euo pipefail

log_info() { printf '[INFO] %s\n' "$*"; }
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

PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR:-${REPO_DIR}/scripts}"
CONFIG_DIR="${CONFIG_DIR:-${REPO_DIR}/configs}"
PATHOGEN_CONFIG_TSV="${PATHOGEN_CONFIG_TSV:-${CONFIG_DIR}/pathogen_panel_3_with_taxid.tsv}"
SAMPLE_FASTQ_PY="${SAMPLE_FASTQ_PY:-${PY_SCRIPTS_DIR}/sample_fastq.py}"
BUILD_MIXED_FASTQ_PY="${BUILD_MIXED_FASTQ_PY:-${PY_SCRIPTS_DIR}/build_mixed_fastq.py}"
COMBINE_NANOSIM_FASTQ_PY="${COMBINE_NANOSIM_FASTQ_PY:-${PY_SCRIPTS_DIR}/combine_nanosim_fastq.py}"
SUMMARISE_METAMAPS_WIMP_PY="${SUMMARISE_METAMAPS_WIMP_PY:-${PY_SCRIPTS_DIR}/summarise_metamaps_wimp.py}"

REAL_FASTQ="${REAL_FASTQ:-/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC1023_AmM008WB.fastq.gz}"
MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna.gz}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna}"
DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/host_plus_plasmo_depletion.fasta}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/MRC1023_AmM008WB.depleted.fastq.gz}"
METAMAPS_DB_DIR="${METAMAPS_DB_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/metamaps_db/custom_latest/custom_metamaps_db}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_multi_metamaps_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-12}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-3}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"
MAX_MEMORY_GB="${MAX_MEMORY_GB:-140}"

for x in python3 minimap2 samtools read_analysis.py simulator.py metamaps; do require_exe "$x"; done
for f in "${SAMPLE_FASTQ_PY}" "${BUILD_MIXED_FASTQ_PY}" "${COMBINE_NANOSIM_FASTQ_PY}" "${SUMMARISE_METAMAPS_WIMP_PY}" "${REAL_FASTQ}" "${PATHOGEN_CONFIG_TSV}"; do require_file "$f"; done
require_dir "${METAMAPS_DB_DIR}"
require_file "${METAMAPS_DB_DIR}/DB.fa"
require_dir "${METAMAPS_DB_DIR}/taxonomy"
mkdir -p "${OUT_DIR}"

mapfile -t PANEL_LINES < <(awk 'BEGIN{FS="\t"} NR>1 && NF>=2 {print $1"\t"$2}' "${PATHOGEN_CONFIG_TSV}")
[[ ${#PANEL_LINES[@]} -gt 0 ]] || { log_error "No pathogen entries found in ${PATHOGEN_CONFIG_TSV}"; exit 1; }

TRAIN_FASTQ="${OUT_DIR}/train_reads.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"
SUMMARY_TSV="${OUT_DIR}/spikein_multi_metamaps_summary.tsv"
WORK_FASTQ="${DEPLETED_FASTQ}"

python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${REAL_FASTQ}" --n_reads "${TRAIN_READS_N}" --seed 1 --out_fastq_gz "${TRAIN_FASTQ}"
if [[ ! -s "${MONKEY_SMALL_FASTA}" ]]; then gzip -cd "${MONKEY_SMALL_GZ}" > "${MONKEY_SMALL_FASTA}"; fi
read_analysis.py genome --read "${TRAIN_FASTQ}" --ref_g "${MONKEY_SMALL_FASTA}" --output "${NS_MODEL_PREFIX}" --num_threads "${THREADS}" --fastq

GENOME_INDEX=0
for line in "${PANEL_LINES[@]}"; do
    GENOME_INDEX=$((GENOME_INDEX + 1))
    pathogen_fasta="$(printf '%s' "${line}" | cut -f1)"
    target_label="$(printf '%s' "${line}" | cut -f2)"
    require_file "${pathogen_fasta}"
    sim_prefix="${OUT_DIR}/simulated_pathogen_${GENOME_INDEX}"
    sim_pool_fastq="${OUT_DIR}/sim_pool_${GENOME_INDEX}.fastq"
    sim_pool_fastq_gz="${OUT_DIR}/sim_pool_${GENOME_INDEX}.fastq.gz"
    simulator.py genome --ref_g "${pathogen_fasta}" --model_prefix "${NS_MODEL_PREFIX}" --output "${sim_prefix}" --number "${SIM_POOL_N}" -dna_type linear --seed "$((40 + GENOME_INDEX))" --num_threads "${THREADS}" --fastq
    python3 "${COMBINE_NANOSIM_FASTQ_PY}" --sim_prefix "${sim_prefix}" --out_fastq "${sim_pool_fastq}"
    gzip -c "${sim_pool_fastq}" > "${sim_pool_fastq_gz}"
done

if [[ -s "${WORK_FASTQ}" ]]; then :; elif [[ "${DO_HOST_DEPLETION}" == "true" ]]; then minimap2 -t "${THREADS}" -a -x map-ont "${DEPLETION_REF_FASTA}" "${REAL_FASTQ}" | samtools view -h -T "${DEPLETION_REF_FASTA}" -b -f 4 -@ "${THREADS}" - | samtools fastq -@ "${THREADS}" - | gzip -c > "${WORK_FASTQ}"; else cp "${REAL_FASTQ}" "${WORK_FASTQ}"; fi

header='replicate\tspike_n\ttotal_spiked_reads\tmixed_fastq_gz\tmetamaps_prefix\tmetamaps_wimp\tmetamaps_total_reads_definedGenomes'
for line in "${PANEL_LINES[@]}"; do target_label="$(printf '%s' "${line}" | cut -f2)"; safe_label="$(printf '%s' "${target_label}" | tr ' ' '_' | tr '/' '_')"; header+="\tmetamaps_${safe_label}_reads"; done
printf '%b\n' "${header}" > "${SUMMARY_TSV}"

unset MALLOC_ARENA_MAX || true

for rep in $(seq 1 "${REPLICATES}"); do
  for spike_n in ${SPIKE_LEVELS}; do
    MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"; mkdir -p "${MIX_DIR}"
    tmp_spikes=(); total_spiked=0; GENOME_INDEX=0
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
    for s in "${tmp_spikes[@]:1}"; do
      mix_args+=(--spike_fastq_gz "${s}")
    done
    mix_args+=(--out_fastq_gz "${combined_spike}")
    python3 "${BUILD_MIXED_FASTQ_PY}" "${mix_args[@]}"
    mix_fastq_gz="${MIX_DIR}/mixed.fastq.gz"; python3 "${BUILD_MIXED_FASTQ_PY}" --background_fastq_gz "${WORK_FASTQ}" --spike_fastq_gz "${combined_spike}" --out_fastq_gz "${mix_fastq_gz}"
    metamaps_prefix="${MIX_DIR}/metamaps"
    metamaps mapDirectly -t "${THREADS}" --all --maxmemory "${MAX_MEMORY_GB}" -r "${METAMAPS_DB_DIR}/DB.fa" -q "${mix_fastq_gz}" -o "${metamaps_prefix}"
    metamaps classify -t "${THREADS}" --mappings "${metamaps_prefix}" --DB "${METAMAPS_DB_DIR}"
    wimp_tsv="${metamaps_prefix}.EM.WIMP"
    require_file "${wimp_tsv}"
    row="${rep}\t${spike_n}\t${total_spiked}\t${mix_fastq_gz}\t${metamaps_prefix}\t${wimp_tsv}"
    wimp_total=''
    for line in "${PANEL_LINES[@]}"; do
      target_label="$(printf '%s' "${line}" | cut -f2)"
      safe_label="$(printf '%s' "${target_label}" | tr ' ' '_' | tr '/' '_')"
      wsum="${MIX_DIR}/metamaps.${safe_label}.summary.tsv"
      python3 "${SUMMARISE_METAMAPS_WIMP_PY}" --wimp_tsv "${wimp_tsv}" --analysis_level definedGenomes --target_label "${target_label}" --out_tsv "${wsum}"
      if [[ -z "${wimp_total}" ]]; then
          wimp_total="$(awk 'NR==2{print $3}' "${wsum}")"
          row+="\t${wimp_total}"
      fi
      wtarget="$(awk 'NR==2{print $4}' "${wsum}")"
      row+="\t${wtarget}"
    done
    printf '%b\n' "${row}" >> "${SUMMARY_TSV}"
  done
done

log_info "Done: ${SUMMARY_TSV}"
