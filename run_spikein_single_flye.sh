#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N spike_single_flye

set -euo pipefail

log_info() { printf '[INFO] %s\n' "$*"; }
log_error() { printf '[ERROR] %s\n' "$*" >&2; }
require_file() { [[ -f "$1" ]] || { log_error "Required file not found: $1"; exit 1; }; }
require_dir() { [[ -d "$1" ]] || { log_error "Required directory not found: $1"; exit 1; }; }
require_exe() { command -v "$1" >/dev/null 2>&1 || { log_error "Required executable not found on PATH: $1"; exit 1; }; }

if [[ -n "${REPO_DIR:-}" ]]; then REPO_DIR="${REPO_DIR}"; elif [[ -n "${SGE_O_WORKDIR:-}" ]]; then REPO_DIR="${SGE_O_WORKDIR}"; else REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection"; fi
PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR:-${REPO_DIR}/scripts}"
SAMPLE_FASTQ_PY="${PY_SCRIPTS_DIR}/sample_fastq.py"
BUILD_MIXED_FASTQ_PY="${PY_SCRIPTS_DIR}/build_mixed_fastq.py"
COMBINE_NANOSIM_FASTQ_PY="${PY_SCRIPTS_DIR}/combine_nanosim_fastq.py"
SUMMARISE_KRAKEN_PY="${PY_SCRIPTS_DIR}/summarise_kraken_report.py"
DEDUP_FASTQ_NAMES_PY="${PY_SCRIPTS_DIR}/dedup_fastq_names.py"
ASSEMBLY_STATS_PY="${PY_SCRIPTS_DIR}/assembly_stats.py"

REAL_FASTQ="${REAL_FASTQ:-/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC1023_AmM008WB.fastq.gz}"
MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna.gz}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna}"
DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/host_plus_plasmo_depletion.fasta}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/MRC1023_AmM008WB.depleted.fastq.gz}"
PATHOGEN_FASTA="${PATHOGEN_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/GCA_014843685.1_ASM1484368v1_genomic.fna}"
TARGET_LABEL="${TARGET_LABEL:-Plasmodium vivax}"
KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"
OUT_DIR="${OUT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_single_flye_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-12}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-3}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"

for x in python3 read_analysis.py simulator.py flye kraken2 gzip; do require_exe "$x"; done
for f in "${SAMPLE_FASTQ_PY}" "${BUILD_MIXED_FASTQ_PY}" "${COMBINE_NANOSIM_FASTQ_PY}" "${SUMMARISE_KRAKEN_PY}" "${DEDUP_FASTQ_NAMES_PY}" "${ASSEMBLY_STATS_PY}" "${REAL_FASTQ}" "${PATHOGEN_FASTA}"; do require_file "$f"; done
require_dir "${KRAKEN_DB_DIR}"
mkdir -p "${OUT_DIR}"

TRAIN_FASTQ="${OUT_DIR}/train_reads.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"
SIM_PREFIX="${OUT_DIR}/simulated_pathogen"
SIM_POOL_FASTQ="${OUT_DIR}/sim_pool.fastq"
SIM_POOL_FASTQ_GZ="${OUT_DIR}/sim_pool.fastq.gz"
SUMMARY_TSV="${OUT_DIR}/spikein_flye_summary.tsv"
WORK_FASTQ="${DEPLETED_FASTQ}"

python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${REAL_FASTQ}" --n_reads "${TRAIN_READS_N}" --seed 1 --out_fastq_gz "${TRAIN_FASTQ}"
if [[ ! -s "${MONKEY_SMALL_FASTA}" ]]; then gzip -cd "${MONKEY_SMALL_GZ}" > "${MONKEY_SMALL_FASTA}"; fi
read_analysis.py genome --read "${TRAIN_FASTQ}" --ref_g "${MONKEY_SMALL_FASTA}" --output "${NS_MODEL_PREFIX}" --num_threads "${THREADS}" --fastq
simulator.py genome --ref_g "${PATHOGEN_FASTA}" --model_prefix "${NS_MODEL_PREFIX}" --output "${SIM_PREFIX}" --number "${SIM_POOL_N}" -dna_type linear --seed 42 --num_threads "${THREADS}" --fastq
python3 "${COMBINE_NANOSIM_FASTQ_PY}" --sim_prefix "${SIM_PREFIX}" --out_fastq "${SIM_POOL_FASTQ}"
gzip -c "${SIM_POOL_FASTQ}" > "${SIM_POOL_FASTQ_GZ}"

if [[ -s "${WORK_FASTQ}" ]]; then :; elif [[ "${DO_HOST_DEPLETION}" == "true" ]]; then minimap2 -t "${THREADS}" -a -x map-ont "${DEPLETION_REF_FASTA}" "${REAL_FASTQ}" | samtools view -h -T "${DEPLETION_REF_FASTA}" -b -f 4 -@ "${THREADS}" - | samtools fastq -@ "${THREADS}" - | gzip -c > "${WORK_FASTQ}"; else cp "${REAL_FASTQ}" "${WORK_FASTQ}"; fi

printf 'replicate\tspike_n\tmixed_fastq_gz\tflye_out_dir\tassembly_fasta\tcontig_count\ttotal_bases\tkraken_report\tkraken_classified_contigs\tkraken_target_contigs\n' > "${SUMMARY_TSV}"

for rep in $(seq 1 "${REPLICATES}"); do
  for spike_n in ${SPIKE_LEVELS}; do
    MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"; mkdir -p "${MIX_DIR}"
    SPIKE_FASTQ_GZ="${MIX_DIR}/spike.fastq.gz"; MIX_FASTQ_GZ="${MIX_DIR}/mixed.fastq.gz"; DEDUP_FASTQ="${MIX_DIR}/mixed.dedup.fastq"; FLYE_OUT_DIR="${MIX_DIR}/flye_out"; ASSEMBLY_FASTA="${FLYE_OUT_DIR}/assembly.fasta"; ASM_STATS_TSV="${MIX_DIR}/assembly.stats.tsv"; KRAKEN_REPORT="${MIX_DIR}/assembly.kraken.report.tsv"; KRAKEN_OUT="${MIX_DIR}/assembly.kraken.classifications.tsv"; KRAKEN_SUM_TSV="${MIX_DIR}/assembly.kraken.summary.tsv"
    SPIKE_SEED=$(( rep * 100000 + spike_n ))
    python3 "${SAMPLE_FASTQ_PY}" --fastq_gz "${SIM_POOL_FASTQ_GZ}" --n_reads "${spike_n}" --seed "${SPIKE_SEED}" --out_fastq_gz "${SPIKE_FASTQ_GZ}"
    python3 "${BUILD_MIXED_FASTQ_PY}" --background_fastq_gz "${WORK_FASTQ}" --spike_fastq_gz "${SPIKE_FASTQ_GZ}" --out_fastq_gz "${MIX_FASTQ_GZ}"
    python3 "${DEDUP_FASTQ_NAMES_PY}" --input_fastq_gz "${MIX_FASTQ_GZ}" --output_fastq "${DEDUP_FASTQ}"
    flye --meta --nano-hq "${DEDUP_FASTQ}" --out-dir "${FLYE_OUT_DIR}" --threads "${THREADS}"
    python3 "${ASSEMBLY_STATS_PY}" --assembly_fasta "${ASSEMBLY_FASTA}" --out_tsv "${ASM_STATS_TSV}"
    kraken2 --db "${KRAKEN_DB_DIR}" --threads "${THREADS}" --report "${KRAKEN_REPORT}" --output "${KRAKEN_OUT}" "${ASSEMBLY_FASTA}" >/dev/null
    python3 "${SUMMARISE_KRAKEN_PY}" --kraken_report_tsv "${KRAKEN_REPORT}" --target_label "${TARGET_LABEL}" --out_tsv "${KRAKEN_SUM_TSV}"
    CONTIG_COUNT="$(awk 'NR==2{print $2}' "${ASM_STATS_TSV}")"; TOTAL_BASES="$(awk 'NR==2{print $3}' "${ASM_STATS_TSV}")"; K_CLASSIFIED="$(awk 'NR==2{print $2}' "${KRAKEN_SUM_TSV}")"; K_TARGET="$(awk 'NR==2{print $3}' "${KRAKEN_SUM_TSV}")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${rep}" "${spike_n}" "${MIX_FASTQ_GZ}" "${FLYE_OUT_DIR}" "${ASSEMBLY_FASTA}" "${CONTIG_COUNT}" "${TOTAL_BASES}" "${KRAKEN_REPORT}" "${K_CLASSIFIED}" "${K_TARGET}" >> "${SUMMARY_TSV}"
  done
done

log_info "Done: ${SUMMARY_TSV}"
