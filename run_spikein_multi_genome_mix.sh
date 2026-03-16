#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N SENS_kraken_multi

set -euo pipefail

# run_spikein_multi_genome_equal.sh
#
# Multi-genome ONT spike-in pipeline with equal contribution from each genome.
#
# This script:
#   1. Builds a depletion/training reference from monkey + plasmo masked FASTA
#   2. Optionally removes reads mapping to that reference from the real sample
#   3. Trains NanoSim once on the cleaned background
#   4. Simulates one ONT-like read pool per pathogen genome
#   5. Spikes an equal number of reads from each genome into the background
#   6. Runs Kraken2 and minimap2
#   7. Writes a TSV summary with per-genome counts
#
# Notes
# -----
# - SPIKE_LEVELS are interpreted as reads per genome.
# - Total reads added = spike_n * number_of_genomes.
# - To add another genome later, add another entry to PATHOGEN_FASTA_ARR
#   and another matching entry to TARGET_LABEL_ARR.

###############################################################################
# User inputs
###############################################################################

REAL_FASTQ="${REAL_FASTQ:-/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC1023_AmM008WB.fastq.gz}"

MONKEY="${MONKEY:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/M.nemestrina_tonkeana_nigra.fasta}"
PLASMO="${PLASMO:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"
HOST_REF_FASTA="${HOST_REF_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/host_plus_plasmo_depletion.fasta}"

MINIMAP_DB_GZ="${MINIMAP_DB_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"
KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"

SCRIPT_DIR="${SCRIPT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection}"
UTILS_PY="${SCRIPT_DIR}/spikein_utils.py"

OUT_DIR="${OUT_DIR:-spikein_multi_genome_equal_out}"

THREADS="${THREADS:-12}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000}"
REPLICATES="${REPLICATES:-3}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
MIN_ALIGN_LEN="${MIN_ALIGN_LEN:-500}"
MIN_MAPQ="${MIN_MAPQ:-15}"

###############################################################################
# Pathogen genomes and labels
#
# Add more genomes here if needed.
# Keep TARGET_LABEL_ARR in the same order as PATHOGEN_FASTA_ARR.
###############################################################################

PATHOGEN_FASTA_ARR=(
  "/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/GCA_014843685.1_ASM1484368v1_genomic.fna"
  "/home/pthorpe001/data/project_back_up_2024/kracken/all_plasmodium_genomes/GCA_001861195.1_383_1_genomic.fna"
  "/home/pthorpe001/data/project_back_up_2024/kracken/all_plasmodium_genomes/GCA_023845545.1_PKCLINC047_genomic.fna"
)

TARGET_LABEL_ARR=(
  "Plasmodium vivax"
  "Plasmodium falciparum"
  "Plasmodium knowlesi"
)

###############################################################################
# Checks
###############################################################################

if [[ -z "${REAL_FASTQ}" ]]; then
  echo "ERROR: REAL_FASTQ is not set."
  exit 1
fi

if [[ "${#PATHOGEN_FASTA_ARR[@]}" -lt 2 ]]; then
  echo "ERROR: define at least two pathogen genomes in PATHOGEN_FASTA_ARR."
  exit 1
fi

if [[ "${#PATHOGEN_FASTA_ARR[@]}" -ne "${#TARGET_LABEL_ARR[@]}" ]]; then
  echo "ERROR: PATHOGEN_FASTA_ARR and TARGET_LABEL_ARR must have the same length."
  exit 1
fi

if [[ "${DO_HOST_DEPLETION}" == "true" && -z "${HOST_REF_FASTA}" ]]; then
  echo "ERROR: DO_HOST_DEPLETION=true but HOST_REF_FASTA is not set."
  exit 1
fi

mkdir -p "${OUT_DIR}"
mkdir -p "$(dirname "${HOST_REF_FASTA}")"

if [[ ! -f "${UTILS_PY}" ]]; then
  echo "ERROR: expected utility script not found: ${UTILS_PY}"
  exit 1
fi

python3 - <<'PY'
import shutil

req = [
    "minimap2",
    "samtools",
    "gzip",
    "kraken2",
    "bamToBed",
    "read_analysis.py",
    "simulator.py",
]
missing = [r for r in req if shutil.which(r) is None]
if missing:
    raise SystemExit(f"ERROR: missing executables on PATH: {missing}")
PY

for f in "${REAL_FASTQ}" "${MONKEY}" "${PLASMO}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: required input file not found: ${f}"
    exit 1
  fi
done

for f in "${PATHOGEN_FASTA_ARR[@]}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: pathogen FASTA not found: ${f}"
    exit 1
  fi
done

###############################################################################
# Build combined depletion/training reference
###############################################################################

echo "[INFO] Building combined depletion reference: ${HOST_REF_FASTA}"
cat "${MONKEY}" > "${HOST_REF_FASTA}"
gzip -cd "${PLASMO}" >> "${HOST_REF_FASTA}"

###############################################################################
# Report configuration
###############################################################################

echo "[INFO] REAL_FASTQ=${REAL_FASTQ}"
echo "[INFO] HOST_REF_FASTA=${HOST_REF_FASTA}"
echo "[INFO] MINIMAP_DB_GZ=${MINIMAP_DB_GZ}"
echo "[INFO] KRAKEN_DB_DIR=${KRAKEN_DB_DIR}"
echo "[INFO] OUT_DIR=${OUT_DIR}"
echo "[INFO] THREADS=${THREADS}"
echo "[INFO] SPIKE_LEVELS=${SPIKE_LEVELS}"
echo "[INFO] REPLICATES=${REPLICATES}"
echo "[INFO] TRAIN_READS_N=${TRAIN_READS_N}"
echo "[INFO] SIM_POOL_N=${SIM_POOL_N}"
echo "[INFO] Number of genomes=${#PATHOGEN_FASTA_ARR[@]}"

for i in "${!PATHOGEN_FASTA_ARR[@]}"; do
  echo "[INFO] Genome $((i + 1)) FASTA: ${PATHOGEN_FASTA_ARR[$i]}"
  echo "[INFO] Genome $((i + 1)) label: ${TARGET_LABEL_ARR[$i]}"
done

###############################################################################
# Step A: optional host/background depletion
###############################################################################

WORK_FASTQ="${OUT_DIR}/background.fastq.gz"

if [[ "${DO_HOST_DEPLETION}" == "true" ]]; then
  echo "[INFO] Host/background depletion enabled"

  minimap2 -t "${THREADS}" -ax map-ont "${HOST_REF_FASTA}" "${REAL_FASTQ}" \
    | samtools view -b -f 4 -@ "${THREADS}" - \
    | samtools fastq -@ "${THREADS}" - \
    | gzip -c > "${WORK_FASTQ}"
else
  echo "[INFO] Host/background depletion disabled; using REAL_FASTQ directly"
  ln -sf "${REAL_FASTQ}" "${WORK_FASTQ}"
fi

###############################################################################
# Step B: subsample for NanoSim training
###############################################################################

TRAIN_FASTQ="${OUT_DIR}/train_${TRAIN_READS_N}.fastq.gz"

echo "[INFO] Creating NanoSim training subsample: ${TRAIN_FASTQ}"

python3 "${UTILS_PY}" sample-fastq \
  --fastq_gz "${WORK_FASTQ}" \
  --n_reads "${TRAIN_READS_N}" \
  --seed 1 \
  --out_fastq_gz "${TRAIN_FASTQ}"

###############################################################################
# Step C: NanoSim training
###############################################################################

NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"

echo "[INFO] Running NanoSim training"
read_analysis.py genome \
  --read "${TRAIN_FASTQ}" \
  --ref_g "${HOST_REF_FASTA}" \
  --output "${NS_MODEL_PREFIX}" \
  --num_threads "${THREADS}" \
  --fastq

###############################################################################
# Step D: simulate one pool per genome
###############################################################################

declare -a SIM_POOL_FASTQ_GZ_ARR=()

for i in "${!PATHOGEN_FASTA_ARR[@]}"; do
  genome_idx=$((i + 1))
  pathogen_fasta="${PATHOGEN_FASTA_ARR[$i]}"
  sim_prefix="${OUT_DIR}/simulated_pathogen_${genome_idx}"
  sim_pool_fastq="${OUT_DIR}/sim_pool_${genome_idx}.fastq"
  sim_pool_fastq_gz="${OUT_DIR}/sim_pool_${genome_idx}.fastq.gz"

  echo "[INFO] Simulating NanoSim pool for genome ${genome_idx}"
  echo "[INFO]   FASTA=${pathogen_fasta}"
  echo "[INFO]   Output prefix=${sim_prefix}"

  simulator.py genome \
    --ref_g "${pathogen_fasta}" \
    --model_prefix "${NS_MODEL_PREFIX}" \
    --output "${sim_prefix}" \
    --number "${SIM_POOL_N}" \
    --dna_type linear \
    --seed "$((42 + genome_idx))" \
    --num_threads "${THREADS}" \
    --fastq

  python3 "${UTILS_PY}" combine-nanosim-fastq \
    --sim_prefix "${sim_prefix}" \
    --out_fastq "${sim_pool_fastq}"

  gzip -c "${sim_pool_fastq}" > "${sim_pool_fastq_gz}"
  SIM_POOL_FASTQ_GZ_ARR+=("${sim_pool_fastq_gz}")
done

###############################################################################
# Step E: prepare minimap reference
###############################################################################

MINIMAP_DB_FASTA="${OUT_DIR}/plas_outgrps_genomes_Hard_MASKED.fasta"

echo "[INFO] Decompressing minimap reference to: ${MINIMAP_DB_FASTA}"
gzip -cd "${MINIMAP_DB_GZ}" > "${MINIMAP_DB_FASTA}"

###############################################################################
# Step F: summary header
###############################################################################

SUMMARY_TSV="${OUT_DIR}/spikein_multi_summary.tsv"

header=(
  "replicate"
  "spike_n_per_genome"
  "n_genomes"
  "total_spike_n"
  "mixed_fastq_gz"
  "kraken_report"
  "kraken_classified_reads"
  "minimap_bam"
)

for i in "${!TARGET_LABEL_ARR[@]}"; do
  genome_idx=$((i + 1))
  header+=("target_label_g${genome_idx}")
  header+=("kraken_target_reads_g${genome_idx}")
  header+=("minimap_target_alignments_g${genome_idx}")
done

{
  IFS=$'\t'
  echo "${header[*]}"
} > "${SUMMARY_TSV}"

echo "[INFO] Summary TSV: ${SUMMARY_TSV}"

###############################################################################
# Step G: spike series
###############################################################################

for rep in $(seq 1 "${REPLICATES}"); do
  for spike_n in ${SPIKE_LEVELS}; do

    MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"
    mkdir -p "${MIX_DIR}"

    MIX_FASTQ_GZ="${MIX_DIR}/mixed.fastq.gz"
    COMBINED_SPIKE_FASTQ="${MIX_DIR}/combined_spike.fastq"
    COMBINED_SPIKE_FASTQ_GZ="${MIX_DIR}/combined_spike.fastq.gz"

    : > "${COMBINED_SPIKE_FASTQ}"

    echo "[INFO] rep=${rep} spike_n_per_genome=${spike_n}"

    ###########################################################################
    # Sample equal reads from each genome pool
    ###########################################################################

    for i in "${!SIM_POOL_FASTQ_GZ_ARR[@]}"; do
      genome_idx=$((i + 1))
      sim_pool_fastq_gz="${SIM_POOL_FASTQ_GZ_ARR[$i]}"
      spike_fastq_gz="${MIX_DIR}/spike_g${genome_idx}.fastq.gz"
      spike_fastq="${MIX_DIR}/spike_g${genome_idx}.fastq"

      spike_seed=$(( rep * 1000000 + genome_idx * 10000 + spike_n ))

      python3 "${UTILS_PY}" sample-fastq \
        --fastq_gz "${sim_pool_fastq_gz}" \
        --n_reads "${spike_n}" \
        --seed "${spike_seed}" \
        --out_fastq_gz "${spike_fastq_gz}"

      gzip -cd "${spike_fastq_gz}" > "${spike_fastq}"
      cat "${spike_fastq}" >> "${COMBINED_SPIKE_FASTQ}"
    done

    gzip -c "${COMBINED_SPIKE_FASTQ}" > "${COMBINED_SPIKE_FASTQ_GZ}"

    ###########################################################################
    # Combine background + all spike reads
    ###########################################################################

    cat "${WORK_FASTQ}" "${COMBINED_SPIKE_FASTQ_GZ}" > "${MIX_FASTQ_GZ}"

    ###########################################################################
    # Kraken2
    ###########################################################################

    KRAKEN_REPORT="${MIX_DIR}/kraken.report.tsv"
    KRAKEN_OUT="${MIX_DIR}/kraken.classifications.tsv"

    kraken2 \
      --db "${KRAKEN_DB_DIR}" \
      --threads "${THREADS}" \
      --report "${KRAKEN_REPORT}" \
      --output "${KRAKEN_OUT}" \
      "${MIX_FASTQ_GZ}" \
      >/dev/null

    ###########################################################################
    # Minimap2
    ###########################################################################

    BAM_OUT="${MIX_DIR}/minimap_sorted.bam"
    BED_OUT="${MIX_DIR}/minimap_sorted.MASKED.bed"

    minimap2 -t "${THREADS}" -ax map-ont "${MINIMAP_DB_FASTA}" "${MIX_FASTQ_GZ}" \
      | samtools view -Sh -q "${MIN_MAPQ}" -e "rlen >= ${MIN_ALIGN_LEN}" - \
      | samtools sort -O bam -@ "${THREADS}" - \
      | tee "${BAM_OUT}" \
      | bamToBed > "${BED_OUT}"

    ###########################################################################
    # Summaries
    ###########################################################################

    TMP_KRAKEN_SUM="${MIX_DIR}/kraken.summary.tmp.tsv"
    FIRST_LABEL="${TARGET_LABEL_ARR[0]}"

    python3 "${UTILS_PY}" summarise-kraken \
      --kraken_report_tsv "${KRAKEN_REPORT}" \
      --target_label "${FIRST_LABEL}" \
      --out_tsv "${TMP_KRAKEN_SUM}"

    K_CLASSIFIED="$(awk 'NR==2{print $2}' "${TMP_KRAKEN_SUM}")"

    row=(
      "${rep}"
      "${spike_n}"
      "${#PATHOGEN_FASTA_ARR[@]}"
      "$(( spike_n * ${#PATHOGEN_FASTA_ARR[@]} ))"
      "${MIX_FASTQ_GZ}"
      "${KRAKEN_REPORT}"
      "${K_CLASSIFIED}"
      "${BAM_OUT}"
    )

    for i in "${!TARGET_LABEL_ARR[@]}"; do
      genome_idx=$((i + 1))
      target_label="${TARGET_LABEL_ARR[$i]}"
      per_label_kraken_tsv="${MIX_DIR}/kraken.summary.g${genome_idx}.tsv"

      python3 "${UTILS_PY}" summarise-kraken \
        --kraken_report_tsv "${KRAKEN_REPORT}" \
        --target_label "${target_label}" \
        --out_tsv "${per_label_kraken_tsv}"

      K_TARGET="$(awk 'NR==2{print $3}' "${per_label_kraken_tsv}")"

      MINIMAP_TARGET_ALN="$(python3 "${UTILS_PY}" count-minimap-target \
        --bed "${BED_OUT}" \
        --target_label "${target_label}")"

      row+=("${target_label}")
      row+=("${K_TARGET}")
      row+=("${MINIMAP_TARGET_ALN}")
    done

    {
      IFS=$'\t'
      echo "${row[*]}"
    } >> "${SUMMARY_TSV}"
  done
done

echo "[INFO] Done. Summary: ${SUMMARY_TSV}"
