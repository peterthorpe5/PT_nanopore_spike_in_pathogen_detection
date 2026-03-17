#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N SENS_multi_kraken_multi


set -euo pipefail

# run_spikein_multi_genome_equal_v2.sh
#
# Multi-genome ONT spike-in pipeline with equal contribution from each genome.
#
# Updated logic:
#   1. Train NanoSim on reads from the original real FASTQ
#      using a smaller monkey reference
#   2. Simulate one ONT-like read pool per pathogen genome
#   3. Perform host/background depletion once and cache the result
#   4. Spike an equal number of reads from each genome into the depleted background
#   5. Run Kraken2 and minimap2
#   6. Write a TSV summary with per-genome counts
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

MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna.gz}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna}"

MONKEY_DEPLETION_FASTA="${MONKEY_DEPLETION_FASTA:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/M.nemestrina_tonkeana_nigra.fasta}"
PLASMO_MASKED_GZ="${PLASMO_MASKED_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"

DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/host_plus_plasmo_depletion.fasta}"
DEPLETED_FASTQ="${DEPLETED_FASTQ:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/MRC1023_AmM008WB.depleted.fastq.gz}"

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

mkdir -p "${OUT_DIR}"
mkdir -p "$(dirname "${MONKEY_SMALL_FASTA}")"
mkdir -p "$(dirname "${DEPLETION_REF_FASTA}")"
mkdir -p "$(dirname "${DEPLETED_FASTQ}")"

if [[ ! -f "${REAL_FASTQ}" ]]; then
  echo "ERROR: REAL_FASTQ not found: ${REAL_FASTQ}"
  exit 1
fi

if [[ ! -f "${MONKEY_SMALL_GZ}" && ! -f "${MONKEY_SMALL_FASTA}" ]]; then
  echo "ERROR: neither MONKEY_SMALL_GZ nor MONKEY_SMALL_FASTA exists."
  exit 1
fi

if [[ ! -f "${MONKEY_DEPLETION_FASTA}" ]]; then
  echo "ERROR: MONKEY_DEPLETION_FASTA not found: ${MONKEY_DEPLETION_FASTA}"
  exit 1
fi

if [[ ! -f "${PLASMO_MASKED_GZ}" ]]; then
  echo "ERROR: PLASMO_MASKED_GZ not found: ${PLASMO_MASKED_GZ}"
  exit 1
fi

if [[ ! -f "${MINIMAP_DB_GZ}" ]]; then
  echo "ERROR: MINIMAP_DB_GZ not found: ${MINIMAP_DB_GZ}"
  exit 1
fi

if [[ ! -d "${KRAKEN_DB_DIR}" ]]; then
  echo "ERROR: KRAKEN_DB_DIR not found: ${KRAKEN_DB_DIR}"
  exit 1
fi

if [[ ! -f "${UTILS_PY}" ]]; then
  echo "ERROR: utility script not found: ${UTILS_PY}"
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

for f in "${PATHOGEN_FASTA_ARR[@]}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: pathogen FASTA not found: ${f}"
    exit 1
  fi
done

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

###############################################################################
# Report configuration
###############################################################################

echo "[INFO] REAL_FASTQ=${REAL_FASTQ}"
echo "[INFO] MONKEY_SMALL_GZ=${MONKEY_SMALL_GZ}"
echo "[INFO] MONKEY_SMALL_FASTA=${MONKEY_SMALL_FASTA}"
echo "[INFO] MONKEY_DEPLETION_FASTA=${MONKEY_DEPLETION_FASTA}"
echo "[INFO] PLASMO_MASKED_GZ=${PLASMO_MASKED_GZ}"
echo "[INFO] DEPLETION_REF_FASTA=${DEPLETION_REF_FASTA}"
echo "[INFO] DEPLETED_FASTQ=${DEPLETED_FASTQ}"
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
# Prepare smaller monkey FASTA for NanoSim training
###############################################################################

if [[ ! -s "${MONKEY_SMALL_FASTA}" ]]; then
  echo "[INFO] Decompressing smaller monkey FASTA for NanoSim training"
  gzip -cd "${MONKEY_SMALL_GZ}" > "${MONKEY_SMALL_FASTA}"
else
  echo "[INFO] Reusing existing smaller monkey FASTA: ${MONKEY_SMALL_FASTA}"
fi

###############################################################################
# Build combined depletion reference
###############################################################################

if [[ ! -s "${DEPLETION_REF_FASTA}" ]]; then
  echo "[INFO] Building depletion reference: ${DEPLETION_REF_FASTA}"
  cat "${MONKEY_DEPLETION_FASTA}" > "${DEPLETION_REF_FASTA}"
  gzip -cd "${PLASMO_MASKED_GZ}" >> "${DEPLETION_REF_FASTA}"
else
  echo "[INFO] Reusing existing depletion reference: ${DEPLETION_REF_FASTA}"
fi

###############################################################################
# Step A: NanoSim training subsample from original FASTQ
###############################################################################

TRAIN_FASTQ="${OUT_DIR}/train_${TRAIN_READS_N}.fastq.gz"

echo "[INFO] Creating NanoSim training subsample from original FASTQ"
python3 "${UTILS_PY}" sample-fastq \
  --fastq_gz "${REAL_FASTQ}" \
  --n_reads "${TRAIN_READS_N}" \
  --seed 1 \
  --out_fastq_gz "${TRAIN_FASTQ}"

###############################################################################
# Step B: NanoSim training on original reads against smaller monkey reference
###############################################################################

NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"

echo "[INFO] Running NanoSim training"
echo "[INFO] Training reads: ${TRAIN_FASTQ}"
echo "[INFO] Training reference: ${MONKEY_SMALL_FASTA}"

read_analysis.py genome \
  --read "${TRAIN_FASTQ}" \
  --ref_g "${MONKEY_SMALL_FASTA}" \
  --output "${NS_MODEL_PREFIX}" \
  --num_threads "${THREADS}" \
  --fastq

###############################################################################
# Step C: simulate one pool per genome
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
# Step D: host/background depletion once, then reuse
###############################################################################

WORK_FASTQ="${DEPLETED_FASTQ}"

if [[ -s "${WORK_FASTQ}" ]]; then
  echo "[INFO] Reusing existing depleted background: ${WORK_FASTQ}"
elif [[ "${DO_HOST_DEPLETION}" == "true" ]]; then
  echo "[INFO] Running host/background depletion"
  minimap2 -t "${THREADS}" -a -x map-ont "${DEPLETION_REF_FASTA}" "${REAL_FASTQ}" \
    | samtools view -h -T "${DEPLETION_REF_FASTA}" -b -f 4 -@ "${THREADS}" - \
    | samtools fastq -@ "${THREADS}" - \
    | gzip -c > "${WORK_FASTQ}"
else
  echo "[INFO] Host/background depletion disabled; using REAL_FASTQ directly"
  ln -sf "${REAL_FASTQ}" "${WORK_FASTQ}"
fi

###############################################################################
# Step E: prepare minimap reference
###############################################################################

MINIMAP_DB_FASTA="${OUT_DIR}/plas_outgrps_genomes_Hard_MASKED.fasta"

if [[ ! -s "${MINIMAP_DB_FASTA}" ]]; then
  echo "[INFO] Decompressing minimap reference to: ${MINIMAP_DB_FASTA}"
  gzip -cd "${MINIMAP_DB_GZ}" > "${MINIMAP_DB_FASTA}"
else
  echo "[INFO] Reusing existing minimap reference: ${MINIMAP_DB_FASTA}"
fi

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