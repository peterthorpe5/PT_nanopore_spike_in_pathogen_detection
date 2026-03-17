#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G

#$ -N SENS_kraken_test_v2

set -euo pipefail

# run_spikein_one_sample_v2.sh
#
# One-sample, one-pathogen ONT spike-in pipeline.
#
# Updated logic:
#   1. Train NanoSim on reads from the original real FASTQ
#      using a smaller monkey reference
#   2. Simulate ONT-like reads from one target pathogen genome
#   3. Perform host/background depletion once and cache the result
#   4. Spike simulated pathogen reads into the depleted background
#   5. Run Kraken2 and minimap2
#   6. Write a TSV summary
#
# Notes
# -----
# - NanoSim training is done on the original FASTQ, not the depleted FASTQ.
# - Host/background depletion is done once and reused if the output exists.
# - Outputs are TSV, not comma-separated.

###############################################################################
# User inputs
###############################################################################

REAL_FASTQ="${REAL_FASTQ:-/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC1023_AmM008WB.fastq.gz}"

MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna.gz}"
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna}"

MONKEY_DEPLETION_FASTA="${MONKEY_DEPLETION_FASTA:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/M.nemestrina_tonkeana_nigra.fasta}"
PLASMO_MASKED_GZ="${PLASMO_MASKED_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"

DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/host_plus_plasmo_depletion.fasta}"

PATHOGEN_FASTA="${PATHOGEN_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/GCA_014843685.1_ASM1484368v1_genomic.fna}"
TARGET_LABEL="${TARGET_LABEL:-Plasmodium vivax}"

MINIMAP_DB_GZ="${MINIMAP_DB_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"
KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"

SCRIPT_DIR="${SCRIPT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection}"
UTILS_PY="${SCRIPT_DIR}/spikein_utils.py"

#OUT_DIR="${OUT_DIR:-spikein_pilot_out}"

OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/runs/spikein_pilot_out_$(date +%Y%m%d_%H%M%S)"


# Cached depleted background. This lets later runs reuse the expensive depletion.
DEPLETED_FASTQ="${DEPLETED_FASTQ:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/MRC1023_AmM008WB.depleted.fastq.gz}"

THREADS="${THREADS:-12}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-3}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
MIN_ALIGN_LEN="${MIN_ALIGN_LEN:-500}"
MIN_MAPQ="${MIN_MAPQ:-15}"

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

if [[ ! -f "${PATHOGEN_FASTA}" ]]; then
  echo "ERROR: PATHOGEN_FASTA not found: ${PATHOGEN_FASTA}"
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

###############################################################################
# Report configuration
###############################################################################

echo "[INFO] REAL_FASTQ=${REAL_FASTQ}"
echo "[INFO] MONKEY_SMALL_GZ=${MONKEY_SMALL_GZ}"
echo "[INFO] MONKEY_SMALL_FASTA=${MONKEY_SMALL_FASTA}"
echo "[INFO] MONKEY_DEPLETION_FASTA=${MONKEY_DEPLETION_FASTA}"
echo "[INFO] PLASMO_MASKED_GZ=${PLASMO_MASKED_GZ}"
echo "[INFO] DEPLETION_REF_FASTA=${DEPLETION_REF_FASTA}"
echo "[INFO] PATHOGEN_FASTA=${PATHOGEN_FASTA}"
echo "[INFO] TARGET_LABEL=${TARGET_LABEL}"
echo "[INFO] DEPLETED_FASTQ=${DEPLETED_FASTQ}"
echo "[INFO] OUT_DIR=${OUT_DIR}"

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
# Step A: NanoSim training subsample from original real FASTQ
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

echo "[INFO] Running NanoSim read_analysis.py genome"
echo "[INFO] Training reads: ${TRAIN_FASTQ}"
echo "[INFO] Training reference: ${MONKEY_SMALL_FASTA}"

read_analysis.py genome \
  --read "${TRAIN_FASTQ}" \
  --ref_g "${MONKEY_SMALL_FASTA}" \
  --output "${NS_MODEL_PREFIX}" \
  --num_threads "${THREADS}" \
  --fastq

###############################################################################
# Step C: simulate pathogen read pool
###############################################################################

SIM_POOL_FASTQ="${OUT_DIR}/sim_pool.fastq"
SIM_POOL_FASTQ_GZ="${OUT_DIR}/sim_pool.fastq.gz"

echo "[INFO] Simulating pathogen pool with NanoSim"
echo "[INFO] Pathogen genome: ${PATHOGEN_FASTA}"
echo "[INFO] Pool size: ${SIM_POOL_N}"

simulator.py genome \
  --ref_g "${PATHOGEN_FASTA}" \
  --model_prefix "${NS_MODEL_PREFIX}" \
  --output "${OUT_DIR}/simulated_pathogen" \
  --number "${SIM_POOL_N}" \
  -dna_type linear \
  --seed 42 \
  --num_threads "${THREADS}" \
  --fastq

python3 "${UTILS_PY}" combine-nanosim-fastq \
  --sim_prefix "${OUT_DIR}/simulated_pathogen" \
  --out_fastq "${SIM_POOL_FASTQ}"

gzip -c "${SIM_POOL_FASTQ}" > "${SIM_POOL_FASTQ_GZ}"

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
  echo "[INFO] Host depletion disabled; using REAL_FASTQ directly"
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
# Step F: spike series
###############################################################################

SUMMARY_TSV="${OUT_DIR}/spikein_summary.tsv"
echo -e "replicate\tspike_n\tmixed_fastq_gz\tkraken_report\tkraken_classified_reads\tkraken_target_reads\tminimap_bam\tminimap_target_alignments" > "${SUMMARY_TSV}"

echo "[INFO] Summary TSV: ${SUMMARY_TSV}"

for rep in $(seq 1 "${REPLICATES}"); do
  for spike_n in ${SPIKE_LEVELS}; do

    MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"
    mkdir -p "${MIX_DIR}"

    MIX_FASTQ_GZ="${MIX_DIR}/mixed.fastq.gz"
    SPIKE_FASTQ_GZ="${MIX_DIR}/spike.fastq.gz"

    echo "[INFO] rep=${rep} spike_n=${spike_n}"

    SPIKE_SEED=$(( rep * 100000 + spike_n ))

    python3 "${UTILS_PY}" sample-fastq \
      --fastq_gz "${SIM_POOL_FASTQ_GZ}" \
      --n_reads "${spike_n}" \
      --seed "${SPIKE_SEED}" \
      --out_fastq_gz "${SPIKE_FASTQ_GZ}"

    cat "${WORK_FASTQ}" "${SPIKE_FASTQ_GZ}" > "${MIX_FASTQ_GZ}"

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

    KRAKEN_SUM_TSV="${MIX_DIR}/kraken.summary.tsv"
    python3 "${UTILS_PY}" summarise-kraken \
      --kraken_report_tsv "${KRAKEN_REPORT}" \
      --target_label "${TARGET_LABEL}" \
      --out_tsv "${KRAKEN_SUM_TSV}"

    K_CLASSIFIED="$(awk 'NR==2{print $2}' "${KRAKEN_SUM_TSV}")"
    K_TARGET="$(awk 'NR==2{print $3}' "${KRAKEN_SUM_TSV}")"

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

    MINIMAP_TARGET_ALN="$(python3 "${UTILS_PY}" count-minimap-target \
      --bed "${BED_OUT}" \
      --target_label "${TARGET_LABEL}")"

    echo -e "${rep}\t${spike_n}\t${MIX_FASTQ_GZ}\t${KRAKEN_REPORT}\t${K_CLASSIFIED}\t${K_TARGET}\t${BAM_OUT}\t${MINIMAP_TARGET_ALN}" \
      >> "${SUMMARY_TSV}"
  done
done

echo "[INFO] Done. Summary: ${SUMMARY_TSV}"