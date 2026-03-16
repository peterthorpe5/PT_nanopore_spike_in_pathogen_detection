#!/usr/bin/env bash
set -euo pipefail

# run_spikein_one_sample.sh
#
# One-sample, one-pathogen spike-in pilot:
# - optional host depletion
# - NanoSim training on a subsample
# - simulate ONT-like reads from pathogen genome
# - spike-in series, Kraken2 + minimap2 detection, TSV summary
#
# Requirements:
#   minimap2, samtools, gzip, python3
#   NanoSim: read_analysis.py and simulator.py on PATH in the active conda env
#   seqkit strongly recommended (for sampling); if not present we fall back to python sampling
#   bedtools (bamToBed) for your minimap method output bed
#
# Notes:
#   - Outputs are TSV, not comma-separated.
#   - Designed for ONT fastq.gz.
#   - Uses NanoSim "genome" mode with --fastq to model qualities. See NanoSim docs. :contentReference[oaicite:2]{index=2}

###############################################################################
# User inputs (edit these or pass via env vars)
###############################################################################

# Background FASTQ (real sample). Example:
# REAL_FASTQ="/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC1023_AmM008WB.fastq.gz"
REAL_FASTQ="${REAL_FASTQ:-}"

# Host reference FASTA (for optional host depletion AND NanoSim profiling).
# You can point this to a single FASTA file (recommended).
# Example:
# HOST_REF_FASTA="/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/host.fasta"

MONKEY=/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/M.nemestrina_tonkeana_nigra.fasta
PLASMO=/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz

HOST_REF_FASTA="${HOST_REF_FASTA:-}"

MONKEY="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/M.nemestrina_tonkeana_nigra.fasta"
PLASMO="/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz"

HOST_REF_FASTA="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/host_plus_plasmo_depletion.fasta"

cat "${MONKEY}" > "${HOST_REF_FASTA}"
gzip -cd "${PLASMO}" >> "${HOST_REF_FASTA}"

# Pathogen genome FASTA (spike source)
PATHOGEN_FASTA="${PATHOGEN_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/GCA_014843685.1_ASM1484368v1_genomic.fna}"

# Minimap detection reference (your hard-masked plas/outgroup/bait database)
MINIMAP_DB_GZ="${MINIMAP_DB_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"

# Kraken2 DB directory
KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"

# Output directory for this run
OUT_DIR="${OUT_DIR:-spikein_pilot_out}"

# Threads
THREADS="${THREADS:-12}"

# Spike levels: absolute number of simulated reads added
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"

# Replicates per spike level (different seeds)
REPLICATES="${REPLICATES:-3}"

# Host depletion toggle: "true" or "false"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"

# NanoSim training: number of reads to sample from background for training
TRAIN_READS_N="${TRAIN_READS_N:-200000}"

# NanoSim simulation: size of the simulated pool from pathogen
SIM_POOL_N="${SIM_POOL_N:-20000}"

# Minimap detection thresholds (matching your earlier pipeline)
MIN_ALIGN_LEN="${MIN_ALIGN_LEN:-500}"
MIN_MAPQ="${MIN_MAPQ:-15}"

###############################################################################
# Checks
###############################################################################

if [[ -z "${REAL_FASTQ}" ]]; then
  echo "ERROR: REAL_FASTQ is not set. Export REAL_FASTQ=/path/to/sample.fastq.gz"
  exit 1
fi

if [[ "${DO_HOST_DEPLETION}" == "true" && -z "${HOST_REF_FASTA}" ]]; then
  echo "ERROR: DO_HOST_DEPLETION=true but HOST_REF_FASTA is not set."
  exit 1
fi

mkdir -p "${OUT_DIR}"

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
# Helper paths
###############################################################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UTILS_PY="${SCRIPT_DIR}/spikein_utils.py"

if [[ ! -f "${UTILS_PY}" ]]; then
  echo "ERROR: expected ${UTILS_PY} next to this script"
  exit 1
fi

###############################################################################
# Step A: (Optional) host depletion
###############################################################################

WORK_FASTQ="${OUT_DIR}/background.fastq.gz"

if [[ "${DO_HOST_DEPLETION}" == "true" ]]; then
  echo "[INFO] Host depletion enabled"
  echo "[INFO] Background input: ${REAL_FASTQ}"
  echo "[INFO] Host reference:   ${HOST_REF_FASTA}"

  # Map to host; keep only unmapped reads (-f 4).
  # For ONT, map-ont preset is appropriate.
  # Produce gzipped fastq.
  minimap2 -t "${THREADS}" -ax map-ont "${HOST_REF_FASTA}" "${REAL_FASTQ}" \
    | samtools view -b -f 4 -@ "${THREADS}" - \
    | samtools fastq -@ "${THREADS}" - \
    | gzip -c > "${WORK_FASTQ}"

else
  echo "[INFO] Host depletion disabled; using REAL_FASTQ directly"
  ln -sf "${REAL_FASTQ}" "${WORK_FASTQ}"
fi

###############################################################################
# Step B: subsample reads for NanoSim training
###############################################################################

TRAIN_FASTQ="${OUT_DIR}/train_${TRAIN_READS_N}.fastq.gz"

echo "[INFO] Creating training subsample: ${TRAIN_FASTQ}"

python3 "${UTILS_PY}" sample-fastq \
  --fastq_gz "${WORK_FASTQ}" \
  --n_reads "${TRAIN_READS_N}" \
  --seed 1 \
  --out_fastq_gz "${TRAIN_FASTQ}"

###############################################################################
# Step C: NanoSim characterisation (genome mode, model base qualities)
###############################################################################

NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"

if [[ -n "${HOST_REF_FASTA}" ]]; then
  echo "[INFO] Running NanoSim read_analysis.py genome"
  echo "[INFO] Model prefix: ${NS_MODEL_PREFIX}"

  # NanoSim genome mode; --fastq models base qualities. :contentReference[oaicite:3]{index=3}
  read_analysis.py genome \
    --read "${TRAIN_FASTQ}" \
    --ref_g "${HOST_REF_FASTA}" \
    --output "${NS_MODEL_PREFIX}" \
    --num_threads "${THREADS}" \
    --fastq
else
  echo "ERROR: HOST_REF_FASTA is required for NanoSim genome-mode training."
  exit 1
fi

###############################################################################
# Step D: Simulate a pathogen read pool using NanoSim simulator.py (fastq output)
###############################################################################

SIM_POOL_FASTQ="${OUT_DIR}/sim_pool.fastq"
SIM_POOL_FASTQ_GZ="${OUT_DIR}/sim_pool.fastq.gz"

echo "[INFO] Simulating pathogen pool with NanoSim"
echo "[INFO] Pathogen genome: ${PATHOGEN_FASTA}"
echo "[INFO] Pool size:       ${SIM_POOL_N}"

# NanoSim simulator genome mode; --fastq emits FASTQ. :contentReference[oaicite:4]{index=4}
simulator.py genome \
  --ref_g "${PATHOGEN_FASTA}" \
  --model_prefix "${NS_MODEL_PREFIX}" \
  --output "${OUT_DIR}/simulated_pathogen" \
  --number "${SIM_POOL_N}" \
  --dna_type linear \
  --seed 42 \
  --num_threads "${THREADS}" \
  --fastq

# NanoSim writes multiple files depending on aligned/unaligned categories in some versions.
# We combine any produced fastq into one pool.
python3 "${UTILS_PY}" combine-nanosim-fastq \
  --sim_prefix "${OUT_DIR}/simulated_pathogen" \
  --out_fastq "${SIM_POOL_FASTQ}"

gzip -c "${SIM_POOL_FASTQ}" > "${SIM_POOL_FASTQ_GZ}"

###############################################################################
# Step E: Prepare minimap reference (decompress once)
###############################################################################

MINIMAP_DB_FASTA="${OUT_DIR}/plas_outgrps_genomes_Hard_MASKED.fasta"

echo "[INFO] Decompressing minimap DB once to: ${MINIMAP_DB_FASTA}"
gzip -cd "${MINIMAP_DB_GZ}" > "${MINIMAP_DB_FASTA}"

###############################################################################
# Step F: Run spike series
###############################################################################

SUMMARY_TSV="${OUT_DIR}/spikein_summary.tsv"
echo -e "replicate\tspike_n\tmixed_fastq_gz\tkraken_report\tkraken_classified_reads\tkraken_target_reads\tminimap_bam\tminimap_target_alignments" > "${SUMMARY_TSV}"


TARGET_LABEL="${TARGET_LABEL:-}"

if [[ -z "${TARGET_LABEL}" ]]; then
  TARGET_LABEL="$(python3 "${UTILS_PY}" guess-label --fasta "${PATHOGEN_FASTA}")"
fi

echo "[INFO] Target label guess (used for parsing): ${TARGET_LABEL}"
echo "[INFO] Summary TSV: ${SUMMARY_TSV}"

for rep in $(seq 1 "${REPLICATES}"); do
  for spike_n in ${SPIKE_LEVELS}; do

    MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"
    mkdir -p "${MIX_DIR}"

    MIX_FASTQ_GZ="${MIX_DIR}/mixed.fastq.gz"
    SPIKE_FASTQ_GZ="${MIX_DIR}/spike.fastq.gz"

    echo "[INFO] rep=${rep} spike_n=${spike_n}"

    # Sample spike_n reads from simulated pool
    SPIKE_SEED=$(( rep * 100000 + spike_n ))

    python3 "${UTILS_PY}" sample-fastq \
      --fastq_gz "${SIM_POOL_FASTQ_GZ}" \
      --n_reads "${spike_n}" \
      --seed "${SPIKE_SEED}" \
      --out_fastq_gz "${SPIKE_FASTQ_GZ}"



    # Concatenate real + simulated
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

    # Summarise Kraken2
    KRAKEN_SUM_TSV="${MIX_DIR}/kraken.summary.tsv"
    python3 "${UTILS_PY}" summarise-kraken \
      --kraken_report_tsv "${KRAKEN_REPORT}" \
      --target_label "${TARGET_LABEL}" \
      --out_tsv "${KRAKEN_SUM_TSV}"

    K_CLASSIFIED="$(awk 'NR==2{print $2}' "${KRAKEN_SUM_TSV}")"
    K_TARGET="$(awk 'NR==2{print $3}' "${KRAKEN_SUM_TSV}")"

    ###########################################################################
    # Minimap detection pipeline (your method, parameterised)
    ###########################################################################
    BAM_OUT="${MIX_DIR}/minimap_sorted_nodup.bam"
    BED_OUT="${MIX_DIR}/minimap_sorted_nodup.MASKED.bed"


    minimap2 -t "${THREADS}" -ax map-ont "${MINIMAP_DB_FASTA}" "${MIX_FASTQ_GZ}" \
      | samtools view -Sh -q "${MIN_MAPQ}" -e "rlen >= ${MIN_ALIGN_LEN}" - \
      | samtools sort -O bam -@ "${THREADS}" - \
      | tee "${BAM_OUT}" \
      | bamToBed > "${BED_OUT}"



    # Count alignments that hit the target contigs/species label
    MINIMAP_TARGET_ALN="$(python3 "${UTILS_PY}" count-minimap-target \
      --bed "${BED_OUT}" \
      --target_label "${TARGET_LABEL}")"

    echo -e "${rep}\t${spike_n}\t${MIX_FASTQ_GZ}\t${KRAKEN_REPORT}\t${K_CLASSIFIED}\t${K_TARGET}\t${BAM_OUT}\t${MINIMAP_TARGET_ALN}" \
      >> "${SUMMARY_TSV}"

  done
done

echo "[INFO] Done. Summary: ${SUMMARY_TSV}"
