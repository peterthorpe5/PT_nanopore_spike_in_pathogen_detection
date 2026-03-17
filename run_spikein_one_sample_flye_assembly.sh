#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 24
#$ -jc long
#$ -mods l_hard mfree 450G
#$ -adds l_hard h_vmem 450G
#$ -N fly_SENS_assembly

set -euo pipefail

# run_spikein_one_sample_flye_assembly.sh
#
# One-sample, one-pathogen ONT spike-in pipeline with assembly-based detection.
#
# Workflow
# --------
# 1. Train NanoSim on the original real FASTQ using a smaller monkey reference
# 2. Simulate ONT-like reads from one target pathogen genome
# 3. Reuse or create a depleted background FASTQ
# 4. Spike simulated reads into the depleted background at multiple levels
# 5. Deduplicate FASTQ read names
# 6. Assemble each mixed dataset with Flye --meta
# 7. Classify the assembly FASTA with Kraken2
#
# Notes
# -----
# - Outputs are TSV, not comma-separated.
# - Each spike level gets its own assembly folder.
# - Kraken2 is run on assembly.fasta rather than the raw reads.

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

PATHOGEN_FASTA="${PATHOGEN_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/GCA_014843685.1_ASM1484368v1_genomic.fna}"
TARGET_LABEL="${TARGET_LABEL:-Plasmodium vivax}"

KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"

SCRIPT_DIR="${SCRIPT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection}"
UTILS_PY="${SCRIPT_DIR}/spikein_utils.py"

OUT_DIR="${OUT_DIR:-spikein_flye_assembly_out}"

THREADS="${THREADS:-24}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000}"
REPLICATES="${REPLICATES:-3}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"

###############################################################################
# Checks
###############################################################################

mkdir -p "${OUT_DIR}"
mkdir -p "$(dirname "${MONKEY_SMALL_FASTA}")"
mkdir -p "$(dirname "${DEPLETION_REF_FASTA}")"
mkdir -p "$(dirname "${DEPLETED_FASTQ}")"

for f in \
  "${REAL_FASTQ}" \
  "${MONKEY_DEPLETION_FASTA}" \
  "${PLASMO_MASKED_GZ}" \
  "${PATHOGEN_FASTA}"
do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: required input not found: ${f}"
    exit 1
  fi
done

if [[ ! -f "${MONKEY_SMALL_GZ}" && ! -f "${MONKEY_SMALL_FASTA}" ]]; then
  echo "ERROR: neither MONKEY_SMALL_GZ nor MONKEY_SMALL_FASTA exists."
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

python3 - <<'PY'
import shutil

req = [
    "minimap2",
    "samtools",
    "gzip",
    "kraken2",
    "read_analysis.py",
    "simulator.py",
    "flye",
]
missing = [r for r in req if shutil.which(r) is None]
if missing:
    raise SystemExit(f"ERROR: missing executables on PATH: {missing}")
PY

###############################################################################
# Report configuration
###############################################################################

echo "[INFO] REAL_FASTQ=${REAL_FASTQ}"
echo "[INFO] PATHOGEN_FASTA=${PATHOGEN_FASTA}"
echo "[INFO] TARGET_LABEL=${TARGET_LABEL}"
echo "[INFO] OUT_DIR=${OUT_DIR}"
echo "[INFO] DEPLETED_FASTQ=${DEPLETED_FASTQ}"
echo "[INFO] THREADS=${THREADS}"
echo "[INFO] SPIKE_LEVELS=${SPIKE_LEVELS}"
echo "[INFO] REPLICATES=${REPLICATES}"

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
# Step B: NanoSim training
###############################################################################

NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"

echo "[INFO] Running NanoSim training"
read_analysis.py genome \
  --read "${TRAIN_FASTQ}" \
  --ref_g "${MONKEY_SMALL_FASTA}" \
  --output "${NS_MODEL_PREFIX}" \
  --num_threads "${THREADS}" \
  --fastq

###############################################################################
# Step C: Simulate pathogen read pool
###############################################################################

SIM_POOL_FASTQ="${OUT_DIR}/sim_pool.fastq"
SIM_POOL_FASTQ_GZ="${OUT_DIR}/sim_pool.fastq.gz"

echo "[INFO] Simulating pathogen pool with NanoSim"
simulator.py genome \
  --ref_g "${PATHOGEN_FASTA}" \
  --model_prefix "${NS_MODEL_PREFIX}" \
  --output "${OUT_DIR}/simulated_pathogen" \
  --number "${SIM_POOL_N}" \
  --dna_type linear \
  --seed 42 \
  --num_threads "${THREADS}" \
  --fastq

python3 "${UTILS_PY}" combine-nanosim-fastq \
  --sim_prefix "${OUT_DIR}/simulated_pathogen" \
  --out_fastq "${SIM_POOL_FASTQ}"

gzip -c "${SIM_POOL_FASTQ}" > "${SIM_POOL_FASTQ_GZ}"

###############################################################################
# Step D: Host/background depletion once, then reuse
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
# Helper: deduplicate FASTQ read names
###############################################################################

dedup_fastq_gz() {
  local input_fastq_gz="$1"
  local output_fastq="$2"

  python3 - "$input_fastq_gz" "$output_fastq" <<'PY'
import gzip
import sys

input_fastq_gz = sys.argv[1]
output_fastq = sys.argv[2]

seen = {}
record_n = 0

with gzip.open(input_fastq_gz, "rt", encoding="utf-8", errors="replace") as infh, \
        open(output_fastq, "wt", encoding="utf-8") as outfh:
    while True:
        h = infh.readline()
        if not h:
            break
        s = infh.readline()
        p = infh.readline()
        q = infh.readline()

        if not q:
            raise SystemExit("ERROR: truncated FASTQ record encountered")

        record_n += 1
        header = h.rstrip("\n")
        seq = s.rstrip("\n")
        plus = p.rstrip("\n")
        qual = q.rstrip("\n")

        if not header.startswith("@"):
            raise SystemExit(f"ERROR: invalid FASTQ header at record {record_n}: {header}")

        core = header[1:].split()[0]
        extra = header[len(core) + 1:] if len(header) > len(core) + 1 else ""

        count = seen.get(core, 0) + 1
        seen[core] = count

        if count == 1:
            new_core = core
        else:
            new_core = f"{core}_dup{count}"

        if extra:
            new_header = f"@{new_core}{extra}"
        else:
            new_header = f"@{new_core}"

        outfh.write(new_header + "\n")
        outfh.write(seq + "\n")
        outfh.write(plus + "\n")
        outfh.write(qual + "\n")
PY
}

###############################################################################
# Helper: count assembly contigs and bases
###############################################################################

assembly_stats() {
  local assembly_fasta="$1"

  python3 - "$assembly_fasta" <<'PY'
import sys

assembly_fasta = sys.argv[1]
n_contigs = 0
n_bases = 0

with open(assembly_fasta, "rt", encoding="utf-8", errors="replace") as fh:
    seq_len = 0
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if seq_len > 0:
                n_contigs += 1
                n_bases += seq_len
            seq_len = 0
        else:
            seq_len += len(line)

    if seq_len > 0:
        n_contigs += 1
        n_bases += seq_len

print(f"{n_contigs}\t{n_bases}")
PY
}

###############################################################################
# Step E: spike, deduplicate, assemble, classify assembly
###############################################################################

SUMMARY_TSV="${OUT_DIR}/spikein_flye_summary.tsv"
echo -e "replicate\tspike_n\tmixed_fastq_gz\tdedup_fastq\tflye_out_dir\tassembly_fasta\tassembly_n_contigs\tassembly_total_bases\tkraken_report\tkraken_classified_contigs\tkraken_target_contigs" > "${SUMMARY_TSV}"

for rep in $(seq 1 "${REPLICATES}"); do
  for spike_n in ${SPIKE_LEVELS}; do

    MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"
    FLYE_DIR="${MIX_DIR}/flye_meta"
    mkdir -p "${MIX_DIR}"

    MIX_FASTQ_GZ="${MIX_DIR}/mixed.fastq.gz"
    SPIKE_FASTQ_GZ="${MIX_DIR}/spike.fastq.gz"
    DEDUP_FASTQ="${MIX_DIR}/mixed.dedup.fastq"

    echo "[INFO] rep=${rep} spike_n=${spike_n}"

    SPIKE_SEED=$(( rep * 100000 + spike_n ))

    python3 "${UTILS_PY}" sample-fastq \
      --fastq_gz "${SIM_POOL_FASTQ_GZ}" \
      --n_reads "${spike_n}" \
      --seed "${SPIKE_SEED}" \
      --out_fastq_gz "${SPIKE_FASTQ_GZ}"

    cat "${WORK_FASTQ}" "${SPIKE_FASTQ_GZ}" > "${MIX_FASTQ_GZ}"

    echo "[INFO] Deduplicating FASTQ read names"
    dedup_fastq_gz "${MIX_FASTQ_GZ}" "${DEDUP_FASTQ}"

    echo "[INFO] Running Flye metagenome assembly"
    flye \
      --meta \
      --nano-hq "${DEDUP_FASTQ}" \
      --out-dir "${FLYE_DIR}" \
      --threads "${THREADS}"

    ASSEMBLY_FASTA="${FLYE_DIR}/assembly.fasta"

    if [[ ! -s "${ASSEMBLY_FASTA}" ]]; then
      echo "ERROR: Flye assembly not found or empty: ${ASSEMBLY_FASTA}"
      exit 1
    fi

    read -r ASSEMBLY_N_CONTIGS ASSEMBLY_TOTAL_BASES < <(assembly_stats "${ASSEMBLY_FASTA}")

    KRAKEN_REPORT="${MIX_DIR}/assembly.kraken.report.tsv"
    KRAKEN_OUT="${MIX_DIR}/assembly.kraken.classifications.tsv"
    KRAKEN_SUM_TSV="${MIX_DIR}/assembly.kraken.summary.tsv"

    echo "[INFO] Running Kraken2 on assembly FASTA"
    kraken2 \
      --db "${KRAKEN_DB_DIR}" \
      --threads "${THREADS}" \
      --report "${KRAKEN_REPORT}" \
      --output "${KRAKEN_OUT}" \
      "${ASSEMBLY_FASTA}" \
      >/dev/null

    python3 "${UTILS_PY}" summarise-kraken \
      --kraken_report_tsv "${KRAKEN_REPORT}" \
      --target_label "${TARGET_LABEL}" \
      --out_tsv "${KRAKEN_SUM_TSV}"

    K_CLASSIFIED="$(awk 'NR==2{print $2}' "${KRAKEN_SUM_TSV}")"
    K_TARGET="$(awk 'NR==2{print $3}' "${KRAKEN_SUM_TSV}")"

    echo -e "${rep}\t${spike_n}\t${MIX_FASTQ_GZ}\t${DEDUP_FASTQ}\t${FLYE_DIR}\t${ASSEMBLY_FASTA}\t${ASSEMBLY_N_CONTIGS}\t${ASSEMBLY_TOTAL_BASES}\t${KRAKEN_REPORT}\t${K_CLASSIFIED}\t${K_TARGET}" \
      >> "${SUMMARY_TSV}"
  done
done

echo "[INFO] Done. Summary: ${SUMMARY_TSV}"
