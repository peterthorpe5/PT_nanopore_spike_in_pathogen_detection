#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 24
#$ -jc long
#$ -mods l_hard mfree 450G
#$ -adds l_hard h_vmem 450G
#$ -N multi_flye_SENS

set -euo pipefail

# Multi-genome ONT spike-in pipeline with equal contribution per genome.
# Detection is done by Flye metagenome assembly followed by Kraken2 on assembly.

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

KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"

SCRIPT_DIR="${SCRIPT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection}"
UTILS_PY="${SCRIPT_DIR}/spikein_utils.py"

OUT_DIR="${OUT_DIR:-spikein_multi_genome_flye_out}"

THREADS="${THREADS:-24}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000}"
REPLICATES="${REPLICATES:-3}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"

###############################################################################
# Pathogen genomes and labels
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

for f in \
  "${REAL_FASTQ}" \
  "${MONKEY_DEPLETION_FASTA}" \
  "${PLASMO_MASKED_GZ}"
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

if [[ "${#PATHOGEN_FASTA_ARR[@]}" -lt 2 ]]; then
  echo "ERROR: define at least two pathogen genomes."
  exit 1
fi

if [[ "${#PATHOGEN_FASTA_ARR[@]}" -ne "${#TARGET_LABEL_ARR[@]}" ]]; then
  echo "ERROR: PATHOGEN_FASTA_ARR and TARGET_LABEL_ARR must have same length."
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
    "read_analysis.py",
    "simulator.py",
    "flye",
]
missing = [r for r in req if shutil.which(r) is None]
if missing:
    raise SystemExit(f"ERROR: missing executables on PATH: {missing}")
PY

###############################################################################
# Prepare references
###############################################################################

if [[ ! -s "${MONKEY_SMALL_FASTA}" ]]; then
  gzip -cd "${MONKEY_SMALL_GZ}" > "${MONKEY_SMALL_FASTA}"
fi

if [[ ! -s "${DEPLETION_REF_FASTA}" ]]; then
  cat "${MONKEY_DEPLETION_FASTA}" > "${DEPLETION_REF_FASTA}"
  gzip -cd "${PLASMO_MASKED_GZ}" >> "${DEPLETION_REF_FASTA}"
fi

###############################################################################
# NanoSim training from original FASTQ
###############################################################################

TRAIN_FASTQ="${OUT_DIR}/train_${TRAIN_READS_N}.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"

python3 "${UTILS_PY}" sample-fastq \
  --fastq_gz "${REAL_FASTQ}" \
  --n_reads "${TRAIN_READS_N}" \
  --seed 1 \
  --out_fastq_gz "${TRAIN_FASTQ}"

read_analysis.py genome \
  --read "${TRAIN_FASTQ}" \
  --ref_g "${MONKEY_SMALL_FASTA}" \
  --output "${NS_MODEL_PREFIX}" \
  --num_threads "${THREADS}" \
  --fastq

###############################################################################
# Simulate one pool per genome
###############################################################################

declare -a SIM_POOL_FASTQ_GZ_ARR=()

for i in "${!PATHOGEN_FASTA_ARR[@]}"; do
  genome_idx=$((i + 1))
  pathogen_fasta="${PATHOGEN_FASTA_ARR[$i]}"
  sim_prefix="${OUT_DIR}/simulated_pathogen_${genome_idx}"
  sim_pool_fastq="${OUT_DIR}/sim_pool_${genome_idx}.fastq"
  sim_pool_fastq_gz="${OUT_DIR}/sim_pool_${genome_idx}.fastq.gz"


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
    --sim_prefix "${sim_prefix}" \
    --out_fastq "${sim_pool_fastq}"

  gzip -c "${sim_pool_fastq}" > "${sim_pool_fastq_gz}"
  SIM_POOL_FASTQ_GZ_ARR+=("${sim_pool_fastq_gz}")
done

###############################################################################
# Host/background depletion once, then reuse
###############################################################################

WORK_FASTQ="${DEPLETED_FASTQ}"

if [[ -s "${WORK_FASTQ}" ]]; then
  echo "[INFO] Reusing existing depleted background: ${WORK_FASTQ}"
elif [[ "${DO_HOST_DEPLETION}" == "true" ]]; then
  minimap2 -t "${THREADS}" -a -x map-ont "${DEPLETION_REF_FASTA}" "${REAL_FASTQ}" \
    | samtools view -h -T "${DEPLETION_REF_FASTA}" -b -f 4 -@ "${THREADS}" - \
    | samtools fastq -@ "${THREADS}" - \
    | gzip -c > "${WORK_FASTQ}"
else
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
# Helper: assembly stats
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
# Summary header
###############################################################################

SUMMARY_TSV="${OUT_DIR}/spikein_multi_flye_summary.tsv"

header=(
  "replicate"
  "spike_n_per_genome"
  "n_genomes"
  "total_spike_n"
  "mixed_fastq_gz"
  "dedup_fastq"
  "flye_out_dir"
  "assembly_fasta"
  "assembly_n_contigs"
  "assembly_total_bases"
  "kraken_report"
  "kraken_classified_contigs"
)

for i in "${!TARGET_LABEL_ARR[@]}"; do
  genome_idx=$((i + 1))
  header+=("target_label_g${genome_idx}")
  header+=("kraken_target_contigs_g${genome_idx}")
done

{
  IFS=$'\t'
  echo "${header[*]}"
} > "${SUMMARY_TSV}"

###############################################################################
# Spike, deduplicate, assemble, classify assembly
###############################################################################

for rep in $(seq 1 "${REPLICATES}"); do
  for spike_n in ${SPIKE_LEVELS}; do
    MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"
    FLYE_DIR="${MIX_DIR}/flye_meta"
    mkdir -p "${MIX_DIR}"

    MIX_FASTQ_GZ="${MIX_DIR}/mixed.fastq.gz"
    COMBINED_SPIKE_FASTQ="${MIX_DIR}/combined_spike.fastq"
    COMBINED_SPIKE_FASTQ_GZ="${MIX_DIR}/combined_spike.fastq.gz"
    DEDUP_FASTQ="${MIX_DIR}/mixed.dedup.fastq"

    : > "${COMBINED_SPIKE_FASTQ}"

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
    cat "${WORK_FASTQ}" "${COMBINED_SPIKE_FASTQ_GZ}" > "${MIX_FASTQ_GZ}"

    dedup_fastq_gz "${MIX_FASTQ_GZ}" "${DEDUP_FASTQ}"

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

    kraken2 \
      --db "${KRAKEN_DB_DIR}" \
      --threads "${THREADS}" \
      --report "${KRAKEN_REPORT}" \
      --output "${KRAKEN_OUT}" \
      "${ASSEMBLY_FASTA}" \
      >/dev/null

    TMP_KRAKEN_SUM="${MIX_DIR}/assembly.kraken.summary.tmp.tsv"
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
      "${DEDUP_FASTQ}"
      "${FLYE_DIR}"
      "${ASSEMBLY_FASTA}"
      "${ASSEMBLY_N_CONTIGS}"
      "${ASSEMBLY_TOTAL_BASES}"
      "${KRAKEN_REPORT}"
      "${K_CLASSIFIED}"
    )

    for i in "${!TARGET_LABEL_ARR[@]}"; do
      genome_idx=$((i + 1))
      target_label="${TARGET_LABEL_ARR[$i]}"
      per_label_kraken_tsv="${MIX_DIR}/assembly.kraken.summary.g${genome_idx}.tsv"

      python3 "${UTILS_PY}" summarise-kraken \
        --kraken_report_tsv "${KRAKEN_REPORT}" \
        --target_label "${target_label}" \
        --out_tsv "${per_label_kraken_tsv}"

      K_TARGET="$(awk 'NR==2{print $3}' "${per_label_kraken_tsv}")"

      row+=("${target_label}")
      row+=("${K_TARGET}")
    done

    {
      IFS=$'\t'
      echo "${row[*]}"
    } >> "${SUMMARY_TSV}"
  done
done

echo "[INFO] Done. Summary: ${SUMMARY_TSV}"
