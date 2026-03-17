#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 300G
#$ -adds l_hard h_vmem 300G
#$ -N SENS_shuffle_ctrl

set -euo pipefail

# run_spikein_shuffled_control.sh
#
# Negative control for the one-sample spike-in pipeline.
#
# Logic:
#   1. Take the original pathogen reference FASTA
#   2. Mononucleotide-shuffle each sequence independently
#   3. Simulate ONT-like reads from the shuffled reference using the same
#      NanoSim model logic as the main pipeline
#   4. Spike those reads into the same depleted background framework
#   5. Run the standard one-sample pipeline unchanged, except PATHOGEN_FASTA
#      points to the shuffled FASTA
#
# Purpose:
#   Test whether the workflow yields spurious Plasmodium detection from
#   composition-matched but biologically meaningless sequence.

###############################################################################
# Inputs aligned to run_spikein_one_sample.sh
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

MINIMAP_DB_GZ="${MINIMAP_DB_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"
KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"

SCRIPT_DIR="${SCRIPT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection}"
MAIN_SCRIPT="${MAIN_SCRIPT:-${SCRIPT_DIR}/run_spikein_one_sample.sh}"

OUT_DIR="${OUT_DIR:-spikein_shuffled_control_out}"
THREADS="${THREADS:-12}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-3}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"

SHUFFLE_SEED="${SHUFFLE_SEED:-123}"
SHUFFLE_MODE="${SHUFFLE_MODE:-mononucleotide}"

###############################################################################
# Checks
###############################################################################

mkdir -p "${OUT_DIR}"

for f in \
    "${REAL_FASTQ}" \
    "${PATHOGEN_FASTA}" \
    "${MAIN_SCRIPT}" \
    "${MONKEY_DEPLETION_FASTA}" \
    "${MINIMAP_DB_GZ}" \
    "${PLASMO_MASKED_GZ}"
do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required input not found: ${f}" >&2
        exit 1
    fi
done

if [[ ! -f "${MONKEY_SMALL_GZ}" && ! -f "${MONKEY_SMALL_FASTA}" ]]; then
    echo "ERROR: neither MONKEY_SMALL_GZ nor MONKEY_SMALL_FASTA exists." >&2
    exit 1
fi

if [[ ! -d "${KRAKEN_DB_DIR}" ]]; then
    echo "ERROR: KRAKEN_DB_DIR not found: ${KRAKEN_DB_DIR}" >&2
    exit 1
fi

python3 - <<'PY'
import shutil

required = ["python3", "gzip", "bash"]
missing = [tool for tool in required if shutil.which(tool) is None]
if missing:
    raise SystemExit(f"ERROR: missing executables on PATH: {missing}")
PY

###############################################################################
# Paths
###############################################################################

SHUFFLED_FASTA="${OUT_DIR}/pathogen.shuffled.fasta"
SHUFFLE_META_TSV="${OUT_DIR}/shuffle_metadata.tsv"
RUN_LOG="${OUT_DIR}/shuffle_control_run.log"

###############################################################################
# Create shuffled pathogen FASTA
###############################################################################

echo "[INFO] Creating shuffled pathogen FASTA: ${SHUFFLED_FASTA}"

python3 - "${PATHOGEN_FASTA}" "${SHUFFLED_FASTA}" "${SHUFFLE_META_TSV}" "${SHUFFLE_SEED}" "${SHUFFLE_MODE}" <<'PY'
import gzip
import random
import sys
from collections import Counter
from pathlib import Path


def open_maybe_gzip(path_str, mode):
    """Open plain-text or gzipped file based on suffix."""
    path = Path(path_str)
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", errors="replace")
    return open(path, mode, encoding="utf-8", errors="replace")


def read_fasta(path_str):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    seq_parts = []
    with open_maybe_gzip(path_str, "rt") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
        if header is not None:
            yield header, "".join(seq_parts)


def wrap_sequence(seq, width=80):
    """Wrap a sequence string to fixed-width FASTA lines."""
    for start in range(0, len(seq), width):
        yield seq[start:start + width]


def shuffle_mono(seq, rng):
    """Shuffle A/C/G/T positions while preserving non-ACGT characters in place."""
    seq_list = list(seq)
    acgt_positions = [
        idx for idx, base in enumerate(seq_list)
        if base in {"A", "C", "G", "T"}
    ]
    acgt_bases = [seq_list[idx] for idx in acgt_positions]
    rng.shuffle(acgt_bases)
    for idx, base in zip(acgt_positions, acgt_bases):
        seq_list[idx] = base
    return "".join(seq_list)


in_fasta = sys.argv[1]
out_fasta = sys.argv[2]
meta_tsv = sys.argv[3]
seed = int(sys.argv[4])
shuffle_mode = sys.argv[5]

if shuffle_mode != "mononucleotide":
    raise SystemExit(f"ERROR: unsupported SHUFFLE_MODE: {shuffle_mode}")

rng = random.Random(seed)

records = list(read_fasta(in_fasta))
if not records:
    raise SystemExit("ERROR: no FASTA records found in PATHOGEN_FASTA")

with open(out_fasta, "wt", encoding="utf-8") as out_handle, \
        open(meta_tsv, "wt", encoding="utf-8") as meta_handle:
    meta_handle.write(
        "record_index\toriginal_header\tshuffled_header\tlength\tA\tC\tG\tT\tN_or_other\n"
    )

    for idx, (header, seq) in enumerate(records, start=1):
        shuffled = shuffle_mono(seq=seq, rng=rng)
        shuffled_header = f"{header}__shuffled"

        out_handle.write(f">{shuffled_header}\n")
        for chunk in wrap_sequence(seq=shuffled, width=80):
            out_handle.write(chunk + "\n")

        counts = Counter(seq)
        n_other = sum(
            count for base, count in counts.items()
            if base not in {"A", "C", "G", "T"}
        )

        meta_handle.write(
            f"{idx}\t{header}\t{shuffled_header}\t{len(seq)}\t"
            f"{counts.get('A', 0)}\t{counts.get('C', 0)}\t"
            f"{counts.get('G', 0)}\t{counts.get('T', 0)}\t{n_other}\n"
        )

print(out_fasta)
PY

###############################################################################
# Record provenance
###############################################################################

{
    echo "parameter\tvalue"
    echo -e "real_fastq\t${REAL_FASTQ}"
    echo -e "original_pathogen_fasta\t${PATHOGEN_FASTA}"
    echo -e "shuffled_fasta\t${SHUFFLED_FASTA}"
    echo -e "target_label\t${TARGET_LABEL}"
    echo -e "shuffle_seed\t${SHUFFLE_SEED}"
    echo -e "shuffle_mode\t${SHUFFLE_MODE}"
    echo -e "out_dir\t${OUT_DIR}"
} > "${OUT_DIR}/run_metadata.tsv"

###############################################################################
# Run the main one-sample pipeline using the shuffled FASTA
###############################################################################

echo "[INFO] Launching main one-sample pipeline with shuffled pathogen FASTA"
echo "[INFO] Output dir: ${OUT_DIR}"
echo "[INFO] Target label still set to: ${TARGET_LABEL}"
echo "[INFO] Log: ${RUN_LOG}"

REAL_FASTQ="${REAL_FASTQ}" \
MONKEY_SMALL_GZ="${MONKEY_SMALL_GZ}" \
MONKEY_SMALL_FASTA="${MONKEY_SMALL_FASTA}" \
MONKEY_DEPLETION_FASTA="${MONKEY_DEPLETION_FASTA}" \
PLASMO_MASKED_GZ="${PLASMO_MASKED_GZ}" \
DEPLETION_REF_FASTA="${DEPLETION_REF_FASTA}" \
DEPLETED_FASTQ="${DEPLETED_FASTQ}" \
PATHOGEN_FASTA="${SHUFFLED_FASTA}" \
TARGET_LABEL="${TARGET_LABEL}" \
MINIMAP_DB_GZ="${MINIMAP_DB_GZ}" \
KRAKEN_DB_DIR="${KRAKEN_DB_DIR}" \
SCRIPT_DIR="${SCRIPT_DIR}" \
OUT_DIR="${OUT_DIR}/main_pipeline_run" \
THREADS="${THREADS}" \
TRAIN_READS_N="${TRAIN_READS_N}" \
SIM_POOL_N="${SIM_POOL_N}" \
SPIKE_LEVELS="${SPIKE_LEVELS}" \
REPLICATES="${REPLICATES}" \
DO_HOST_DEPLETION="${DO_HOST_DEPLETION}" \
bash "${MAIN_SCRIPT}" 2>&1 | tee "${RUN_LOG}"

echo "[INFO] Shuffled control completed successfully."