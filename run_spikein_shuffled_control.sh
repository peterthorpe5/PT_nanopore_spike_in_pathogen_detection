#!/usr/bin/env bash
set -euo pipefail

# run_spikein_shuffled_control.sh
#
# Negative control:
#   shuffle the pathogen genome sequence (mononucleotide shuffle),
#   simulate reads from the shuffled genome with the same NanoSim model,
#   spike into the same background, run Kraken2 + minimap2.
#
# This tests whether you get spurious Plas detection.

REAL_FASTQ="${REAL_FASTQ:-}"
HOST_REF_FASTA="${HOST_REF_FASTA:-}"
PATHOGEN_FASTA="${PATHOGEN_FASTA:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/genomes/GCF_000524495.1_Plas_inui_San_Antonio_1_V1_genomic.fasta}"
MINIMAP_DB_GZ="${MINIMAP_DB_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"
KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-/home/pthorpe001/data/project_back_up_2024/kracken}"

OUT_DIR="${OUT_DIR:-spikein_shuffled_control_out}"
THREADS="${THREADS:-12}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-3}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UTILS_PY="${SCRIPT_DIR}/spikein_utils.py"

if [[ -z "${REAL_FASTQ}" ]]; then
  echo "ERROR: REAL_FASTQ is not set."
  exit 1
fi
if [[ "${DO_HOST_DEPLETION}" == "true" && -z "${HOST_REF_FASTA}" ]]; then
  echo "ERROR: DO_HOST_DEPLETION=true but HOST_REF_FASTA is not set."
  exit 1
fi

mkdir -p "${OUT_DIR}"

SHUFFLED_FASTA="${OUT_DIR}/pathogen.shuffled.fasta"
python3 - <<PY
import random
from pathlib import Path

in_fa = Path("${PATHOGEN_FASTA}")
out_fa = Path("${SHUFFLED_FASTA}")

random.seed(123)

def read_fasta(path):
    name = None
    seqs = []
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seqs)
                name = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line.upper())
        if name is not None:
            yield name, "".join(seqs)

def shuffle_seq(seq):
    bases = list(seq)
    random.shuffle(bases)
    return "".join(bases)

with open(out_fa, "wt", encoding="utf-8") as out:
    for name, seq in read_fasta(in_fa):
        out.write(f">{name}_shuffled\n")
        sh = shuffle_seq(seq)
        for i in range(0, len(sh), 80):
            out.write(sh[i:i+80] + "\n")

print(str(out_fa))
PY

# Now re-use the main pipeline script by overriding PATHOGEN_FASTA
export PATHOGEN_FASTA="${SHUFFLED_FASTA}"
export OUT_DIR="${OUT_DIR}"

# You can call the main script (ensure it's executable and in same folder):
bash "${SCRIPT_DIR}/run_spikein_one_sample.sh"
