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
# Self-contained shuffled negative-control ONT spike-in pipeline.
#
# Workflow:
#   1. Sample training reads from the original real FASTQ
#   2. Train NanoSim on those reads against the smaller monkey reference
#   3. Create a mononucleotide-shuffled pathogen FASTA
#   4. Simulate ONT-like reads from the shuffled FASTA
#   5. Reuse or create depleted background FASTQ
#   6. Spike shuffled reads into background across a series of spike levels
#   7. Run Kraken2 and minimap2 for each mixture
#   8. Write a TSV summary
#
# Notes:
#   - This script does not call any other shell script.
#   - It uses spikein_utils.py for helper functions already used elsewhere.
#   - It builds mixed FASTQ safely by decompressing and recompressing once.
#   - It writes short shuffled contig names to reduce the risk of oversized
#     downstream read headers.

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

MINIMAP_DB_GZ="${MINIMAP_DB_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"
KRAKEN_DB_DIR="${KRAKEN_DB_DIR:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"

SCRIPT_DIR="${SCRIPT_DIR:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection}"
UTILS_PY="${SCRIPT_DIR}/spikein_utils.py"

OUT_DIR="${OUT_DIR:-spikein_shuffled_control_out}"
THREADS="${THREADS:-12}"
TRAIN_READS_N="${TRAIN_READS_N:-200000}"
SIM_POOL_N="${SIM_POOL_N:-20000}"
SPIKE_LEVELS="${SPIKE_LEVELS:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES="${REPLICATES:-3}"
DO_HOST_DEPLETION="${DO_HOST_DEPLETION:-true}"

MIN_MAPQ="${MIN_MAPQ:-15}"
MIN_ALIGN_LEN="${MIN_ALIGN_LEN:-500}"

SHUFFLE_SEED="${SHUFFLE_SEED:-123}"
SHUFFLE_MODE="${SHUFFLE_MODE:-mononucleotide}"

###############################################################################
# Helper functions
###############################################################################

log_info() {
    printf '[INFO] %s\n' "$*"
}

log_warn() {
    printf '[WARN] %s\n' "$*" >&2
}

log_error() {
    printf '[ERROR] %s\n' "$*" >&2
}

require_file() {
    local file_path="$1"
    if [[ ! -f "${file_path}" ]]; then
        log_error "Required file not found: ${file_path}"
        exit 1
    fi
}

require_dir() {
    local dir_path="$1"
    if [[ ! -d "${dir_path}" ]]; then
        log_error "Required directory not found: ${dir_path}"
        exit 1
    fi
}

require_exe() {
    local exe_name="$1"
    if ! command -v "${exe_name}" >/dev/null 2>&1; then
        log_error "Required executable not found on PATH: ${exe_name}"
        exit 1
    fi
}

gzip_test_if_gz() {
    local file_path="$1"
    if [[ "${file_path}" == *.gz ]]; then
        log_info "Testing gzip integrity: ${file_path}"
        gzip -t "${file_path}"
    fi
}

make_mixed_fastq_gz() {
    local background_fastq_gz="$1"
    local spike_fastq_gz="$2"
    local out_fastq_gz="$3"
    local tmp_fastq="${out_fastq_gz%.gz}"

    gzip -t "${background_fastq_gz}"
    gzip -t "${spike_fastq_gz}"

    gzip -cd "${background_fastq_gz}" > "${tmp_fastq}"
    gzip -cd "${spike_fastq_gz}" >> "${tmp_fastq}"
    gzip -c "${tmp_fastq}" > "${out_fastq_gz}"
    gzip -t "${out_fastq_gz}"
    rm -f "${tmp_fastq}"
}

###############################################################################
# Checks
###############################################################################

mkdir -p "${OUT_DIR}"

require_file "${REAL_FASTQ}"
require_file "${PATHOGEN_FASTA}"
require_file "${MONKEY_DEPLETION_FASTA}"
require_file "${PLASMO_MASKED_GZ}"
require_file "${MINIMAP_DB_GZ}"
require_file "${UTILS_PY}"
require_dir "${KRAKEN_DB_DIR}"

if [[ ! -f "${MONKEY_SMALL_GZ}" && ! -f "${MONKEY_SMALL_FASTA}" ]]; then
    log_error "Neither MONKEY_SMALL_GZ nor MONKEY_SMALL_FASTA exists."
    exit 1
fi

require_exe python3
require_exe gzip
require_exe minimap2
require_exe samtools
require_exe kraken2
require_exe bamToBed
require_exe read_analysis.py
require_exe simulator.py

gzip_test_if_gz "${REAL_FASTQ}"
gzip_test_if_gz "${PATHOGEN_FASTA}"
gzip_test_if_gz "${PLASMO_MASKED_GZ}"
gzip_test_if_gz "${MINIMAP_DB_GZ}"
if [[ -f "${MONKEY_SMALL_GZ}" ]]; then
    gzip_test_if_gz "${MONKEY_SMALL_GZ}"
fi

###############################################################################
# Paths
###############################################################################

TRAIN_FASTQ="${OUT_DIR}/train_reads.fastq.gz"
NS_MODEL_PREFIX="${OUT_DIR}/nanosim_training"

SHUFFLED_FASTA="${OUT_DIR}/pathogen.shuffled.fasta"
SHUFFLE_META_TSV="${OUT_DIR}/shuffle_metadata.tsv"

SIM_PREFIX="${OUT_DIR}/simulated_pathogen"
SIM_POOL_FASTQ="${OUT_DIR}/sim_pool.fastq"
SIM_POOL_FASTQ_GZ="${OUT_DIR}/sim_pool.fastq.gz"

MINIMAP_DB_FASTA="${OUT_DIR}/plas_outgrps_genomes_Hard_MASKED.fasta"

SUMMARY_TSV="${OUT_DIR}/spikein_summary.tsv"
RUN_META_TSV="${OUT_DIR}/run_metadata.tsv"

WORK_FASTQ="${DEPLETED_FASTQ}"

###############################################################################
# Step A: create NanoSim training subsample
###############################################################################

log_info "Creating NanoSim training subsample from original FASTQ"
python3 "${UTILS_PY}" sample-fastq \
    --fastq_gz "${REAL_FASTQ}" \
    --n_reads "${TRAIN_READS_N}" \
    --seed 1 \
    --out_fastq_gz "${TRAIN_FASTQ}"

gzip -t "${TRAIN_FASTQ}"

###############################################################################
# Step B: prepare smaller monkey reference if needed
###############################################################################

if [[ ! -s "${MONKEY_SMALL_FASTA}" ]]; then
    if [[ -s "${MONKEY_SMALL_GZ}" ]]; then
        log_info "Decompressing smaller monkey reference"
        gzip -cd "${MONKEY_SMALL_GZ}" > "${MONKEY_SMALL_FASTA}"
    else
        log_error "Smaller monkey FASTA not found and no gzipped source available"
        exit 1
    fi
else
    log_info "Reusing existing smaller monkey FASTA: ${MONKEY_SMALL_FASTA}"
fi

###############################################################################
# Step C: NanoSim training
###############################################################################

log_info "Running NanoSim read_analysis.py genome"
log_info "Training reads: ${TRAIN_FASTQ}"
log_info "Training reference: ${MONKEY_SMALL_FASTA}"

read_analysis.py genome \
    --read "${TRAIN_FASTQ}" \
    --ref_g "${MONKEY_SMALL_FASTA}" \
    --output "${NS_MODEL_PREFIX}" \
    --num_threads "${THREADS}" \
    --fastq

###############################################################################
# Step D: create shuffled FASTA
###############################################################################

log_info "Creating shuffled pathogen FASTA"

python3 - "${PATHOGEN_FASTA}" "${SHUFFLED_FASTA}" "${SHUFFLE_META_TSV}" "${SHUFFLE_SEED}" "${SHUFFLE_MODE}" <<'PY'
"""Create a mononucleotide-shuffled FASTA with short headers."""

import gzip
import random
import sys
from collections import Counter
from pathlib import Path


def open_maybe_gzip(path_str, mode):
    """Open plain-text or gzipped file depending on suffix."""
    path = Path(path_str)
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", errors="replace")
    return open(path, mode, encoding="utf-8", errors="replace")


def read_fasta(path_str):
    """Yield FASTA records as (header, sequence) tuples."""
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
    """Yield fixed-width sequence chunks."""
    for start in range(0, len(seq), width):
        yield seq[start:start + width]


def shuffle_mono(seq, rng):
    """Shuffle A/C/G/T positions while preserving non-ACGT characters."""
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


input_fasta = sys.argv[1]
output_fasta = sys.argv[2]
metadata_tsv = sys.argv[3]
seed = int(sys.argv[4])
shuffle_mode = sys.argv[5]

if shuffle_mode != "mononucleotide":
    raise SystemExit(f"Unsupported SHUFFLE_MODE: {shuffle_mode}")

rng = random.Random(seed)
records = list(read_fasta(input_fasta))

if not records:
    raise SystemExit("No FASTA records found in PATHOGEN_FASTA.")

with open(output_fasta, "wt", encoding="utf-8") as out_handle, \
        open(metadata_tsv, "wt", encoding="utf-8") as meta_handle:
    meta_handle.write(
        "record_index\toriginal_header\tshuffled_header\tlength\t"
        "A\tC\tG\tT\tN_or_other\n"
    )

    for idx, (header, seq) in enumerate(records, start=1):
        shuffled_seq = shuffle_mono(seq=seq, rng=rng)
        shuffled_header = f"shuffle_contig_{idx}"

        out_handle.write(f">{shuffled_header}\n")
        for chunk in wrap_sequence(seq=shuffled_seq, width=80):
            out_handle.write(chunk + "\n")

        counts = Counter(seq)
        n_other = sum(
            value for key, value in counts.items()
            if key not in {"A", "C", "G", "T"}
        )

        meta_handle.write(
            f"{idx}\t{header}\t{shuffled_header}\t{len(seq)}\t"
            f"{counts.get('A', 0)}\t{counts.get('C', 0)}\t"
            f"{counts.get('G', 0)}\t{counts.get('T', 0)}\t{n_other}\n"
        )
PY

require_file "${SHUFFLED_FASTA}"
require_file "${SHUFFLE_META_TSV}"

###############################################################################
# Step E: simulate shuffled pathogen read pool
###############################################################################

log_info "Simulating shuffled pathogen pool with NanoSim"
log_info "Shuffled genome: ${SHUFFLED_FASTA}"
log_info "Pool size: ${SIM_POOL_N}"

simulator.py genome \
    --ref_g "${SHUFFLED_FASTA}" \
    --model_prefix "${NS_MODEL_PREFIX}" \
    --output "${SIM_PREFIX}" \
    --number "${SIM_POOL_N}" \
    --dna_type linear \
    --seed 42 \
    --num_threads "${THREADS}" \
    --fastq

python3 "${UTILS_PY}" combine-nanosim-fastq \
    --sim_prefix "${SIM_PREFIX}" \
    --out_fastq "${SIM_POOL_FASTQ}"

gzip -c "${SIM_POOL_FASTQ}" > "${SIM_POOL_FASTQ_GZ}"
gzip -t "${SIM_POOL_FASTQ_GZ}"

###############################################################################
# Step F: host/background depletion once, then reuse
###############################################################################

if [[ -s "${WORK_FASTQ}" ]]; then
    log_info "Reusing existing depleted background: ${WORK_FASTQ}"
elif [[ "${DO_HOST_DEPLETION}" == "true" ]]; then
    log_info "Running host/background depletion"
    minimap2 -t "${THREADS}" -a -x map-ont "${DEPLETION_REF_FASTA}" "${REAL_FASTQ}" \
        | samtools view -h -T "${DEPLETION_REF_FASTA}" -b -f 4 -@ "${THREADS}" - \
        | samtools fastq -@ "${THREADS}" - \
        | gzip -c > "${WORK_FASTQ}"
    gzip -t "${WORK_FASTQ}"
else
    log_info "Host depletion disabled; copying REAL_FASTQ"
    cp "${REAL_FASTQ}" "${WORK_FASTQ}"
    gzip -t "${WORK_FASTQ}"
fi

###############################################################################
# Step G: prepare minimap reference
###############################################################################

if [[ ! -s "${MINIMAP_DB_FASTA}" ]]; then
    log_info "Decompressing minimap reference to: ${MINIMAP_DB_FASTA}"
    gzip -cd "${MINIMAP_DB_GZ}" > "${MINIMAP_DB_FASTA}"
else
    log_info "Reusing existing minimap reference: ${MINIMAP_DB_FASTA}"
fi

###############################################################################
# Step H: write provenance
###############################################################################

{
    printf 'parameter\tvalue\n'
    printf 'real_fastq\t%s\n' "${REAL_FASTQ}"
    printf 'pathogen_fasta\t%s\n' "${PATHOGEN_FASTA}"
    printf 'shuffled_fasta\t%s\n' "${SHUFFLED_FASTA}"
    printf 'target_label\t%s\n' "${TARGET_LABEL}"
    printf 'threads\t%s\n' "${THREADS}"
    printf 'train_reads_n\t%s\n' "${TRAIN_READS_N}"
    printf 'sim_pool_n\t%s\n' "${SIM_POOL_N}"
    printf 'spike_levels\t%s\n' "${SPIKE_LEVELS}"
    printf 'replicates\t%s\n' "${REPLICATES}"
    printf 'min_mapq\t%s\n' "${MIN_MAPQ}"
    printf 'min_align_len\t%s\n' "${MIN_ALIGN_LEN}"
    printf 'do_host_depletion\t%s\n' "${DO_HOST_DEPLETION}"
    printf 'shuffle_seed\t%s\n' "${SHUFFLE_SEED}"
    printf 'shuffle_mode\t%s\n' "${SHUFFLE_MODE}"
} > "${RUN_META_TSV}"

###############################################################################
# Step I: spike series
###############################################################################

printf 'replicate\tspike_n\tmixed_fastq_gz\tkraken_report\tkraken_classified_reads\tkraken_target_reads\tminimap_bam\tminimap_target_alignments\n' \
    > "${SUMMARY_TSV}"

log_info "Summary TSV: ${SUMMARY_TSV}"

for rep in $(seq 1 "${REPLICATES}"); do
    for spike_n in ${SPIKE_LEVELS}; do

        MIX_DIR="${OUT_DIR}/mix_rep${rep}_n${spike_n}"
        mkdir -p "${MIX_DIR}"

        SPIKE_FASTQ_GZ="${MIX_DIR}/spike.fastq.gz"
        MIX_FASTQ_GZ="${MIX_DIR}/mixed.fastq.gz"

        KRAKEN_REPORT="${MIX_DIR}/kraken.report.tsv"
        KRAKEN_OUT="${MIX_DIR}/kraken.classifications.tsv"
        KRAKEN_SUM_TSV="${MIX_DIR}/kraken.summary.tsv"

        BAM_OUT="${MIX_DIR}/minimap_sorted.bam"
        BED_OUT="${MIX_DIR}/minimap_sorted.MASKED.bed"

        log_info "rep=${rep} spike_n=${spike_n}"

        SPIKE_SEED=$(( rep * 100000 + spike_n ))

        python3 "${UTILS_PY}" sample-fastq \
            --fastq_gz "${SIM_POOL_FASTQ_GZ}" \
            --n_reads "${spike_n}" \
            --seed "${SPIKE_SEED}" \
            --out_fastq_gz "${SPIKE_FASTQ_GZ}"

        gzip -t "${SPIKE_FASTQ_GZ}"

        make_mixed_fastq_gz \
            "${WORK_FASTQ}" \
            "${SPIKE_FASTQ_GZ}" \
            "${MIX_FASTQ_GZ}"

        #######################################################################
        # Kraken2
        #######################################################################

        kraken2 \
            --db "${KRAKEN_DB_DIR}" \
            --threads "${THREADS}" \
            --report "${KRAKEN_REPORT}" \
            --output "${KRAKEN_OUT}" \
            "${MIX_FASTQ_GZ}" \
            >/dev/null

        python3 "${UTILS_PY}" summarise-kraken \
            --kraken_report_tsv "${KRAKEN_REPORT}" \
            --target_label "${TARGET_LABEL}" \
            --out_tsv "${KRAKEN_SUM_TSV}"

        K_CLASSIFIED="$(awk 'NR==2{print $2}' "${KRAKEN_SUM_TSV}")"
        K_TARGET="$(awk 'NR==2{print $3}' "${KRAKEN_SUM_TSV}")"

        #######################################################################
        # Minimap2
        #######################################################################

        minimap2 -t "${THREADS}" -ax map-ont "${MINIMAP_DB_FASTA}" "${MIX_FASTQ_GZ}" \
            | samtools view -Sh -q "${MIN_MAPQ}" -e "rlen >= ${MIN_ALIGN_LEN}" - \
            | samtools sort -O bam -@ "${THREADS}" - \
            | tee "${BAM_OUT}" \
            | bamToBed > "${BED_OUT}"

        MINIMAP_TARGET_ALN="$(python3 "${UTILS_PY}" count-minimap-target \
            --bed "${BED_OUT}" \
            --target_label "${TARGET_LABEL}")"

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${rep}" \
            "${spike_n}" \
            "${MIX_FASTQ_GZ}" \
            "${KRAKEN_REPORT}" \
            "${K_CLASSIFIED}" \
            "${K_TARGET}" \
            "${BAM_OUT}" \
            "${MINIMAP_TARGET_ALN}" \
            >> "${SUMMARY_TSV}"
    done
done

log_info "Done. Summary: ${SUMMARY_TSV}"