#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 100G
#$ -adds l_hard h_vmem 100G
#$ -N strict_deplete

set -euo pipefail

log_info() {
    printf '[INFO] %s\n' "$*"
}

log_error() {
    printf '[ERROR] %s\n' "$*" >&2
}

require_file() {
    [[ -f "$1" ]] || {
        log_error "Required file not found: $1"
        exit 1
    }
}

require_exe() {
    command -v "$1" >/dev/null 2>&1 || {
        log_error "Required executable not found on PATH: $1"
        exit 1
    }
}

THREADS="${THREADS:-12}"

WORK_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity"
REAL_FASTQ="${WORK_DIR}/data/../project_back_up_2024/jcs_blood_samples/MRC1023_AmM008WB.fastq.gz"
REAL_FASTQ="/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC1023_AmM008WB.fastq.gz"
DEPLETION_REF="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/host_plus_plasmo_depletion.fasta"

OUT_FASTQ="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/MRC1023_AmM008WB.depleted.fastq.gz"
OUT_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/strict_depletion_$(date +%Y%m%d_%H%M%S)"

mkdir -p "${OUT_DIR}"

SAM_BAM="${OUT_DIR}/all_vs_depletion.bam"
MAPPED_IDS="${OUT_DIR}/mapped_read_ids.txt"
SUMMARY_TSV="${OUT_DIR}/strict_depletion_summary.tsv"

require_exe minimap2
require_exe samtools
require_exe python3
require_file "${REAL_FASTQ}"
require_file "${DEPLETION_REF}"

log_info "Input FASTQ: ${REAL_FASTQ}"
log_info "Depletion reference: ${DEPLETION_REF}"
log_info "Output FASTQ: ${OUT_FASTQ}"
log_info "Work dir: ${OUT_DIR}"

log_info "Aligning reads to depletion reference"
minimap2 \
    -t "${THREADS}" \
    -a \
    -x map-ont \
    -N 20 \
    "${DEPLETION_REF}" \
    "${REAL_FASTQ}" \
    | samtools view \
        -@ "${THREADS}" \
        -T "${DEPLETION_REF}" \
        -b \
        -o "${SAM_BAM}" \
        -

log_info "Collecting all read IDs with any mapped alignment"
samtools view \
    -@ "${THREADS}" \
    -F 4 \
    "${SAM_BAM}" \
    | cut -f 1 \
    | sort -u \
    > "${MAPPED_IDS}"

log_info "Filtering original FASTQ by read name blacklist"


python3 - <<PY
import gzip
from pathlib import Path

real_fastq = Path("${REAL_FASTQ}")
mapped_ids = Path("${MAPPED_IDS}")
out_fastq = Path("${OUT_FASTQ}")
summary_tsv = Path("${SUMMARY_TSV}")

with mapped_ids.open("r", encoding="utf-8") as handle:
    blacklist = {line.strip() for line in handle if line.strip()}

total_reads = 0
kept_reads = 0
removed_reads = 0

with gzip.open(real_fastq, "rt", encoding="utf-8") as infile, \
     gzip.open(out_fastq, "wt", encoding="utf-8") as outfile:
    while True:
        header = infile.readline()
        if not header:
            break
        seq = infile.readline()
        plus = infile.readline()
        qual = infile.readline()

        if not seq or not plus or not qual:
            raise ValueError(f"Truncated FASTQ record in {real_fastq}")

        total_reads += 1
        read_id = header.split()[0].lstrip("@")

        if read_id in blacklist:
            removed_reads += 1
            continue

        outfile.write(header)
        outfile.write(seq)
        outfile.write(plus)
        outfile.write(qual)
        kept_reads += 1

with summary_tsv.open("w", encoding="utf-8") as handle:
    handle.write("metric\tvalue\n")
    handle.write(f"total_reads\t{total_reads}\n")
    handle.write(f"removed_reads\t{removed_reads}\n")
    handle.write(f"kept_reads\t{kept_reads}\n")
    handle.write(f"blacklist_size\t{len(blacklist)}\n")
PY


log_info "Done"
log_info "Depleted FASTQ: ${OUT_FASTQ}"
log_info "Mapped read IDs: ${MAPPED_IDS}"
log_info "Summary: ${SUMMARY_TSV}"
