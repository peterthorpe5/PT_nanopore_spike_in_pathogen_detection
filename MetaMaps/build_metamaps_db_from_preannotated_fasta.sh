#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 8
#$ -jc long
#$ -mods l_hard mfree 32G
#$ -adds l_hard h_vmem 32G

set -euo pipefail

###############################################################################
# Build a MetaMaps custom database from a single pre-annotated FASTA.
#
# Expected input FASTA headers already contain taxon IDs, for example:
# >contig_name some description, kraken:taxid|508771
#
# Required environment variables (can be overridden with export VAR=...):
#   METAMAPS_DB_DIR   Output MetaMaps database directory
#   INPUT_FASTA_GZ    Input FASTA (gzipped or plain text)
#   TAXONOMY_DIR      Directory containing NCBI taxonomy files for MetaMaps
#   THREADS           Number of threads
#
# Example:
#   export METAMAPS_DB_DIR=/path/to/metamaps_plas_db
#   export INPUT_FASTA_GZ=/path/to/plas_outgrps_genomes_Hard_MASKED.fasta.gz
#   export TAXONOMY_DIR=/path/to/ncbi_taxonomy
#   qsub build_metamaps_db_from_preannotated_fasta.sh
###############################################################################

THREADS="${THREADS:-1}"
INPUT_FASTA_GZ="${INPUT_FASTA_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"
METAMAPS_DB_DIR="${METAMAPS_DB_DIR:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/metamaps_plas_outgroups_db}"
TAXONOMY_DIR="${TAXONOMY_DIR:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal/taxonomy}"

if [[ -z "${TAXONOMY_DIR}" ]]; then
    echo "[ERROR] TAXONOMY_DIR is not set. Please export TAXONOMY_DIR to the NCBI taxonomy directory." >&2
    exit 1
fi

if [[ ! -f "${INPUT_FASTA_GZ}" ]]; then
    echo "[ERROR] Input FASTA not found: ${INPUT_FASTA_GZ}" >&2
    exit 1
fi

if [[ ! -d "${TAXONOMY_DIR}" ]]; then
    echo "[ERROR] Taxonomy directory not found: ${TAXONOMY_DIR}" >&2
    exit 1
fi

mkdir -p "${METAMAPS_DB_DIR}"
mkdir -p "${METAMAPS_DB_DIR}/tmp"

WORK_FASTA="${METAMAPS_DB_DIR}/DB.fa"
LOG_FILE="${METAMAPS_DB_DIR}/build_metamaps_db.log"

exec > >(tee -a "${LOG_FILE}") 2>&1

echo "[INFO] Starting MetaMaps database build"
echo "[INFO] Input FASTA: ${INPUT_FASTA_GZ}"
echo "[INFO] Output DB dir: ${METAMAPS_DB_DIR}"
echo "[INFO] Taxonomy dir: ${TAXONOMY_DIR}"
echo "[INFO] Threads: ${THREADS}"

action_decompress() {
    if [[ "${INPUT_FASTA_GZ}" == *.gz ]]; then
        zcat "${INPUT_FASTA_GZ}" > "${WORK_FASTA}"
    else
        cp "${INPUT_FASTA_GZ}" "${WORK_FASTA}"
    fi
}

echo "[INFO] Copying/decompressing FASTA to ${WORK_FASTA}"
action_decompress

echo "[INFO] Checking that FASTA headers contain kraken taxids"
HEADER_TEST_COUNT=$(grep -c '^>' "${WORK_FASTA}" || true)
HEADER_TAXID_COUNT=$(grep '^>' "${WORK_FASTA}" | grep -c 'kraken:taxid|' || true)

echo "[INFO] Total headers: ${HEADER_TEST_COUNT}"
echo "[INFO] Headers with kraken:taxid|: ${HEADER_TAXID_COUNT}"

if [[ "${HEADER_TEST_COUNT}" -eq 0 ]]; then
    echo "[ERROR] No FASTA headers found in ${WORK_FASTA}" >&2
    exit 1
fi

if [[ "${HEADER_TAXID_COUNT}" -ne "${HEADER_TEST_COUNT}" ]]; then
    echo "[ERROR] Not all FASTA headers contain kraken:taxid| entries" >&2
    echo "[ERROR] Please ensure every sequence header is taxon-annotated before building the DB" >&2
    exit 1
fi

echo "[INFO] Running MetaMaps buildDB.pl"
buildDB.pl \
    --DB "${METAMAPS_DB_DIR}" \
    --FASTAs "${WORK_FASTA}" \
    --taxonomy "${TAXONOMY_DIR}"

echo "[INFO] MetaMaps database build completed"
echo "[INFO] Database FASTA: ${METAMAPS_DB_DIR}/DB.fa"
echo "[INFO] Database taxonomy dir: ${METAMAPS_DB_DIR}/taxonomy"
