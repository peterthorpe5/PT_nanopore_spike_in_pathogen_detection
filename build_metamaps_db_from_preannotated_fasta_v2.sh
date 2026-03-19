#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -jc long
#$ -mods l_hard mfree 32G
#$ -adds l_hard h_vmem 32G

set -euo pipefail

THREADS="${THREADS:-1}"
INPUT_FASTA_GZ="${INPUT_FASTA_GZ:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz}"
METAMAPS_DB_DIR="${METAMAPS_DB_DIR:-/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/metamaps_plas_outgroups_db}"
TAXONOMY_DIR="${TAXONOMY_DIR:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal/taxonomy}"

WORK_FASTA="${METAMAPS_DB_DIR}/DB.fa"
LOG_FILE="${METAMAPS_DB_DIR}/build_metamaps_db.log"

mkdir -p "${METAMAPS_DB_DIR}"

exec > >(tee -a "${LOG_FILE}") 2>&1

echo "[INFO] Starting MetaMaps database build"
echo "[INFO] INPUT_FASTA_GZ=${INPUT_FASTA_GZ}"
echo "[INFO] METAMAPS_DB_DIR=${METAMAPS_DB_DIR}"
echo "[INFO] TAXONOMY_DIR=${TAXONOMY_DIR}"
echo "[INFO] WORK_FASTA=${WORK_FASTA}"
echo "[INFO] THREADS=${THREADS}"

if [[ ! -f "${INPUT_FASTA_GZ}" ]]; then
    echo "[ERROR] Input FASTA not found: ${INPUT_FASTA_GZ}" >&2
    exit 1
fi

if [[ ! -d "${TAXONOMY_DIR}" ]]; then
    echo "[ERROR] Taxonomy directory not found: ${TAXONOMY_DIR}" >&2
    exit 1
fi

echo "[INFO] Decompressing FASTA to ${WORK_FASTA}"
if [[ "${INPUT_FASTA_GZ}" == *.gz ]]; then
    zcat "${INPUT_FASTA_GZ}" > "${WORK_FASTA}"
else
    cp "${INPUT_FASTA_GZ}" "${WORK_FASTA}"
fi

if [[ ! -s "${WORK_FASTA}" ]]; then
    echo "[ERROR] WORK_FASTA was not created correctly: ${WORK_FASTA}" >&2
    exit 1
fi

echo "[INFO] Checking FASTA headers"
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
    exit 1
fi

echo "[INFO] Running buildDB.pl"
buildDB.pl \
    --DB "${METAMAPS_DB_DIR}" \
    --FASTAs "${WORK_FASTA}" \
    --taxonomy "${TAXONOMY_DIR}"

echo "[INFO] MetaMaps database build completed"
echo "[INFO] Database FASTA: ${METAMAPS_DB_DIR}/DB.fa"
echo "[INFO] Database taxonomy dir: ${METAMAPS_DB_DIR}/taxonomy"
