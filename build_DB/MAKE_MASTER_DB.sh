#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 8
#$ -jc long
#$ -mods l_hard mfree 32G
#$ -adds l_hard h_vmem 32G

set -euo pipefail
shopt -s nullglob

###############################################################################
# Configuration
###############################################################################

LIBRARY_DIR="/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal/library"
KRAKEN_DB_DIR="/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal"
OUT_DIR="/home/pthorpe001/data/project_back_up_2024/shared_pathogen_reference"
METAMAPS_DIR="${HOME}/apps/MetaMaps"
THREADS="8"

USE_EXISTING_KRAKEN_LIBRARY="yes"

MASTER_FASTA="${OUT_DIR}/master_pathogen_reference.fna"
MASTER_MINIMAP_INDEX="${OUT_DIR}/master_pathogen_reference.mmi"

METAMAPS_WORK_DIR="${OUT_DIR}/metamaps_build"
METAMAPS_SPLIT_DIR="${METAMAPS_WORK_DIR}/split_records"
METAMAPS_INPUT_LIST="${METAMAPS_WORK_DIR}/metamaps_input_filelist.tsv"
METAMAPS_COMBINED_FASTA="${METAMAPS_WORK_DIR}/MetaMaps_combined.fa"
METAMAPS_TAXONOMY_OUT="${METAMAPS_WORK_DIR}/taxonomy"

LOG_FILE="${OUT_DIR}/build_master_reference.log"

###############################################################################
# Setup
###############################################################################

mkdir -p "${OUT_DIR}"
mkdir -p "${METAMAPS_WORK_DIR}"
mkdir -p "${METAMAPS_SPLIT_DIR}"

exec > >(tee "${LOG_FILE}") 2>&1

echo "[INFO] Started: $(date)"
echo "[INFO] LIBRARY_DIR=${LIBRARY_DIR}"
echo "[INFO] OUT_DIR=${OUT_DIR}"

###############################################################################
# Step 1: Build or reuse master FASTA
###############################################################################

if [[ "${USE_EXISTING_KRAKEN_LIBRARY}" == "yes" ]]; then
    if [[ ! -s "${LIBRARY_DIR}/library.fna" ]]; then
        echo "[ERROR] Expected existing Kraken library not found: ${LIBRARY_DIR}/library.fna"
        exit 1
    fi

    echo "[INFO] Reusing existing Kraken combined library:"
    echo "[INFO]   ${LIBRARY_DIR}/library.fna"

    cp -f "${LIBRARY_DIR}/library.fna" "${MASTER_FASTA}"
else
    echo "[INFO] Rebuilding master FASTA from components"

    SOURCE_FASTAS=()

    if [[ -s "${LIBRARY_DIR}/bacteria/library.fna" ]]; then
        SOURCE_FASTAS+=("${LIBRARY_DIR}/bacteria/library.fna")
    fi

    if [[ -s "${LIBRARY_DIR}/viral/library.fna" ]]; then
        SOURCE_FASTAS+=("${LIBRARY_DIR}/viral/library.fna")
    fi

    for FASTA in "${LIBRARY_DIR}"/added/*.fna; do
        [[ "${FASTA}" == *.masked ]] && continue
        [[ "$(basename "${FASTA}")" == prelim_map* ]] && continue
        [[ -s "${FASTA}" ]] || continue
        SOURCE_FASTAS+=("${FASTA}")
    done

    if [[ "${#SOURCE_FASTAS[@]}" -eq 0 ]]; then
        echo "[ERROR] No source FASTA files were found."
        exit 1
    fi

    : > "${MASTER_FASTA}"
    for FASTA in "${SOURCE_FASTAS[@]}"; do
        echo "[INFO] Adding ${FASTA}"
        cat "${FASTA}" >> "${MASTER_FASTA}"
    done
fi

echo "[INFO] Master FASTA written to: ${MASTER_FASTA}"
echo "[INFO] Master FASTA size:"
ls -lh "${MASTER_FASTA}"

###############################################################################
# Step 2: Build minimap2 index
###############################################################################

echo "[INFO] Building minimap2 index"
minimap2 -d "${MASTER_MINIMAP_INDEX}" "${MASTER_FASTA}"

echo "[INFO] Minimap2 index written to: ${MASTER_MINIMAP_INDEX}"

###############################################################################
# Step 3: Split master FASTA into one file per record for MetaMaps
###############################################################################

echo "[INFO] Cleaning previous MetaMaps split files"
rm -f "${METAMAPS_SPLIT_DIR}"/*.fna
rm -f "${METAMAPS_INPUT_LIST}"

echo "[INFO] Splitting master FASTA into per-record FASTA files"

awk \
    -v outdir="${METAMAPS_SPLIT_DIR}" \
    -v filelist="${METAMAPS_INPUT_LIST}" \
    '
    BEGIN {
        RS=">";
        FS="\n";
        OFS="\t";
    }
    NR > 1 {
        header = $1;
        taxid = "";

        if (match(header, /kraken:taxid\|[0-9]+/)) {
            token = substr(header, RSTART, RLENGTH);
            gsub(/kraken:taxid\|/, "", token);
            taxid = token;
        } else {
            next;
        }

        record_count++;
        outfile = sprintf("%s/record_%09d_taxid_%s.fna", outdir, record_count, taxid);

        print ">" header > outfile;
        for (i = 2; i <= NF; i++) {
            print $i > outfile;
        }
        close(outfile);

        print taxid, outfile >> filelist;
        close(filelist);
    }
    ' "${MASTER_FASTA}"

if [[ ! -s "${METAMAPS_INPUT_LIST}" ]]; then
    echo "[ERROR] MetaMaps input list was not created."
    exit 1
fi

echo "[INFO] MetaMaps input list written to: ${METAMAPS_INPUT_LIST}"
echo "[INFO] Number of split FASTA records:"
wc -l "${METAMAPS_INPUT_LIST}"

###############################################################################
# Step 4: Prepare MetaMaps-compatible combined FASTA and taxonomy
###############################################################################

if [[ ! -d "${KRAKEN_DB_DIR}/taxonomy" ]]; then
    echo "[ERROR] Kraken taxonomy directory not found: ${KRAKEN_DB_DIR}/taxonomy"
    exit 1
fi

if [[ ! -f "${METAMAPS_DIR}/combineAndAnnotateReferences.pl" ]]; then
    echo "[ERROR] MetaMaps helper script not found: ${METAMAPS_DIR}/combineAndAnnotateReferences.pl"
    exit 1
fi

rm -f "${METAMAPS_COMBINED_FASTA}"
rm -rf "${METAMAPS_TAXONOMY_OUT}"

echo "[INFO] Running MetaMaps combineAndAnnotateReferences.pl"
perl "${METAMAPS_DIR}/combineAndAnnotateReferences.pl" \
    --inputFileList "${METAMAPS_INPUT_LIST}" \
    --outputFile "${METAMAPS_COMBINED_FASTA}" \
    --taxonomyInDirectory "${KRAKEN_DB_DIR}/taxonomy" \
    --taxonomyOutDirectory "${METAMAPS_TAXONOMY_OUT}"

echo "[INFO] MetaMaps combined FASTA written to: ${METAMAPS_COMBINED_FASTA}"
echo "[INFO] MetaMaps taxonomy written to: ${METAMAPS_TAXONOMY_OUT}"

###############################################################################
# Step 5: Summary
###############################################################################

echo "[INFO] Finished: $(date)"
echo "[INFO] Outputs:"
echo "[INFO]   Master FASTA: ${MASTER_FASTA}"
echo "[INFO]   Minimap2 index: ${MASTER_MINIMAP_INDEX}"
echo "[INFO]   MetaMaps file list: ${METAMAPS_INPUT_LIST}"
echo "[INFO]   MetaMaps combined FASTA: ${METAMAPS_COMBINED_FASTA}"
echo "[INFO]   MetaMaps taxonomy: ${METAMAPS_TAXONOMY_OUT}"
