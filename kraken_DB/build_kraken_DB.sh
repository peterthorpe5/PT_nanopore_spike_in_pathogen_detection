#!/usr/bin/env bash
#$ -cwd
#$ -V
#$ -pe smp 16
#$ -l mfree=300G
#$ -l h_vmem=300G
#$ -N kraken_db_build

set -euo pipefail

###############################################################################
# User settings
###############################################################################

WORKDIR="/home/pthorpe001/data/project_back_up_2024/kracken"
DB="kraken_bact_virus_plasmo_fungal"

THREADS="${NSLOTS:-16}"

GENOMES_DIR="${WORKDIR}/genomes"
ALL_PLAS_DIR="${WORKDIR}/all_plasmodium_genomes"

META_DIR="${WORKDIR}/ncbi_metadata"

GENBANK_SUMMARY="${META_DIR}/assembly_summary_genbank.txt"
REFSEQ_SUMMARY="${META_DIR}/assembly_summary_refseq.txt"

# Optional cap to keep DB smaller (and reduce RAM for classification)
# Leave empty to disable.
MAX_DB_SIZE="${MAX_DB_SIZE:-30G}"

# Manual taxid assignments for genomes without GCA/GCF in filename or not found.
MANUAL_TSV="${WORKDIR}/manual_taxid.tsv"

###############################################################################
# Activate environment
###############################################################################

cd "${WORKDIR}"

# Adjust to your cluster's conda initialisation if needed
if [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
  # Common install location
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "${HOME}/conda/etc/profile.d/conda.sh" ]]; then
  source "${HOME}/conda/etc/profile.d/conda.sh"
fi

conda activate kraken2

###############################################################################
# Prepare directories
###############################################################################

mkdir -p "${DB}"
mkdir -p "${META_DIR}"

echo "[INFO] Working directory: ${WORKDIR}"
echo "[INFO] Kraken DB:         ${DB}"
echo "[INFO] Threads:           ${THREADS}"

###############################################################################
# Download taxonomy (use ftp to avoid rsync issues)
###############################################################################

echo "[INFO] Downloading taxonomy"
kraken2-build --download-taxonomy --db "${DB}" --use-ftp

###############################################################################
# Download standard libraries
###############################################################################

echo "[INFO] Downloading standard libraries"
kraken2-build --download-library bacteria     --db "${DB}" --threads "${THREADS}" --use-ftp
kraken2-build --download-library viral       --db "${DB}" --threads "${THREADS}" --use-ftp
kraken2-build --download-library fungi       --db "${DB}" --threads "${THREADS}" --use-ftp
kraken2-build --download-library protozoa    --db "${DB}" --threads "${THREADS}" --use-ftp
kraken2-build --download-library UniVec_Core --db "${DB}" --threads "${THREADS}" --use-ftp

###############################################################################
# Fetch assembly summary files (for mapping GCA/GCF -> taxid)
###############################################################################

echo "[INFO] Downloading NCBI assembly summaries"
curl -L -A "Mozilla/5.0" -o "${REFSEQ_SUMMARY}" \
  https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt

curl -L -A "Mozilla/5.0" -o "${GENBANK_SUMMARY}" \
  https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt

echo "[INFO] Checking assembly summary headers"
head -n 2 "${GENBANK_SUMMARY}" | sed 's/\t/ /g' | head -n 2
head -n 2 "${REFSEQ_SUMMARY}" | sed 's/\t/ /g' | head -n 2

###############################################################################
# Build taxid maps for custom genomes
###############################################################################

echo "[INFO] Mapping taxids for: ${ALL_PLAS_DIR}"
python make_kraken_taxid_map.py \
  --genomes_dir "${ALL_PLAS_DIR}" \
  --assembly_summary_genbank "${GENBANK_SUMMARY}" \
  --assembly_summary_refseq "${REFSEQ_SUMMARY}" \
  --out_tsv "${WORKDIR}/all_plasmodium_genome_to_taxid.tsv"

echo "[INFO] Mapping taxids for: ${GENOMES_DIR}"
python make_kraken_taxid_map.py \
  --genomes_dir "${GENOMES_DIR}" \
  --assembly_summary_genbank "${GENBANK_SUMMARY}" \
  --assembly_summary_refseq "${REFSEQ_SUMMARY}" \
  --out_tsv "${WORKDIR}/genome_to_taxid.tsv"

###############################################################################
# Extract genome_file + taxid only, and merge
###############################################################################

awk -F $'\t' 'NR>1 && $3!="NA" {print $1 "\t" $3}' "${WORKDIR}/all_plasmodium_genome_to_taxid.tsv" \
  > "${WORKDIR}/all_plasmodium_genome_file_taxid.tsv"

awk -F $'\t' 'NR>1 && $3!="NA" {print $1 "\t" $3}' "${WORKDIR}/genome_to_taxid.tsv" \
  > "${WORKDIR}/genome_file_taxid.tsv"

cat "${WORKDIR}/all_plasmodium_genome_file_taxid.tsv" "${WORKDIR}/genome_file_taxid.tsv" \
  | awk 'NF>=2' \
  | sort -u \
  > "${WORKDIR}/custom_genomes_file_taxid.tsv"

###############################################################################
# Report unmapped genomes (so you can fix the 18 not_found cases cleanly)
###############################################################################

echo "[INFO] Unmapped (all_plasmodium_genomes):"
grep -P "\tnot_found$|\tmissing_assembly$" "${WORKDIR}/all_plasmodium_genome_to_taxid.tsv" || true

echo "[INFO] Unmapped (genomes):"
grep -P "\tnot_found$|\tmissing_assembly$" "${WORKDIR}/genome_to_taxid.tsv" || true

echo "[INFO] Counts:"
echo -n "  mapped custom genomes: " ; wc -l < "${WORKDIR}/custom_genomes_file_taxid.tsv"
echo -n "  not_found total:       " ; cat "${WORKDIR}/all_plasmodium_genome_to_taxid.tsv" "${WORKDIR}/genome_to_taxid.tsv" | grep -c "not_found" || true
echo -n "  missing_assembly total:" ; cat "${WORKDIR}/all_plasmodium_genome_to_taxid.tsv" "${WORKDIR}/genome_to_taxid.tsv" | grep -c "missing_assembly" || true

###############################################################################
# Manual taxids (optional but recommended for the 18 not_found / odd filenames)
###############################################################################

# Create file if it doesn't exist
if [[ ! -f "${MANUAL_TSV}" ]]; then
  cat > "${MANUAL_TSV}" <<'EOF'
genome_file	taxid	species
# Example:
# /home/pthorpe001/data/project_back_up_2024/kracken/genomes/HepatitusA.fasta	12092	Hepatovirus A
EOF
fi

# If user has added manual entries, append them (genome_file + taxid only)
if [[ "$(grep -v '^#' "${MANUAL_TSV}" | grep -c -E '^[^[:space:]]+\t[0-9]+' || true)" -gt 0 ]]; then
  awk -F $'\t' 'NR>1 && $1!~/^#/ && $2~/^[0-9]+$/ {print $1 "\t" $2}' "${MANUAL_TSV}" \
    >> "${WORKDIR}/custom_genomes_file_taxid.tsv"
fi

# De-duplicate again
sort -u "${WORKDIR}/custom_genomes_file_taxid.tsv" -o "${WORKDIR}/custom_genomes_file_taxid.tsv"

###############################################################################
# Add custom genomes to Kraken DB
###############################################################################

echo "[INFO] Adding custom genomes to library with explicit taxid"
while IFS=$'\t' read -r genome_file taxid
do
  if [[ ! -f "${genome_file}" ]]; then
    echo "[WARN] Missing file, skipping: ${genome_file}"
    continue
  fi
  echo "[INFO] add-to-library: ${genome_file} taxid=${taxid}"
  kraken2-build --add-to-library "${genome_file}" --db "${DB}" --taxid "${taxid}"
done < "${WORKDIR}/custom_genomes_file_taxid.tsv"

###############################################################################
# Build DB
###############################################################################

echo "[INFO] Building Kraken DB"
if [[ -n "${MAX_DB_SIZE}" ]]; then
  kraken2-build --threads "${THREADS}" --build --db "${DB}" --max-db-size "${MAX_DB_SIZE}"
else
  kraken2-build --threads "${THREADS}" --build --db "${DB}"
fi

echo "[INFO] Done"