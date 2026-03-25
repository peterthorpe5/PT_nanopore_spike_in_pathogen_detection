#!/usr/bin/env bash
#SBATCH -J kraken  #jobname
#SBATCH --partition=himem 
#SBATCH --cpus-per-task=16
#SBATCH --mem=700GB
set -e
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

WORKDIR="/home/pthorpe/scratch/kracken"
DB="kraken_bact_virus_plasmo_fungal"

THREADS="${NSLOTS:-16}"

GENOMES_DIR="${WORKDIR}/genomes"
ALL_PLAS_DIR="${WORKDIR}/all_plasmodium_genomes"

META_DIR="${WORKDIR}/ncbi_metadata"

GENBANK_SUMMARY="${META_DIR}/assembly_summary_genbank.txt"
REFSEQ_SUMMARY="${META_DIR}/assembly_summary_refseq.txt"

# Optional cap to keep DB smaller (and reduce RAM for classification)
# Leave empty to disable.
#MAX_DB_SIZE="${MAX_DB_SIZE:-80G}"

# Manual taxid assignments for genomes without GCA/GCF in filename or not found.
MANUAL_TSV="${WORKDIR}/manual_taxid.tsv"

###############################################################################
# Activate environment
###############################################################################

cd "${WORKDIR}"



#conda activate kraken2

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
kraken2-build --download-taxonomy --db "${DB}"  --use-ftp

###############################################################################
# Download standard libraries
###############################################################################

echo "[INFO] Downloading standard libraries"
kraken2-build --download-library bacteria    --db "${DB}" --threads "${THREADS}" --use-ftp
kraken2-build --download-library viral       --db "${DB}" --threads "${THREADS}" --use-ftp
kraken2-build --download-library fungi       --db "${DB}" --threads "${THREADS}" --use-ftp
kraken2-build --download-library protozoa    --db "${DB}" --threads "${THREADS}" --use-ftp
kraken2-build --download-library UniVec_Core --db "${DB}" --threads "${THREADS}"  --use-ftp



###############################################################################
# Fetch assembly summary files (for mapping GCA/GCF -> taxid)
###############################################################################



DB="kraken_bact_virus_plasmo_fungal"

for f in */*.fna.gz
do
  kraken2-build --add-to-library "${f}" --db "${DB}"
done




###############################################################################
# Build taxid maps for custom genomes
###############################################################################



##############################################################################
# Extract genome_file + taxid only, and merge
###############################################################################


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
genomes/HepatitusA.fasta	12092	Hepatovirus A
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