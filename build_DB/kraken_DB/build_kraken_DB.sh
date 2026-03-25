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
KRAKEN2_DB_NAME="kraken_bact_virus_plasmo_fungal"

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

# Adjust to your cluster's conda initialisation if needed
#if [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
#  # Common install location
#  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
#elif [[ -f "${HOME}/conda/etc/profile.d/conda.sh" ]]; then
#  source "${HOME}/conda/etc/profile.d/conda.sh"
#fi

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


#rsync -avP kraken_bact_virus_plasmo_fungal/ \
#   pthorpe001@login.compute.dundee.ac.uk:/home/pthorpe001/data/project_back_up_2024/kracken/https_genomes/

###############################################################################
# Fetch assembly summary files (for mapping GCA/GCF -> taxid)
###############################################################################

#echo "[INFO] Downloading NCBI assembly summaries"
#curl -L -A "Mozilla/5.0" -o "${REFSEQ_SUMMARY}" \
#  https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt

#curl -L -A "Mozilla/5.0" -o "${GENBANK_SUMMARY}" \
#  https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt

  # get viral genomes from RefSeq (since many are missing from GenBank)
  awk -F $'\t' '
    BEGIN{OFS="\t"}
    $0 ~ /^#/ {next}
    $25 == "viral" && $20 != "na" {print $20}
  ' ncbi_metadata/assembly_summary_refseq.txt > viral_ftp_paths.txt


DB="kraken_bact_virus_plasmo_fungal"

for f in viral_genomes/*.fna.gz
do
  kraken2-build --add-to-library "${f}" --db "${DB}"
done


# for bacteria
awk -F $'\t' '
  $0 ~ /^#/ {next}
  $25 == "bacteria" && $20 != "na" {print $20}
' ncbi_metadata/assembly_summary_refseq.txt > bacteria_ftp_paths.txt

wc -l bacteria_ftp_paths.txt
head -n 5 bacteria_ftp_paths.txt

awk -F $'\t' '
  $0 ~ /^#/ {next}
  $25 == "bacteria" && $12 == "Complete Genome" && $20 != "na" {print $20}
' ncbi_metadata/assembly_summary_refseq.txt > bacteria_complete_ftp_paths.txt

wc -l bacteria_complete_ftp_paths.txt
head -n 5 bacteria_complete_ftp_paths.txt


mkdir -p bacteria_genomes

cat bacteria_complete_ftp_paths.txt \
  | awk '{base=$0; sub(".*/","",base); print $0 "/" base "_genomic.fna.gz\tbacteria_genomes/" base "_genomic.fna.gz"}' \
  | xargs -P 6 -n 1 bash -c '
      url=$(echo "$0" | cut -f1)
      out=$(echo "$0" | cut -f2)
      if [[ -s "$out" ]]; then exit 0; fi
      echo "[INFO] $url"
      curl -L --retry 5 --retry-delay 2 --connect-timeout 20 -o "${out}.tmp" "$url" && mv "${out}.tmp" "$out"
    '

DB="kraken_bact_virus_plasmo_fungal"

for f in bacteria_genomes/*.fna.gz
do
  kraken2-build --add-to-library "${f}" --db "${DB}"
done




awk -F $'\t' '
  $0 ~ /^#/ {next}
  $25 == "protozoa" && $12 == "Complete Genome" && $8 ~ /^Plasmodium / && $20 != "na" {print $20}
' ncbi_metadata/assembly_summary_refseq.txt > plasmodium_complete_paths.txt

wc -l plasmodium_complete_paths.txt
head -n 10 plasmodium_complete_paths.txt

awk -F $'\t' '
  $0 ~ /^#/ {next}
  $25 == "protozoa" && $8 ~ /^Plasmodium / && $20 != "na" {print $20}
' ncbi_metadata/assembly_summary_refseq.txt > plasmodium_paths.txt


mkdir -p plasmodium_genomes

while IFS= read -r base_url
do
  base=$(basename "${base_url}")
  url="${base_url}/${base}_genomic.fna.gz"
  out="plasmodium_genomes/${base}_genomic.fna.gz"

  if [[ -s "${out}" ]]; then
    continue
  fi

  echo "[INFO] ${url}"
  curl -L --retry 5 --retry-delay 2 --connect-timeout 20 \
    -o "${out}.tmp" "${url}" && mv "${out}.tmp" "${out}"
done < plasmodium_complete_paths.txt

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
/home/pthorpe001/data/project_back_up_2024/kracken/genomes/HepatitusA.fasta	12092	Hepatovirus A
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

# De-duplicate exact duplicate lines first
sort -u "${WORKDIR}/custom_genomes_file_taxid.tsv" \
  -o "${WORKDIR}/custom_genomes_file_taxid.tsv"

# Keep a full copy
cp "${WORKDIR}/custom_genomes_file_taxid.tsv" \
  "${WORKDIR}/custom_genomes_file_taxid_all.tsv"

# Keep one representative per taxid
awk -F $'\t' '!seen[$2]++' "${WORKDIR}/custom_genomes_file_taxid.tsv" \
  > "${WORKDIR}/custom_genomes_file_taxid_representatives.tsv"

# Replace original with representative-only version
mv "${WORKDIR}/custom_genomes_file_taxid_representatives.tsv" \
   "${WORKDIR}/custom_genomes_file_taxid.tsv"

###############################################################################
# Add custom genomes to Kraken DB
###############################################################################
echo "[INFO] DB=${DB}"
echo "[INFO] WORKDIR=${WORKDIR}"
echo "[INFO] pwd=$(pwd)"
DB=/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal
WORKDIR=/home/pthorpe001/data/project_back_up_2024/kracken
pwd=/home/pthorpe001/data/project_back_up_2024/kracken


set -euo pipefail

input_tsv="${WORKDIR}/custom_genomes_file_taxid.tsv"
clean_tsv="${WORKDIR}/custom_genomes_file_taxid.clean.tsv"
tagged_dir="${WORKDIR}/kraken_tagged_fastas"

mkdir -p "${tagged_dir}"

awk -F $'\t' '
  NF == 2 && $2 ~ /^[0-9]+$/ && $1 !~ /\/\._/ {
    print $1 "\t" $2
  }
' "${input_tsv}" > "${clean_tsv}"

echo "[INFO] Cleaned TSV written to ${clean_tsv}"



while IFS=$'\t' read -r genome_file taxid
do
  if [[ "${genome_file}" != /* ]]; then
    genome_path="${WORKDIR}/${genome_file}"
  else
    genome_path="${genome_file}"
  fi

  if [[ ! -f "${genome_path}" ]]; then
    echo "[WARN] Missing file, skipping: ${genome_path}"
    continue
  fi

  base_name="$(basename "${genome_path}")"
  base_name="${base_name%.gz}"
  output_fasta="${tagged_dir}/${base_name}"

  echo "[INFO] Rewriting headers: ${genome_path} taxid=${taxid}"

  if [[ "${genome_path}" == *.gz ]]; then
    gzip -cd "${genome_path}" | \
    awk -v taxid="${taxid}" '
      /^>/ {
        sub(/^>/, "")
        print ">kraken:taxid|" taxid "|" $0
        next
      }
      {
        print
      }
    ' > "${output_fasta}"
  else
    awk -v taxid="${taxid}" '
      /^>/ {
        sub(/^>/, "")
        print ">kraken:taxid|" taxid "|" $0
        next
      }
      {
        print
      }
    ' "${genome_path}" > "${output_fasta}"
  fi

  echo "[INFO] add-to-library: ${output_fasta}"
  kraken2-build --add-to-library "${output_fasta}" --db "${DB}"
done < "${clean_tsv}"




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