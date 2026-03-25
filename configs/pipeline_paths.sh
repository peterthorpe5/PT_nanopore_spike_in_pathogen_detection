#!/usr/bin/env bash

# Shared path and database defaults for the spike-in pipeline.
# Source this file after REPO_DIR has been resolved.

REPO_DIR="${REPO_DIR:?REPO_DIR must be set before sourcing pipeline_paths.sh}"

SHELLS_DIR="${SHELLS_DIR:-${REPO_DIR}/shells}"
PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR:-${REPO_DIR}/scripts}"
CONFIG_DIR="${CONFIG_DIR:-${REPO_DIR}/configs}"
BUILD_DB_DIR="${BUILD_DB_DIR:-${REPO_DIR}/build_DB}"
SUMMARY_DIR="${SUMMARY_DIR:-${REPO_DIR}/summary}"

DB_BASE="${DB_BASE:-/home/pthorpe001/data/databases}"
KRAKEN_DB_DIR_DEFAULT="${KRAKEN_DB_DIR_DEFAULT:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"
MINIMAP_DB_FASTA_DEFAULT="${MINIMAP_DB_FASTA_DEFAULT:-${DB_BASE}/minimap2/shared_bact_viral_plasmo_refs.fa}"
MINIMAP_DB_INDEX_DEFAULT="${MINIMAP_DB_INDEX_DEFAULT:-${DB_BASE}/minimap2/shared_bact_viral_plasmo_refs.mmi}"


# Core database paths
DB_BASE="${DB_BASE:-/home/pthorpe001/data/databases}"
KRAKEN_DB_DIR_DEFAULT="${KRAKEN_DB_DIR_DEFAULT:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"
MINIMAP_DB_FASTA_DEFAULT="${MINIMAP_DB_FASTA_DEFAULT:-${DB_BASE}/minimap2/shared_bact_viral_plasmo_refs.fa}"
MINIMAP_DB_INDEX_DEFAULT="${MINIMAP_DB_INDEX_DEFAULT:-${DB_BASE}/minimap2/shared_bact_viral_plasmo_refs.mmi}"


# Core input data paths
REAL_FASTQ_DEFAULT="${REAL_FASTQ_DEFAULT:-/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC1023_AmM008WB.fastq.gz}"
DEPLETED_FASTQ_DEFAULT="${DEPLETED_FASTQ_DEFAULT:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/MRC1023_AmM008WB.depleted.fastq.gz}"

# Host / depletion references
SMALL_MONKEY_REF_DEFAULT="${SMALL_MONKEY_REF_DEFAULT:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/Mnem_1.0.fna.gz}"
DEPLETION_REF_DEFAULT="${DEPLETION_REF_DEFAULT:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/host_plus_plasmo_depletion.fasta}"

# Shared defaults
THREADS_DEFAULT="${THREADS_DEFAULT:-12}"
SPIKE_LEVELS_DEFAULT="${SPIKE_LEVELS_DEFAULT:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES_DEFAULT="${REPLICATES_DEFAULT:-3}"


# Common external inputs
REAL_FASTQ_DEFAULT="${REAL_FASTQ_DEFAULT:-/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC1023_AmM008WB.fastq.gz}"
MONKEY_SMALL_GZ_DEFAULT="${MONKEY_SMALL_GZ_DEFAULT:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna.gz}"
MONKEY_SMALL_FASTA_DEFAULT="${MONKEY_SMALL_FASTA_DEFAULT:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna}"
DEPLETION_REF_FASTA_DEFAULT="${DEPLETION_REF_FASTA_DEFAULT:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/host_plus_plasmo_depletion.fasta}"
DEPLETED_FASTQ_DEFAULT="${DEPLETED_FASTQ_DEFAULT:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/MRC1023_AmM008WB.depleted.fastq.gz}"
DEFAULT_PATHOGEN_FASTA="${DEFAULT_PATHOGEN_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/GCA_014843685.1_ASM1484368v1_genomic.fna}"
DEFAULT_TARGET_LABEL="${DEFAULT_TARGET_LABEL:-Plasmodium vivax}"
DEFAULT_PATHOGEN_PANEL_2="${DEFAULT_PATHOGEN_PANEL_2:-${CONFIG_DIR}/pathogen_panel_2.tsv}"
DEFAULT_PATHOGEN_PANEL_3="${DEFAULT_PATHOGEN_PANEL_3:-${CONFIG_DIR}/pathogen_panel_3.tsv}"
DEFAULT_PATHOGEN_PANEL_3_WITH_TAXID="${DEFAULT_PATHOGEN_PANEL_3_WITH_TAXID:-${CONFIG_DIR}/pathogen_panel_3.tsv}"



# Shared path and database defaults for the spike-in pipeline.
# Source this file after REPO_DIR has been resolved.

REPO_DIR="${REPO_DIR:?REPO_DIR must be set before sourcing pipeline_paths.sh}"

SHELLS_DIR="${SHELLS_DIR:-${REPO_DIR}/shells}"
PY_SCRIPTS_DIR="${PY_SCRIPTS_DIR:-${REPO_DIR}/scripts}"
CONFIG_DIR="${CONFIG_DIR:-${REPO_DIR}/configs}"
BUILD_DB_DIR="${BUILD_DB_DIR:-${REPO_DIR}/build_DB}"
SUMMARY_DIR="${SUMMARY_DIR:-${REPO_DIR}/summary}"

# Database locations
DB_BASE="${DB_BASE:-/home/pthorpe001/data/databases}"
KRAKEN_DB_DIR_DEFAULT="${KRAKEN_DB_DIR_DEFAULT:-/home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal}"
MINIMAP_DB_FASTA_DEFAULT="${MINIMAP_DB_FASTA_DEFAULT:-${DB_BASE}/minimap2/shared_bact_viral_plasmo_refs.fa}"
MINIMAP_DB_INDEX_DEFAULT="${MINIMAP_DB_INDEX_DEFAULT:-${DB_BASE}/minimap2/shared_bact_viral_plasmo_refs.mmi}"
METAMAPS_DB_DIR_DEFAULT="${METAMAPS_DB_DIR_DEFAULT:-/home/pthorpe001/data/databases/metamaps_from_kraken_shared/custom_metamaps_db}"


# Backwards-compatibility alias for older shells still expecting a gz variable
MINIMAP_DB_GZ="${MINIMAP_DB_GZ:-${MINIMAP_DB_FASTA_DEFAULT}}"

# Core input data
REAL_FASTQ_DEFAULT="${REAL_FASTQ_DEFAULT:-/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC1023_AmM008WB.fastq.gz}"
DEPLETED_FASTQ_DEFAULT="${DEPLETED_FASTQ_DEFAULT:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/MRC1023_AmM008WB.depleted.fastq.gz}"

# Host and depletion references
MONKEY_SMALL_GZ_DEFAULT="${MONKEY_SMALL_GZ_DEFAULT:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna.gz}"
MONKEY_SMALL_FASTA_DEFAULT="${MONKEY_SMALL_FASTA_DEFAULT:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/genome/GCF_000956065.1_Mnem_1.0_genomic.fna}"
DEPLETION_REF_FASTA_DEFAULT="${DEPLETION_REF_FASTA_DEFAULT:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/host_plus_plasmo_depletion.fasta}"

# Shared workflow defaults
THREADS_DEFAULT="${THREADS_DEFAULT:-12}"
SPIKE_LEVELS_DEFAULT="${SPIKE_LEVELS_DEFAULT:-0 1 5 10 25 50 100 250 500 1000 2500 5000}"
REPLICATES_DEFAULT="${REPLICATES_DEFAULT:-3}"

# Single-genome defaults
DEFAULT_PATHOGEN_FASTA="${DEFAULT_PATHOGEN_FASTA:-/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/GCA_014843685.1_ASM1484368v1_genomic.fna}"
DEFAULT_TARGET_LABEL="${DEFAULT_TARGET_LABEL:-Plasmodium vivax}"

# Multi-genome config defaults
DEFAULT_PATHOGEN_PANEL_2="${DEFAULT_PATHOGEN_PANEL_2:-${CONFIG_DIR}/pathogen_panel_2.tsv}"
DEFAULT_PATHOGEN_PANEL_3="${DEFAULT_PATHOGEN_PANEL_3:-${CONFIG_DIR}/pathogen_panel_3.tsv}"
DEFAULT_PATHOGEN_PANEL_3_WITH_TAXID="${DEFAULT_PATHOGEN_PANEL_3_WITH_TAXID:-${CONFIG_DIR}/pathogen_panel_with_taxid_template.tsv}"