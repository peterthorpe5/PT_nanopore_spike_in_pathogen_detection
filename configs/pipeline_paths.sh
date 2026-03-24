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
METAMAPS_DB_DIR_DEFAULT="${METAMAPS_DB_DIR_DEFAULT:-${DB_BASE}/metamaps/current/custom_metamaps_db/custom_metamaps_db}"

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
