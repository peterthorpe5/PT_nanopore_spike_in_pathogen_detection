#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 2
#$ -jc long
#$ -mods l_hard mfree 64G
#$ -adds l_hard h_vmem 64G
#$ -N build_metabuili_db

cd /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity



python /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection/scripts/prepare_kraken_libraries_for_metabuli.py \
  --input_fastas \
    /home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal/library/library.fna \
    /home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal/library/viral/library.fna \
   /home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal/library/bacteria/library.fna \
  --out_dir /home/pthorpe001/data/databases/metabuli



metabuli build \
  /home/pthorpe001/data/databases/metabuli/custom_metabuli_db \
  /home/pthorpe001/data/databases/metabuli/fasta_list.txt \
  /home/pthorpe001/data/databases/metabuli/accession2taxid.tsv \
  --taxonomy-path /home/pthorpe001/data/databases/metamaps_from_kraken_shared/taxonomy \
  --make-library 1 \
  --threads 12 \
  --max-ram 180 \
  --validate-input 1 \
  --validate-db 1