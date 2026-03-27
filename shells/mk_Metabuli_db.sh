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



python /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection/scripts/prepare_kraken_library_for_metamaps.py \
  --input_fastas \
    /home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal/library/library.fna \
    /home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal/library/viral/library.fna \
   /home/pthorpe001/data/project_back_up_2024/kraken_bact_virus_plasmo_fungal/library/bacteria/library.fna \
  --out_dir /home/pthorpe001/data/databases/metabuli
