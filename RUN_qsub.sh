

cd /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity

qsub -v REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection" \
  /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection/run_spikein_single_readlevel.sh

qsub -v REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection" \
  /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection/run_spikein_single_flye.sh

qsub -v REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection",PATHOGEN_CONFIG_TSV="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection/configs/pathogen_panel_2.tsv" \
  /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection/run_spikein_multi_readlevel.sh

qsub -v REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection",PATHOGEN_CONFIG_TSV="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection/configs/pathogen_panel_3.tsv" \
  /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection/run_spikein_multi_flye.sh

qsub -v REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection" \
  /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection/run_spikein_shuffled_readlevel.sh

qsub -v REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection" \
  /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection/run_spikein_shuffled_flye.sh

  
