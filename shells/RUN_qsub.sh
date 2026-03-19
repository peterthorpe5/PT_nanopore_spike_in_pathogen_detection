# 1. Move to your working directory
cd /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity

# 2. Define and export the base repository path
export REPO_DIR="/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection"

# 3. Submit the jobs using -V to pass the exported environment variables
# Note: For jobs requiring specific configs, we export the config just before the qsub call.

# --- Single Read Level ---
qsub -V "$REPO_DIR/run_spikein_single_readlevel.sh"

qsub -V "$REPO_DIR/run_spikein_single_flye_medaka.sh"

# --- Multi Read Level (with specific configs) ---
export PATHOGEN_CONFIG_TSV="$REPO_DIR/configs/pathogen_panel_2.tsv"
qsub -V "$REPO_DIR/run_spikein_multi_readlevel.sh"

export PATHOGEN_CONFIG_TSV="$REPO_DIR/configs/pathogen_panel_3.tsv"
qsub -V "$REPO_DIR/run_spikein_multi_flye_medaka.sh"

# --- Shuffled Read Level ---
qsub -V "$REPO_DIR/run_spikein_shuffled_readlevel.sh"

qsub -V "$REPO_DIR/run_spikein_shuffled_flye_medaka.sh"



# --- Melon Single ---
qsub -V "$REPO_DIR/melon/run_melon_single.sh"

# --- Melon Multi ---
export PATHOGEN_CONFIG_TSV="$REPO_DIR/configs/pathogen_panel_2.tsv"
qsub -V "$REPO_DIR/melon/run_melon_mixed.sh"



# --- Melon Multi ---
export PATHOGEN_CONFIG_TSV="$REPO_DIR/configs/pathogen_panel_3.tsv"
qsub -V "$REPO_DIR/melon/run_melon_mixed.sh"