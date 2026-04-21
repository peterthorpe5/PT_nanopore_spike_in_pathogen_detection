#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 2
#$ -jc long
#$ -N summary

### This script is used to summarise the results of multiple spike-in runs and generate a report.


cd /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity

mv runs_* runs/  # move all run output dirs into a single 'runs' dir for easier summarisation


python PT_nanopore_spike_in_pathogen_detection/summary/run_spikein_summary_pipeline.py \
  --runs_dir runs \
  --out_dir spikein_summary_report \
  --threshold_mode fixed \
  --min_detect_value 1 \
  --target_fpr 0.05 \
  --make_krona_inputs \
  --verbose

