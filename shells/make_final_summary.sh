

### This script is used to summarise the results of multiple spike-in runs and generate a report.

mv runs_* runs/  # move all run output dirs into a single 'runs' dir for easier summarisation


  python PT_nanopore_spike_in_pathogen_detection/summary/summarise_spikein_runs_v4_reported_taxa_full.py \
    --input_dirs runs \
    --out_dir spikein_summary_report \
    --verbose


# get performance stats

# or less conservative:

  python PT_nanopore_spike_in_pathogen_detection/summary/build_method_performance_table_with_reported_taxa.py \
    --combined_long_tsv spikein_summary_report/combined_long.tsv \
    --out_dir spikein_summary_report \
    --reported_taxa_long_tsv spikein_summary_report/reported_taxa_long.tsv \
    --threshold_mode fixed \
    --min_detect_value 1

########################
## After running the above, you can then generate a report with:

PT_nanopore_spike_in_pathogen_detection/shells/make_krona_from_kraken_outputs.sh



 python PT_nanopore_spike_in_pathogen_detection/summary/make_spikein_report_v5.py \
    --summary_dir spikein_summary_report \
    --title "ONT spike-in summary report"



python PT_nanopore_spike_in_pathogen_detection/summary/make_spikein_replicate_report_v2.py \
  --summary_dir spikein_summary_report \
  --out_dir spikein_summary_report/replicate_resolved_report \
  --title "ONT spike-in replicate-resolved report" \
  --threshold_mode fixed \
  --min_detect_value 1



python PT_nanopore_spike_in_pathogen_detection/summary/make_spikein_threshold_calibration_report_v3.py \
  --summary_dir spikein_summary_report \
  --out_dir spikein_summary_report/threshold_calibration_report_v3 \
  --title "ONT spike-in threshold calibration report" \
  --min_detect_value 1 \
  --target_fpr 0.05




# final real world report:
python3 PT_nanopore_spike_in_pathogen_detection/summary/build_combined_real_world_report.py \
  --method_performance_xlsx spikein_summary_report/method_performance.xlsx \
  --replicate_report_xlsx spikein_summary_report/replicate_resolved_report/replicate_resolved_report.xlsx \
  --threshold_report_xlsx spikein_summary_report/threshold_calibration_report_v3/threshold_calibration_report.xlsx \
  --out_dir spikein_summary_report/combined_real_world_report \
  --report_title "Combined ONT spike-in benchmark report with real-world taxonomic burden"





