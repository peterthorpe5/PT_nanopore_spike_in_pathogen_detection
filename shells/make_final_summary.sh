

### This script is used to summarise the results of multiple spike-in runs and generate a report.


  python PT_nanopore_spike_in_pathogen_detection/summary/summarise_spikein_runs_v4.py \
  --input_dirs runs \
  --out_dir spikein_summary_report \
  --verbose


# get performance stats

python PT_nanopore_spike_in_pathogen_detection/summary/build_method_performance_table.py \
  --combined_long_tsv spikein_summary_report/combined_long.tsv \
  --out_dir spikein_summary_report


# or less conservative:

python PT_nanopore_spike_in_pathogen_detection/summary/build_method_performance_table.py \
  --combined_long_tsv spikein_summary_report/combined_long.tsv \
  --out_dir spikein_summary_report \
  --threshold_mode fixed \
  --min_detect_value 1

########################
## After running the above, you can then generate a report with:

PT_nanopore_spike_in_pathogen_detection/shells/make_krona_from_kraken_outputs.sh

 python PT_nanopore_spike_in_pathogen_detection/summary/make_spikein_report_v5.py \
  --summary_dir spikein_summary_report \
  --title "ONT spike-in summary report"



########################
## You can also generate Krona inputs for visualisation of taxonomic classifications:
python3 make_krona_inputs.py \
  --input_dirs /path/to/runs \
  --out_dir /path/to/krona_inputs

 python PT_nanopore_spike_in_pathogen_detection/summary/make_krona_inputs.py \
  --input_dirs spikein_summary_report \
  --out_dir krona_inputs


