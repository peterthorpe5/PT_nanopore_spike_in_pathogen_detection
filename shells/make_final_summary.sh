

### This script is used to summarise the results of multiple spike-in runs and generate a report.

python3 summarise_spikein_runs_v3.py \
  --input_dirs /path/to/runs \
  --out_dir /path/to/spikein_summary_report \
  --verbose

 python PT_nanopore_spike_in_pathogen_detection/summary/summarise_spikein_runs_v3.py \
  --input_dirs runs --out_dir spikein_summary_report --verbose


  python PT_nanopore_spike_in_pathogen_detection/summary/summarise_spikein_runs_v4_with_multi_metamaps.py \
  --input_dirs runs \
  --out_dir spikein_summary_report \
  --verbose

########################
## After running the above, you can then generate a report with:

python3 make_spikein_report_v5.py \
  --summary_dir /path/to/spikein_summary_report \
  --title "ONT spike-in summary report"


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


