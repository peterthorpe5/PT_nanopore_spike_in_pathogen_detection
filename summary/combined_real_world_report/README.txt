Combined real-world spike-in report bundle

Contents
--------
- build_combined_real_world_report.py
  Main Python script.

- example_output/
  Example HTML, Excel, and TSV outputs generated from the current benchmark
  workbooks in this conversation.

How to run
----------
python build_combined_real_world_report.py \
  --method_performance_xlsx /path/to/method_performance.xlsx \
  --replicate_report_xlsx /path/to/replicate_resolved_report.xlsx \
  --threshold_report_xlsx /path/to/threshold_calibration_report.xlsx \
  --out_dir /path/to/output_dir \
  --report_title "Combined ONT spike-in benchmark report"

Alternative
-----------
If the three workbooks are already in one directory with their standard names:

python build_combined_real_world_report.py \
  --input_dir /path/to/summary_dir \
  --out_dir /path/to/output_dir

Important
---------
The script reports two layers of metrics:
1. Existing tracked-target metrics from the current benchmark workbooks.
2. Stricter taxonomic-clean metrics that penalise off-target taxa and require
   taxonomically clean negative observations.