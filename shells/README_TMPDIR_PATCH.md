# TMPDIR-patched shell scripts

These files keep the original script filenames. Heavy scripts run in `$TMPDIR` and copy selected outputs back to the permanent `OUT_DIR`.

Set `COPY_WORKING_FILES=true` to copy FASTQ, BAM, FASTA, Flye, and Medaka intermediates as well. Set `USE_TMPDIR=false` only for local debugging.

Patched scripts:
- `run_spikein_multi_flye.sh`
- `run_spikein_multi_flye_clean.sh`
- `run_spikein_multi_flye_medaka.sh`
- `run_spikein_multi_flye_medaka_clean.sh`
- `run_spikein_multi_minimap_only.sh`
- `run_spikein_multi_readlevel.sh`
- `run_spikein_multi_readlevel_kraken_confidence.sh`
- `run_spikein_multi_readlevel_metamaps.sh`
- `run_spikein_readlevel_metabuli.sh`
- `run_spikein_single_flye.sh`
- `run_spikein_single_flye_medaka.sh`
- `run_spikein_single_flye_medaka2.sh`
- `run_spikein_single_minimap_only.sh`
- `run_spikein_single_readlevel.sh`
- `run_spikein_single_readlevel_kraken_confidence.sh`
- `run_spikein_single_readlevel_metamaps.sh`

Skipped:
- None
