# ONT in silico spike-in framework for pathogen detection

This repository contains a reproducible framework for benchmarking low-abundance pathogen detection in Oxford Nanopore Technologies (ONT) sequencing data using in silico spike-in experiments. Synthetic ONT-like reads are generated from one or more pathogen genomes using NanoSim, added to a real empirical ONT background dataset, and then analysed using either direct read-level classification or assembly-based detection.

The framework supports:

- single-genome spike-in experiments
- multi-genome equal-contribution spike-in experiments
- read-level detection using Kraken2 and minimap2
- assembly-based detection using Flye followed by Kraken2
- shuffled-sequence negative controls
- automated summary collation and plotting across runs

All summary tables are written as tab-separated files.

---

## Overview

The workflow is designed to test how sensitively low-abundance pathogen reads can be detected against a realistic host-derived ONT background. The same empirical FASTQ is used both to train the NanoSim error model and to provide the background into which simulated reads are spiked.

The general strategy is:

1. sample empirical reads from a real ONT FASTQ
2. train NanoSim on those reads using a host reference
3. simulate ONT-like reads from one or more pathogen genomes
4. optionally shuffle pathogen sequence for negative controls
5. deplete the real background FASTQ against a combined host plus masked parasite reference
6. spike simulated reads into the depleted background at defined abundances
7. analyse mixed datasets using:
   - Kraken2 directly on reads
   - minimap2 against a masked Plasmodium/outgroup reference
   - Flye assembly followed by Kraken2 on assembled contigs
8. summarise all runs across replicates, spike levels, and workflows

---

## Repository structure

```text
.
├── configs
│   ├── pathogen_panel_2.tsv
│   └── pathogen_panel_3.tsv
├── kraken_DB
│   ├── build_kraken_DB.sh
│   ├── build_kraken_DB_slurm.sh
│   ├── create_viral_list.sh
│   └── make_kraken_taxid_map.py
├── old
├── run_spikein_multi_flye.sh
├── run_spikein_multi_readlevel.sh
├── run_spikein_shuffled_flye.sh
├── run_spikein_shuffled_readlevel.sh
├── run_spikein_single_flye.sh
├── run_spikein_single_readlevel.sh
├── scripts
│   ├── assembly_stats.py
│   ├── build_mixed_fastq.py
│   ├── combine_nanosim_fastq.py
│   ├── count_bed_hits.py
│   ├── dedup_fastq_names.py
│   ├── sample_fastq.py
│   ├── shuffle_fasta.py
│   ├── summarise_kraken_report.py
│   └── summarise_spikein_runs.py
├── spikein_utils.py
└── README.md


Notes:

scripts/ contains the current standalone Python helper scripts used by the newer shell workflows.

configs/ contains pathogen panels for multi-genome experiments.

old/ contains earlier monolithic shell scripts retained for reference.

spikein_utils.py may still be present for backwards compatibility with older runs, but the current recommended workflow is based on the standalone scripts in scripts/.

Core workflows
1. Single-genome read-level spike-in

run_spikein_single_readlevel.sh

This workflow:

samples empirical reads from the real ONT FASTQ

trains NanoSim on the sampled reads

simulates ONT-like reads from a single pathogen genome

creates spike-in mixtures across a series of spike levels and replicates

classifies mixed reads with Kraken2

aligns mixed reads to the masked Plasmodium/outgroup reference with minimap2

writes a per-run spikein_summary.tsv

Primary outputs:

mixed FASTQ files per replicate and spike level

Kraken2 report and classification files

minimap2 BAM and BED files

spikein_summary.tsv

2. Single-genome assembly-based spike-in

run_spikein_single_flye.sh

This workflow is similar to the single-genome read-level workflow, but after spike-in generation:

read identifiers are deduplicated if required

each mixed FASTQ is assembled independently using Flye in metagenome mode

assembly.fasta is classified with Kraken2

assembly statistics are computed

a per-run assembly summary table is written

Primary outputs:

deduplicated FASTQ files

Flye assembly directories

Kraken2 reports for assemblies

assembly summary tables

3. Multi-genome equal-contribution read-level spike-in

run_spikein_multi_readlevel.sh

This workflow supports two or more pathogen genomes at equal contribution per spike level. Pathogens are specified in a TSV config file.

At each spike level:

the same number of reads is sampled independently from each pathogen pool

all pathogen spike reads are combined

the combined spike is added to the depleted background

Kraken2 and minimap2 analyses are run as for the single-genome read-level workflow

Primary outputs:

combined spike FASTQ per replicate and spike level

mixed FASTQ files

Kraken2 report and classification files

minimap2 BAM and BED files

spikein_multi_summary.tsv

4. Multi-genome equal-contribution assembly-based spike-in

run_spikein_multi_flye.sh

This workflow extends the multi-genome spike-in design to assembly-based detection:

equal-contribution multi-pathogen mixtures are generated

mixed FASTQ files are deduplicated if required

Flye assemblies are produced for each replicate and spike level

assembled contigs are classified with Kraken2

assembly summaries are written

Primary outputs:

Flye assembly directories

Kraken2 assembly reports

spikein_multi_flye_summary.tsv

5. Shuffled negative controls

run_spikein_shuffled_readlevel.sh
run_spikein_shuffled_flye.sh

These workflows generate negative-control pathogen references by shuffling each pathogen contig independently before simulation.

The purpose is to test whether apparent detection can arise from:

sequence composition alone

non-specific Kraken classification

non-specific minimap2 alignment

assembly artefacts

The shuffled controls then follow the same read-level or assembly-based workflow as the true pathogen runs.

Primary outputs:

shuffled FASTA file

shuffle metadata table

standard read-level or assembly-based outputs

Python helper scripts

The current workflows rely on standalone helper scripts in scripts/.

sample_fastq.py

Reservoir-samples reads from a gzipped FASTQ file without loading the full file into memory.

Typical use:

python3 scripts/sample_fastq.py \
  --fastq_gz input.fastq.gz \
  --n_reads 200000 \
  --seed 1 \
  --out_fastq_gz sampled.fastq.gz
combine_nanosim_fastq.py

Combines NanoSim FASTQ outputs into a single FASTQ file for downstream sampling.

build_mixed_fastq.py

Builds a mixed gzipped FASTQ file from a background FASTQ and one or more spike FASTQ files. This is preferred over raw gzip concatenation because it provides a safer and more transparent construction step.

shuffle_fasta.py

Creates a shuffled FASTA for negative-control experiments. Standard A/C/G/T positions are shuffled, while non-ACGT characters are preserved in place.

summarise_kraken_report.py

Parses Kraken2 reports and extracts:

classified read count

unclassified read count

target-specific read or contig count

This script replaces earlier logic that could mis-handle classified read totals.

count_bed_hits.py

Counts BED entries matching a given target label or alias set. This is used for minimap2-derived target counting and can also be used to inspect which reference names are being hit.

dedup_fastq_names.py

Ensures that FASTQ read identifiers are unique before Flye assembly.

assembly_stats.py

Computes basic assembly metrics from an assembly FASTA, including contig count and total assembled bases.

summarise_spikein_runs.py

Recursively scans run output directories, collates summary tables across workflows, writes harmonised combined tables, logs missing or malformed outputs, and generates overview plots.

Configuration files for multi-genome runs

Multi-genome workflows use TSV config files in configs/.

Examples:

configs/pathogen_panel_2.tsv

configs/pathogen_panel_3.tsv

The expected structure is:

pathogen_fasta	target_label
/path/to/genome1.fna	Plasmodium vivax
/path/to/genome2.fna	Plasmodium falciparum

For three-genome runs:

pathogen_fasta	target_label
/path/to/genome1.fna	Plasmodium vivax
/path/to/genome2.fna	Plasmodium falciparum
/path/to/genome3.fna	Plasmodium knowlesi

The first column should contain the reference FASTA path.
The second column should contain the label used for downstream matching and reporting.

Required software

The shell workflows assume the following tools are available on PATH:

python3

gzip

minimap2

samtools

kraken2

bamToBed

read_analysis.py from NanoSim

simulator.py from NanoSim

flye for assembly workflows

Additional Python packages may be required for the helper scripts, depending on your local environment.

Typical inputs

The pipeline requires:

a real ONT background FASTQ

a smaller host reference for NanoSim training

a larger depletion reference for host/background filtering

a masked Plasmodium/outgroup reference for minimap2

one or more pathogen genome FASTA files

a Kraken2 database containing the desired taxonomic scope

The real ONT FASTQ is used in two distinct ways:

as the empirical basis for NanoSim model training

as the background into which simulated pathogen reads are spiked

Output structure

Each run creates a separate output directory. Within each run, replicate and spike-level combinations are typically placed in separate subdirectories.

For read-level workflows, outputs commonly include:

mixed.fastq.gz

kraken.report.tsv

kraken.classifications.tsv

kraken.summary.tsv

minimap_sorted.bam

minimap_sorted.MASKED.bed

spikein_summary.tsv or spikein_multi_summary.tsv

For assembly-based workflows, outputs commonly include:

deduplicated FASTQ

Flye output directories

assembly.fasta

Kraken2 report for assemblies

assembly metrics

spikein_flye_summary.tsv or spikein_multi_flye_summary.tsv

Shuffled-control runs additionally include:

pathogen.shuffled.fasta

shuffle_metadata.tsv

Example job submission

The newer shell scripts are designed to work on an SGE cluster. To avoid issues with spool directories, it is recommended to pass the repository root explicitly via REPO_DIR.

Move to the repository root first:

cd /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/PT_nanopore_spike_in_pathogen_detection
Single-genome read-level
qsub -v REPO_DIR="$(pwd)" run_spikein_single_readlevel.sh
Single-genome assembly-based
qsub -v REPO_DIR="$(pwd)" run_spikein_single_flye.sh
Two-genome read-level
qsub -v REPO_DIR="$(pwd)",PATHOGEN_CONFIG_TSV="$(pwd)/configs/pathogen_panel_2.tsv" run_spikein_multi_readlevel.sh
Three-genome assembly-based
qsub -v REPO_DIR="$(pwd)",PATHOGEN_CONFIG_TSV="$(pwd)/configs/pathogen_panel_3.tsv" run_spikein_multi_flye.sh
Shuffled controls
qsub -v REPO_DIR="$(pwd)" run_spikein_shuffled_readlevel.sh
qsub -v REPO_DIR="$(pwd)" run_spikein_shuffled_flye.sh

If preferred, absolute paths can be used directly for PATHOGEN_CONFIG_TSV and other variables.

Summarising runs

The scripts/summarise_spikein_runs.py script can be used to collate and visualise results across multiple run directories.

Typical use:

python3 scripts/summarise_spikein_runs.py \
  --input_dirs /path/to/run_parent_1 /path/to/run_parent_2 \
  --out_dir spikein_summary_report \
  --verbose

This script will:

recursively find recognised spike-in summary files

merge outputs across workflows

write harmonised wide and long tabular summaries

record missing or malformed files

generate overview plots in PNG and PDF format

Expected summary outputs include:

run_manifest.tsv

combined_wide.tsv

combined_long.tsv

missing_or_problematic_files.tsv

plot_manifest.tsv

plots/

Interpretation notes
Kraken2

Kraken2 is used both at the read level and on assembled contigs. Parsed summary tables should be interpreted carefully, especially when comparing baseline background signal against low spike-in levels.

minimap2

Minimap2 alignments are counted using BED output and target-label matching. The quality of minimap-based target summaries depends on how well target labels correspond to the actual reference names in the masked alignment reference. For complex references containing multiple related taxa, it may be necessary to use explicit alias mapping rather than simple free-text matching.

Shuffled controls

Detection in shuffled-sequence controls should be interpreted as evidence of non-specific classification, alignment, or assembly behaviour rather than genuine biological signal.

Historical scripts

The old/ directory contains earlier shell workflows that were used during development. These may still be useful for provenance or comparison, but they are not the recommended current entry points for new analyses.

The current recommended workflows are:

run_spikein_single_readlevel.sh

run_spikein_single_flye.sh

run_spikein_multi_readlevel.sh

run_spikein_multi_flye.sh

run_spikein_shuffled_readlevel.sh

run_spikein_shuffled_flye.sh

with helper scripts from scripts/.

Notes on maintenance

The current codebase separates orchestration from parsing and manipulation:

shell scripts handle looping, environment setup, and tool execution

Python helper scripts handle FASTQ sampling, FASTA shuffling, Kraken parsing, BED hit counting, mixed FASTQ construction, assembly summaries, and run collation

This separation is intended to make the pipeline easier to debug, maintain, and extend.
