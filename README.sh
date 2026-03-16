# Nanopore spike-in detection threshold pipeline (Kraken2 + minimap2)

This workflow performs in silico spike-in experiments for Oxford Nanopore FASTQ datasets. It takes a real background sample (`.fastq.gz`), optionally removes host reads, simulates ONT-like reads from one or more pathogen genomes using NanoSim, spikes simulated reads into the real dataset at varying amounts, and then runs two classification/detection approaches:

1. Kraken2 taxonomic classification
2. A minimap2 alignment-based method against a hard-masked Plasmodium/outgroup/bait reference database

The primary output is a tab-separated summary table reporting, per spike level and replicate, how many reads were detected by each method.

## Overview

For a single sample and a single pathogen genome:

1. Host depletion (optional but recommended)
2. Subsample real reads for NanoSim training
3. Train a NanoSim error/length profile on the real reads
4. Simulate a pool of ONT-like reads from the pathogen genome
5. Create spike-in mixtures across a series of spike levels (absolute simulated read counts)
6. Run Kraken2 and minimap2 per mixture
7. Write a single TSV summary

A negative-control run is also supported where the pathogen genome is shuffled prior to simulation. This tests whether shuffled sequence can generate off-target hits in Kraken2 or minimap2.

## Repository contents

- `run_spikein_one_sample.sh`
  - end-to-end run for one background sample and one pathogen
  - optional host depletion
  - NanoSim profiling and simulation
  - spike-in series
  - Kraken2 + minimap2 runs
  - writes `spikein_summary.tsv`

- `spikein_utils.py`
  - Python utilities for reservoir sampling FASTQ.GZ, combining NanoSim outputs, parsing Kraken2 reports, and counting minimap2 target alignments
  - all outputs are TSV (no comma-separated outputs)

- `run_spikein_shuffled_control.sh`
  - shuffles the pathogen genome, then calls `run_spikein_one_sample.sh`
  - used as a negative control for off-target detection

## Requirements

### Software

The pipeline assumes Oxford Nanopore reads (FASTQ.GZ). You will need:

- minimap2
- samtools
- kraken2
- bedtools (for `bamToBed`)
- python3
- gzip

NanoSim must be installed and available on your `PATH` in the active environment:

- `read_analysis.py`
- `simulator.py`

Optional but helpful:

- `seqkit` (the scripts do not require it, but it is useful for independent QC)

### Conda (example)

Create a dedicated environment (example only; adjust to your local setup):

```bash
conda create --name spikein_ont \
  python=3.11 \
  minimap2 samtools bedtools kraken2 -c bioconda -c conda-forge

conda activate spikein_ont

# Install NanoSim according to its documentation
# https://github.com/bcgsc/NanoSim



Inputs

You will need:

A real background FASTQ.GZ sample

A host reference genome FASTA (single FASTA file) for host depletion and NanoSim profiling

A pathogen genome FASTA to simulate spike-in reads from

A minimap2 reference FASTA database for the alignment-based detection method

in this workflow, the reference is a hard-masked Plasmodium/outgroup/bait database

A Kraken2 database directory

Example paths used in this project:

Background FASTQ folder:

/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples

Host reference genome folder (you must point to a single FASTA file for HOST_REF_FASTA):

/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use

Pathogen genome example:

/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/genomes/GCF_000524495.1_Plas_inui_San_Antonio_1_V1_genomic.fasta

Hard-masked minimap2 detection database (gzipped FASTA):

/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz

Kraken2 DB folder:

/home/pthorpe001/data/project_back_up_2024/kracken

Configuration

The main script reads configuration from environment variables (with sensible defaults for many settings). At minimum, you must set:

REAL_FASTQ (path to a single .fastq.gz background sample)

HOST_REF_FASTA (path to a single host reference .fasta file)

KRAKEN_DB_DIR (path to Kraken2 DB directory)

Optional overrides include:

PATHOGEN_FASTA (pathogen genome fasta)

MINIMAP_DB_GZ (hard-masked minimap database fasta.gz)

OUT_DIR (output folder)

THREADS (threads for minimap2/kraken2/nanosim)

DO_HOST_DEPLETION (true or false)

TRAIN_READS_N (reads used for NanoSim training; default 200000)

SIM_POOL_N (simulated pathogen read pool size; default 20000)

SPIKE_LEVELS (absolute read counts to spike)

REPLICATES (number of replicates per spike level)

MIN_ALIGN_LEN (min alignment length filter used in minimap2 pipeline; default 500)

MIN_MAPQ (min MAPQ filter used in minimap2 pipeline; default 15)

Running the one-pathogen pilot

Activate your environment and export the required paths:




conda activate spikein_ont

export REAL_FASTQ="/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC0123_AmM001WB.fastq.gz"
export HOST_REF_FASTA="/path/to/your/host_reference.fasta"
export PATHOGEN_FASTA="/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/genomes/GCF_000524495.1_Plas_inui_San_Antonio_1_V1_genomic.fasta"
export MINIMAP_DB_GZ="/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz"
export KRAKEN_DB_DIR="/home/pthorpe001/data/project_back_up_2024/kracken"

bash run_spikein_one_sample.sh



Output structure

The OUT_DIR (default spikein_pilot_out) contains:

background.fastq.gz

host-depleted background reads if host depletion is enabled, otherwise a symlink to REAL_FASTQ

train_<N>.fastq.gz

subsample used for NanoSim profiling

nanosim_training*

NanoSim model/profile outputs

sim_pool.fastq.gz

large pool of simulated pathogen reads

mix_rep<rep>_n<spike_n>/

per spike level and replicate subfolder containing:

mixed.fastq.gz

kraken.report.tsv

kraken.classifications.tsv

kraken.summary.tsv

minimap_sorted_nodup.bam

minimap_sorted_nodup.MASKED.bed

spikein_summary.tsv

a single TSV summary table across all spike levels and replicates

The summary TSV columns are:

replicate

spike_n

mixed_fastq_gz

kraken_report

kraken_classified_reads

kraken_target_reads

minimap_bam

minimap_target_alignments

Negative control: shuffled genome spike-in

This run shuffles the pathogen genome sequence, simulates reads from the shuffled genome using the same NanoSim model approach, and repeats the spike series. It is intended to test for off-target detection artefacts.


conda activate spikein_ont

export REAL_FASTQ="/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC0123_AmM001WB.fastq.gz"
export HOST_REF_FASTA="/path/to/your/host_reference.fasta"
export PATHOGEN_FASTA="/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/genomes/GCF_000524495.1_Plas_inui_San_Antonio_1_V1_genomic.fasta"
export MINIMAP_DB_GZ="/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz"
export KRAKEN_DB_DIR="/home/pthorpe001/data/project_back_up_2024/kracken"

export OUT_DIR="spikein_shuffled_control_out"

bash run_spikein_shuffled_control.sh