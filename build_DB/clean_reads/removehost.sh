#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V
#$ -pe smp 12
#$ -jc long
#$ -mods l_hard mfree 100G
#$ -adds l_hard h_vmem 100G


cd /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity

minimap2 -t 12 -a -x map-ont \
  /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/host_plus_plasmo_depletion.fasta \
  /home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC1023_AmM008WB.fastq.gz \
  | samtools view -h \
      -T /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/host_plus_plasmo_depletion.fasta \
      -b -f 4 -@ 12 - \
  | samtools fastq -@ 12 - \
  | gzip -c \
  > /home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/data/MRC1023_AmM008WB.depleted.fastq.gz


  # Also
python remove_reads_by_kraken.py --kraken_classifications kraken.classifications.tsv \
   --kraken_report kraken.report.tsv --input_fastq  MRC1023_AmM008WB.depleted.strict.fastq \
   --output_fastq background_depleted.no_plasmodium.fastq --read_ids_out \
     plasmodium_read_ids.txt --summary_tsv removed_plasmodium_summary.tsv




