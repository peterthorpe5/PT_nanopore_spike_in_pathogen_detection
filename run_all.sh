

conda activate minimap   # or whichever env has minimap2/samtools/kraken2/nanosim

export REAL_FASTQ="/home/pthorpe001/data/project_back_up_2024/jcs_blood_samples/MRC0123_AmM001WB.fastq.gz"
export HOST_REF_FASTA="/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/<YOUR_HOST_FASTA>.fasta"
export PATHOGEN_FASTA="/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/genomes/GCF_000524495.1_Plas_inui_San_Antonio_1_V1_genomic.fasta"

export MINIMAP_DB_GZ="/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta.gz"
export KRAKEN_DB_DIR="/home/pthorpe001/data/project_back_up_2024/kracken"

bash run_spikein_one_sample.sh
