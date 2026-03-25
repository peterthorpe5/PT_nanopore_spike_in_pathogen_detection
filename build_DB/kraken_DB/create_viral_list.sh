


awk -F $'\t' '
NR>1 && $12=="Complete Genome" && $8!="na" && $5!="na" && $6!="na" && $20!="na" && $8 ~ /virus|Virus/ {
    print $20
}' ncbi_metadata/assembly_summary_refseq.txt > viral_ftp_paths.txt
