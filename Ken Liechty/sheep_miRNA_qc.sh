#!/bin/bash
# Move to base folder
cd ~/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Vigers/sheep_miRNA
# QC
for file in data_clean/FASTQ/*; do
    fastqc ${file} -t 32 -o data_clean/QC
done
# Combine
multiqc data_clean/QC -o data_clean/QC
# Trim adapters
cd data_clean/FASTQ
for fn in $(find . -name "*_R1_001.fastq.gz" -type f | sed 's/_R._001\..*//' | sort | uniq); do
    trim_galore --length 0 -cores 8 --no_report_file --paired ${fn}_R1_001.fastq.gz ${fn}_R2_001.fastq.gz
done
# Move to new folder
cd ~/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Vigers/sheep_miRNA
find data_clean/FASTQ -name '*trimm*' -exec mv {} data_clean/FASTQ_trimmed \;
find data_clean/FASTQ -name '*_val_*' -exec mv {} data_clean/FASTQ_trimmed \;
# QC again
for file in $(find data_clean/FASTQ_trimmed -name "*gz" -type f); do
    fastqc ${file} -t 32 -o data_clean/QC_trimmed
done
# Combine
multiqc data_clean/QC_trimmed -o data_clean/QC_trimmed
