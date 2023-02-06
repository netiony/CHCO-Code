#!/bin/bash
cd /mnt/HD2/Tim/scRNA
# Run FastQC on all fastq files
for file in $(find . -name "*.fastq.gz"); do
    SAMPLE=$(basename $file)
    fastqc -t 32 Data_Raw/${SAMPLE} -o /mnt/HD2/Tim/scRNA/Data_Clean/FastQC
done
# MultiQC to put everything together
python3 -m multiqc Data_Clean/FastQC
