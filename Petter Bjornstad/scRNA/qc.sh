#!/bin/bash
cd /mnt/HD2/Tim/scRNA
# Run FastQC on all fastq files
for file in $(find . -name "*.fastq.gz"); do
    SAMPLE=$(basename $file)
    echo Data_Raw/${SAMPLE}
    # fastqc -t 16 Data_Raw/${SAMPLE} -o /mnt/HD2/Tim/scRNA/Data_Clean/FastQC
done
# MultiQC to put everything together
multiqc /Data_Clean/FastQC
