#!/bin/bash
# Move to base folder
cd ~/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Vigers/sheep_miRNA
# QC
for file in data_clean/FASTQ/*; do
    fastqc ${file} -o data_clean/QC
done
