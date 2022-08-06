#!/bin/bash
# Move to base folder
cd ~/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Vigers/sheep_miRNA
# Use bowtie to index reference genome
bowtie-build data_clean/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna data_clean/indexed_genome.fa
# Collapse trimmed files
for file in $(find data_clean/FASTQ_trimmed -name "*gz" -type f); do
    mapper.pl file -e -h -i -j -l 18 -m -p refdb.fa -s reads_collapsed.fa -t reads_vs_refdb.arf -v -o 4
done
