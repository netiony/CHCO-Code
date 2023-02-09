#!/bin/bash
# This code needs to be run in an environment where STAR is installed.
# The easiest way to do this is with conda (conda activate scRNA-seq).
cd /mnt/HD2/Tim/scRNA
# # Download reference data
# wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz -P Miscellaneous
# tar -xvzf Miscellaneous/refdata-gex-GRCh38-2020-A.tar.gz -C Miscellaneous
# # Use STAR to generate genome index
# STAR \
#     --runMode genomeGenerate \
#     --runThreadN 36 \
#     --genomeSAsparseD 3 \
#     --genomeDir Miscellaneous/indexed_reference \
#     --genomeFastaFiles Miscellaneous/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
#     --sjdbGTFfile Miscellaneous/refdata-gex-GRCh38-2020-A/genes/genes.gtf
# STARsolo for alignment (per Phil's recommendation)
# See instructions at https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
# List files in STARsolo format
R2=$(find Data_Raw/ -type f -name "*_R2_*" | tr '\n' ',')
R1=$(find Data_Raw/ -type f -name "*_R1_*" | tr '\n' ',')
# Run everything with STAR
STAR \
    --genomeDir Miscellaneous/indexed_reference \
    --outFileNamePrefix Data_Clean/Aligned/ \
    --readFilesIn $R2 $R1 \
    --readFilesCommand zcat \
    --soloUMIlen 12 \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist Miscellaneous/3M-february-2018.txt \
    --runThreadN 36
# Sync everything to the shared drive
rsync -av --delete /mnt/HD2/Tim/scRNA/ /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/scRNA