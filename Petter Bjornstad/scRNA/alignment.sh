#!/bin/bash
# This code needs to be run in an environment where STAR is installed.
# The easiest way to do this is with conda (conda activate scRNA-seq).
cd /Volumes/Work/scRNA
# Download reference data
# wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz -P Miscellaneous
# Use STAR to generate genome index
STAR  --runMode genomeGenerate --runThreadN 4 \
    --genomeDir Miscellaneous/indexed_reference \
    --genomeFastaFiles Miscellaneous/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
    --sjdbGTFfile Miscellaneous/refdata-gex-GRCh38-2020-A/genes/genes.gtf
# STARsolo for alignment (per Phil's recommendation)
# See instructions at https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
STAR --genomeDir /path/to/genome/dir/ --readFilesIn ...  --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist ~/GitHub/cellranger/lib/python/cellranger/barcodes/3M-february-2018.txt