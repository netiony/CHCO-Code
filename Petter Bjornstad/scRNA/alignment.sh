#!/bin/bash
cd /Volumes/Work/scRNA
# STARsolo for alignment (per Phil's recommendation)
STAR --genomeDir /path/to/genome/dir/ --readFilesIn ...  --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist ~/GitHub/cellranger/lib/python/cellranger/barcodes/3M-february-2018.txt