#!/bin/bash
# This code needs to be run in an environment where STAR is installed.
# The easiest way to do this is with conda.
cd /Volumes/Work/scRNA
# STARsolo for alignment (per Phil's recommendation)
# See instrcutions at https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
STAR --genomeDir /path/to/genome/dir/ --readFilesIn ...  --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist ~/GitHub/cellranger/lib/python/cellranger/barcodes/3M-february-2018.txt