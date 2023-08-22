#!/bin/bash
cd /mnt/HD2/Bjornstad/scRNA/data_clean
cellranger count --id=cellranger_count \
   --fastqs=/mnt/HD2/Bjornstad/scRNA/data_raw/fastqs \
   --sample=2327-EO-1_TGTCCCAA-TGGACATC \
   --transcriptome=/home/tim/cellranger/refdata-gex-GRCh38-2020-A \
   --no-bam \
   --localcores 16
# Sync to shared drive as a backup
rsync -auv /mnt/HD2/Bjornstad/scRNA /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad --dry-run
rsync -auv /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/scRNA /mnt/HD2/Bjornstad --dry-run
