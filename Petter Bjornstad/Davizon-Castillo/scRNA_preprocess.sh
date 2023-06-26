#!/bin/bash
# cd to top project folder
cd /mnt/HD2/Davizon-Castillo
# Run FastQC on all fastq files
for file in $(find Data_Raw/fastq -name "*.fastq.gz"); do
    SAMPLE=$(basename $file)
    fastqc -q -t 16 Data_Raw/fastq/${SAMPLE} -o /mnt/HD2/Davizon-Castillo/Data_Clean/QC/FastQC
done
# MultiQC to put everything together
python3 -m multiqc Data_Clean/QC/FastQC -o Data_Clean/QC --no-data-dir
# Prepare STAR - no need to run every time, just want the code saved somewhere
# STAR  --runMode genomeGenerate \
#     --runThreadN 16 \
#     --genomeDir /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/RNAseq/STAR \
#     --genomeFastaFiles /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/RNAseq/GRCh38.p13.genome.fa \
#     --sjdbGTFfile /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/RNAseq/gencode.v43.primary_assembly.annotation.gtf
# Alignment with STARsolo
STAR --genomeDir ./Miscellaneous/STAR \
    --readFilesIn $(find Data_Raw/fastqs/ -name "*_R1_*.fastq.gz" | tr '\n' ',') $(find Data_Raw/fastqs/ -name "*_R2_*.fastq.gz" | tr '\n' ',') \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix ./Data_Clean/Mapped/\
    --runThreadN 16