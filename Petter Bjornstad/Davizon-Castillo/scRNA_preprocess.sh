#!/bin/bash
# cd to top project folder
cd /mnt/HD2/Davizon-Castillo

# # Run FastQC on all fastq files
# for file in $(find data_raw/fastqs -name "*.fastq.gz"); do
#     SAMPLE=$(basename $file)
#     fastqc -q -t 16 data_raw/fastqs/${SAMPLE} -o /mnt/HD2/Davizon-Castillo/data_clean/qc/FastQC
# done

# # MultiQC to put everything together
# python3 -m multiqc data_clean/QC/FastQC -o data_clean/QC --no-data-dir

# # Prepare STAR - no need to run every time, just want the code saved somewhere
# STAR  --runMode genomeGenerate \
#     --runThreadN 16 \
#     --genomeDir /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/RNAseq/STAR \
#     --genomeFastaFiles /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/RNAseq/GRCh38.p13.genome.fa \
#     --sjdbGTFfile /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/RNAseq/gencode.v43.primary_assembly.annotation.gtf

# Alignment with STARsolo
for file in $(find data_raw/fastqs -name "*_R1_*.fastq.gz"); do
    echo $file
    STAR --genomeDir ./miscellaneous/STAR \
        --readFilesIn $file ${file/R1/"R2"} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix ./data_clean/mapped/\
        --runThreadN 16
done