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

# # Prepare STAR genome to align to - no need to run every time, just want the code saved somewhere
# STAR  --runMode genomeGenerate \
#     --runThreadN 20 \
#     --genomeDir /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/RNAseq/STAR \
#     --genomeFastaFiles /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/RNAseq/GRCh38.p13.genome.fa \
#     --sjdbGTFfile /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/RNAseq/gencode.v43.primary_assembly.annotation.gtf

# Load index into memory for faster processing
STAR --genomeLoad LoadAndExit --genomeDir ./miscellaneous/STAR
# Alignment with STAR
for file in $(find data_raw/fastqs -name "*_R1_*.fastq.gz"); do
    name="$(basename -- $file)"
    name=${name/_R1*/""}
    STAR --genomeDir ./miscellaneous/STAR \
        --genomeLoad LoadAndKeep \
        --readFilesIn $file ${file/R1/"R2"} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix "./data_clean/aligned/${name}"\
        --quantMode GeneCounts \
        --runThreadN 20
done
# Remove genome from memory
STAR --genomeLoad Remove --genomeDir ./miscellaneous/STAR