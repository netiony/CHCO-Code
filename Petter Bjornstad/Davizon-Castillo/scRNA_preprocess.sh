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
# Prepare STARsolo to mimic Cell Ranger as close as possible
STAR  --runMode genomeGenerate \
    --runThreadN 16 \
    --genomeDir ./Miscellaneous/STAR\ Genome \
    --genomeSAindexNbases 10 \
    --genomeFastaFiles /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/scRNA/Miscellaneous/cellranger-7.1.0/external/cellranger_tiny_ref/fasta/genome.fa \
    --sjdbGTFfile /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/scRNA/Miscellaneous/cellranger-7.1.0/external/cellranger_tiny_ref/genes/genes.gtf
# Alignment with STARsolo
STAR --genomeDir ./Miscellaneous/STAR\ Genome \
    --readFilesIn $(find Data_Raw/fastq/ -name "*_R1_*.fastq.gz" | tr '\n' ',') $(find Data_Raw/fastq/ -name "*_R2_*.fastq.gz" | tr '\n' ',') \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist /home/tim/GitHub/cellranger/lib/python/cellranger/barcodes/737K-august-2016.txt \
    --soloCBstart 1 \
    --soloCBlen 16 \
    --soloUMIstart 17 \
    --soloUMIlen 10 \
    --soloBarcodeMate 1 \
    --clip5pNbases 39 0 \
    --readFilesCommand gunzip -c \
    --genomeSAsparseD 3 \
    --soloMultiMappers EM \
    --outFileNamePrefix ./Data_Clean/Mapped/\
    --runThreadN 16