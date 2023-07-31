#!/bin/bash
# cd to top project folder
cd /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/scRNA
# Run FastQC on all fastq files
for file in $(find data_raw/fastqs -name "*.fastq.gz"); do
    SAMPLE=$(basename $file)
    fastqc -t 16 data_raw/fastqs/${SAMPLE} -o data_clean/qc/fastqc
done
# MultiQC to put everything together
python3 -m multiqc data_clean/qc/fastqc -o data_clean/qc --no-data-dir --interactive
# # Prepare STARsolo to mimic Cell Ranger as close as possible - don't need to run every time
# STAR  --runMode genomeGenerate \
#     --runThreadN 16 \
#     --genomeDir ./Miscellaneous/STAR\ Genome \
#     --genomeSAindexNbases 10 \
#     --genomeFastaFiles /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/scRNA/Miscellaneous/cellranger-7.1.0/external/cellranger_tiny_ref/fasta/genome.fa \
#     --sjdbGTFfile /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/scRNA/Miscellaneous/cellranger-7.1.0/external/cellranger_tiny_ref/genes/genes.gtf
# Alignment with STARsolo
STAR --genomeDir ./Miscellaneous/STAR\ Genome \
    --readFilesIn $(find data_raw/fastqs/ -name "*_R1_*.fastq.gz" | tr '\n' ',') $(find data_raw/fastqs/ -name "*_R2_*.fastq.gz" | tr '\n' ',') \
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
    --outFileNamePrefix ./data_clean/mapped/ \
    --outTmpDir /mnt/HD2/Bjornstad/scRNA/STAR_temp \
    --runThreadN 16