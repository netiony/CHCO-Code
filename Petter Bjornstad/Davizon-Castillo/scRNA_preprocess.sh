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

# Prepare RSEM reference genome for alignment - no need to run every time, just want the code saved somewhere
# cd /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Peds\ Endo/Petter\ Bjornstad/RNAseq
# rsem-prepare-reference --star --star-path /home/tim/Tools/STAR/source -p 20 \
#     --gtf gencode.v43.primary_assembly.annotation.gtf \
#     GRCh38.p13.genome.fa \
#     /mnt/HD2/Davizon-Castillo/reference_genome/ref
# cd /mnt/HD2/Davizon-Castillo

# Calculate expression levels
rsem-calculate-expression -p 20 \
    --star-gzipped-read-file \
    --paired-end \
	--star \
    --star-path /home/tim/Tools/STAR/source \
	$(find data_raw/fastqs/ -name "*_R1_*.fastq.gz" | sort | tr '\n' ',') \
    $(find data_raw/fastqs/ -name "*_R2_*.fastq.gz" | sort | tr '\n' ',') \
	reference_genome/ref \
    data_clean/rsem_out