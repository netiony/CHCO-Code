#!/bin/bash
# Move to base folder
cd /mnt/HD1/Tim/sheep_miRNA
# Copy fastq files to data folder
find data_raw -name \*.fastq.gz -exec cp -n {} data_clean/fastq \;
# Move to clean data folder
cd data_clean
# Index the sheep transcriptome
salmon index -t Ovis_aries_rambouillet.Oar_rambouillet_v1.0.cdna.all.fa.gz -i ovis_index
# Quantify reads
for fn in $(find fastq/* -type f | sed 's/_R._001\..*//' | sort | uniq);
do
samp=`basename ${fn}`
samp=${samp%_R1_001.fastq.gz}
samp=${samp%_R2_001.fastq.gz}
echo "Processing sample ${samp}"
salmon quant -i ovis_index -l A \
         -1 ${fn}/${samp}_R1_001.fastq.gz \
         -2 ${fn}/${samp}_R2_001.fastq.gz \
         -p 12 --validateMappings --gcBias -o quants/${samp}_quant
done 