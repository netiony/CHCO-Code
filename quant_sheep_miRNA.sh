#!/bin/bash
# Move to base folder
cd /mnt/HD1/Tim/sheep_miRNA
# Index the sheep transcriptome
salmon index -t data_clean/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.cdna.all.fa.gz -i data_clean/ovis_index
# Quantify reads
for fn in $(find data_raw/tanner_all_new -name "*_R1_001.fastq.gz" -type f | sed 's/_R._001\..*//' | sort | uniq);
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i data_clean/ovis_index -l A \
         -1 ${fn}_R1_001.fastq.gz \
         -2 ${fn}_R2_001.fastq.gz \
         -p 24 --validateMappings --gcBias -o data_clean/quants/${samp}_quant
done