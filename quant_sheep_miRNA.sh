#!/bin/bash
# Move to folder
cd /mnt/HD1/Tim/sheep_miRNA/
# Index the sheep transcriptome
salmon index -t Ovis_aries_rambouillet.Oar_rambouillet_v1.0.cdna.all.fa.gz -i ovis_index
# Quantify reads
for fn in data/DRR0161{25..40};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i athal_index -l A \
         -1 ${fn}/${samp}_1.fastq.gz \
         -2 ${fn}/${samp}_2.fastq.gz \
         -p 12 --validateMappings --gcBias -o data_clean/${samp}_quant
done 