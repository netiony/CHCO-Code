#!/bin/bash
# Move to base folder
cd /home/tim/UCD/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Vigers/sheep_miRNA
# List of genome targets (decoys)
grep "^>" <(gunzip -c data_clean/Ovis_aries.Oar_v3.1.dna.toplevel.fa.gz) | cut -d " " -f 1 >data_clean/decoys.txt
sed -i.bak -e 's/>//g' data_clean/decoys.txt
cat data_clean/Ovis_aries.Oar_v3.1.cdna.all.fa.gz data_clean/Ovis_aries.Oar_v3.1.dna.toplevel.fa.gz >data_clean/gentrome.fa.gz
# Index the sheep transcriptome
salmon index -t data_clean/gentrome.fa.gz -d data_clean/decoys.txt -p 24 -i data_clean/ovis_index
# Quantify reads
for fn in $(find data_raw/tanner_all_new -name "*_R1_001.fastq.gz" -type f | sed 's/_R._001\..*//' | sort | uniq); do
    samp=$(basename ${fn})
    echo "Processing sample ${samp}"
    salmon quant -i data_clean/ovis_index -l A \
        -1 ${fn}_R1_001.fastq.gz \
        -2 ${fn}_R2_001.fastq.gz \
        -p 24 --validateMappings --gcBias -o data_clean/quants/${samp}_quant
done
