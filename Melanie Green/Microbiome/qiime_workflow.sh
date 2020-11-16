#!/bin/bash
# Start QIIME
conda activate
source activate qiime2-2019.4
# Change directory
cd /Volumes/som/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Laura\ Tim\ projects/Melanie\ Green/Microbiome/Data_Cleaned

# Import sample data to QIIME format
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.txt \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path demux.qza

# Summarize quality
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

# De-noise with DADA2 - tried several different lengths, 240 was best
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 240 \
  --o-table denoised/table240.qza \
  --o-representative-sequences denoised/rep-seqs240.qza \
  --o-denoising-stats denoised/denoising-stats240.qza \
  --p-n-threads 0 \
  --verbose

  # Trimming at 220 is too short. Many of the samples have a feature count of 0.
  # 260 looks better, but it might be good to be slightly more conservative.
  # 240 looks pretty good to me (Tim)!

# Visualize
qiime metadata tabulate \
  --m-input-file denoised/denoising-stats240.qza \
  --o-visualization denoised/stats-dada240.qzv

# Remove samples with no metadata
qiime feature-table filter-samples \
  --i-table denoised/table240.qza \
  --m-metadata-file clinical.txt \
  --o-filtered-table filtered/table240.qza

# Feature tables
qiime feature-table summarize \
  --i-table filtered/table240.qza \
  --o-visualization feature_tables/table240.qzv \
  --m-sample-metadata-file clinical.txt

qiime feature-table tabulate-seqs \
  --i-data denoised/rep-seqs240.qza \
  --o-visualization denoised/rep-seqs240.qzv

# Phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences denoised/rep-seqs240.qza \
  --o-alignment alignment/aligned-rep-seqs.qza \
  --o-masked-alignment alignment/masked-aligned-rep-seqs.qza \
  --o-tree trees/unrooted-tree.qza \
  --o-rooted-tree trees/rooted-tree.qza

# Diversity analysis - not normalized
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny trees/rooted-tree.qza \
  --i-table filtered/table240.qza \
  --p-sampling-depth 3059 \
  --m-metadata-file clinical.txt \
  --output-dir core-metrics-results-original

# A sampling depth of 3059 was chosen based on the number of sequences in the
# PCOS63 sample because it’s close to the number of sequences in the next few samples
# that have higher sequence counts, and because it is considerably higher (relatively)
# than the number of sequences in the samples that have fewer sequences

# Alpha diversity
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-original/faith_pd_vector.qza \
  --m-metadata-file clinical.txt \
  --o-visualization core-metrics-results-original/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-original/evenness_vector.qza \
  --m-metadata-file clinical.txt \
  --o-visualization core-metrics-results-original/evenness-group-significance.qzv

# Beta Diversity
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-original/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file clinical.txt \
  --m-metadata-column SEQRun \
  --o-visualization core-metrics-results-original/unweighted-unifrac-seq-run-significance.qzv \
  --p-pairwise

# Feature classifier
qiime feature-classifier classify-sklearn \
    --i-reads denoised/rep-seqs240.qza \
    --i-classifier /Volumes/som/PEDS/RI\ Biostatistics\ Core/Shared/Shared\ Projects/Laura/Laura\ Tim\ projects/General\ Microbiome/silva-132-99-nb-classifier.qza \
    --o-classification taxa/taxa.qza

# Taxa barplot
qiime taxa barplot \
    --i-table filtered/table240.qza \
    --i-taxonomy taxa/taxa.qza \
    --m-metadata-file clinical.txt \
    --o-visualization taxa/taxa_barplot.qzv

# Relative frequency table
qiime feature-table relative-frequency \
    --i-table filtered-table.qza \
    --o-relative-frequency-table relative.qza

# Compositional analysis
# Correlation-clustering
qiime gneiss correlation-clustering \
  --i-table filtered-table.qza \
  --o-clustering hierarchy.qza
# ILR transform
qiime gneiss ilr-hierarchical \
  --i-table filtered-table.qza \
  --i-tree hierarchy.qza \
  --o-balances balances.qza
# Linear model
qiime gneiss ols-regression \
  --p-formula "Group+SEQRun" \
  --i-table balances.qza \
  --i-tree hierarchy.qza \
  --m-metadata-file clinical.txt \
  --o-visualization regression_summary.qzv

# Not normalized
# Correlation-clustering
qiime gneiss correlation-clustering \
  --i-table pnorm_freq.qza \
  --o-clustering hierarchy_norm.qza
# ILR transform
qiime gneiss ilr-hierarchical \
  --i-table pnorm_freq.qza \
  --i-tree hierarchy_norm.qza \
  --o-balances balances_norm.qza
# Linear model
qiime gneiss ols-regression \
  --p-formula "Group+SEQRun" \
  --i-table balances_norm.qza \
  --i-tree hierarchy_norm.qza \
  --m-metadata-file clinical.txt \
  --o-visualization regression_summary_norm.qzv
# 10-fold validation is used to check overfitting. If prediction accuracy (pred_mse) < model error (mse) overfitting is unlikely.
# Heatmap is of coefficient p-values for each balance. Columns are balances, rows are covariates.
# For prediction and residual plots, only the top two balances are shown.
# In the SS tree, the longest branch corresponds to the most informative balance.
# If any balances appear to be significant, pass them to qiime gneiss balance-taxonomy

# Normalized
# Show balances in a heatmap
qiime gneiss dendrogram-heatmap \
  --i-table pnorm_freq.qza \
  --i-tree hierarchy_norm.qza \
  --m-metadata-file clinical.txt \
  --m-metadata-column Group \
  --p-color-map seismic \
  --o-visualization heatmap_norm.qzv
# Which taxa might explain group-wise differences?
qiime gneiss balance-taxonomy \
  --i-table pnorm_freq.qza \
  --i-tree hierarchy_norm.qza \
  --i-taxonomy taxa.qza \
  --p-taxa-level 2 \
  --p-balance-name 'y0' \
  --m-metadata-file clinical.txt \
  --m-metadata-column Group \
  --o-visualization taxa_summary_norm.qzv

  # Not normalized
  # Show balances in a heatmap
qiime gneiss dendrogram-heatmap \
  --i-table filtered-table.qza \
  --i-tree hierarchy.qza \
  --m-metadata-file clinical.txt \
  --m-metadata-column Group \
  --p-color-map seismic \
  --o-visualization heatmap.qzv
# Which taxa might explain group-wise differences?
qiime gneiss balance-taxonomy \
  --i-table filtered-table.qza \
  --i-tree hierarchy.qza \
  --i-taxonomy taxa.qza \
  --p-taxa-level 2 \
  --p-balance-name 'y0' \
  --m-metadata-file clinical.txt \
  --m-metadata-column Group \
  --o-visualization taxa_summary.qzv





# To do
# - In the regression model from normalized data, why does group explain more variability than in the model using non-normalized?
