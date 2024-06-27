library(Seurat)
library(future)
library(tidyverse)
library(ggdendro)
library(ggpubr)
library(rstatix)
library(GSVA)
library(Biobase)
library(ReactomeGSA)
library(GSEABase)
library(msigdbr)

# Parallel processing
plan(multicore, workers = 16)
options(future.globals.maxSize=6e9)
# Import
so <- readRDS("/run/user/1023/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/scRNA/data_clean/seurat_data_no_computations.RDS")
# Exclude control with IgA
so <- subset(so, T2D_HC_Phil != "HC_igA")
so$Group <- so$T2D_HC_Phil
# Normalize and scale
so <- NormalizeData(so)
so <- ScaleData(so)
# PCA
so <- RunPCA(so, features = VariableFeatures(object = so))
# ElbowPlot(so)
# Cluster cells
so <- FindNeighbors(so)
so <- FindClusters(so)
# Perform UMAP and tSNE
so <- RunUMAP(so, dims = 1:30)
# so = RunTSNE(so,dim.embed = 3)
# General cell types as identifiers
so$generaltype <- sub("_.*", "", so$LR_clusters)
Idents(so) <- so$LR_clusters
# Tubular cells
tubular_cells <- c("PT", "TAL", "PC", "ATL", "IC", "DCT", "CNT")
so$celltype_tubular <- ifelse(so$LR_clusters %in% tubular_cells, "tubular", "non-tubular")
so$tubular_id <- paste0(so$celltype_tubular, "_", so$michigan_id)

saveRDS(so, "/home/yejichoi2/Documents/seurat_data_compiled.RDS")