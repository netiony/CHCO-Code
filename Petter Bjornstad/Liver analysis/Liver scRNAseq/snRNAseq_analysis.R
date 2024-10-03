#snRNAseq Analysis ----
#1. Load Libraries ----
library(tidyverse)
library(BiocManager)
library(arsenal)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(future)
library(colorspace)
library(patchwork)
library(ggdendro)
library(cowplot)
library(ggpubr)
library(rstatix)
library(table1)
library(Biobase)
library(ReactomeGSA)
library(GSEABase)
library(msigdbr)
library(kableExtra)
library(knitr)

#2. Set up Directories ----
dir.data <- c("/Users/hhampson/Dropbox/Bjornstad data/Liver")
dir.code <- c("/Users/hhampson/Documents/CHCO-Code/Petter Bjornstad/Liver analysis/Liver scRNAseq")

#3. Format Saurat Object for Single Nuc analysis ----
#Load liver sn saurat object
so_liver_sn <- readRDS(fs::path(dir.data,"NoRef_PetterLiver_ClinData_Labels_Char_041924.RDS"))

#Load missing meta data
meta_liver_raw <- read.csv(fs::path(dir.data,"liver_biopsy_metadata_PN.csv"))

#Add missing liver meta data to so
meta_liver_sn <-  so_liver_sn@meta.data[,1:11] %>%
  dplyr::mutate(RNAlater_ID = SampleID) %>%
  left_join(meta_liver_raw)

rownames(meta_liver_sn) <- rownames(so_liver_sn@meta.data)
so_liver_sc <- AddMetaData(so_liver_sc, meta_liver_sc)

#4. Run Analyses ----
# Set a higher memory limit
# mem.maxVSize(64000000000)

#Normalize & Scale the data
so_liver_sn_rna <- NormalizeData(so_liver_sn@assays$RNA)
so_liver_sn_rna <- ScaleData(so_liver_sn@assays$RNA)

so_liver_sn_int <- NormalizeData(so_liver_sn@assays$integrated)
so_liver_sn_int <- ScaleData(so_liver_sn@assays$integrated)

# PCA
so_liver_sn_int <- RunPCA(so_liver_sn_int, features = VariableFeatures(object = so_liver_sn_int))
ElbowPlot(so_liver_sn_int)
# Cluster cells
so_liver_sn_int <- FindNeighbors(so_liver_sn_int)
so_liver_sn_int <- FindClusters(so_liver_sn_int)
# Perform UMAP and tSNE
so_liver_sn <- RunUMAP(so_liver_sn, dims = 1:15)
DimPlot(so_liver_sn, reduction = "umap") 





#ssh hampsonh@10.45.122.105
#Yo8EJh2iSHhyT8hyXjN7
