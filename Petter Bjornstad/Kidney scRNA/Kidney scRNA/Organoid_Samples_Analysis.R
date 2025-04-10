#DEG Analysis for Organoid Samples 
#1.Set Up ----
###a. Load Packages ----
library(reprex)
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
library(venn)
library(rstatix)
library(table1)
library(Biobase)
library(ReactomeGSA)
library(GSEABase)
library(msigdbr)
library(kableExtra)
library(knitr)
library(SingleCellExperiment)
library(fgsea)
library(EnhancedVolcano)
library(openxlsx)
library(BiocManager)
library(MAST)
library(ggrepel)
# library(qpcR)
library(ggpubr)
library(openxlsx)
library(ggplot2)
library(GGally)
library(GSEABase)
library(limma)
library(reshape2)
library(data.table)
library(knitr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
#library(NMF)
library(rsvd)
library(RColorBrewer)
library(MAST)
library(devtools)
# install_github("Sun-lab/ideas",force=T)
#library(ideas)
library(foreach)
library(parallel)
library(doRNG)
library(doParallel)
library(fs)
# registerDoParallel(cores = 6)
library(VennDiagram)
library(janitor)
# devtools::install_github('immunogenomics/presto')
library(presto)
library(knitr)
library(lme4)
library(lmerTest)
#install.packages("glmmTMB")
# Reinstall glmmTMB from source
# install.packages("glmmTMB", type = "source")
library(glmmTMB)
# Install DoubletFinder (if not already installed)
# devtools::install_github("chris-mcginnis-ucsf/DoubletFinder",force=T)
# Load the package
# Install DoubletFinder from GitHub (use devtools to install)
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("chris-mcginnis-ucsf/DoubletFinder",force=T)
# library(DoubletFinder)
# install.packages("emmeans")
library(emmeans)

##b. Set up directories ----
dir.dat <- c("//Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive")
dir.dat2 <- c("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/scRNA")
dir.code <- c("/Users/hhampson/Documents/CHCO-Code/Petter Bjornstad/Liver analysis/Liver scRNAseq")
dir.results <- c("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Kidney scRNAseq Project/Organoid Results")

##c. Load functions ----
source("Kidney_functions_sc.R")

#Set ids for organoid samples
ids <- c("CRC-10","CRC-11","CRC-03","RH-50-T","RH-72-T","RH-62-T")

##d. Load Data ----
so_kpmp_sc <- readRDS(fs::path(dir.dat2,"data_raw","PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds"))
so_kpmp_sc <- subset(so_kpmp_sc,record_id %in% ids)

# #Fix Typos in kit ids in PB90
# so_kpmp_sc$kit_id[which(so_kpmp_sc$kit_id=="KI-0014643")] <- "KL-0014643"
# so_kpmp_sc$kit_id[which(so_kpmp_sc$kit_id=="kl-0023998")] <- "KL-0023998"

#Load harmonized data that has been filtered from 90 to the 83 participants that have baseline single cell data
harm_meta_data <- read.csv(fs::path(dir.dat,"Kidney scRNAseq Project","Data","harmonized_data_kidney_sc_all_metadata2.csv"))

##e. Filter to 6 organoid samples ----
meta <- harm_meta_data %>% 
dplyr::select(-X) %>% 
  filter(record_id %in% ids)
rm(harm_meta_data)
#Select metadata from seurat object to facilitate merge of new metadata into seurat object
meta_kidney_sc <-  so_kpmp_sc@meta.data
rownames(meta_kidney_sc) <- rownames(so_kpmp_sc@meta.data)

#Merge metadata from 83 participants at baseline into seurat object metadata
meta_kidney_sc <- meta_kidney_sc %>%
  left_join(meta,by="kit_id")
rownames(meta_kidney_sc) <- rownames(so_kpmp_sc@meta.data)

#Pull ids from 83 participants at baseline to filter seurat object to these participants only 
ids <- meta$kit_id

#Merge metadata back into seurat object
so_kpmp_sc <- AddMetaData(so_kpmp_sc, meta_kidney_sc)

#Check number of unique ids
length(unique(so_kpmp_sc$kit_id)) #should be 6

#Remove metadatasets
rm(meta_kidney_sc,meta)

##f. Merge metadata into filtered seurat object ----
#Load in most up to date medication data to update medication information
med <- read.xlsx(fs::path(dir.dat,"Kidney scRNAseq Project/Data/Biopsies_w_mrn_Oct3.xlsx"))
#Select Metformin, RASSI, Insulin data
med <- med %>%
  dplyr::select(all_of(c("record_id","mrn","raasi_1","insulin_1","mfm_1")))
#Pull seurat object metadata to help harmoinize in new metadata
meta_kidney_sc <-  so_kpmp_sc@meta.data
#Filter to only those with a unique identifier id in the seurat object metadata
med <- med %>%
  filter(mrn %in% as.character(meta_kidney_sc$mrn)) 
length(unique(med$mrn)) #6 participants
#Filter to only those that have a unique record id in the seurat object
med <- med %>%
  filter(record_id %in% meta_kidney_sc$record_id) 
length(unique(med$mrn)) #6 remain
length(unique(med$record_id)) #6
rownames(meta_kidney_sc) <- rownames(so_kpmp_sc@meta.data)
med$mrn <- as.numeric(med$mrn) #Make numeric to merge
#Merge med data with seurat metadata
meta_kidney_sc <- meta_kidney_sc %>%
  left_join(med,by=c("mrn","record_id"))
rownames(meta_kidney_sc) <- rownames(so_kpmp_sc@meta.data)
length(unique(meta_kidney_sc$mrn)) #6 remain
length(unique(meta_kidney_sc$record_id)) #6

#Add Med Meta Data to Seurat object
so_kpmp_sc <- AddMetaData(so_kpmp_sc, meta_kidney_sc)
#Remove med metadatset
rm(med,meta_kidney_sc)

#Create medication & disease status groups of interest
so_kpmp_sc@meta.data <- so_kpmp_sc@meta.data %>%
  mutate(glp1_sglt2=ifelse(epic_glp1ra_1=="Yes" & epic_sglti2_1=="Yes","Yes","No")) %>%
  mutate(sglt2=ifelse(epic_sglti2_1=="Yes" & epic_glp1ra_1=="No","Yes","No")) %>%
  mutate(glp1=ifelse(epic_sglti2_1=="No" & epic_glp1ra_1=="Yes","Yes","No")) %>%
  mutate(no_med=ifelse(epic_sglti2_1=="No" & epic_glp1ra_1=="No","Yes","No"))

#Define 4 exposure groups:
#SGLT2i(+)/GLP1RA(+), SGLT2i(+)/GLP1RA(-), SGLT2i(-)/GLPRA(+), SGLT2i(-)/GLP1RA(-)
invisible(gc())
so_kpmp_sc@meta.data <- so_kpmp_sc@meta.data %>%
  mutate(medication = case_when(glp1_sglt2 == "Yes" ~ "glp1_sglt2",
                                sglt2 == "Yes" ~ "sglt2",
                                glp1 == "Yes" ~ "glp1",
                                no_med == "Yes" ~ "no_med"))
so_kpmp_sc@meta.data$medication <- factor(so_kpmp_sc@meta.data$medication, levels = c("no_med", "sglt2", "glp1","glp1_sglt2"))

#Ensure default assay in seurat object to RNA
DefaultAssay(so_kpmp_sc) <- "RNA"
invisible(gc())

# #Perform Quality Control & Preprocessing Steps
ncol(so_kpmp_sc) #11386 cells
nrow(so_kpmp_sc) #31332 genes
#YE JI's filtering code for percent expression 
#Filter out rare genes expressed in less than "gene_pct" of cells
# expr_matrix <- as.matrix(GetAssayData(so_kpmp_sc, layer = "counts"))
expr_matrix <- as.matrix(GetAssayData(so_kpmp_sc, assay = "RNA", layer = "counts"))
# expr_matrix <- as.matrix(so_kpmp_sc@assays$RNA@layers$counts)
# Calculate the proportion of cells expressing each gene
num_cells_per_gene <- rowSums(expr_matrix > 0)  # Count nonzero cells per gene
total_cells <- ncol(expr_matrix)  # Total number of cells
gene_proportion <- num_cells_per_gene / total_cells  # Fraction of cells expressing each gene
remove(expr_matrix)
# Keep genes expressed in at least "gene_pct" of cells
genes_to_keep <- names(gene_proportion[gene_proportion >= 0.05])
so_kpmp_sc <- subset(so_kpmp_sc, features = genes_to_keep)
ncol(so_kpmp_sc) #11386 nuclei
nrow(so_kpmp_sc) # 9209 genes

# Step 2: Remove mitochondrial genes (those starting with "MT")
mito_genes <- grep("^MT-", rownames(so_kpmp_sc), value = TRUE)
#keep_ids <- unique(rownames(so_kpmp_sc)[which(!rownames(so_kpmp_sc) %in% mito_genes)])
# so_kpmp_sc <- subset(so_kpmp_sc, features = setdiff(rownames(so_kpmp_sc), mito_genes))
# so_kpmp_sc <- subset(so_kpmp_sc,kit_id!="KL-0029535")
#so_kpmp_sc$Gene <- rownames(so_kpmp_sc)
#so_kpmp_sc <- subset(so_kpmp_sc, Gene %in% keep_ids)
so_kpmp_sc <- subset(so_kpmp_sc, features = setdiff(rownames(so_kpmp_sc), mito_genes))
#so_kpmp_sc <- subset(so_kpmp_sc, features = setdiff(rownames(so_kpmp_sc@assays$RNA@counts), mito_genes))
# so_kpmp_sc <- subset(so_kpmp_sc, !rownames(so_kpmp_sc) %in% mito_genes)
# grep("^MT-", rownames(so_kpmp_sc@assays$RNA@counts), value = TRUE)
# dim(so_kpmp_sc@assays$RNA@counts) #9276 186125
# dim(so_kpmp_sc@assays$RNA@data) #9276 186125
# dim(so_kpmp_sc@assays$RNA)#9276 186125
sum(grepl("^MT-", rownames(so_kpmp_sc))) #0
ncol(so_kpmp_sc) #11386 cells
nrow(so_kpmp_sc) #9196 genes

#Renormalize & Scale after filtering
so_kpmp_sc <- NormalizeData(so_kpmp_sc)
so_kpmp_sc <- ScaleData(so_kpmp_sc, features = VariableFeatures(so_kpmp_sc))

#Create general hepatocyte cell type variable
#Create PT and TAL pseudobulk cell type variable
so_kpmp_sc$celltype1 <- case_when(grepl("PT-",so_kpmp_sc$celltype_rpca)~"PT",
                                  grepl("TAL-",so_kpmp_sc$celltype_rpca)~"TAL",
                                  grepl("EC-",so_kpmp_sc$celltype_rpca)~"EC",
                                  grepl("POD",so_kpmp_sc$celltype_rpca)~"POD",
                                  grepl("MAC",so_kpmp_sc$celltype_rpca)~"MAC",
                                  grepl("MON",so_kpmp_sc$celltype_rpca)~"MON",
                                  grepl("PC-",so_kpmp_sc$celltype_rpca)~"PC",
                                  grepl("FIB",so_kpmp_sc$celltype_rpca)~"FIB_MC_VSMC",
                                  grepl("DTL",so_kpmp_sc$celltype_rpca)~"DTL",
                                  so_kpmp_sc$celltype_rpca=="DCT"~"DCT",
                                  so_kpmp_sc$celltype_rpca=="ATL"~"ATL",
                                  so_kpmp_sc$celltype_rpca=="B"~"B",
                                  so_kpmp_sc$celltype_rpca=="T"~"T")
so_kpmp_sc$celltype1 <- as.character(so_kpmp_sc$celltype1)

so_kpmp_sc$KPMP_celltype2 <- as.character(so_kpmp_sc$KPMP_celltype)
so_kpmp_sc$celltype2 <- ifelse(so_kpmp_sc$KPMP_celltype=="aPT" | 
                                 so_kpmp_sc$KPMP_celltype=="PT-S1/S2" | 
                                 so_kpmp_sc$KPMP_celltype == "PT-S3","PT",
                               ifelse(grepl("TAL",so_kpmp_sc$KPMP_celltype),"TAL",
                                      ifelse(grepl("EC-",so_kpmp_sc$KPMP_celltype),"EC",so_kpmp_sc$KPMP_celltype2)))

#2. Visualize ----
# PCA
so_kpmp_sc <- FindVariableFeatures(object = so_kpmp_sc)
so_kpmp_sc <- RunPCA(so_kpmp_sc, features = VariableFeatures(object = so_kpmp_sc),assay="RNA")
ElbowPlot(so_kpmp_sc)

# # Find neighbors and clusters (again using integrated data)
so_kpmp_sc <- FindNeighbors(so_kpmp_sc, assay = "RNA", dims = 1:20)
so_kpmp_sc <- FindClusters(so_kpmp_sc, resolution = 0.5)

# Perform UMAP and tSNE
# so_liver_sc <- RunUMAP(so_liver_sc, dims = 1:15)
so_kpmp_sc@reductions
DimPlot(so_kpmp_sc, reduction = "umap.harmony", raster = F)
DimPlot(so_kpmp_sc, reduction = "umap.harmony", group.by = "KPMP_celltype", raster = FALSE)
plot1 <- DimPlot(so_kpmp_sc, reduction = "umap.harmony", group.by = "KPMP_celltype",raster = FALSE, label = TRUE, label.size = 4)+
  ggtitle("Harmony UMAP Plot with Cell Type Labels")
# pdf(plot1,fs::path(dir.results,"UMAP_Celltypes_Organoid_Samples.pdf"),width=10,height=10)
# plot(plot1)
# dev.off()

so_kpmp_sc$group <- factor(so_kpmp_sc$group)
DimPlot(so_kpmp_sc, reduction = "umap.harmony",group.by = "group",label=F,raster=F) +
  ggtitle(paste0("UMAP by Disease Category"))

DimPlot(so_kpmp_sc, reduction = "umap.rpca", raster = F)
DimPlot(so_kpmp_sc, reduction = "umap.rpca", group.by = "KPMP_celltype", raster = FALSE)
DimPlot(so_kpmp_sc, reduction = "umap.rpca", group.by = "KPMP_celltype",raster = FALSE, label = TRUE, label.size = 4)+
  ggtitle("rpca UMAP Plot with Cell Type Labels")

so_kpmp_sc$group <- factor(so_kpmp_sc$group)
DimPlot(so_kpmp_sc, reduction = "umap.rpca",group.by = "group",label=F,raster=F) +
  ggtitle(paste0("rpca UMAP by Disease Category"))

#Barcharts of proportions
# By PT subtypes
cellcount <- so_kpmp_sc@meta.data %>% 
  filter(celltype2 == "PT")

label(cellcount$group) <- "Disease Status"
# label(cellcount$steatosis_cat) <- "Steatosis Category"
# cellcount$group <- ifelse(cellcount$group=="0+1", "Low Steatosis (0+1)","High Steatosis (2+3)")
# cellcount$group <- factor(cellcount$group,levels=c("Lean_Control","Obese_Control","Type_1_Diabetes","Type_2_Diabetes"),labels=c("Lean Controls","Obese Controls","Type 1 Diabetes","Type 2 Diabetes"))
cellcount$KPMP_celltype2 <- factor(cellcount$KPMP_celltype2,levels=c("PT-S1/S2","PT-S3","aPT"),labels=c("PT-S1/S2","PT-S3","aPT"))
cellcount <- cellcount %>% 
  filter(group!="Obese Controls")

prop_plot1 <- ggplot(data=cellcount,aes(group, fill = KPMP_celltype)) + 
  geom_bar(stat = "count", position = "fill") +
  theme_classic() +
  labs(x = NULL,
       y = "Proportion of Cells",
       fill = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        text = element_text(size = 20)) +
  ggtitle("Disease Status") +
  scale_fill_manual(values = c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51","darkred"))

prop_plot1 

#Barcharts of proportions
# By PT subtypes
cellcount <- so_kpmp_sc@meta.data %>% 
  filter(celltype2 == "TAL")

label(cellcount$group) <- "Disease Status"
# label(cellcount$steatosis_cat) <- "Steatosis Category"
# cellcount$group <- ifelse(cellcount$group=="0+1", "Low Steatosis (0+1)","High Steatosis (2+3)")
# cellcount$group <- factor(cellcount$group,levels=c("Lean_Control","Obese_Control","Type_1_Diabetes","Type_2_Diabetes"),labels=c("Lean Controls","Obese Controls","Type 1 Diabetes","Type 2 Diabetes"))
# cellcount$KPMP_celltype2 <- factor(cellcount$KPMP_celltype2,levels=c("PT-S1/S2","PT-S3","aPT"),labels=c("PT-S1/S2","PT-S3","aPT"))
cellcount <- cellcount %>%
  filter(group!="Obese Controls")

prop_plot1 <- ggplot(data=cellcount,aes(group, fill = KPMP_celltype)) + 
  geom_bar(stat = "count", position = "fill") +
  theme_classic() +
  labs(x = NULL,
       y = "Proportion of Cells",
       fill = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        text = element_text(size = 20)) +
  ggtitle("Disease Status") +
  scale_fill_manual(values = c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51","darkred"))

prop_plot1 

#Barcharts of proportions
# By PT subtypes
cellcount <- so_kpmp_sc@meta.data %>% 
  filter(grepl("EC-",KPMP_celltype) | grepl("EC/VSMC",KPMP_celltype))

label(cellcount$group) <- "Disease Status"
# label(cellcount$steatosis_cat) <- "Steatosis Category"
# cellcount$group <- ifelse(cellcount$group=="0+1", "Low Steatosis (0+1)","High Steatosis (2+3)")
# cellcount$group <- factor(cellcount$group,levels=c("Lean_Control","Obese_Control","Type_1_Diabetes","Type_2_Diabetes"),labels=c("Lean Controls","Obese Controls","Type 1 Diabetes","Type 2 Diabetes"))
# cellcount$KPMP_celltype2 <- factor(cellcount$KPMP_celltype2,levels=c("PT-S1/S2","PT-S3","aPT"),labels=c("PT-S1/S2","PT-S3","aPT"))
cellcount <- cellcount %>%
  filter(group!="Obese Controls")

prop_plot1 <- ggplot(data=cellcount,aes(group, fill = KPMP_celltype)) + 
  geom_bar(stat = "count", position = "fill") +
  theme_classic() +
  labs(x = NULL,
       y = "Proportion of Cells",
       fill = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        text = element_text(size = 20)) +
  ggtitle("Disease Status") +
  scale_fill_manual(values = c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51","darkred"))

prop_plot1


#3. Analysis ----
#ROBO and FN1 genes DEGs and Pathways (IPA & GSEA)

