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

#2. Visualize & Descriptive Stats ----
##a. UMAPS & Barcharts ----
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
rm(cellcount,plot1,prop_plot1)



##b. Descriptive Statistics ----
#Get metadata for everyone
dat <- so_kpmp_sc@meta.data %>%
  group_by(record_id) %>%
  summarise(across(everything(), first)) %>%
  ungroup() 

dat$hba1c <- as.numeric(dat$hba1c)
dat$eGFR_CKD_epi <- as.numeric(dat$eGFR_CKD_epi)

#Table 1. 
table1(~ age + sex + race_ethnicity  + bmi + triglycerides + hba1c + medication| study, data=dat)
table1(~pah_clear_bsa + eGFR_CKD_epi +acr_u +mfm_1+insulin_1| study, data=dat)

#Covariates to adjust for: 
#Acru, metformin/insulin,bmi,age,tg
#Although check on metformin/insulin variables...

#3. Analysis ----
#ROBO and FN1 genes DEGs and Pathways (IPA & GSEA)
##a. PT Cells ----
#PT
so_celltype <- subset(so_kpmp_sc,celltype1=="PT")
ncol(so_celltype) #2246 cells
nrow(so_celltype) #9196 genes
# Extract the gene expression data for all genes
gene_expression <- as.data.frame(GetAssayData(so_celltype, layer = "data"))

#Assign gene names as rownames, cell names as colnames
# rownames(gene_expression) <- rownames(so_celltype) #Gene Names
# colnames(gene_expression) <- colnames(so_celltype) #Cell Names
#Transpose gene expression dataset to merge with metadata
gene_expression <- t(gene_expression)
gene_expression <- data.frame(gene_expression) #Make a dataframe again after transposing
#Set gene list
gene_list_total <- colnames(gene_expression)
gene_expression$cellname <- rownames(gene_expression) 
rownames(gene_expression) <- NULL

# Extract the metadata
metadata <- so_celltype@meta.data
# metadata <- metadata %>%
#   mutate(across(everything(),~ifelse(.==".",NA,.)))
metadata$cellname <- rownames(metadata)
rownames(metadata) <- NULL

# Combine the gene expression data and metadata
data <- tidylog::full_join(metadata,gene_expression,by="cellname")
rm(metadata,gene_expression)

#Compare Health Controls to Type 2 Diabetes on no meds
data_subset <- data
  # filter(group=="Type_2_Diabetes" | group=="Lean_Control") %>%
  # filter(medication=="no_med")
rm(data)
data_subset$group <- factor(data_subset$medication)
data_subset$group <- relevel(data_subset$group,ref="Lean_Control")
#Select covariates to keep in the model: Acru, metformin/insulin,bmi,age,tg
data_subset <- data_subset %>% 
  dplyr::select(c("kit_id","group","medication","age","bmi","triglycerides","insulin_1","mfm_1","acr_u",all_of(gene_list_total)))

#Try as batch loop
# Set batch size
batch_size <- 2000
total_cores <- 50

# Deduplicate gene list to avoid double-processing
gene_list_total <- unique(gene_list_total)

# Split the gene list into batches
gene_batches <- split(gene_list_total, ceiling(seq_along(gene_list_total) / batch_size))

# Define the function that processes a single gene
process_gene <- function(gene) {
  if (sum(data_subset[[gene]]) > 0) {
    m1 <- as.formula(paste0(gene, " ~ group + age + bmi + acr_u + triglycerides + insulin_1 + mfm_1 + (1 | kit_id)"))
    
    model1 <- tryCatch({
      glmmTMB(m1, data = data_subset, family = gaussian)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(model1)) {
      Beta <- summary(model1)$coef$cond[2,1]
      PValue <- summary(model1)$coef$cond[2,4]
    } else {
      Beta <- NA
      PValue <- NA
    }
  } else {
    Beta <- NA
    PValue <- NA
  }
  
  return(data.frame(Gene = gene, Beta = Beta, PValue = PValue))
}

# Set number of cores
total_cores <- parallel::detectCores() - 1

# Initialize and register cluster
cl <- makeCluster(total_cores)
registerDoParallel(cl)

# Run analysis per batch
all_results <- lapply(seq_along(gene_batches), function(batch_idx) {
  batch_genes <- gene_batches[[batch_idx]]
  
  batch_results <- foreach(
    gene = batch_genes,
    .combine = rbind,
    .packages = c("glmmTMB", "lme4"),
    .export = c("process_gene", "data_subset")
  ) %dopar% {
    process_gene(gene)
  }
  
  cat("Processed batch", batch_idx, "with", length(batch_genes), "genes\n")
  return(batch_results)
})

# Combine all batch results into a single data frame
final_results <- do.call(rbind, all_results)
# Stop cluster
stopCluster(cl)

# #Make volcano plot of all gene results for group
full_results <- final_results %>%
  mutate(fdr=p.adjust(PValue,method="fdr"))  
# mutate(fdr3=p.adjust(PValue3,method="fdr"))
full_results$PValue10 <- -log10(pmax(full_results$PValue, 1e-10))  # Avoid log(0)
# full_results$PValue10_3 <- -log10(pmax(full_results$PValue3, 1e-10))  # Avoid log(0)
Nonconvergence_Rate <- paste0(round(((length(which(is.na(full_results$PValue)))+length(which(full_results$PValue=="NaN")))/length(full_results$Gene))*100,0),"%")
# view(full_results)

write.csv(full_results,fs::path(dir.results,"PT_Cells_No_Med_T2D_LC_adj_cov_lmm_no_zi.csv"))
# full_results <- read.csv(fs::path(dir.results,"Hep_5_Fibrosis_Steatosis.csv"))

full_results$color <- ifelse(full_results$fdr < 0.05 & full_results$Beta > 0, "lightcoral",
                             ifelse(full_results$fdr < 0.05 & full_results$Beta < 0, "lightblue", "gray"))

# Identify significant points (fdr < 0.05)
significant_df <- full_results[full_results$fdr < 0.05, ]
# 
# full_results$color3 <- ifelse(full_results$fdr3 < 0.2 & full_results$Beta3 > 0, "lightcoral",
#                               ifelse(full_results$fdr3 < 0.2 & full_results$Beta3 < 0, "lightblue", "gray"))
# 
# # Identify significant points (fdr < 0.05)
# significant_df3 <- full_results[full_results$fdr3 < 0.2, ]

Genes <- length(unique(full_results$Gene))
Cells <- ncol(so_celltype)

#Figure Range
max <- max(significant_df$Beta,na.rm=T)
min <- min(significant_df$Beta,na.rm=T)


# full_results$Log2FC <- 2^(full_results$Beta)
# Create the volcano plot using ggplot
volcano_plot <- ggplot(full_results, aes(x = Beta, y = PValue10, color = color)) +
  geom_point(alpha = 0.7) +  # Plot points with transparency
  scale_color_identity() +  # Use the color column directly
  theme_minimal() +  # Minimal theme
  labs(
    title = "Type 2 Diabetes vs. Lean Controls (All No Medication)",
    subtitle = "PT Cells, Adjusted for Age, Sex & BMI",
    x = "Fold Change",
    y = "-log10(P-Value)",
    color = "FC Direction Direction",
    caption = paste0("FDR < 0.05, Genes = ",Genes,", Cells = ",Cells,", Non-Convergence Rate: ",Nonconvergence_Rate)
  ) +
  xlim(min,max)+
  theme(
    plot.title = element_text(hjust = 0),
    axis.text.x = element_text(angle = 0, hjust = 1)
  )+
  # # Add labels for significant points
  geom_text(data = significant_df, aes(label = Gene),
            vjust = 1, hjust = 1, size = 3, check_overlap = TRUE, color = "black")
# Add labels for significant points with ggrepel
# geom_text_repel(data = significant_df, aes(label = Gene),
#                 size = 3, color = "black", box.padding = 0.5, max.overlaps = 20)

volcano_plot
pdf(fs::path(dir.results,"Plot_PT_Cells_T2D_LC_NoMed.pdf"),width=10,height=7)
print(volcano_plot)
dev.off()
