# ```{r libraries, echo=F, include = F}
# library(tidyverse)
# library(BiocManager)
# library(arsenal)
# library(dplyr)
# library(ggplot2)
# library(ggrepel)
# library(Seurat)
# library(future)
# library(colorspace)
# library(patchwork)
# library(ggdendro)
# library(cowplot)
# library(ggpubr)
# library(rstatix)
# library(table1)
# library(Biobase)
# library(ReactomeGSA)
# library(GSEABase)
# library(msigdbr)
# library(kableExtra)
# library(knitr)
# library(EnhancedVolcano)
# library(MAST)
# library(future)
# #library(slingshot)
# library(SingleCellExperiment)
# library(RColorBrewer)
# library(scales)
# library(viridis)
# #library(UpSetR)
# #library(pheatmap)
# #library(fgsea)
# #library(tradeSeq)
# #library(DescTools)
# # remotes::install_github("dynverse/dynfeature")
# # remotes::install_github("dynverse/dynplot")
# # remotes::install_github("elolab/Totem",force=T)
# #library(Totem)
# #library(dyndimred)
# #library(pushoverr)
# library(future)
# 
# #Increase Memory
# mem.maxVSize(64000000000)
# 
# # Parallel processing
# # plan(multicore, workers = 16)
# # options(future.globals.maxSize=2e9)
# # options(future.globals.maxSize = NULL)   # 24 GB in bytes
# plan()
# future::plan("sequential")
# options(future.globals.maxSize = 3e9)
# 
# #Set up directories
# dir.dat <- c("/Users/hhampson/Dropbox/Bjornstad data")
# dir.home <- c("/Users/hhampson/Documents/CHCO-Code/Petter Bjornstad/OATs and PAH/OATs and PAH Clearance")
# dir.dat2 <- c("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean")
# 
# #Load functions
# # source(fs::path(dir.home,"OAT_PAH_Functions.R"))
# ```

```{r echo = F}
##a. Kidney scRNA----
so_kidney_sc <- readRDS(fs::path(dir.dat,"Kidney scRNA","PB_90samples_Harmony_rpca_Fadhl_PhilApproved_091024.RDS"))

##b. Metadata ----
meta_raw <- read.csv(fs::path(dir.dat,"Clinical Data","renal_clearance_biopsy.csv"))
harm_data <- read.csv("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv") %>% 
  dplyr::select(record_id,cryostor_id,visit,pah_clear_bsa,pah_clear_abs,group)

##c. Process Data----
#Filter out meta data participants without sc 
ids <- unique(so_kidney_sc@meta.data$cryostor_id)
meta_raw <- meta_raw %>% 
  filter(cryostor_id %in% ids) 
ids2 <- unique(so_kidney_sc@meta.data$record_id)
harm_data <- harm_data %>% 
  filter(record_id %in% ids2) %>% 
  filter(!is.na(pah_clear_bsa)) %>% 
  # filter(visit=="baseline") %>% 
  mutate(diabetes=ifelse(group=="Type 1 Diabetes"|group=="Type 2 Diabetes","Yes","No")) %>% 
  dplyr::select(record_id,visit,pah_clear_bsa,diabetes,pah_clear_abs,pah_bsa,
                pah_bsa_plasma_urine)
meta_raw <- tidylog::left_join(meta_raw,harm_data,by=c("record_id","visit"))

#Merge new metadata with 
meta_kidney_sc <-  tidylog::left_join(so_kidney_sc@meta.data, meta_raw,by="cryostor_id")

#Check there are no duplicates
length(which(duplicated(meta_kidney_sc$barcode)))

#Add rownames back to metadata to merge back into saurat object 
rownames(meta_kidney_sc) <- rownames(so_kidney_sc@meta.data)

#Create study visit variale to filter to RH, CROCCODILE and IMPROVE
so_kidney_sc@meta.data$study= ifelse(grepl("RH",so_kidney_sc@meta.data$record_id) | grepl("CRC",so_kidney_sc@meta.data$record_id) | grepl("IT",so_kidney_sc@meta.data$record_id),"subset","non-subset")

# Add new meta data to so 
so_kidney_sc <- AddMetaData(so_kidney_sc, meta_kidney_sc)

gc()
#Filter to baseline visits only
so_kidney_sc <- subset(so_kidney_sc, visit != "12_months_post_surgery")

gc()
#Filter to RH/RH2, CROCODILE, and IMPROVE
so_kidney_sc <- subset(so_kidney_sc, study == "subset")


# #Normalize & Scale Data
# so_kidney_sc <- NormalizeData(so_kidney_sc)
# so_kidney_sc <- ScaleData(so_kidney_sc)

```

#3. Visualize Data 
```{r Visualize Data}
#Perform PCA
so_kidney_sc <- RunPCA(so_kidney_sc, features = VariableFeatures(object = so_kidney_sc))
ElbowPlot(so_kidney_sc)
# Cluster cells
so_kidney_sc <- FindNeighbors(so_kidney_sc)
so_kidney_sc <- FindClusters(so_kidney_sc)
# Perform UMAP and tSNE
so_kidney_sc <- RunUMAP(so_kidney_sc, dims = 1:15)
DimPlot(so_kidney_sc, reduction = "umap")
DimPlot(so_kidney_sc, reduction = "umap", group.by = "sglt2i_ever")
DimPlot(so_kidney_sc, reduction = "umap", group.by = "diabetes")
DimPlot(so_kidney_sc, reduction = "umap", group.by = "group")
```