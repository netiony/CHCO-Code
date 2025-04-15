############################################################
#                                                          #
#   Author: Long Yuan                                      #
#   Email : lyuan13@jhmi.edu                               #
#                                                          #
#           __                                             #
#       (___()'`;                                          #
#       /,    /`                                           #
#       \\"--\\                                            #
#                                                          #
#   Script Description:                                    #
#   [Mixed Model Analysis To Identify DEGs]                #
#                                                          #
############################################################


library(lme4)
library(lmerTest) 
library(dplyr)
library(purrr)
library(tidyr)
library(broom.mixed)
library(Seurat)
library(ggplot2)

scrna <- readRDS("/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/scRNA/data_raw/PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds")
scrna <- subset(scrna, subset = percent.mt < 50 & nFeature_RNA < 5000 & nFeature_RNA > 500) 
# why subset to features <5000? doublets, etc.

#############
## IMPROVE ##
#############
improve <- subset(scrna, subset = cohort == "IMPROVE"| record_id %in% c("RH-59-T", "RH-60-T", "RH-65-T", "RH-66-T"))
improve@meta.data[improve@meta.data$record_id == "RH-65-T", ]$pre_post <- 'pre'
improve@meta.data[improve@meta.data$record_id == "RH-59-T", ]$record_id <- 'IT_07'
improve@meta.data[improve@meta.data$record_id == "RH-60-T", ]$record_id <- 'IT_08'
improve@meta.data[improve@meta.data$record_id == "RH-65-T", ]$record_id <- 'IT_09'
improve@meta.data[improve@meta.data$record_id == "RH-66-T", ]$record_id <- 'IT_10'


##########
### RH ###
##########
rh <- subset(scrna, subset = cohort == c('RENAL HEIR', 'RENAL HEIRITAGE') & 
               record_id %in% c('RH-23-T', 'RH2-14-T', 'RH-67-T', 'RH2-19-T'))
rh@meta.data[rh@meta.data$record_id == "RH-23-T", ]$pre_post <- 'pre'
rh@meta.data[rh@meta.data$record_id == "RH2-14-T", ]$pre_post <- 'post'
rh@meta.data[rh@meta.data$record_id == "RH-67-T", ]$pre_post <- 'pre'
rh@meta.data[rh@meta.data$record_id == "RH2-19-T", ]$pre_post <- 'post'

improve@meta.data$treatment <- "VSG"
rh@meta.data$treatment <- "Standard"
combined <- merge(improve, rh)

combined@meta.data[combined@meta.data$record_id == "RH2-14-T", ]$record_id <- 'RH-23-T'
combined@meta.data[combined@meta.data$record_id == "RH2-19-T", ]$record_id <- 'RH-67-T'

#############
# Get Genes #
#############
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

# Extract HVGs
hvg_layer1 <- VariableFeatures(combined, assay = "RNA", layer = "counts.1")
hvg_layer2 <- VariableFeatures(combined, assay = "RNA", layer = "counts.2")

# Take the union which yields 2643 genes
top2000_union <- union(hvg_layer1, hvg_layer2)

# Extract metadata
meta <- combined@meta.data
meta$cell <- rownames(meta)
meta$pre_post <- factor(meta$pre_post, levels = c("pre", "post"))
meta$treatment <- factor(meta$treatment, levels = c("Standard", "VSG"))
meta$record_id <- as.factor(meta$record_id)

data1 <- GetAssayData(improve, slot = "data", layer = "data") #dim 31332 33170
data2 <- GetAssayData(rh, slot = "data", layer = "data") #dim 31332 10903
expr_matrix <- cbind(data1, data2) #31332 44073
expr_matrix <- expr_matrix[top2000_union, ]


################
# Run Pipeline #
################
# Initialize list
results <- list()

# Loop through all genes and all KPMP_celltypes
for (gene in top2000_union) {
  for (cell_type in unique(meta$KPMP_celltype)) {
    
    # Subset metadata
    sub_meta <- meta %>% filter(KPMP_celltype == cell_type)
    if (nrow(sub_meta) < 50) next  # Skip small cell types
    
    # Extract expression
    if (!gene %in% rownames(expr_matrix)) next
    sub_expr <- expr_matrix[gene, sub_meta$cell, drop = FALSE]
    if (all(sub_expr == 0)) next  # Skip genes not expressed
    
    # Add expression
    df <- sub_meta %>% mutate(expr = as.numeric(sub_expr))
    
    #
    try({
      model <- lmerTest::lmer(expr ~ pre_post * treatment + (1 | record_id), data = df)
      tidy_mod <- tidy(model, effects = "fixed", conf.int = TRUE)
      
      # Extract info
      interaction <- tidy_mod %>% filter(term == "pre_postpost:treatmentVSG")
      if (nrow(interaction) == 0) next
      
      results[[length(results) + 1]] <- tibble(
        gene = gene,
        cell_type = cell_type,
        estimate = interaction$estimate,
        conf.low = interaction$conf.low,
        conf.high = interaction$conf.high,
        p.value = interaction$p.value,
        converged = TRUE,
        is_singular = isSingular(model, tol = 1e-4)
      )
    }, silent = TRUE)
  }
}

# Combine
final_df <- bind_rows(results)

# Adjust FDR (per cell type, not global)
final_df2 <- final_df %>%
  group_by(cell_type) %>%
  mutate(FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()

# Save
write.csv(final_df2, "mixed_model_results_summary.csv", row.names = FALSE)