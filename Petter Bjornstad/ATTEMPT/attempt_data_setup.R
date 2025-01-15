library(Seurat)
library(future)
library(tidyverse)
library(colorspace)
library(patchwork)
library(ggdendro)
library(cowplot)
library(ggpubr)
library(rstatix)
library(arsenal)
library(Biobase)
library(ReactomeGSA)
library(GSEABase)
library(msigdbr)
library(kableExtra)
library(knitr)
library(slingshot); library(SingleCellExperiment)
library(REDCapR)
library(edgeR)
library(data.table)
library(MAST)
library(limma)      # For linear modeling of microarray and RNA-seq data
library(muscat)

# clinical data clean & save
# Load dictionary function and file
source("/home/choiyej/Documents/CHCO-Code/Petter Bjornstad/Resources/YC/R Functions/label_harmonized_function.R")
source("/home/choiyej/Documents/CHCO-Code/Petter Bjornstad/CROCODILE/crocodile_functions.R")

token_dat <- read.csv("/run/user/778527649/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/api_tokens.csv")
attempt_token <- token_dat$Token[token_dat$Study == "ATTEMPT"]
attempt_dat <- redcap_read(redcap_uri = "https://redcap.ucdenver.edu/api/",
                           token = attempt_token)
attempt_grp <- read.csv("/run/user/778527649/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/ATTEMPT/Data Raw/ATTEMPT_unblinding.csv")

race_names <- c(
  "American Indian or Alaskan Native",
  "Asian",
  "Hawaiian or Pacific Islander",
  "Black or African American",
  "White",
  "Unknown",
  "Other"
)
ethnicity_names <- c(
  "Hispanic",
  "Non-Hispanic",
  "Unknown")

dat <- attempt_dat$data %>%
  group_by(subject_id) %>%
  fill(everything(), .direction = "down") %>%
  ungroup() %>% rowwise() %>%
  dplyr::mutate(visit = case_when(redcap_event_name == "visit_2_arm_1" ~ "PRE",
                                  redcap_event_name == "visit_3_arm_1" ~ "POST",
                                  T ~ "Screening"),
                age = case_when(is.na(bx_date) ~ age_consent,
                                T ~ as.integer((bx_date - dob) / dyears(1))),
                race = case_when(sum(c_across(starts_with("race___")), na.rm = TRUE) == 0 ~ "Unknown",
                                 sum(c_across(starts_with("race___")), na.rm = TRUE) > 1 ~ "More than one race",
                                 race___1 == 1 ~ race_names[1],
                                 race___2 == 1 ~ race_names[2],
                                 race___3 == 1 ~ race_names[3],
                                 race___4 == 1 ~ race_names[4],
                                 race___5 == 1 ~ race_names[5],
                                 race___6 == 1 ~ race_names[6],
                                 race___7 == 1 ~ race_names[7]),
                ethnicity = case_when(sum(c_across(starts_with("ethnicity___")), na.rm = TRUE) == 0 ~ "Unknown",
                                      sum(c_across(starts_with("ethnicity___")), na.rm = TRUE) > 1 ~ "More than one ethnicity",
                                      ethnicity___1 == 1 ~ ethnicity_names[1],
                                      ethnicity___2 == 1 ~ ethnicity_names[2],
                                      ethnicity___3 == 1 ~ ethnicity_names[3]),
                diabetes_duration = case_when(visit == "Screening" ~ as.integer((consent_date - t1d_date) / dyears(1)),
                                              T ~ as.integer((bx_date - t1d_date) / dyears(1)))) %>%
  filter(visit != "Screening") %>%
  left_join(attempt_grp)

# manually adding three participants missing in REDCap for now...
dat[nrow(dat) + 1, "subject_id"] <- 30058
dat[nrow(dat), "visit"] <- "PRE"
dat[nrow(dat), "treatment"] <- "Placebo"
dat[nrow(dat) + 1, "subject_id"] <- 30058
dat[nrow(dat), "visit"] <- "POST"
dat[nrow(dat), "treatment"] <- "Placebo"
dat[nrow(dat) + 1, "subject_id"] <- 30173
dat[nrow(dat), "visit"] <- "PRE"
dat[nrow(dat), "treatment"] <- "Placebo"
dat[nrow(dat) + 1, "subject_id"] <- 30173
dat[nrow(dat), "visit"] <- "POST"
dat[nrow(dat), "treatment"] <- "Placebo"
dat[nrow(dat) + 1, "subject_id"] <- 30186
dat[nrow(dat), "visit"] <- "PRE"
dat[nrow(dat), "treatment"] <- "Placebo"
dat[nrow(dat) + 1, "subject_id"] <- 30186
dat[nrow(dat), "visit"] <- "POST"
dat[nrow(dat), "treatment"] <- "Placebo"

dat$visit <- factor(dat$visit, levels = c("PRE", "POST"))

saveRDS(dat, file = "/run/user/778527649/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/ATTEMPT/Data Clean/attempt_clinical_data.RDS")

# scRNA
plan(multisession, workers = 16)
options(future.globals.maxSize=2e9)
so_attempt <- readRDS("/home/choiyej/Documents/Local data/PB_attempt_harmony_rpca_Sept2024.RDS")
so_attempt_meta <- so_attempt@meta.data %>%
  mutate(subject_id = Subject.ID,
         visit = case_when(Visit == "BL" ~ "PRE", 
                           Visit == "4M" ~ "POST"))
so_attempt_meta <- left_join(so_attempt_meta, dat, by = c("subject_id", "visit"))
rownames(so_attempt_meta) <- so_attempt_meta$barcode
so_attempt <- AddMetaData(so_attempt, so_attempt_meta)
saveRDS(so_attempt, file = "/run/user/778527649/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/ATTEMPT/Data Clean/so_attempt.RDS")

so_attempt$celltype_pt <- ifelse(grepl("PT-", so_attempt$celltype),
                                 "PT", as.character(so_attempt$celltype))
so_attempt_pt <- subset(so_attempt, celltype_pt == "PT" & celltype != "PT_lowQuality")
options(future.globals.maxSize = 3000 * 1024^3)
so_attempt_pt <- NormalizeData(so_attempt_pt)
so_attempt_pt <- ScaleData(so_attempt_pt)
ElbowPlot(so_attempt_pt)
so_attempt_pt <- RunPCA(so_attempt_pt, ncomponents = 10, features = VariableFeatures(object = so_attempt_pt))
so_attempt_pt <- FindNeighbors(so_attempt_pt)
so_attempt_pt <- FindClusters(so_attempt_pt)
so_attempt_pt <- RunUMAP(so_attempt_pt, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_pt$visit <- factor(so_attempt_pt$visit, levels = c("PRE", "POST"))
so_attempt_pt$treatment <- factor(so_attempt_pt$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                  labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_pt, file = "/run/user/778527649/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/ATTEMPT/Data Clean/attempt_so_pt.RDS")

# Fetch all data upfront
gene_names <- rownames(so_attempt_pt)
all_data <- FetchData(so_attempt_pt, vars = c("subject_id", "treatment", "visit", "ident", gene_names)) %>%
  group_by(subject_id, treatment, visit) %>%
  dplyr::summarise(across(where(is.numeric), ~ mean(.x > 0, na.rm = TRUE)), .groups = 'drop')

all_data <- all_data %>%
  group_by(subject_id) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  mutate(treatment = case_when(treatment == "Placebo" ~ "Placebo",
                               T ~ "Dapagliflozin"))
gc()

all_data$treatment <- factor(all_data$treatment)
all_data$treatment <- relevel(all_data$treatment, ref = "Placebo")
saveRDS(all_data, file = "/run/user/778527649/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/ATTEMPT/Data Clean/attempt_pb.RDS")