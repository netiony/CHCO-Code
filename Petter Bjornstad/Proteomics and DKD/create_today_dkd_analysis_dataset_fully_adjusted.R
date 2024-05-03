library(tidyverse)
library(readxl)
library(Seurat)
library(pedbp)
# Import top proteins for MIC, MAC, etc. at baseline
top_mic_df <- read_excel("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted new covariates.xlsx", sheet = "MIC CPH")
de_genes_mic <- top_mic_df[top_mic_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_mic <- setNames(de_genes_mic$estimate, de_genes_mic$EntrezGeneID)
top_mic <- top_mic_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_mac_df <- read_excel("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted new covariates.xlsx", sheet = "MAC CPH")
de_genes_mac <- top_mac_df[top_mac_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_mac <- setNames(de_genes_mac$estimate, de_genes_mac$EntrezGeneID)
top_mac <- top_mac_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_mic.or.mac_df <- read_excel("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted new covariates.xlsx", sheet = "MIC.OR.MAC CPH")
de_genes_mic.or.mac <- top_mic.or.mac_df[top_mic.or.mac_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_mic.or.mac <- setNames(de_genes_mic.or.mac$estimate, de_genes_mic.or.mac$EntrezGeneID)
top_mic.or.mac <- top_mic.or.mac_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_hyp_df <- read_excel("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted new covariates.xlsx", sheet = "HYP CPH")
de_genes_hyp <- top_hyp_df[top_hyp_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_hyp <- setNames(de_genes_hyp$estimate, de_genes_hyp$EntrezGeneID)
top_hyp <- top_hyp_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_rapid_df <- read_excel("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted new covariates.xlsx", sheet = "RAPID CPH")
de_genes_rapid <- top_rapid_df[top_rapid_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_rapid <- setNames(de_genes_rapid$estimate, de_genes_rapid$EntrezGeneID)
top_rapid <- top_rapid_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_htn_sbp_df <- read_excel("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted new covariates.xlsx", sheet = "HTN with SBP CPH")
de_genes_htn_sbp <- top_htn_sbp_df[top_htn_sbp_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_htn_sbp <- setNames(de_genes_htn_sbp$estimate, de_genes_htn_sbp$EntrezGeneID)
top_htn_sbp <- top_htn_sbp_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_htn_df <- read_excel("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted new covariates.xlsx", sheet = "HTN CPH")
de_genes_htn <- top_htn_df[top_htn_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_htn <- setNames(de_genes_htn$estimate, de_genes_htn$EntrezGeneID)
top_htn <- top_htn_df %>%
  filter(Target %in% c("WFKN2","SEZ6L","SCG3","PSA","LSAMP","H6ST3","T132B","Nr-CAM","PEDF","IGLO5",
                       "PSB3","Myosin light chain 1","PCD10:ECD","UNC5H4","TLR5","SLIK1","PSPC1",
                       "STA10","Secretoglobin family 3A member 1","sICAM-5")) %>%
  slice_max(abs(log(estimate)), n = 20) %>%
  pull(AptName)
top_htn_uacr_df <- read_excel("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted new covariates.xlsx", sheet = "HTN with UACR CPH")
de_genes_htn_uacr <- top_htn_uacr_df[top_htn_uacr_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_htn_uacr <- setNames(de_genes_htn_uacr$estimate, de_genes_htn_uacr$EntrezGeneID)
top_htn_uacr <- top_htn_uacr_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_htn_egfr_df <- read_excel("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted new covariates.xlsx", sheet = "HTN with eGFR CPH")
de_genes_htn_egfr <- top_htn_egfr_df[top_htn_egfr_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_htn_egfr <- setNames(de_genes_htn_egfr$estimate, de_genes_htn_egfr$EntrezGeneID)
top_htn_egfr <- top_htn_egfr_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_neuro_df <- read_excel("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted new covariates.xlsx", sheet = "NEURO CPH")
de_genes_neuro <- top_neuro_df[top_neuro_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_neuro <- setNames(de_genes_neuro$estimate, de_genes_neuro$EntrezGeneID)
top_neuro <- top_neuro_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_retino_df <- read_excel("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted new covariates.xlsx", sheet = "RETINO CPH")
de_genes_retino <- top_retino_df[top_retino_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_retino <- setNames(de_genes_retino$estimate, de_genes_retino$EntrezGeneID)
top_retino <- top_retino_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_glyc_df <- read_excel("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted new covariates.xlsx", sheet = "GLYC CPH")
de_genes_glyc <- top_glyc_df[top_glyc_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_glyc <- setNames(de_genes_glyc$estimate, de_genes_glyc$EntrezGeneID)
top_glyc <- top_glyc_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_glyc_a1c_df <- read_excel("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted new covariates.xlsx", sheet = "GLYC with A1c CPH")
de_genes_glyc_a1c <- top_glyc_a1c_df[top_glyc_a1c_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_glyc_a1c <- setNames(de_genes_glyc_a1c$estimate, de_genes_glyc_a1c$EntrezGeneID)
top_glyc_a1c <- top_glyc_a1c_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
# Save
save(top_mac, top_mic, top_mic.or.mac, top_hyp, top_rapid, top_htn, top_htn_sbp, top_htn_uacr, top_htn_egfr, 
  top_mac_df, top_mic_df, top_mic.or.mac_df, top_hyp_df, top_rapid_df, top_htn_df, top_htn_sbp_df, top_htn_uacr_df, top_htn_egfr_df,
  de_genes_mac, de_genes_mic, de_genes_mic.or.mac, de_genes_hyp, de_genes_rapid, de_genes_htn, de_genes_htn_sbp, de_genes_htn_uacr, de_genes_htn_egfr,
  file = "/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Proteomics and DKD/Data_Cleaned/TODAY_top_proteins_new_covariates.RData"
)
