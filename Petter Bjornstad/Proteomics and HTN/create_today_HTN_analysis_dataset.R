library(tidyverse)
library(readxl)
library(Seurat)
library(pedbp)
top_htn_df <- read_excel("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted HTN model selection.xlsx", sheet = "HTN CPH base")
de_genes_htn <- top_htn_df[top_htn_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_htn <- setNames(de_genes_htn$estimate, de_genes_htn$EntrezGeneID)
top_htn <- top_htn_df %>%
 # filter(Target %in% c("WFKN2","SEZ6L","SCG3","PSA","LSAMP","H6ST3","T132B","Nr-CAM","PEDF","IGLO5",
#                       "PSB3","Myosin light chain 1","PCD10:ECD","UNC5H4","TLR5","SLIK1","PSPC1",
#                       "STA10","Secretoglobin family 3A member 1","sICAM-5")) %>%
    filter(adj.p.value <= 0.05) %>%
    slice_max(abs(log(estimate)), n = 21) %>%
  pull(AptName)
top_htn_df_10 <- read_excel("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic limma yr10 adjusted.xlsx", sheet = "htn_moderated_FDR")
de_genes_htn_10 <- top_htn_df_10[top_htn_df_10$P.Value <= 0.05, c("EntrezGeneID", "logFC")]
de_genes_htn_10 <- setNames(de_genes_htn_10$logFC, de_genes_htn_10$EntrezGeneID)

##########################################################

# OLD CODE PRIOR TO NEW HARMONIZED SOMA DATASET

# Import proteomics data for RH and IMPROVE
 load("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Proteomics HTN/Copy of old soma dataset for HTN response/analytes.Rdata")
 load("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Proteomics HTN/Copy of old soma dataset for HTN response/rh_soma.Rdata")
 load("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Proteomics HTN/Copy of old soma dataset for HTN response/improve_soma.Rdata")
# # Format and combine
 improve_soma <- improve_soma %>% dplyr::select(SampleDescription, TimePoint, contains("seq."))
 rh_soma <- rh_soma %>% dplyr::select(SampleDescription, TimePoint, contains("seq."))
 soma <- rbind(improve_soma, rh_soma)
# Transform
 soma[, 3:ncol(soma)] <- lapply(soma[, 3:ncol(soma)], log)
# Merge clinical and SOMA data
 soma <- soma %>%
   rename(record_id = "SampleDescription", visit = "TimePoint") %>%
   mutate(
     record_id = sub("IT2D-", "IT_", record_id),
     visit = case_when(
       visit == "BL" ~ "baseline",
       visit == "Baseline" ~ "baseline",
       visit == "3M" ~ "3_months_post_surgery",
       visit == "12M" ~ "12_months_post_surgery"
     )
   )

 # Import and clean harmonized data
 df <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))
 df <- df %>%
   dplyr::select(record_id, co_enroll_id, study, visit, date, group, age, sex, race, ethnicity,
                 diabetes_duration, sglti_timepoint, sglt2i_ever, elevated_albuminuria,
                 bmi, hba1c, gfr_bsa_plasma, gfr_raw_plasma, gfr_bsa_plasma_urine,
                 gfr_raw_plasma_urine, acr_u, map, sbp, dbp, height, eGFR_fas_cr, participation_status)
 coenroll_id <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Renal HERITAGE/Data_Cleaned/coenrolled_ids.csv") %>%
   pivot_longer(cols = 'improve_id':'rh2_id',
                values_to = "record_id") %>%
   dplyr::select(merged_id, record_id) %>%
   filter(record_id != "")
 df <- df %>%
   filter(study == "RENAL-HEIRitage"|study == "RENAL-HEIR"|study == "IMPROVE") %>%
   dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                    .by = c(record_id)) %>%
   filter(participation_status!="Removed"|is.na(participation_status)) %>%
   mutate(race_ethnicity_condensed = case_when(race == "White" &
                                                 ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                               race == "Black or African American" &
                                                 ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                               ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                               T ~ "Not Hispanic or Latino Other")) %>%
   mutate(elevated_albuminuria = case_when(elevated_albuminuria == "Yes" ~ "Elevated albuminuria",
                                           elevated_albuminuria == "No" ~ "Normoalbuminuria")) %>%
   mutate_at(vars(starts_with("fsoc")), function(x) case_when(x < 0 ~ 0, T~x)) %>%
   left_join(coenroll_id)
  # BP percentiles and HTN
 df$bp_age <- df$age * 12
 df$bp_age[df$bp_age >= 19 * 12] <- 19 * 12 - 1e-8
 bps <- p_bp(
   q_sbp = df$sbp, q_dbp = df$dbp, age = df$bp_age,
   male = df$sex == "Male", height = df$height
 )
 df$sbp_perc <- bps$sbp_percentile
 df$sbp_perc[df$age >= 19] <- NA
 df$dbp_perc <- bps$dbp_percentile
 df$dbp_perc[df$age >= 19] <- NA
 df$htn <- df$sbp >= 130 | df$dbp >= 80 | df$sbp_perc >= 0.95 | df$dbp_perc >= 0.95
 df$htn[is.na(df$htn)] <- F
 df$htn[is.na(df$sbp) & is.na(df$dbp) & is.na(df$sbp_perc) & is.na(df$dbp_perc)] <- NA
 df$htn <- factor(df$htn, levels = c(F, T), labels = c("HTN-", "HTN+"))
 # Hyperfiltration
 df$hyp <- factor(df$eGFR_fas_cr >= 135, levels = c(F, T), labels = c("eGFR < 135", "eGFR >= 135"))
 # UACR
 df$elevated_uacr <- factor(df$acr_u >= 30, labels = c("UACR < 30", "UACR >= 30"))
 # Merge
 df <- merge(df, soma, by = c("record_id", "visit"), all.x = F, all.y = T)
 # Baseline visit only for IMPROVE
 df <- df %>% filter(visit == "baseline")

 coenroll <- df %>% filter(!is.na(co_enroll_id))
 coenroll$merged_id <- str_replace(coenroll$merged_id, "IT2D", "IT_")
 coenroll <- coenroll %>%
   mutate(
     merged_id = case_when(
       (is.na(merged_id) & str_starts("IT", record_id)) ~ paste0(record_id, "-", co_enroll_id),
       (is.na(merged_id) & !str_starts("IT", record_id)) ~ paste0(co_enroll_id, "-", record_id),
       !is.na(merged_id) ~ merged_id
     )
   )
 temp <- coenroll %>% group_by(merged_id) %>% arrange(date) %>% filter(row_number()==1)
 out <- temp$co_enroll_id
 df <- df %>% filter(!co_enroll_id %in% out)

#######################################################

# Import and clean data
# df <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/soma_harmonized_dataset.csv", na.strings = "")
# df <- df %>%
#   filter(
#     study %in% c("IMPROVE", "RENAL-HEIR","RENAL-HEIRitage"),
#     !grepl("IT2D", co_enroll_id), participation_status == "Participated"
#   ) %>%
  #%>%
  #  select(
  #    record_id, co_enroll_id, visit, group, age, sex, race, ethnicity,
  #    diabetes_duration, sglti_timepoint, sglt2i_ever, elevated_albuminuria,
  #    bmi, hba1c, gfr_bsa_plasma, gfr_raw_plasma, gfr_bsa_plasma_urine,
  #    gfr_raw_plasma_urine, acr_u, map, sbp, dbp, height, eGFR_fas_cr
  #  ) 
#   group_by(record_id, visit) %>%
#   summarise(across(where(is.character), ~ last(na.omit(.x))),
#             across(where(is.factor), ~ last(na.omit(.x))),
#             across(where(is.numeric), ~ mean(.x, na.rm = T)),
#             .groups = "drop"
#   ) %>%
#   mutate_all(~ ifelse(is.nan(.), NA, .))
# # make sure we are only using data from the first batch of Soma results
# check <- read_delim("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Local cohort Somalogic data/WUS-22-002/samplesheet_WUS-22-002_nonQC.samples.txt")
# test <- df %>% filter(record_id %in% check$SampleDescription)
# # Import proteomics data for RH and IMPROVE
# load("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Combined SomaScan/analytes.Rdata")
# # below not needed now that we have combined soma and harmonized data
# #load("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Renal HERITAGE/Somalogic data/rh_soma.Rdata")
# #load("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/IMPROVE T2D/Somalogic data/improve_soma.Rdata")
# # Format and combine
# #improve_soma <- improve_soma %>% select(SampleDescription, TimePoint, contains("seq."))
# #rh_soma <- rh_soma %>% select(SampleDescription, TimePoint, contains("seq."))
# soma <- df %>% select(record_id, visit, contains("seq."))
# # Transform
# soma[, 3:ncol(soma)] <- lapply(soma[, 3:ncol(soma)], log)
# # Merge clinical and SOMA data
# #soma <- soma %>%
# # rename(record_id = "SampleDescription", visit = "TimePoint") %>%
# #  mutate(
# #    record_id = sub("IT2D-", "IT_", record_id),
# #    visit = case_when(
# #      visit == "BL" ~ "baseline",
# #      visit == "3M" ~ "3_months_post_surgery",
# #      visit == "12M" ~ "12_months_post_surgery"
# #    )
# #  )



# Save
save(df, analytes, 
  #top_mac, top_mic, top_hyp, top_rapid, 
  top_htn, 
  #top_htn_sbp,
  #top_mac_df, top_mic_df, top_hyp_df, top_rapid_df, 
  top_htn_df, 
  #top_htn_sbp_df,
  #de_genes_mac, de_genes_mic, de_genes_hyp, de_genes_rapid, 
  de_genes_htn, 
  #de_genes_htn_sbp,de_genes_mac_10, de_genes_mic_10, de_genes_hyp_10, de_genes_rapid_10, de_genes_htn_10, de_genes_htn_sbp_10,
  file = "/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Proteomics HTN/copy_of_old_analysis_dataset_for_HTN_response.RData"
)

