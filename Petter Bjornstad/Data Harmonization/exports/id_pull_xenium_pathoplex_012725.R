# load libraries
library(tidyverse)
library(dplyr)
library(REDCapR)
library(readxl)
# load harmonized libraries
dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))

# read in biopsy master spreadsheet (updated real time)
biopsy_master <- read_excel("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Lab/Biopsy master tracker.xlsx", skip = 2) 
colnames(biopsy_master) <- gsub(" ", "_", colnames(biopsy_master))  # Replace spaces with underscores
colnames(biopsy_master) <- gsub("/", "_", colnames(biopsy_master))  # Replace slashes with underscores
colnames(biopsy_master) <- gsub("\\(", "_", colnames(biopsy_master))  # Replace opening parentheses with underscores
colnames(biopsy_master) <- gsub("\\)", "", colnames(biopsy_master))  # Remove closing parentheses
biopsy_master <- biopsy_master %>%
  dplyr::rename(record_id =Study_ID ,
         visit = Visit_ID) %>%
  mutate(visit = case_when(Study == "ATTEMPT" ~ visit,
                           Study == "REMODEL" ~ visit, 
                           Study == "IMPROVE" & visit == "4M" ~ "3_months_post_surgery",
                           visit == "12M" ~ "12_months_post_surgery",
                           visit == "Follow-up" ~ "follow-up",
                           T ~ "baseline"),
         record_id = case_when(startsWith(record_id, "RPC") & !is.na(Coenroll_ID__Same_visit_only) ~ Coenroll_ID__Same_visit_only,
                               T ~ record_id)) # RH2-60-T and RH2-48-T coenrolled into RPC2 as of 02/03/25 data pull

# RH/RH2/IMPROVE sequential biopsies
seq_dat <- dat %>%
  filter(study %in% c("RENAL-HEIR", "RENAL-HEIRitage", "IMPROVE")) %>%
  distinct(record_id, .keep_all = T) %>%
  dplyr::select(1:12)

seq_biopsies <- biopsy_master %>%
  filter(Study %in% c("RENAL HEIR", "RENAL HEIRITAGE", "IMPROVE")) %>%
  left_join(seq_dat) %>%
  filter(!is.na(Biopsy_date)) %>%
  distinct(paste0(Kit_ID, Biopsy_date), .keep_all=T) %>%
  arrange(mrn, Biopsy_date) %>%
  group_by(mrn) %>%
  filter(n() > 1) %>%
  mutate(date_diff_days = c(NA, diff(Biopsy_date)),
         date_diff_months = date_diff_days/30.41667,
         date_diff_years = date_diff_days/365) %>%
  ungroup() %>%
  dplyr::select(rh_id, rh2_id, improve_id, 
                ends_with("_ID", ignore.case = F), -Pathology_report_ID, Biopsy_date, 
                date_diff_days, date_diff_months, date_diff_years)
# write.csv(seq_biopsies, "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Exports/sequential_biopsies.csv",
#           row.names = F, na = "")


# ID request from Petter of kidney biopsies we have from T2D, OB and HC (RH, RH2, IT2D, CRC) to prioritize for Xenium 5000 & PathoPlex
xenium_patho_ids <- dat %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  filter(study %in% c("RENAL-HEIR", "RENAL-HEIRitage", "IMPROVE", "CROCODILE") & group != "Type 1 Diabetes") %>%
  left_join(biopsy_master) %>%
  filter(Shipped_Y_N =="Yes") %>%
  dplyr::select(record_id, visit, ends_with("_ID", ignore.case = F), -Pathology_report_ID)

# write.csv(xenium_patho_ids, "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Exports/Xenium_Pathoplex_IDs.csv",
#           row.names = F, na = "")

# ID request from Petter of all ATTEMPT + the coenrolled PANDAs with kidney biopsies & scRNA data to prioritize for Xenium 5000 & PathoPlex
# ATTEMPT is not in the harmonized dataset as of Jan 27 2025, need to pull directly from REDCap
# api_tok <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/api_tokens.csv")
# attempt_tok <- api_tok[api_tok$Study == "ATTEMPT",]$Token
# ATTEMPT_rc <- redcap_read(redcap_uri = "https://redcap.ucdenver.edu/api/",
#                           token = attempt_tok)$data
# ATTEMPT_mrn <- ATTEMPT_rc %>%
#   dplyr::select(subject_id, mrn, name_last, name_first) %>% 
#   distinct(subject_id, .keep_all = T) 

# list pulled from seurat object of ATTEMPT "PB_attempt_harmony_rpca_Sept2024.RDS"
attempt_biopsy_ids <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/scRNA/data_clean/ATTEMPT_record_ids.csv")

attempt_biopsy <- attempt_biopsy_ids %>%
  mutate(visit = Visit) %>%
  dplyr::select(subject_id, visit, Cryostor_ID) %>%
  left_join(biopsy_master) %>%
  filter(!is.na(record_id)) %>% #filter out 30058 and 30173 who had biopsies as part of PANDA only, not as part of ATTEMPT (email thread with MH and PB)
  dplyr::rename(panda_id= Coenroll_ID__Same_visit_only ) %>%
  arrange(subject_id) %>%
  dplyr::select(record_id, panda_id, visit, ends_with("_ID", ignore.case = F), -Pathology_report_ID)

# attempt_biopsy_combined <- left_join(attempt_biopsy_ids, ATTEMPT_mrn)

# PANDA is in the harmonized dataset but there is already code for merging ATTEMPT data with PANDA coenrollment from the RPPR reports with REDCap pulls
# panda_tok <- api_tok[api_tok$Study == "PANDA",]$Token
# PANDA_rc <- redcap_read(redcap_uri = "https://redcap.ucdenver.edu/api/",
#                         token = panda_tok)$data
# PANDA_mrn <- PANDA_rc %>% dplyr::select(record_id, mrn) %>% distinct(record_id, .keep_all = T) 
# 
# combined_subset <- left_join(attempt_biopsy_combined, PANDA_mrn)

# write.csv(attempt_biopsy, "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Exports/ATTEMPT_biopsy_IDs.csv",
#           row.names = F, na = "")
