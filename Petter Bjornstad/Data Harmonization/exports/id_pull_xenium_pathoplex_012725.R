# load libraries
library(tidyverse)
library(dplyr)
library(REDCapR)

# load harmonized libraries
dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))

# ID request from Petter of kidney biopsies we have from T2D, OB and HC (RH, RH2, IT2D, CRC) to prioritize for Xenium 5000 & PathoPlex
xenium_patho_ids <- dat %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  filter(study %in% c("RENAL-HEIR", "RENAL-HEIRitage", "IMPROVE", "CROCODILE") & group != "Type 1 Diabetes") %>%
  filter(!is.na(kit_id)) %>%
  dplyr::select(record_id)
write.csv(xenium_patho_ids, "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Exports/Xenium_Pathoplex_IDs.csv",
          row.names = F)


# ID request from Petter of all ATTEMPT + the coenrolled PANDAs with kidney biopsies & scRNA data to prioritize for Xenium 5000 & PathoPlex
# ATTEMPT is not in the harmonized dataset as of Jan 27 2025, need to pull directly from REDCap
api_tok <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/api_tokens.csv")
attempt_tok <- api_tok[api_tok$Study == "ATTEMPT",]$Token
ATTEMPT_rc <- redcap_read(redcap_uri = "https://redcap.ucdenver.edu/api/",
                          token = attempt_tok)$data
ATTEMPT_mrn <- ATTEMPT_rc %>%
  dplyr::select(subject_id, mrn, name_last, name_first) %>% 
  distinct(subject_id, .keep_all = T) 
ATTEMPT_mrn$subject_id <- as.character(ATTEMPT_mrn$subject_id)

# PANDA is in the harmonized dataset but there is already code for merging ATTEMPT data with PANDA coenrollment from the RPPR reports with REDCap pulls
panda_tok <- api_tok[api_tok$Study == "PANDA",]$Token
PANDA_rc <- redcap_read(redcap_uri = "https://redcap.ucdenver.edu/api/",
                        token = panda_tok)$data
PANDA_mrn <- PANDA_rc %>% select(record_id, mrn, name_last, name_first) %>% distinct(record_id, .keep_all = T) 

combined_mrn <- left_join(ATTEMPT_mrn, PANDA_mrn)

# list pulled from seurat object of ATTEMPT "PB_attempt_harmony_rpca_Sept2024.RDS"
attempt_biopsy_ids <- c(30001, 30013, 30051, 30020, 30132, 30184, 30058, 30073, 30191, 30098, 30245, 30230, 30403, 30472, 30193, 30173)

combined_mrn_subset <- combined_mrn %>%
  filter(subject_id %in% attempt_biopsy_ids) %>%
  select(subject_id, record_id) %>%
  rename(ATTEMPT_ID = subject_id,
         PANDA_ID = record_id)
write.csv(combined_mrn_subset, "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Exports/ATTEMPT_biopsy_IDs.csv",
          row.names = F, na = "")
