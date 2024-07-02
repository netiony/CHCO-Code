library(purrr)
library(dplyr)

dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "")
ids <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/AstraZeneca/ACA_UM_COLORADO_URINE_TargetedMetabolomicsResultstoRPC2_240809.csv", na.strings = "") %>%
  dplyr::mutate(record_id = gsub("IT2D-", "IT_", E_CODE_SUBJECT_ID),
                record_id = gsub("_BL", "", record_id),
                record_id = gsub("_12M", "", record_id),
                record_id = gsub("PNDA-", "PNDA-1", record_id),
                record_id = gsub("RH2-38-O", "RH2-38-T", record_id),
                record_id = gsub("RH2-39-T", "RH2-39-O", record_id)) %>%
  dplyr::select(record_id, E_CODE_SUBJECT_ID)

filtered_dat <- dat %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, last(na.omit(.x)))),
                   .by = c(record_id)) %>%
  right_join(ids) %>%
  dplyr::select(E_CODE_SUBJECT_ID, group)


write.csv(filtered_dat, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/AZ_groups.csv",
         row.names = F, na = "")
