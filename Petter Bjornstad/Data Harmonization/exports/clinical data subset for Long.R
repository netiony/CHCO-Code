library(dplyr)
library(tidyr)
library(Hmisc)
library(purrr)

dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))
st_ids <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Spatial transcriptomics/spatial_transcriptomics_ids.csv")

dat_sub <- dat %>%
  filter(paste0(record_id, visit) %in% paste0(st_ids$record_id, st_ids$visit)) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  mutate(race_ethnicity_condensed = case_when(race == "White" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                              race == "Black or African American" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                              ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                              T ~ "Not Hispanic or Latino Other"),
         sglt2i_ever = case_when(group == "Lean Control" ~ "No", T~ sglt2i_ever)) %>% 
  dplyr::select(record_id, visit, kit_id, cryostor_id, sex, group, diabetes_duration, 
                age, race, ethnicity, weight, height, bmi, acr_u, hba1c, sglt2i_ever, sglti_timepoint,
                eGFR_CKD_epi, sbp, dbp, map)

write.csv(dat_sub, "/Volumes/Peds Endo/Petter Bjornstad/Spatial transcriptomics/spatial_transcriptomics_clin.csv", row.names = F)
