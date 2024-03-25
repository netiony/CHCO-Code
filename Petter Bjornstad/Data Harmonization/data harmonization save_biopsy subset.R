library(dplyr)
library(tidyr)
library(table1)
library(arsenal)
library(Hmisc)
library(purrr)

dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))

biopsy_dat <- dat %>%
  filter(procedure == "kidney_biopsy") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  filter(participation_status != "Removed"|is.na(participation_status)) %>%
  mutate(race_ethnicity_condensed = case_when(race == "White" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                              race == "Black or African American" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                              ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                              T ~ "Not Hispanic or Latino Other"),
         obtained_lm_em = if_else(rowSums(across(ends_with("core_spec___1"))) > 0, "Yes", "No"),
         obtained_cryo_hypotherm = if_else(rowSums(across(ends_with("core_spec___3"))) > 0, "Yes", "No"),
         obtained_oct = if_else(rowSums(across(ends_with("core_spec___2"))) > 0, "Yes", "No"),
         obtained_rnalater = if_else(rowSums(across(ends_with("core_spec___4"))) > 0, "Yes", "No"),
         biopsy_yn = if_else(rowSums(!is.na(select(., date, kit_id, glut_id, form_id, rnalater_id, cryomold_id, cryostor_id, ln2_id))) > 0, "Yes", "No")) %>%
  dplyr::select(record_id, co_enroll_id, biopsy_yn, date, group, visit, sglt2i_ever, sglti_timepoint, 
                kit_id,glut_id, form_id, rnalater_id, cryomold_id, cryostor_id, ln2_id) %>%
  arrange(desc(biopsy_yn), record_id)

write.csv(biopsy_dat, 
          "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/harmonized_dataset_biopsy_subset.csv", 
          row.names = F,
          na = "")
