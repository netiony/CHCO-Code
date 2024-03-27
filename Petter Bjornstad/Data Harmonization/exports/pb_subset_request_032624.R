# run a report through the harmonized dataset / redcap providing a list or patient names / MRNs / IDs for participants that meet these criteria:
# Age >= 18 years of age
# BMI >= 27 kg/m2
# UACR >= 30 mg/g
# eGFR <75 ml/min per 1.73m2 OR serum creatinine >= 1.0 
# It is OK if they meet 3 and/or 4, i.e., don't necessarily need to meet both 

dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "")

filtered_dat <- dat %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(mrn, screen_date)) %>%
  dplyr::mutate(ind = case_when(age >= 18 & bmi >= 27 & (acr_u >= 30 | eGFR_CKD_epi < 75 | creatinine_s >= 1) ~ 1,
                         T ~ 0)) %>%
  filter(ind == 1) %>%
  dplyr::select(casper_id, coffee_id, croc_id, improve_id, penguin_id, rh_id, rh2_id, panther_id, panda_id,
                mrn, age, bmi, acr_u, eGFR_CKD_epi, creatinine_s)

write.csv(filtered_dat, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/petter_subset_032624.csv",
          row.names = F, na = "")
