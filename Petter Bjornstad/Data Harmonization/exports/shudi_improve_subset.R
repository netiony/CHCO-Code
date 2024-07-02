library(purrr)
library(dplyr)
library(lubridate)

# Shudi Pan data subset from IMPROVE for AHA fellowship

dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "")
blood_met <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/IMPROVE T2D/Data_Raw/blood_met_id.csv", na.strings = "")

mround <- function(x,base){
  base*round(x/base)
}

filtered_dat <- dat %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, last(na.omit(.x)))),
                   .by = c(mrn, visit)) %>%
  dplyr::group_by(mrn) %>%
  dplyr::mutate(date_diff = round(as.period(interval(min(date), date)) / months(1),1)) %>%
  filter(!is.na(improve_id)) %>% ungroup() %>%
  dplyr::mutate(date_diff = mround(date_diff, 6)) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, last(na.omit(.x)))),
                   .by = c(mrn, date_diff)) %>%
  dplyr::group_by(improve_id) %>%
  dplyr::mutate(age_baseline = floor(min(age, na.rm = T)),
                surgery_year = min(year(date)),
                ace_inhibitor = case_when(hypertension == "No" & is.na(ace_inhibitor) ~ "No",
                                          T ~ ace_inhibitor),
                beta_blocker = case_when(hypertension == "No" & is.na(beta_blocker) ~ "No",
                                          T ~ beta_blocker),
                diuretic = case_when(hypertension == "No" & is.na(diuretic) ~ "No",
                                          T ~ diuretic)) %>%
  dplyr::select(improve_id, race, ethnicity, hypertension, ace_inhibitor, beta_blocker, diuretic, age_baseline, date_diff, surgery_year, sex, sbp, dbp, bmi, tolower(blood_met$compound_id)) %>%
  dplyr::arrange(date_diff, .by_group = T) %>%
  fill(c(ace_inhibitor, beta_blocker, diuretic, hypertension), .direction = "updown") %>%
  dplyr::rename(visit = date_diff)

write.csv(filtered_dat, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/improve_subset_shubi_051524.csv",
         row.names = F, na = "")
