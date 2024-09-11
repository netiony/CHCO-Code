library(dplyr)
library(purrr)

# Subsetting for Hailey and Jesse for PAH and OATS analysis

harm_dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "")

dat <- harm_dat %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  dplyr::mutate(pah_coerce = (coalesce(pah_raw, pah_bsa, pah_clear_abs, pah_clear_bsa))) %>%
  filter(!is.na(kit_id)) %>%
  filter(!is.na(pah_coerce)) %>%
  dplyr::select(record_id, group, visit, sex, age, race, ethnicity, weight, height, bmi, hba1c, acr_u)

write.csv(dat, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/renal_clearance_biopsy.csv",
          row.names = F, na = "")
