library(dplyr)
library(purrr)
dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))

di_dat <- dat %>%
  filter(study == "IMPROVE") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c("record_id", "visit")) %>%
  select(record_id, visit, group, age, race, ethnicity, bmi, raw_m, steady_state_insulin, ffa_suppression, airg, acprg, di)

write.csv(di_dat, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/improve_clamp.csv", row.names = F)
