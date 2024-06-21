library(purrr)
library(dplyr)
library(lubridate)

# Phoom data subset from RH, RH2, IT2D, and CRC

dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "")

filtered_dat <- dat %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, last(na.omit(.x)))),
                   .by = c(record_id)) %>%
  filter(study %in% c("RENAL-HEIR", "RENAL-HEIRitage", "IMPROVE", "CROCODILE")) %>%
  dplyr::select(record_id, mrn)


write.csv(filtered_dat, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/rh_rh2_imp_crc_mrn_phoom.csv",
          row.names = F, na = "")
