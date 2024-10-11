library(dplyr)
library(purrr)
library(Hmisc)
harm_dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/soma_harmonized_dataset.csv", na.strings = "")

dat <- harm_dat %>% 
  filter(study %nin% c("RENAL-HEIR", "RENAL-HEIRitage", "IMPROVE")) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  filter(!is.na(kit_id)) %>%
  filter(!is.na(seq.10000.28))

write.csv(dat, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/biopsies_proteomics.csv",
          row.names = F, na = "")
