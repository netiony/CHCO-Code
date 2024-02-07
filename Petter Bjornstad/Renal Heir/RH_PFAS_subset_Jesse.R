library(readxl)
library(dplyr)
library(purrr)
dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings=c(""," ","NA"))

vars <- readxl::read_excel("/Volumes/Peds Endo/Petter Bjornstad/Renal HEIR/Data_Raw/Copy of Renal_HEIR variables PFAS DKD pilot.xlsx",
                   range = "C1:C100", col_names = TRUE) %>%
  filter(!is.na(...1))
vars <- as.vector(vars)[[1]]

pfas_subset <- dat %>%
  filter(study == "RENAL-HEIR") %>%
  dplyr::select(c(record_id, group, vars)) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id)) 

write.csv(pfas_subset, "/Volumes/Peds Endo/Petter Bjornstad/Renal HEIR/Data_Cleaned/RH PFAS subset.csv", row.names = F)
