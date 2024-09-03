library(dplyr)
library(purrr)
harm_dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "")

ips_ids <- c("CRC-03", "RH-62-T", 
"CRC-01", "CRC-06",
"CRC-07", "CRC-10",
"CRC-11", "CRC-12",
"CRC-13", "CRC-14", 
"IT_15", "IT_19", 
"RH-50-T", "RH-68-T", 
"RH-72-T", "RH-77-T",
"RH-81-T", "RH-87-T")

dat <- harm_dat %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id)) %>%
  dplyr::mutate(age = floor(age)) %>%
  filter(sex == "Male") %>%
  filter(group == "Type 2 Diabetes" | group == "Lean Control") %>% 
  filter(record_id %in% ips_ids) %>%
  dplyr::select(record_id, group, age, sex, race_ethnicity)

write.csv(dat, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/ipsc_male_subset.csv",
          row.names = F, na = "")
