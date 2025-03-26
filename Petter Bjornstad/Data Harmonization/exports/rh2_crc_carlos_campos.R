# Clinical data for Carlos (example dataset for MRI for DECODE)
# Load clinical data
harm_dat <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "")

rh2_crc <- harm_dat %>% filter(study == "CROCODILE" | study == "RENAL-HEIRitage") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id)) %>%
  dplyr::mutate(race_ethnicity_condensed = case_when(race == "White" & 
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                                     race == "Black or African American" & 
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                                     ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                                     T ~ "Not Hispanic or Latino Other")) %>%
  arrange(record_id) %>%
  select(-mrn)

write.csv(rh2_crc, "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Exports/harmonized_dataset_rh2_crc.csv", row.names = F, na = "")
