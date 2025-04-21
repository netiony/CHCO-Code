library(reticulate)
library(dplyr)
library(tidyr)
library(purrr)

# Import python harmonization function & run
source_python('/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization/data_harmonization.py')
clean <- harmonize_data()
clean <- data.frame(lapply(clean, as.character))
clean[clean == "NaN"] <- NA # Replace NaN from Python to NA
clean[clean == ""] <- NA

# ATTEMPT from Antoine
load(file = "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Clean/ATTEMPT_AC.RData")
names(merged_data)[names(merged_data) %in% names(clean)]
attempt <- merged_data %>%
  mutate(group = "Type 1 Diabetes",
         study = "ATTEMPT",
         diabetes_duration = diabetes_dx_duration,
         procedure = "attempt_visit") %>%
  select(names(merged_data)[names(merged_data) %in% names(clean)], group, study, PWV, treatment_arm, diabetes_duration, procedure, -diabetes_dx_duration)
attempt[] <- lapply(attempt, function(x) {
  if (is.numeric(x)) as.character(x) else x
})

# Create screen_date (screening date for each participant or earliest date available)
clean <- clean %>%
  full_join(attempt) %>%
  dplyr::group_by(record_id, visit) %>% 
  dplyr::mutate(across(where(is.character), ~ na_if(., ""))) %>%
  dplyr::mutate(screen_date = case_when(procedure == "screening" | visit == "screening" ~ date)) %>%
  fill(screen_date, .direction = "updown") %>%
  dplyr::mutate(screen_date = case_when(is.na(screen_date) ~ min(date, na.rm = T), 
                                 T ~ screen_date)) %>%
  fill(screen_date, .direction = "updown") %>% ungroup() %>%
  dplyr::select(record_id, attempt_id, casper_id, coffee_id, croc_id, improve_id, penguin_id, 
                rh_id, rh2_id, panther_id, panda_id, co_enroll_id,
                mrn, date, screen_date, everything())

# Save clinical harmonized dataset
write.csv(clean, 
          "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", row.names = F,
          na = "")

##########################################################################################

# Combine with SOMA
load("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Combined SomaScan/soma_combined_anml_2.RData")

soma <- soma_combined %>%
  as.matrix() %>% as.data.frame() %>%
  mutate_at(vars(starts_with("seq")), as.numeric) %>%
  dplyr::rename("record_id" = SampleDescription) %>%
  dplyr::mutate(record_id = gsub("IT2D-", "IT_", record_id),
         record_id = gsub("PAN_", "PAN-", record_id), # PAN-44-O entered as PAN_44-O in soma
         record_id = gsub("PAN-14-O", "PAN-14-C", record_id), # PAN-14-O to PAN-14-C (group changed, name should have been changed but changed later)
         record_id = gsub("PAN-103-0", "PAN-103-O", record_id), # PAN-103-O entered as PAN-103-0 in soma
         visit = case_when(grepl("PAN-", record_id) ~ SampleGroup,
                           T ~ TimePoint),
         visit = case_when(visit == "Baseline" ~ "baseline",
                           visit == "BL" ~ "baseline",
                           visit == "Year 1" ~ "year_1",
                           visit == "Year 2" ~ "year_2",
                           visit == "Year 3" ~ "year_3",
                           visit == "3M" ~ "3_months_post_surgery",
                           visit == "12M" ~ "12_months_post_surgery",
                           visit == "V1" ~ "baseline",
                           visit == "V4" ~ "4_months_post"),
         visit = case_when(grepl("PAN-", record_id) & visit == "year_1" ~ "baseline", # PANTHER dates skips baseline in SOMA (shifted)
                           grepl("PAN-", record_id) & visit == "year_2" ~ "year_1",
                           grepl("PAN-", record_id) & visit == "year_3" ~ "year_2",
                           T ~ visit)) %>%
  dplyr::select(record_id, visit, starts_with("seq")) %>% 
  mutate_at(vars(starts_with("seq")), as.numeric)

# Check for mismatches between SOMA harmonized and data harmonized (have proteomics but no clinical data)
mismatches <- soma$record_id[!soma$record_id %in% test$record_id]
mismatches

# Combine SOMA & data harmonized
soma_clean <- clean %>% full_join(soma, by = c("record_id", "visit")) %>%
  group_by(mrn, screen_date) %>%
  fill(starts_with("seq"), .direction = "downup") %>%
  ungroup()

# Save harmonized data with SOMA
write.csv(soma_clean, 
          "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/soma_harmonized_dataset.csv", row.names = F,
          na = "")

##########################################################################################

# Combine with Olink
# Olink data
# olink_map <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Olink Data/Data_Clean/olink_id_map.csv")

## Plasma
olink_plasma <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Olink Data/Data_Clean/plasma_cleaned.csv")

olink_plasma <- olink_plasma %>%
  mutate(visit = case_when(endsWith(record_id, "BL") ~ "baseline",
                           endsWith(record_id, "12M") ~ "12_months_post_surgery",
                           T ~ "baseline"),
         record_id = gsub("_BL|_12M", "", record_id)) %>%
  rename_with(~ paste0(.x, "_p"), starts_with("OID"))

olink_plasma_clean <- clean %>% left_join(olink_plasma, by = c("record_id", "visit")) %>%
  group_by(mrn, screen_date) %>%
  fill(starts_with("OID"), .direction = "downup") %>%
  ungroup() 

# Save harmonized data with Olink Plasma
write.csv(olink_plasma_clean, 
          "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/olink_plasma_harmonized_dataset.csv", row.names = F,
          na = "")

## Urine
olink_urine <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Olink Data/Data_Clean/urine_cleaned.csv")
olink_urine <- olink_urine %>%
  mutate(visit = case_when(endsWith(record_id, "BL") ~ "baseline",
                           endsWith(record_id, "12M") ~ "12_months_post_surgery",
                           T ~ "baseline"),
         record_id = gsub("_BL|_12M", "", record_id)) %>%
  rename_with(~ paste0(.x, "_u"), starts_with("OID"))

olink_urine_clean <- clean %>% left_join(olink_urine, by = c("record_id", "visit")) %>%
  group_by(mrn, screen_date) %>%
  fill(starts_with("OID"), .direction = "downup") %>%
  ungroup() 

# Save harmonized data with Olink Urine
write.csv(olink_urine_clean, 
          "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/olink_urine_harmonized_dataset.csv", row.names = F,
          na = "")

# Save harmonized data with Olink (both plasma and urine) and SOMA
clean_proteomics_comb <- clean %>% 
  left_join(olink_plasma, by = c("record_id", "visit")) %>%
  left_join(olink_urine, by = c("record_id", "visit")) %>%
  left_join(soma, by = c("record_id", "visit")) %>%
  group_by(mrn, screen_date) %>% 
  fill(starts_with("OID"), .direction = "downup") %>%
  fill(starts_with("seq"), .direction = "downup")
write.csv(clean_proteomics_comb, 
          "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/soma_olink_harmonized_dataset.csv", row.names = F,
          na = "")

#### INDICATORS FOR TABLEAU ####
clean_ind <- clean %>% 
  left_join(soma, by = c("record_id", "visit")) %>%
  dplyr::mutate(soma_ind = if_else(rowSums(!is.na(dplyr::select(., starts_with("seq")))) > 0, 
                                   1,0),
                ivgtt_ind = if_else(grepl("ivgtt", procedure) & !is.na(date), 
                                   1,0),
                minmod_ind = if_else(rowSums(!is.na(dplyr::select(., starts_with("mm_")))) > 0, 
                                     1,0),
                clamp_ind = if_else(procedure == "clamp" & rowSums(!is.na(dplyr::select(., starts_with("glucose_1")))) > 0, 
                                    1,0),
                mmtt_ind = if_else(procedure == "mmtt" & rowSums(!is.na(dplyr::select(., starts_with("glucose_1")))) > 0, 
                                   1,0),
                kidney_biopsy_ind = if_else(procedure == "kidney_biopsy" & !is.na(date), 
                                    1,0),
                rc_ind = if_else(procedure == "renal_clearance_testing" & !is.na(date), 
                                1,0),
                hemo_ind = if_else(rowSums(!is.na(dplyr::select(., starts_with("gfr_raw")))) > 0, 
                                   1,0),
                dxa_ind = if_else(procedure == "dxa" & !is.na(date), 
                                  1,0),
                mri_ind = if_else(procedure == "bold_mri" & !is.na(date), 
                                  1,0),
                pet_ind = if_else(procedure == "pet_scan" & !is.na(date), 
                                  1,0),
                brain_bio_ind = if_else(rowSums(!is.na(dplyr::select(., ends_with("_avg_conc")))) > 0, 
                                        1,0),
                morpho_ind = if_else(rowSums(!is.na(dplyr::select(., starts_with("mes_")))) > 0, 
                                     1,0),
                metab_tissue_ind = if_else(rowSums(!is.na(dplyr::select(., ends_with("_tissue")))) > 0, 
                                           1,0),
                metab_blood_ind = if_else(!is.na(alanine), 
                                          1,0),
                hispanic_ind = if_else(grepl("^Hispanic", ethnicity), 
                                       1,0),
                non_hispanic_ind = if_else(grepl("^Not", ethnicity), 
                                           1,0)) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, max(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, max(.x, na.rm = TRUE))),
                   .by = c(mrn, screen_date)) %>%
  ungroup() %>%
  dplyr::select(record_id, ends_with("ind")) %>%
  distinct(record_id, .keep_all = T)

# Save indicators
write.csv(clean_ind, 
          "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/procedure_indicator.csv", row.names = F,
          na = "")
