harm_dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "")
meta_dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Biopsies/Copy of Colorado Meta Data_90.csv") %>%
  dplyr::mutate(record_id_3 = paste0(gsub("-T|-O", "", record_id), visit))

dat <- harm_dat %>% 
  dplyr::mutate(record_id_3 = paste0(gsub("-T|-O", "", record_id), visit)) %>%
  filter(record_id_3 %in% meta_dat$record_id_3) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id_3, visit)) %>%
  dplyr::mutate(race_ethnicity_condensed = case_when(race == "White" & 
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                                     race == "Black or African American" & 
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                                     ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino")) %>%
  dplyr::select(record_id_3, hba1c, eGFR_CKD_epi)

join_dat <- left_join(meta_dat, dat) 

write.csv(join_dat, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Biopsies/Copy of Colorado Meta Data_90_YC.csv",
          row.names = F, na = "")
