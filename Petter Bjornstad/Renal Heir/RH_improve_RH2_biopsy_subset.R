dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))
coenroll_id <- read.csv("/Users/choiyej/GitHub/YC_CHCO/RH2/coenrolled_ids_LR.csv")
dat <- dat %>%
  filter(study %in% c("RENAL-HEIR", "IMPROVE", "RENAL-HEIRitage")) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = record_id) %>%
  mutate(race_ethnicity_condensed = case_when(race == "White" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                              race == "Black or African American" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                              ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                              T ~ "Not Hispanic or Latino Other")) %>%
  mutate(elevated_albuminuria = case_when(elevated_albuminuria == "Yes" ~ "Elevated albuminuria",
                                          elevated_albuminuria == "No" ~ "Normoalbuminuria")) %>% 
  mutate_at(vars(starts_with("fsoc")), function(x) case_when(x < 0 ~ 0, T~x)) 

biopsy_dat <- dat %>%
  dplyr::select(record_id, co_enroll_id, study, group, kidney_side, kidney_location) %>%
  filter(!is.na(kidney_side))

write.csv(biopsy_dat, "/Volumes/Peds Endo/Petter Bjornstad/Renal HEIR/Data_Cleaned/biopsy_available.csv")
