library(dplyr)
library(REDCapR)

dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))

dat <- dat %>%
  filter(study == "RENAL-HEIR") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = record_id) %>%
  filter(participation_status=="Participated") %>%
  mutate(race_ethnicity_condensed = case_when(race == "White" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                              race == "Black or African American" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                              ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                              T ~ "Not Hispanic or Latino Other")) %>%
  mutate(elevated_albuminuria = case_when(elevated_albuminuria == "Yes" ~ "Elevated albuminuria",
                                          elevated_albuminuria == "No" ~ "Normoalbuminuria")) %>% 
  mutate_at(vars(starts_with("fsoc")), function(x) case_when(x < 0 ~ 0, T~x)) 

kvolume_missing <- dat %>%
  filter(is.na(left_kidney_volume_ml|right_kidney_volume_ml)) %>%
  dplyr::select(record_id, group, date)

redcap_dat <- redcap_project$new(redcap_uri = "https://redcap.ucdenver.edu/api/", 
                  token = "476F5830A52A4E79672E8A47A94C869F")$read()$data

mrn_subset <- redcap_dat %>%
  dplyr::select(subject_id, mr_number) %>%
  distinct(subject_id, .keep_all = T) %>%
  rename(record_id = subject_id)

kvolume_missing_mrn <- kvolume_missing %>%
  left_join(mrn_subset) %>%
  dplyr::select(mr_number, record_id, group, date)

write.csv(kvolume_missing_mrn,
          "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/kidney_volume_missing_RH_Phoom.csv" ,row.names = F)
