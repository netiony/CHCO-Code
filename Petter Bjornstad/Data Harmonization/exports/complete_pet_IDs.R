# Import data and data dictionary
dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))
dict <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/data_dictionary_master.csv", na.strings = c(" ", "", "-9999",-9999)) %>%
  dplyr::select(variable_name, label)

# collapse to 1 row per participant
dat <- dat %>%
  filter(study != "PENGUIN") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id)) %>%
  filter(participation_status!="Removed"|is.na(participation_status)) %>%
  mutate(race_ethnicity_condensed = case_when(race == "White" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                              race == "Black or African American" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                              ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                              T ~ "Not Hispanic or Latino Other")) %>%
  mutate(elevated_albuminuria = case_when(elevated_albuminuria == "Yes" ~ "Elevated albuminuria",
                                          elevated_albuminuria == "No" ~ "Normoalbuminuria")) %>% 
  mutate_at(vars(starts_with("fsoc")), function(x) case_when(x < 0 ~ 0, T~x)) 

pet_dat <- dat %>%
  filter(if_any(c(lc_k1, rc_k1, lc_f, rc_f, lc_k2, rc_k2), Negate(is.na))) %>%
  dplyr::select(record_id) 

write.csv(pet_dat,
          "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/complete_pet_id.csv", row.names = F)
