dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))

phil_dat <- dat %>%
  filter(study == "CROCODILE" | study == "RENAL-HEIR") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = record_id) %>%
  filter(participation_status != "Removed (screen failed, withdrawn, dropped out, etc.)") %>% 
  select(record_id, group, study, triglycerides, hdl, ldl, 
         cholesterol, starts_with("glom"), fia, starts_with("mes"), starts_with("pod"))

write.csv(phil_dat, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/crc_rh_subset_for_phil.csv", row.names = F)