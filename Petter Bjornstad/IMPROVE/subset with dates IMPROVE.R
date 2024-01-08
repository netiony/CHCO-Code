dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))

imp_dat <- dat %>%
  filter(study == "IMPROVE") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by =  c("record_id", "visit")) %>%
  select(record_id, visit, date) %>%
  arrange(record_id, date)

write.csv(imp_dat, "/Volumes/Peds Endo/Petter Bjornstad/IMPROVE T2D/Data_Cleaned/IMPROVE_dates.csv", row.names = F)
