dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))

pen_dat <- dat %>%
  filter(study == "PENGUIN") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by =  c("record_id", "visit")) %>%
  select(record_id, hematocrit_avg)

write.csv(pen_dat, "/Volumes/Peds Endo/Petter Bjornstad/PENGUIN/Data_Cleaned/PEN_hematocrit.csv", row.names = F)
