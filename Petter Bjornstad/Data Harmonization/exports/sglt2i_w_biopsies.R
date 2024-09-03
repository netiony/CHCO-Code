library(dplyr)
library(purrr)
harm_dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "")

dat <- harm_dat %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  filter(!is.na(kit_id)) %>%
  dplyr::mutate(sglt2i_ever = case_when(study == "PANDA" | study == "CROCODILE" ~ "No", T ~ sglt2i_ever),
                exclude = case_when(study != "CROCODILE" & !is.na(croc_id) ~ T)) %>%
  filter(is.na(exclude)) %>%
  dplyr::select(record_id, group, visit, sglt2i_ever, sglti_timepoint)

write.csv(dat, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/sglt2i_status.csv",
          row.names = F, na = "")
