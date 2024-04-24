library(dplyr)
library(tidyr)
library(Hmisc)
library(purrr)
library(REDCapR)

dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))
requested <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/Raji_sc_metadata_request.csv")
tokens <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/api_tokens.csv")
attempt_token <- subset(tokens, Study == "ATTEMPT")$Token
attempt <- redcap_read(redcap_uri = "https://redcap.ucdenver.edu/api/",
                       token = attempt_token)$data %>%
  filter(subject_id == 30051)

dat <- dat %>%
  mutate(record_id_2 = gsub("-O$|-T$|-L$|-C$", "", record_id))

requested_dat <- requested %>%
  mutate(record_id_2 = gsub("PNDA-", "PNDA-1", record_id_2),
         record_id_2 = gsub("-O$|-T$|-L$|-C$", "", record_id_2),
         visit = case_when(visit == "4M" ~ "post_treatment", T ~ visit)) %>%
  left_join(dat, by = c("record_id_2", "visit")) %>% 
  mutate(record_id = case_when(is.na(record_id) ~ record_id_2, T~ record_id)) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  dplyr::select(croyster_id, kit_id, record_id, record_id_2, visit, group, sex, age, sglti_timepoint)

write.csv(requested_dat, "/Users/choiyej/Raji_sc_metadata_request_filled.csv")
