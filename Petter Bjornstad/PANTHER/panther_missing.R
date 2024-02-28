library(REDCapR)
library(dplyr)

panther_meta <- REDCapR::redcap_metadata_read(redcap_uri = "https://redcap.ucdenver.edu/api/",
                     token = "1EC857A51835CF17893236DB1262E77E")$data
form_names <- unique(panther_meta$form_name)
screening <- form_names[1:4]
baseline <- form_names[4:15]
year1_2 <- form_names[5:15]

panther_screening <- REDCapR::redcap_read(redcap_uri = "https://redcap.ucdenver.edu/api/",
                     token = "1EC857A51835CF17893236DB1262E77E",
                     raw_or_label = "label", export_checkbox_label = T,
                     events = "screening_arm_1",
                     forms = screening)$data %>%
  dplyr::mutate(race = coalesce(!!!select(., starts_with("race"))),
                ethnicity = coalesce(!!!select(., starts_with("ethnicity"))),
                steroid_type = coalesce(!!!select(., starts_with("steroid_type"))),
                antipsychotic_type = coalesce(!!!select(., starts_with("antipsychotic_type"))),
                htn_med = coalesce(!!!select(., starts_with("htn_med"))),
                hx_cv_positive = coalesce(!!!select(., starts_with("hx_cv_positive"))),
                hx_met_positive = coalesce(!!!select(., starts_with("hx_met_positive")))) %>%
  relocate(race, .before = race___1) %>%
  relocate(ethnicity, .before = ethnicity___1) %>%
  relocate(steroid_type, .before = steroid_type___1) %>%
  relocate(antipsychotic_type, .before = antipsychotic_type___1) %>%
  relocate(htn_med, .before = htn_med___1) %>%
  relocate(hx_cv_positive, .before = hx_cv_positive___1) %>%
  relocate(hx_met_positive, .before = hx_met_positive___1) %>%
  select(-matches("^race_|^ethnicity_|^steroid_type_|^antipsychotic_type_|^htn_med_|^hx_cv_positive_|^hx_met_positive_")) %>%
  mutate_at(vars(-1,-2), ~ ifelse(is.na(.), "MISSING", "")) %>%
  dplyr::select(where(~any(. != "")))


panther_baseline <- REDCapR::redcap_read(redcap_uri = "https://redcap.ucdenver.edu/api/",
                                         token = "1EC857A51835CF17893236DB1262E77E",
                                         raw_or_label = "label", export_checkbox_label = T,
                                         events = "baseline_arm_1",
                                         fields = "record_id",
                                         forms =  baseline)$data %>%
  mutate_at(vars(-1,-2), ~ ifelse(is.na(.), "MISSING", "")) %>%
  dplyr::select(where(~any(. != "")))

panther_yr1_2 <- REDCapR::redcap_read(redcap_uri = "https://redcap.ucdenver.edu/api/",
                                      token = "1EC857A51835CF17893236DB1262E77E",
                                      raw_or_label = "label", export_checkbox_label = T,
                                      events = c("year_1_arm_1", "year_2_arm_1"),
                                      fields = "record_id",
                                      forms =  year1_2)$data %>%
  mutate_at(vars(-1,-2), ~ ifelse(is.na(.), "MISSING", "")) %>%
  dplyr::select(where(~any(. != "")))

write.csv(panther_screening, "/Volumes/Peds Endo/Petter Bjornstad/PANTHER/Data_Cleaned/panther_screening.csv", row.names = F, na = "")
write.csv(panther_baseline, "/Volumes/Peds Endo/Petter Bjornstad/PANTHER/Data_Cleaned/panther_baseline.csv", row.names = F, na = "")
write.csv(panther_yr1_2, "/Volumes/Peds Endo/Petter Bjornstad/PANTHER/Data_Cleaned/panther_year1_2.csv", row.names = F, na = "")
