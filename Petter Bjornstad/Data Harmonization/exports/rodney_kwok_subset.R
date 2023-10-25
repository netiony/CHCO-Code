dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))
coenroll_id <- read.csv("/Users/choiyej/GitHub/YC_CHCO/RH2/coenrolled_ids_LR.csv")

dat_subset <- dat %>%
  filter(study %in% c("RENAL-HEIR", "RENAL-HEIRitage", "CROCODILE",
                      "PANDA", "ULTRA", "PENGUIN")) %>%
  dplyr::mutate(merged_id = ifelse(test = !is.na(match(record_id, coenroll_id$rh_id)), 
                                   yes = coenroll_id[match(record_id, coenroll_id$rh_id), 1],
                                   no = ifelse(test = !is.na(match(record_id, coenroll_id$rh2_id)), 
                                               yes = coenroll_id[match(record_id, coenroll_id$rh2_id), 1],
                                               no = ifelse(test = !is.na(match(record_id, coenroll_id$improve_id)), 
                                                           yes = coenroll_id[match(record_id, coenroll_id$improve_id), 1],
                                                           no = record_id)))) %>%
  select(record_id,
         merged_id,
         age,
         sex,
         gfr_bsa_plasma,
         gfr_raw_plasma,
         visit,
         date,
         creatinine_s,
         acr_u,
         height,
         weight,
         sbp,
         dbp,
         group,
         diabetes_dx_date,
         hba1c,
         cholesterol,
         hdl,
         ldl,
         triglycerides,
         ethnicity) %>%
  group_by(merged_id, date) %>%
  fill(colnames(dat_subset), .direction = "downup")%>%
  filter(!is.na(age)) %>%
  distinct(merged_id, date, .keep_all = T) %>% ungroup()
%>%
  select(-merged_id) %>%
  filter(rowSums(!is.na(.)) >= 9)

write.csv(dat_subset, 
          "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/rodney_kwok.csv", row.names = F)

