dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv")
coenroll_id <- read.csv("/Users/choiyej/GitHub/YC_CHCO/RH2/coenrolled_ids_LR.csv")

sub_dat <- dat %>%
  filter(group != "Type 2 Diabetes") %>%
  dplyr::mutate(merged_id = ifelse(test = !is.na(match(record_id, coenroll_id$rh_id)), 
                                   yes = coenroll_id[match(record_id, coenroll_id$rh_id), 1],
                                   no = ifelse(test = !is.na(match(record_id, coenroll_id$rh2_id)), 
                                               yes = coenroll_id[match(record_id, coenroll_id$rh2_id), 1],
                                               no = ifelse(test = !is.na(match(record_id, coenroll_id$improve_id)), 
                                                           yes = coenroll_id[match(record_id, coenroll_id$improve_id), 1],
                                                           no = record_id)))) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = merged_id) %>%
  select(record_id, co_enroll_id, merged_id, group, age, sex, race, ethnicity,
          gfr_bsa_plasma, gfr_raw_plasma, pah_clear_bsa, 
          acr_u, lc_k1, rc_k1, lm_k1, rm_k1,
          lc_k2, rc_k2, lm_k2, rm_k2) %>%
  filter(!is.na(lc_k1|rc_k1|lm_k1|rm_k1)) 

write.csv(sub_dat, file = "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/MG_petgrant.csv", row.names = F)
