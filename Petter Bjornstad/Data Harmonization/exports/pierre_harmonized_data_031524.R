dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))
requested <- readxl::read_xls("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/PJS data request list.xls")

requested <- requested %>%
  filter(!is.na(requested))

pierre_dat <- dat %>%
  filter(study == "CROCODILE" | study == "RENAL-HEIR" | study == "CASPER") %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = record_id) %>%
  filter(participation_status != "Removed (screen failed, withdrawn, dropped out, etc.)") %>% 
  dplyr::mutate(raw_m = p2_raw_m,
                avg_c_f = rowMeans(dplyr::select(., lc_f, rc_f), na.rm = T),
                avg_c_k1 = rowMeans(dplyr::select(., lc_k1, rc_k1), na.rm = T),
                avg_c_k2 = rowMeans(dplyr::select(., lc_k2, rc_k2), na.rm = T),
                avg_c_vb = rowMeans(dplyr::select(., lc_vb, rc_vb), na.rm = T),
                avg_m_f = rowMeans(dplyr::select(., lm_f, rm_f), na.rm = T),
                avg_m_k1 = rowMeans(dplyr::select(., lm_k1, rm_k1), na.rm = T),
                avg_m_k2 = rowMeans(dplyr::select(., lm_k2, rm_k2), na.rm = T),
                avg_m_vb = rowMeans(dplyr::select(., lm_vb, rm_vb), na.rm = T)) %>%
  dplyr::select(requested$variable_name) %>%
  dplyr::select(where(~ any(!is.na(.))))
  # dplyr::mutate(ind1 = if_else(rowSums(!is.na(dplyr::select(., ends_with("_plasma")))) > 0, 1, 0),
  #               ind2 = if_else(rowSums(!is.na(dplyr::select(., ends_with("_ml")))) > 0, 1, 0),
  #               ind3 = ind1+ind2)

# number of participants with both hemodynamic and MRI data
write.csv(pierre_dat, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/pierre_harmonized_dataset_090524.csv", row.names = F)
