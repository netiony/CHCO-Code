library(dplyr)
library(purrr)

source("~/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization/sto2_calculation_function.R")
harm_dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "")

dat <- harm_dat %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  dplyr::mutate(avg_c_vb = rowMeans(dplyr::select(., lc_vb, rc_vb), na.rm = T),
                avg_m_vb = rowMeans(dplyr::select(., lm_vb, rm_vb), na.rm = T),
                hct = coalesce(hematocrit_210, hct),
                sto2_c = sto2_cor(vb_cor = avg_c_vb,
                                    hct = hct,
                                    r2star_cor = avg_c_r2),
                sto2_m = sto2_med(vb_med = avg_m_vb,
                                    hct = hct,
                                    r2star_med = avg_m_r2)) %>%
  filter(!is.na(sto2_c)) %>%
  dplyr::select(record_id, visit, sto2_c, sto2_m)

write.csv(dat, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/sto2_calculations.csv",
          row.names = F, na = "")
