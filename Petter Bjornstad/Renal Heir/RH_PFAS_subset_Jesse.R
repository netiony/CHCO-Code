library(readxl)
library(dplyr)
library(purrr)
dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings=c(""," ","NA"))

vars <- readxl::read_excel("/Volumes/Peds Endo/Petter Bjornstad/Renal HEIR/Data_Raw/Copy of Renal_HEIR variables PFAS DKD pilot.xlsx",
                   range = "C1:C100", col_names = TRUE) %>%
  filter(!is.na(...1))
vars <- as.vector(vars)[[1]]

ids <- read_xlsx("/Volumes/Peds Endo/Petter Bjornstad/RENAL HEIR/Data_Cleaned/Copy of phenotype_t_pfas_u_prot.xlsx",
                 range = "B1:B150")

pfas_subset <- dat %>%
  filter(study == "RENAL-HEIR" | study == "RENAL-HEIRitage" | study == "IMPROVE" | study == "CROCODILE") %>%
  dplyr::select(c(record_id, visit, study, group, vars, eGFR_fas_cr,
                  di, gir_190)) %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  filter(visit != "3_months_post_surgery") %>%
  dplyr::mutate(record_id = case_when(visit == "12_months_post_surgery" ~ paste0(record_id, "_12M"),
                                      visit == "baseline" & study == "IMPROVE" ~ paste0(record_id, "_BL"),
                                      T ~ record_id),
                m_i = gir_190/steady_state_insulin) %>%
  dplyr::select(-c(study, visit)) %>%
  filter(record_id %in% ids$`Phenotype IDs`)

write.csv(pfas_subset, "/Volumes/Peds Endo/Petter Bjornstad/Renal HEIR/Data_Cleaned/RH PFAS subset.csv", row.names = F)
