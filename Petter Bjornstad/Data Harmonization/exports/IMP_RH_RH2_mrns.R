library(dplyr)
library(REDCapR)
library(magrittr)

redcap_rh <- redcap_project$new(redcap_uri = "https://redcap.ucdenver.edu/api/", 
                                 token = "476F5830A52A4E79672E8A47A94C869F")$read()$data

redcap_rh2 <- redcap_project$new(redcap_uri = "https://redcap.ucdenver.edu/api/", 
                                token = "757A974510ABE1B2CEF124847818B959")$read()$data

redcap_improve <- redcap_project$new(redcap_uri = "https://redcap.ucdenver.edu/api/", 
                                token = "57A4C81389BCF44694CE4170FB534690")$read()$data

redcap_rh %<>%
  dplyr::select(subject_id, mr_number)

redcap_rh2 %<>%
  dplyr::select(record_id, mrn) %>%
  rename(subject_id = record_id,
         mr_number = mrn)

redcap_improve %<>%
  dplyr::select(subject_id, mr_number)

combined_mrn <- rbind(redcap_rh, redcap_rh2, redcap_improve) %>%
  distinct(mr_number, .keep_all = T)

write.csv(combined_mrn, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/imp_rh_rh2_mrn.csv", row.names = F)
