library(dplyr)

# PANTHER raw lab data cleaning
panther_dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "") %>% filter(study == "PANTHER") %>%
  dplyr::select(mrn, record_id) %>% distinct(record_id, .keep_all = T)
panther_lab <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/PANTHER/Data_Raw/PANTHER glucose_insulin labs 030124.csv", na.strings = "")

panther_lab_tst <- panther_lab %>%
  filter(grepl("gluc|ins", TestName, ignore.case=T)) %>%
  mutate(TestName = gsub("\\s+", "", TestName)) %>%
  mutate(test = gsub("glucose", "gluc_", TestName, ignore.case = T),
         test = gsub("insulin", "ins_", test, ignore.case = T),
         test = gsub("-", "minus_", test, ignore.case = T),
         test = gsub("ins_minus_endo", "ins_", test, ignore.case = T),
         test = gsub("min$", "", test, ignore.case = T)) %>%
  rename("mrn" = MRN, "date" = Collection.Date, "value" = Result.Numeric) %>%
  dplyr::select(mrn, date, value, test) %>% arrange(mrn) %>%
  mutate(date = as.Date(date, format = "%m/%d/%y")) %>% 
  group_by(mrn) %>%
  mutate(visit = case_when(date == min(date) ~ "baseline",
                           date == max(date) ~ "year_1",
                           T ~ "")) %>% ungroup()

panther_lab_clean <- spread(panther_lab_tst, key = test, value = value) %>%
  left_join(panther_dat, by = "mrn") %>%
  dplyr::select(-c(gluc_, date)) %>%
  dplyr::select(record_id, mrn, visit, everything()) %>%
  filter(visit != "")
panther_lab_clean[panther_lab_clean==-999.99] <- NA
write.csv(panther_lab_clean, "/Volumes/Peds Endo/Petter Bjornstad/PANTHER/Data_Cleaned/PANTHER glucose_insulin labs 030124_clean.csv", row.names = F, na = "")

## Import from RedCap
library(REDCapR)
ivgtt_dat <- redcap_read(redcap_uri = "https://redcap.ucdenver.edu/api/", 
                         token = "1EC857A51835CF17893236DB1262E77E", 
                         fields = c("record_id"), 
                         forms = c("ivgtt_labs"), 
                         verbose = T)$data %>%
  filter(redcap_event_name == "baseline_arm_1")
ivgtt_dat <- ivgtt_dat[rowSums(!is.na(ivgtt_dat[, grepl("^gluc_|^ins_", names(ivgtt_dat))])) > 0, ]

ivgtt_time <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/PANTHER/fsIVGTT/PANTHER_IVGTT_time_discrepency030624.csv")

ivgtt_time <- ivgtt_time %>% filter(time!=150)

ivgtt_missing <- ivgtt_dat %>%
  rowwise() %>%
  dplyr::mutate(glucose_missing = dplyr::case_when(sum(c_across(starts_with("gluc")), na.rm = T) == 0 ~ "Missing", T ~ "Complete"),
         insulin_missing = dplyr::case_when(sum(c_across(starts_with("ins")), na.rm = T) == 0 ~ "Missing", T ~ "Complete")) %>%
  dplyr::select(record_id, glucose_missing, insulin_missing)

# write.csv(ivgtt_missing, "/Volumes/Peds Endo/Petter Bjornstad/PANTHER/fsIVGTT/ivgtt_missing.csv", row.names = F)

ivgtt_dat <- ivgtt_dat %>% dplyr::select(-redcap_event_name, -ivgtt_labs_complete)
ivgtt_dat <- ivgtt_dat %>%
  filter(if_any(-record_id, Negate(is.na)))

time_column <- c(0, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 19, 22, 23, 24, 25, 27, 30, 35, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180) 

ivgtt_list <- list()
subject_list <- ivgtt_dat$record_id

for (i in 1:length(ivgtt_dat$record_id)) {
  subject_level_df <- ivgtt_dat[i,]
  record_id_temp <- gsub("-", ".", subject_level_df$record_id)
  subject_time <- ivgtt_time[record_id_temp]
  subject_time_col <- round(time_column + (subject_time/60), 2)
  names(subject_time_col) = "time"
  glucose_column = subject_level_df %>% dplyr::select(starts_with("gluc")) %>%
    as.character()
  insulin_column = subject_level_df %>% dplyr::select(starts_with("ins")) %>%
    as.character()
  result_dat <- data.frame(time = subject_time_col,
                           glucose = glucose_column,
                           insulin = insulin_column)
  df_name <- subject_level_df$record_id
  ivgtt_list[[df_name]] <- result_dat
  write.csv(result_dat, paste0("/Volumes/Peds Endo/Petter Bjornstad/PANTHER/fsIVGTT/subject_level_data/", df_name, ".csv"), quote = FALSE, row.names = FALSE)
}


