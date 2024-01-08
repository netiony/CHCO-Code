library(dplyr)
ivgtt_dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/PANTHER/fsIVGTT/PANTHER_DATA_2023-10-23_1044_fsIVGTT_labs_filled.csv")
ivgtt_time <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/PANTHER/fsIVGTT/PANTHER_IVGTT_time_discrepency.csv")

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
