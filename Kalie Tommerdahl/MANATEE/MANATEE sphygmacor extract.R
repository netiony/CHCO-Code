sphygo <- read.csv("/Volumes/Peds Endo/Kalie Tommerdahl/MANATEE/Data clean/MANATEE sphygmacor.csv", na.strings = "")
# Participant arms
control <- c(6, 16, 20, 21, 22)
treatment <- c(2, 3, 5, 8, 9, 10, 11, 12, 14, 15, 17, 19, 23, 24, 25, 26, 27, 28, 29)

is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

sphygo <- sphygo %>%
  mutate(sc_pwa_datetime = paste0(sc_pwa_datetime, " ", time),
         sc_pwv_datetime_1 = paste0(sc_pwv_datetime_1, " ", time.1),
         sc_pwv_datetime_2 = paste0(sc_pwv_datetime_2, " ", time.2),
         sc_pwv_datetime_3 = paste0(sc_pwv_datetime_3, " ", time.3),
         redcap_event_name = 
           case_when(visit == 1 & record_id %in% treatment ~ "study_visit_baseli_arm_1",
                     visit == 2 & record_id %in% treatment ~ "study_visit_follow_arm_1",
                     visit == 1 & record_id %in% control ~ "study_visit_baseli_arm_2",
                     visit == 2 & record_id %in% control ~ "study_visit_follow_arm_2"),
         sc_height = case_when(is.wholenumber((sc_height - floor(sc_height))*10) ~ 
           (floor(sc_height) * 12) + ((sc_height - floor(sc_height))*10),
           T ~ (floor(sc_height) * 12) + ((sc_height - floor(sc_height))*100))) %>%
  dplyr::select(-starts_with("time"), -visit)

write.csv(sphygo, 
          "/Volumes/Peds Endo/Kalie Tommerdahl/MANATEE/Data clean/MANATEE-sphygo spreadsheet.csv", 
          row.names = F, na = "")
