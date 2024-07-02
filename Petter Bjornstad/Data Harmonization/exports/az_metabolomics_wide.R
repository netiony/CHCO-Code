az_dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/AstraZeneca/az_targeted_normalized.csv")

az_dat_wide <-
  reshape(az_dat %>%
            dplyr::select(record_id, metabolite, peak_area, time), 
          direction = "wide",  
          idvar = c("record_id", "time"), 
          timevar = "metabolite") %>%
  rename_with(~ gsub("^peak_area\\.", "", .))

write.csv(az_dat_wide, "/Volumes/Peds Endo/Petter Bjornstad/AstraZeneca/az_targeted_normalized_wide.csv", row.names = F, 
          na = "")
