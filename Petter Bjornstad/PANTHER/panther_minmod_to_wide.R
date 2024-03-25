minmod <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/PANTHER/fsIVGTT/PANTHER_IVGTT_MINMOD_following protocol 030524.csv")

minmod_wide <- reshape(minmod, timevar = "variable", idvar = "record_id", direction = "wide", sep = "_")
names(minmod_wide) <- sub("^value_", "", names(minmod_wide))
minmod_wide$redcap_event_name <- "baseline_arm_1"

write.csv(minmod_wide, "/Volumes/Peds Endo/Petter Bjornstad/PANTHER/fsIVGTT/PANTHER_minmod_wide.csv", row.names = F)
