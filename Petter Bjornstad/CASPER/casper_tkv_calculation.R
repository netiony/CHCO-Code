casper_tkv <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/CASPER/Data Cleaned/casper_tkv_calculation.csv")

casper_tkv <- casper_tkv %>% dplyr::mutate(skv_left = as.numeric(length_left/10) * as.numeric(width_left/10) * as.numeric(depth_left/10) * (pi/6),
              skv_right = as.numeric(length_right/10) * as.numeric(width_right/10) * as.numeric(depth_right/10) * (pi/6),
              tkv = coalesce(skv_left, 0) + coalesce(skv_right, 0))

write.csv(casper_tkv, "/Volumes/Peds Endo/Petter Bjornstad/CASPER/Data Cleaned/casper_tkv_calculation.csv", row.names = F)
