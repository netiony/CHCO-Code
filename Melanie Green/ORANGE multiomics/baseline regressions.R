# read in data
data <- read.csv("/Volumes/Peds Endo/Melanie Green/ORANGE/Data clean/merged multiomics data.csv")
data$X <- NULL

# create file for storing results
dataOut <- "/Volumes/Peds Endo/Melanie Green/ORANGE/Data clean/baseline_regressions.csv"
write(paste("analyte1","analyte2",""))