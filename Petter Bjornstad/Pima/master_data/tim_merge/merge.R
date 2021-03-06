# Load Hmisc here, not the R files from REDCap
library(Hmisc)
# Set working directory based on operating system
home_dir = ifelse(.Platform$OS.type == "unix",
                  "/Volumes/som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Pima/Master data",
                  "B:/Peds Endo/Petter Bjornstad/Pima/Master data")
setwd(home_dir)
