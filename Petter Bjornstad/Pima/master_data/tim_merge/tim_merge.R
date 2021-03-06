# Load Hmisc here, not the R files from REDCap
library(Hmisc)
library(dplyr)
# Set github directory based on operating system
git_dir = ifelse(.Platform$OS.type == "unix",
                  "/Users/timvigers/GitHub/CHCO-Code",
                  "path/to/Laura/or/Cameron's/GitHub/folder")
# Set working directory based on operating system
home_dir = ifelse(.Platform$OS.type == "unix",
                  "/Volumes/som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Pima/Master data/Raw data/",
                  "B:/Peds Endo/Petter Bjornstad/Pima/Master data/Raw data/")
setwd(home_dir)
# Go through in the order of Rob's email:
# 1. Vital status: this contains the death and dialysis dates for participants from all protocols.
source(paste0(git_dir,"/Petter Bjornstad/Pima/master_data/tim_merge/vital_status.R"))
vital_status = data
# 2. Group 4 – this is the kidney study protocol conducted in a small group of participants with advanced kidney disease. 
# It includes clearance studies as well as routine home visits.
source(paste0(git_dir,"/Petter Bjornstad/Pima/master_data/tim_merge/group4.R"))
# Combine vital status and group 4
final_merge = full_join(vital_status,data)
# 3. Ficoll – this is the kidney study protocol involving participants with early kidney disease. 
# Clearance studies were performed at intervals of every 6 months and there were no home visits.
source(paste0(git_dir,"/Petter Bjornstad/Pima/master_data/tim_merge/ficoll.R"))
# Add new data to merged
final_merge = full_join(final_merge,data)
# 4. Losartan – this is the clinical trial that included both annual clearance visits and routine quarterly home visits.  
# It is by far our longest running study with the most data.  It is called “PECRBRenoprote” in the R and .csv files.
source(paste0(git_dir,"/Petter Bjornstad/Pima/master_data/tim_merge/losartan.R"))
# Fix f10spturalbu - warnings are for Questions for "****.*" treat as missing for now
final_merge$f10spturalbu = as.numeric(final_merge$f10spturalbu)
# Add new data to merged
final_merge = full_join(final_merge,data)
# 5. DDN – this is our most recent study—referred to as “determinants of diabetic nephropathy”. 
# These were annual visits and only clearance visits were done.
source(paste0(git_dir,"/Petter Bjornstad/Pima/master_data/tim_merge/ddn.R"))
# Add new data to merged
final_merge = full_join(final_merge,data)
# Check for duplicate columns - none!
# grep("/.x",colnames(final_merge))
# grep("/.y",colnames(final_merge))
