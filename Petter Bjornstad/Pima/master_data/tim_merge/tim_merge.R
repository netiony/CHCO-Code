# Load Hmisc here, not the R files from REDCap
library(Hmisc)
library(tidyverse)
# Set github directory based on operating system
git_dir = ifelse(.Platform$OS.type == "unix",
                  "~/GitHub/CHCO-Code",
                  "path/to/Laura/or/Cameron's/GitHub/folder")
# Set working directory based on operating system
home_dir = ifelse(.Platform$OS.type == "unix",
                  "/run/user/1000/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Pima/Master data/Raw data/",
                  "B:/Peds Endo/Petter Bjornstad/Pima/Master data/Raw data/") # Laura and Cameron, you may need to fix this as well as the GitHub path above
setwd(home_dir)
# Go through in the order of Rob's email:
# 1. Vital status: this contains the death and dialysis dates for participants from all protocols.
source(paste0(git_dir,"/Petter Bjornstad/Pima/master_data/tim_merge/vital_status.R"))
# Get relevant information for each person. Per Rob:
# the vital status database tracks date of renal replacement therapy, type of 
# renal replacement therapy, whether the participant received a transplant, 
# date of death, and cause of death for each participant ever enrolled in any of 
# the kidney studies.  These can be collapsed into a single record for each person.
vital_status = vital_status %>% group_by(record_id) %>%
  summarise(across(esrd_start_date:death_notice_form_complete.factor,
                   ~ .x[max(which(!is.na(.x)))]), # Get the last non-missing row for each patient and each variable
            .groups = "drop")
# 2. Group 4 – this is the kidney study protocol conducted in a small group of participants with advanced kidney disease. 
# It includes clearance studies as well as routine home visits.
source(paste0(git_dir,"/Petter Bjornstad/Pima/master_data/tim_merge/group4.R"))
# Combine visit date columns - there are a lot of repeats within a row
dates = colnames(group4)[grep("f\\d{,2}visitdate$",colnames(group4))]
dates = c(dates,"f12expdatevis")
group4$visit_date = apply(group4[,dates],1,function(r){
  d = unique(r[!is.na(r)])
  ifelse(length(d)>0,d,NA)
  })
# Combine vital status and group 4
final_merge = full_join(vital_status,group4)
rm(vital_status,group4)
# 3. Ficoll – this is the kidney study protocol involving participants with early kidney disease. 
# Clearance studies were performed at intervals of every 6 months and there were no home visits.
source(paste0(git_dir,"/Petter Bjornstad/Pima/master_data/tim_merge/ficoll.R"))
# Rename visit date column
colnames(ficoll)[which(colnames(ficoll)=="clrncvisitdate")] = "visit_date"
# Add new data to merged
final_merge = full_join(final_merge,ficoll)
rm(ficoll)
# 4. Losartan – this is the clinical trial that included both annual clearance visits and routine quarterly home visits.  
# It is by far our longest running study with the most data.  It is called “PECRBRenoprote” in the R and .csv files.
source(paste0(git_dir,"/Petter Bjornstad/Pima/master_data/tim_merge/losartan.R"))
# Dates
dates = colnames(losartan)[grep("f\\d{,2}visitdate$",colnames(losartan))]
dates = c(dates,"f12expdatevis")
losartan$visit_date = apply(losartan[,dates],1,function(r){
  d = unique(r[!is.na(r)])
  ifelse(length(d)>0,d,NA)
})
# Fix f10spturalbu - warnings are for Questions for "****.*" treat as missing for now
final_merge$f10spturalbu = as.numeric(final_merge$f10spturalbu)
# Add new data to merged
final_merge = full_join(final_merge,losartan)
rm(losartan)
# 5. DDN – this is our most recent study—referred to as “determinants of diabetic nephropathy”. 
# These were annual visits and only clearance visits were done.
source(paste0(git_dir,"/Petter Bjornstad/Pima/master_data/tim_merge/ddn.R"))
# Visit dates
dates = colnames(ddn)[grep("f\\d{,2}visitdate$",colnames(ddn))]
dates = c(dates,"f12expdatevis")
ddn$visit_date = apply(ddn[,dates],1,function(r){
  d = unique(r[!is.na(r)])
  ifelse(length(d)>0,d,NA)
})
# Add new data to merged
final_merge = full_join(final_merge,ddn)
rm(ddn)
# Check for duplicate columns - none!
# grep("/.x",colnames(final_merge))
# grep("/.y",colnames(final_merge))
# Format columns
final_merge$visit_date = parse_date(final_merge$visit_date)
# Remove columns with all missing
final_merge[,colSums(is.na(final_merge))==nrow(final_merge)] = NULL
# Sort by id then visit name
final_merge = final_merge %>% select(record_id,redcap_event_name,visit_date,everything())
# Calculate new variables grouped by person
final_merge = final_merge %>% group_by(record_id) %>%
  mutate(days_from_earliest = difftime(visit_date,min(visit_date,na.rm = T),units = "days"))
# Clean up
rm(home_dir,git_dir)
