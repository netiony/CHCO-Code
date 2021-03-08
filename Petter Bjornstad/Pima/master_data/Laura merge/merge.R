library(tidyverse)

code_dir = ifelse(.Platform$OS.type == "unix",
                  "",
                  "C:/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/Pima/master_data/Laura merge")
setwd(code_dir)

# Questions for Tim
# 1. what did he do with all the factor variables created?
# 2. Do we need to convert dates to date variables?
# 3. Why multiple records per pt in vital status - should we collapse these to 1 record?  double check that ESRD, death, and transplant data are not
#    on the same row?
# 4. What should the structure of the final dataset be?  1 row per visit, if visits in different protocols have the same date, combine.  
#    Then vital status fields get appended to every visit?
# 5. Records with missing visit dates in group4 and ficoll?
# 6. There are variables that show up in multiple dataframes (e.g., sacaton number), and therefore get .x added, etc.

# read in datasets using provided R code in order of Rob's email
# use copies of code moved to github directory so they can be edited

# Vital status
source("NelsonVitalStatusDEC_R_2021-02-04_1616.r")
vital_status <- data
rm(data)
# check that death and esrd are not in the same row
vital_status$flag <- ifelse(!is.na(vital_status$esrd_start_date) & !is.na(vital_status$dod),1,0)
# break VS into "forms" and then recombine
esrd <- vital_status[c("record_id","redcap_event_name","esrd_start_date","esrd_start_date_est","esrd_type","esrd_is_transplanted","esrd_diagnosis",
                       "esrd_complete_date","esrd_initials","end_stage_renal_disease_form_complete","esrd_start_date_est.factor",
                       "esrd_type.factor","esrd_is_transplanted.factor","esrd_diagnosis.factor")]
esrd <- esrd[esrd$redcap_event_name=="esrd_arm_1",]
transplant <- vital_status[,c("record_id","redcap_event_name","transplant_date","transplant_complete_date","transplant_initials","transplant_form_complete")]
transplant <- transplant[transplant$redcap_event_name=="transplant_arm_1",]
death <- vital_status[,c("record_id","redcap_event_name","dod","codunerlying","coddrunerlying","physiciancr","initials","death_comments",
                         "death_notice_form_complete","physiciancr.factor")]
death <- death[death$redcap_event_name=="death_arm_1",]
comment <- vital_status[,c("record_id","comment","comment_form_complete")]
comment <- comment[!is.na(comment$comment),]
vital_status_onerow <- merge(esrd,transplant,by="record_id",all.x = T, all.y = T)
vital_status_onerow <- merge(vital_status_onerow,death,by="record_id",all.x=T,all.y=T)
vital_status_onerow <- merge(vital_status_onerow,comment,by="record_id",all.x=T,all.y=T)
vital_status_onerow$redcap_event_name <- NULL
vital_status_onerow$redcap_event_name.x <- NULL
vital_status_onerow$redcap_event_name.y <- NULL

# Group 4
source("Group4UofMRemodel112_R_2021-02-04_1617.r")
group4 <- data
rm(data)
#write.table(names(group4),"variable names group 4.txt")
# there are lots of visit dates per row
dates = colnames(group4)[grep("f\\d{,2}visitdate$",colnames(group4))]
# are all of the visit dates the same?
# they are all the same
group4$visit_date <- NA
for (i in 1:nrow(group4)) {
  temp <- group4[i,dates]
  temp <- temp[!is.na(temp)]
  #group4$datecount <- length(unique(temp))
  group4[i,"visit_date"] <- ifelse(length(unique(temp))>0,unique(temp),NA)
}
dates <- NULL

# Ficoll
source("FicollUniversityOfMi_R_2021-02-04_1618.r")
ficoll <- data
rm(data)
dates = colnames(ficoll)[grep("f\\d{,2}visitdate$",colnames(ficoll))]
ficoll$visit_date <- NA
for (i in 1:nrow(ficoll)) {
  temp <- ficoll[i,dates]
  temp <- temp[!is.na(temp)]
  #ficoll$datecount <- length(unique(temp))
  ficoll[i,"visit_date"] <- ifelse(length(unique(temp))>0,unique(temp),NA)
}
dates <- NULL

# merge Group 4 and Ficoll
final_merge <- merge(group4,ficoll,by=c("record_id","visit_date"),all.x=T,all.y=T)

# Losartan
source("NelsonPECRBRenoprote_R_2021-02-04_1619.r")
losartan <- data
rm(data)

# DDN
source("Nelson13DKN151Determ_R_2021-02-04_1610.r")
ddn <- data
rm(data)

# finally, merge in vital status data