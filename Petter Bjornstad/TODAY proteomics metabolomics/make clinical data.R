library(dplyr)

if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad/TODAY subaward"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = "/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward"
}

setwd(home_dir)

# COMORB dataset
comorb <- read.csv("./Clinical data/COMORB.csv")
comorb$MIC.OR.MAC <- ifelse(comorb$MAC==1 | comorb$MIC==1,1,
                            ifelse(is.na(comorb$MAC) & is.na(comorb$MIC), NA, 0))
comorb$releaseid <- comorb$RELEASEID
comorb$RELEASEID <- NULL
# Save
save(comorb,file = "./Clinical data/comorb.Rdata")

# CBL
CBL <- read.csv("./Clinical data/TODAY/CBL.csv")

# ADDCBL
ADDCBL <- read.csv("./Clinical data/TODAY/ADDCBL.csv")

# BASELINE
BASELINE <- read.csv("./Clinical data/TODAY/BASELINE.csv")

# PAT
PAT <- read.csv("./Clinical data/TODAY/PAT.csv")
keepPAT <- PAT %>% select(releaseid,age,sex)

# AGEBASE - uncollapsed age at baseline
AGEBASE <- read.csv("./Clinical data/TODAY/AGEBASE.csv")
AGEBASE$releaseid <- AGEBASE$RELEASEID
AGEBASE$RELEASEID <- NULL

# create new dataset of baseline risk factors
basecbl <- CBL %>% filter(mvisit=="M00")
baseaddcbl <- ADDCBL %>% filter(mvisit=="M00")
baserisk <- merge(BASELINE, basecbl, by="releaseid", all.x=T, ally=T)
baserisk <- merge(baserisk, baseaddcbl, by="releaseid", all.x=T, ally=T)
baserisk$si_1_ins0 <- 1/baserisk$ins0min
baserisk$log_trig <- log(baserisk$Trig)
baserisk <- baserisk %>% select(releaseid, HbA1c, log_trig, sbp, uacid, si_1_ins0, UAlbCreat, bmi)
baserisk <- merge(baserisk,keepPAT,by="releaseid",all.x = T,all.y = F)
baserisk$age <- NULL
baserisk <- merge(baserisk,AGEBASE,by="releaseid",all.x = T,all.y = F)

# Save
save(baserisk,file = "./Clinical data/TODAY/baserisk.Rdata")

###############################################################
# create dataset of risk factors at 10 year visit from TODAY2 #
###############################################################

# CBL
CBL_TODAY2 <- read.csv("./Clinical data/TODAY2/CBL.csv")
CBL_TODAY2_KEEP <- CBL_TODAY2 %>% filter(pvisit=="P120") %>% select(releaseid, hba1c, trig, ins)

# ADDCBL
#ADDCBL_TODAY2 <- read.csv("./Clinical data/TODAY2/ADDCBL.csv")

# VISIT
VISIT_TODAY2 <- read.csv("./Clinical data/TODAY2/VISIT.csv")
VISIT_TODAY2$releaseid <- VISIT_TODAY2$RELEASEID
VISIT_TODAY2$sbp <- VISIT_TODAY2$SBP
VISIT_TODAY2_KEEP <- VISIT_TODAY2 %>% filter(PVISIT=="P120") %>% select(releaseid, sbp)

# insulin was measured at 9 year visit not 10 year
