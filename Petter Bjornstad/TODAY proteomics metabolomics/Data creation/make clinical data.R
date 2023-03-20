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
# Merge in bariatric surgery info
# read in TME dataset to find those who had bariatric surgery
tme <- read.csv("./Clinical data/TODAY2/TME.csv")
tme <- tme %>% filter(TMETYPE==1)
tme <- tme %>% select(RELEASEID,TMETYPE,DAYSTOTME)
colnames(tme) <- c("releaseid","TMETYPE","DAYSTOTME")
comorb <- merge(comorb, tme, by="releaseid", all.x = T, all.y=F)
# now censor all outcomes at time of TME for those who had bariatric surgery
# need to reset the indicator variable as well as the time variable
# HTN
comorb <- comorb %>% mutate(DAYSTOHTN = case_when(
  HTN==0 ~ DAYSTOHTN,
  HTN==1 & is.na(DAYSTOTME) ~ DAYSTOHTN,
  HTN==1 & !is.na(DAYSTOTME) & DAYSTOHTN>DAYSTOTME ~ DAYSTOTME,
  HTN==1 & !is.na(DAYSTOTME) & DAYSTOHTN<=DAYSTOTME ~ DAYSTOHTN
))  
comorb <- comorb %>% mutate(HTN = case_when(
  HTN==0 ~ HTN,
  HTN==1 & is.na(DAYSTOTME) ~ HTN,
  HTN==1 & !is.na(DAYSTOTME) & DAYSTOHTN>DAYSTOTME ~ as.integer(0),
  HTN==1 & !is.na(DAYSTOTME) & DAYSTOHTN<=DAYSTOTME ~ HTN
))  
# LDL
comorb <- comorb %>% mutate(DAYSTOLDL = case_when(
  LDLDLP==0 ~ DAYSTOLDL,
  LDLDLP==1 & is.na(DAYSTOTME) ~ DAYSTOLDL,
  LDLDLP==1 & !is.na(DAYSTOTME) & DAYSTOLDL>DAYSTOTME ~ DAYSTOTME,
  LDLDLP==1 & !is.na(DAYSTOTME) & DAYSTOLDL<=DAYSTOTME ~ DAYSTOLDL
))  
comorb <- comorb %>% mutate(LDLDLP = case_when(
  LDLDLP==0 ~ LDLDLP,
  LDLDLP==1 & is.na(DAYSTOTME) ~ LDLDLP,
  LDLDLP==1 & !is.na(DAYSTOTME) & DAYSTOLDL>DAYSTOTME ~ as.integer(0),
  LDLDLP==1 & !is.na(DAYSTOTME) & DAYSTOLDL<=DAYSTOTME ~ LDLDLP
))  
# TG
comorb <- comorb %>% mutate(DAYSTOTG = case_when(
  TGDLP==0 ~ DAYSTOTG,
  TGDLP==1 & is.na(DAYSTOTME) ~ DAYSTOTG,
  TGDLP==1 & !is.na(DAYSTOTME) & DAYSTOTG>DAYSTOTME ~ DAYSTOTME,
  TGDLP==1 & !is.na(DAYSTOTME) & DAYSTOTG<=DAYSTOTME ~ DAYSTOTG
))  
comorb <- comorb %>% mutate(TGDLP = case_when(
  TGDLP==0 ~ TGDLP,
  TGDLP==1 & is.na(DAYSTOTME) ~ TGDLP,
  TGDLP==1 & !is.na(DAYSTOTME) & DAYSTOTG>DAYSTOTME ~ as.integer(0),
  TGDLP==1 & !is.na(DAYSTOTME) & DAYSTOTG<=DAYSTOTME ~ TGDLP
))  
# ANYDLP
comorb <- comorb %>% mutate(DAYSTOANYDLP = case_when(
  ANYDLP==0 ~ DAYSTOANYDLP,
  ANYDLP==1 & is.na(DAYSTOTME) ~ DAYSTOANYDLP,
  ANYDLP==1 & !is.na(DAYSTOTME) & DAYSTOANYDLP>DAYSTOTME ~ DAYSTOTME,
  ANYDLP==1 & !is.na(DAYSTOTME) & DAYSTOANYDLP<=DAYSTOTME ~ DAYSTOANYDLP
))  
comorb <- comorb %>% mutate(ANYDLP = case_when(
  ANYDLP==0 ~ ANYDLP,
  ANYDLP==1 & is.na(DAYSTOTME) ~ ANYDLP,
  ANYDLP==1 & !is.na(DAYSTOTME) & DAYSTOANYDLP>DAYSTOTME ~ as.integer(0),
  ANYDLP==1 & !is.na(DAYSTOTME) & DAYSTOANYDLP<=DAYSTOTME ~ ANYDLP
))  
# MIC
comorb <- comorb %>% mutate(DAYSTOMIC = case_when(
  MIC==0 ~ DAYSTOMIC,
  MIC==1 & is.na(DAYSTOTME) ~ DAYSTOMIC,
  MIC==1 & !is.na(DAYSTOTME) & DAYSTOMIC>DAYSTOTME ~ DAYSTOTME,
  MIC==1 & !is.na(DAYSTOTME) & DAYSTOMIC<=DAYSTOTME ~ DAYSTOMIC
))  
comorb <- comorb %>% mutate(MIC = case_when(
  MIC==0 ~ MIC,
  MIC==1 & is.na(DAYSTOTME) ~ MIC,
  MIC==1 & !is.na(DAYSTOTME) & DAYSTOMIC>DAYSTOTME ~ as.integer(0),
  MIC==1 & !is.na(DAYSTOTME) & DAYSTOMIC<=DAYSTOTME ~ MIC
))  
# MAC
comorb <- comorb %>% mutate(DAYSTOMAC = case_when(
  MAC==0 ~ DAYSTOMAC,
  MAC==1 & is.na(DAYSTOTME) ~ DAYSTOMAC,
  MAC==1 & !is.na(DAYSTOTME) & DAYSTOMAC>DAYSTOTME ~ DAYSTOTME,
  MAC==1 & !is.na(DAYSTOTME) & DAYSTOMAC<=DAYSTOTME ~ DAYSTOMAC
))  
comorb <- comorb %>% mutate(MAC = case_when(
  MAC==0 ~ MAC,
  MAC==1 & is.na(DAYSTOTME) ~ MAC,
  MAC==1 & !is.na(DAYSTOTME) & DAYSTOMAC>DAYSTOTME ~ as.integer(0),
  MAC==1 & !is.na(DAYSTOTME) & DAYSTOMAC<=DAYSTOTME ~ MAC
))  
# NEPHRO
comorb <- comorb %>% mutate(DAYSTONEPHRO = case_when(
  NEPHRO==0 ~ DAYSTONEPHRO,
  NEPHRO==1 & is.na(DAYSTOTME) ~ DAYSTONEPHRO,
  NEPHRO==1 & !is.na(DAYSTOTME) & DAYSTONEPHRO>DAYSTOTME ~ DAYSTOTME,
  NEPHRO==1 & !is.na(DAYSTOTME) & DAYSTONEPHRO<=DAYSTOTME ~ DAYSTONEPHRO
))  
comorb <- comorb %>% mutate(NEPHRO = case_when(
  NEPHRO==0 ~ NEPHRO,
  NEPHRO==1 & is.na(DAYSTOTME) ~ NEPHRO,
  NEPHRO==1 & !is.na(DAYSTOTME) & DAYSTONEPHRO>DAYSTOTME ~ as.integer(0),
  NEPHRO==1 & !is.na(DAYSTOTME) & DAYSTONEPHRO<=DAYSTOTME ~ NEPHRO
))  
# HYP
comorb <- comorb %>% mutate(DAYSTOHYP = case_when(
  HYP==0 ~ DAYSTOHYP,
  HYP==1 & is.na(DAYSTOTME) ~ DAYSTOHYP,
  HYP==1 & !is.na(DAYSTOTME) & DAYSTOHYP>DAYSTOTME ~ DAYSTOTME,
  HYP==1 & !is.na(DAYSTOTME) & DAYSTOHYP<=DAYSTOTME ~ DAYSTOHYP
))  
comorb <- comorb %>% mutate(HYP = case_when(
  HYP==0 ~ HYP,
  HYP==1 & is.na(DAYSTOTME) ~ HYP,
  HYP==1 & !is.na(DAYSTOTME) & DAYSTOHYP>DAYSTOTME ~ as.integer(0),
  HYP==1 & !is.na(DAYSTOTME) & DAYSTOHYP<=DAYSTOTME ~ HYP
))  
# FILAM
comorb <- comorb %>% mutate(DAYSTOFILAM = case_when(
  FILAM==0 ~ DAYSTOFILAM,
  FILAM==1 & is.na(DAYSTOTME) ~ DAYSTOFILAM,
  FILAM==1 & !is.na(DAYSTOTME) & DAYSTOFILAM>DAYSTOTME ~ DAYSTOTME,
  FILAM==1 & !is.na(DAYSTOTME) & DAYSTOFILAM<=DAYSTOTME ~ DAYSTOFILAM
))  
comorb <- comorb %>% mutate(FILAM = case_when(
  FILAM==0 ~ FILAM,
  FILAM==1 & is.na(DAYSTOTME) ~ FILAM,
  FILAM==1 & !is.na(DAYSTOTME) & DAYSTOFILAM>DAYSTOTME ~ as.integer(0),
  FILAM==1 & !is.na(DAYSTOTME) & DAYSTOFILAM<=DAYSTOTME ~ FILAM
))  
# NEURO
comorb <- comorb %>% mutate(DAYSTONEURO = case_when(
  NEURO==0 ~ DAYSTONEURO,
  NEURO==1 & is.na(DAYSTOTME) ~ DAYSTONEURO,
  NEURO==1 & !is.na(DAYSTOTME) & DAYSTONEURO>DAYSTOTME ~ DAYSTOTME,
  NEURO==1 & !is.na(DAYSTOTME) & DAYSTONEURO<=DAYSTOTME ~ DAYSTONEURO
))  
comorb <- comorb %>% mutate(NEURO = case_when(
  NEURO==0 ~ NEURO,
  NEURO==1 & is.na(DAYSTOTME) ~ NEURO,
  NEURO==1 & !is.na(DAYSTOTME) & DAYSTONEURO>DAYSTOTME ~ as.integer(0),
  NEURO==1 & !is.na(DAYSTOTME) & DAYSTONEURO<=DAYSTOTME ~ NEURO
))  
# RETINO
comorb <- comorb %>% mutate(DAYSTORETINO = case_when(
  RETINO==0 ~ DAYSTORETINO,
  RETINO==1 & is.na(DAYSTOTME) ~ DAYSTORETINO,
  RETINO==1 & !is.na(DAYSTOTME) & DAYSTORETINO>DAYSTOTME ~ DAYSTOTME,
  RETINO==1 & !is.na(DAYSTOTME) & DAYSTORETINO<=DAYSTOTME ~ DAYSTORETINO
))  
comorb <- comorb %>% mutate(RETINO = case_when(
  RETINO==0 ~ RETINO,
  RETINO==1 & is.na(DAYSTOTME) ~ RETINO,
  RETINO==1 & !is.na(DAYSTOTME) & DAYSTORETINO>DAYSTOTME ~ as.integer(0),
  RETINO==1 & !is.na(DAYSTOTME) & DAYSTORETINO<=DAYSTOTME ~ RETINO
))  
# MVD
comorb <- comorb %>% mutate(DAYSTOMVD = case_when(
  MVD==0 ~ DAYSTOMVD,
  MVD==1 & is.na(DAYSTOTME) ~ DAYSTOMVD,
  MVD==1 & !is.na(DAYSTOTME) & DAYSTOMVD>DAYSTOTME ~ DAYSTOTME,
  MVD==1 & !is.na(DAYSTOTME) & DAYSTOMVD<=DAYSTOTME ~ DAYSTOMVD
))  
comorb <- comorb %>% mutate(MVD = case_when(
  MVD==0 ~ MVD,
  MVD==1 & is.na(DAYSTOTME) ~ MVD,
  MVD==1 & !is.na(DAYSTOTME) & DAYSTOMVD>DAYSTOTME ~ as.integer(0),
  MVD==1 & !is.na(DAYSTOTME) & DAYSTOMVD<=DAYSTOTME ~ MVD
))  
# GLYC
comorb <- comorb %>% mutate(DAYSTOGLYC = case_when(
  GLYC==0 ~ DAYSTOGLYC,
  GLYC==1 & is.na(DAYSTOTME) ~ DAYSTOGLYC,
  GLYC==1 & !is.na(DAYSTOTME) & DAYSTOGLYC>DAYSTOTME ~ DAYSTOTME,
  GLYC==1 & !is.na(DAYSTOTME) & DAYSTOGLYC<=DAYSTOTME ~ DAYSTOGLYC
))  
comorb <- comorb %>% mutate(GLYC = case_when(
  GLYC==0 ~ GLYC,
  GLYC==1 & is.na(DAYSTOTME) ~ GLYC,
  GLYC==1 & !is.na(DAYSTOTME) & DAYSTOGLYC>DAYSTOTME ~ as.integer(0),
  GLYC==1 & !is.na(DAYSTOTME) & DAYSTOGLYC<=DAYSTOTME ~ GLYC
))  
# drop bariatric surgery variables
comorb <- comorb %>% select(-c(TMETYPE,DAYSTOTME))
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
PAT$racedesc[PAT$race==1]<- "Non-Hispanic Black"
PAT$racedesc[PAT$race==2]<- "Hispanic"
PAT$racedesc[PAT$race==3]<- "Non-Hispanic White"
PAT$racedesc[PAT$race==4]<- "Other"
PAT$houseincdesc[PAT$houseinc==1]<- "<$24,999"
PAT$houseincdesc[PAT$houseinc==2]<- "$25,000-$49,999"
PAT$houseincdesc[PAT$houseinc==3]<- ">$50,000"
keepPAT <- PAT %>% select(releaseid,age,sex,dxtime,race,houseinc,racedesc,houseincdesc)

# AGEBASE - uncollapsed age at baseline
AGEBASE <- read.csv("./Clinical data/TODAY/AGEBASE.csv")
AGEBASE$releaseid <- AGEBASE$RELEASEID
AGEBASE$RELEASEID <- NULL

# PRIMOUT - treatment group assignment
PRIMOUT <- read.csv("./Clinical data/TODAY/PRIMOUT.csv")
PRIMOUT$txdesc[PRIMOUT$tx==1]<- "Metformin only"
PRIMOUT$txdesc[PRIMOUT$tx==2]<- "Metformin + rosiglitazone"
PRIMOUT$txdesc[PRIMOUT$tx==3]<- "Metformin + lifestyle"
keepPRIMOUT <- PRIMOUT %>% select(releaseid,tx,txdesc)

# create new dataset of baseline risk factors
basecbl <- CBL %>% filter(mvisit=="M00")
baseaddcbl <- ADDCBL %>% filter(mvisit=="M00")
baserisk <- merge(BASELINE, basecbl, by="releaseid", all.x=T, ally=T)
baserisk <- merge(baserisk, baseaddcbl, by="releaseid", all.x=T, ally=T)
baserisk$si_1_ins0 <- 1/baserisk$ins0min
baserisk$log_trig <- log(baserisk$Trig)
baserisk <- baserisk %>% select(releaseid, HbA1c, log_trig, sbp, uacid, si_1_ins0, UAlbCreat, bmi, HDL, codi)
baserisk <- merge(baserisk,keepPAT,by="releaseid",all.x = T,all.y = F)
baserisk$age <- NULL
baserisk <- merge(baserisk,AGEBASE,by="releaseid",all.x = T,all.y = F)
baserisk <- merge(baserisk, keepPRIMOUT,by="releaseid",all.x = T,all.y = F)

# Save
save(baserisk,file = "./Clinical data/TODAY/baserisk.Rdata")

###############################################################
# create dataset of risk factors at 10 year visit from TODAY2 #
###############################################################

# CBL
CBL_TODAY2 <- read.csv("./Clinical data/TODAY2/CBL.csv")
CBL_TODAY2_KEEP <- CBL_TODAY2 %>% filter(pvisit=="P120") %>% select(releaseid, hba1c, trig)

# INS 
# insulin was measured at 9 year visit not 10 year
INS_TODAY2 <- read.csv("./Clinical data/TODAY2/CBL.csv")
INS_TODAY2_KEEP <- INS_TODAY2 %>% filter(pvisit=="P108") %>% select(releaseid, ins)

# ADDCBL
#ADDCBL_TODAY2 <- read.csv("./Clinical data/TODAY2/ADDCBL.csv")

# VISIT
VISIT_TODAY2 <- read.csv("./Clinical data/TODAY2/VISIT.csv")
VISIT_TODAY2$releaseid <- VISIT_TODAY2$RELEASEID
VISIT_TODAY2$sbp <- VISIT_TODAY2$SBP
VISIT_TODAY2$bmi <- VISIT_TODAY2$BMI
VISIT_TODAY2_KEEP <- VISIT_TODAY2 %>% filter(PVISIT=="P120") %>% select(releaseid, sbp, bmi)

yr10risk <- merge(CBL_TODAY2_KEEP, INS_TODAY2_KEEP, by="releaseid", all.x = T, all.y = T)
yr10risk <- merge(yr10risk, VISIT_TODAY2_KEEP, by="releaseid", all.x = T, all.y = T)
yr10risk$si_1_ins0 <- 1/yr10risk$ins
yr10risk$log_trig <- log(yr10risk$trig)

# Save
save(yr10risk,file = "./Clinical data/TODAY/yr10risk.Rdata")

###############################################################
# create longitudinal datasets from TODAY and TODAY2          #
###############################################################

# TODAY BASELINE - baseline BMI 
BASELINE <- read.csv("./Clinical data/TODAY/BASELINE.csv")
BASELINE$visit <- "M00"
BASELINE_keep <- BASELINE %>% select(releaseid,visit,bmi)

# TODAY VISIT - BMI 
VISIT <- read.csv("./Clinical data/TODAY/VISIT.csv")
VISIT$visit <- VISIT$mvisit
VISIT_keep <- VISIT %>% select(releaseid,visit,bmi)

# TODAY CBL - eIS
CBL <- read.csv("./Clinical data/TODAY/CBL.csv")
CBL$visit <- CBL$mvisit
CBL_keep <- CBL %>% select(releaseid,visit,ins0min)

# TODAY ADDCBL - coDI
ADDCBL <- read.csv("./Clinical data/TODAY/ADDCBL.csv")
ADDCBL$visit <- ADDCBL$mvisit
ADDCBL_keep <- ADDCBL %>% select(releaseid,visit,codi)

# TODAY2 VISIT - BMI and waist circ
VISIT_TODAY2 <- read.csv("./Clinical data/TODAY2/VISIT.csv")
VISIT_TODAY2$releaseid <- VISIT_TODAY2$RELEASEID
VISIT_TODAY2$bmi <- VISIT_TODAY2$BMI
VISIT_TODAY2$visit <- VISIT_TODAY2$PVISIT
VISIT_TODAY2_KEEP <- VISIT_TODAY2 %>% select(releaseid, visit, bmi)

# TODAY2 CBL - eIS and coDI
CBL_TODAY2 <- read.csv("./Clinical data/TODAY2/CBL.csv")
CBL_TODAY2$visit <- CBL_TODAY2$pvisit
CBL_TODAY2$ins0min <- CBL_TODAY2$ins
CBL_TODAY2_KEEP <- CBL_TODAY2 %>% select(releaseid, visit, ins0min, codi)

# merge all data

# calculated variables - eIS

# write file
