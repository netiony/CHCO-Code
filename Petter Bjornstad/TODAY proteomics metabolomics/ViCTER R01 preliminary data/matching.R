library(dplyr)

if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad/TODAY subaward"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = "/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward"
}

setwd(home_dir)

# find those with normal kidney function at baseline
comorb <- read.csv("./Clinical data/COMORB.csv")
normal_kidney_function <- comorb %>% filter(MIC0==0 & MAC0==0 & NEPHRO0==0 & RAPID0==0)
comorb$releaseid <- comorb$RELEASEID
comorb$RELEASEID <- NULL

# AGEBASE - uncollapsed age at baseline
AGEBASE <- read.csv("./Clinical data/TODAY/AGEBASE.csv")
AGEBASE$releaseid <- AGEBASE$RELEASEID
AGEBASE$RELEASEID <- NULL
for_matching <- merge(comorb, AGEBASE, by="releaseid",all.x = T,all.y = F)

# PAT
PAT <- read.csv("./Clinical data/TODAY/PAT.csv")
keepPAT <- PAT %>% select(releaseid,sex)
for_matching <- merge(for_matching,PAT,by="releaseid",all.x = T,all.y = F)

# get treatment arm - PRIMOUT
PRIMOUT <- read.csv("./Clinical data/TODAY/PRIMOUT.csv")
PRIMOUT$releaseid <- PRIMOUT$PTID
#for_matching <- merge(for_matching,PRIMOUT,by="releaseid",all.x = T,all.y = F)
# IDs in PRIMOUT don't seem to match the release IDs

# write file to read into SAS
write.csv(for_matching,"E:/Petter Bjornstad/TODAY subaward/ViCTER matching/for_matching.csv", row.names = F)


