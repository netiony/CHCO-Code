library(dplyr)
library(stringr)
library(CGEN)

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
comorb$releaseid <- comorb$RELEASEID
comorb$RELEASEID <- NULL
normal_kidney_function <- comorb %>% filter(MIC0==0 & MAC0==0 & NEPHRO0==0 & RAPID0==0)

# AGEBASE - uncollapsed age at baseline
AGEBASE <- read.csv("./Clinical data/TODAY/AGEBASE.csv")
AGEBASE$releaseid <- AGEBASE$RELEASEID
AGEBASE$RELEASEID <- NULL
for_matching <- merge(normal_kidney_function, AGEBASE, by="releaseid",all.x = T,all.y = F)

# PAT
PAT <- read.csv("./Clinical data/TODAY/PAT.csv")
keepPAT <- PAT %>% select(releaseid,sex)
for_matching <- merge(for_matching,PAT,by="releaseid",all.x = T,all.y = F)

# PE
PE <- read.csv("./Clinical data/TODAY/PE.csv")
PE_baseline <- PE %>% filter(mvisit=="M00")
PE_baseline <- PE_baseline %>% select(releaseid,tanner)
for_matching <- merge(for_matching,PE_baseline,by="releaseid",all.x = T,all.y = F)

# get treatment arm - PRIMOUT
PRIMOUT <- read.csv("./Clinical data/TODAY/PRIMOUT.csv")
for_matching <- merge(for_matching,PRIMOUT,by="releaseid",all.x = T,all.y = F)
# IDs in PRIMOUT don't seem to match the release IDs

# read in soma data so we only pick people in the ancillary study
load(file = "./Somalogic data raw/soma.Rdata")
keep_soma <- soma %>% select(releaseid)
keep_soma <- unique(keep_soma)
for_matching <- merge(for_matching,keep_soma,by="releaseid",all.x = F, all.y = T)

# write file to read into SAS
for_matching <- for_matching %>% select(releaseid, AGEBASE, sex, tanner, tx, GLYC, DAYSTOGLYC)
write.csv(for_matching,"./ViCTER matching/for_matching.csv", row.names = F, na=".")

# create index variable
for_matching$index <- str_c(for_matching$AGEBASE,for_matching$sex,for_matching$tanner,for_matching$tx)
for_matching <- for_matching %>% filter(!is.na(GLYC))
for_matching <- for_matching %>% filter(!is.na(index))

# find 20 with shortest time to glycemic failure - CASES
f <- for_matching %>% filter(GLYC==1)
f <- f %>% arrange(DAYSTOGLYC)
cases <- f %>% slice_head(n=22)
cases <- cases %>% filter(!releaseid=="65-96152")
cases <- cases %>% filter(!releaseid=="65-92022")
cases$case1_control0 <- 1

# those who did not reach glycemic failure - CONTROLS
controls <- for_matching %>% filter(GLYC==0)
controls$case1_control0 <- 0

final <- rbind(cases,controls)
matched <- getMatchedSets(final, CC=TRUE, NN=FALSE, ccs.var = "case1_control0",dist.vars = c("AGEBASE","sex"), strata.var = "tx")                                                                                 
#matched <- getMatchedSets(final, CC=TRUE, NN=FALSE, ccs.var = "case1_control0",dist.vars = "index")                                                                                 

summary(matched)
final <- cbind(final, matched$CC)
final <- final %>% filter(!`matched$CC`==21)
table(final$`matched$CC`,final$tx)
table(final$`matched$CC`,final$AGEBASE)
table(final$`matched$CC`,final$sex)
table(final$case1_control0,final$tx)
table(final$case1_control0,final$AGEBASE)
write.csv(final,"./ViCTER matching/matched_pairs.csv", row.names = F, na=".")

# merge sample info with matched pairs
# need to select the baseline sample info
# pairs <- read.csv("./ViCTER matching/matched_pairs.csv")
# base <- soma %>% arrange(releaseid,Date.Drawn) %>% group_by(releaseid) %>% filter(row_number()==1)
# sampleinfo <- base[,1:37]
# ids <- base$releaseid
# sampleinfo <- cbind(sampleinfo,ids)
# colnames(sampleinfo) <- c(colnames(sampleinfo)[1:37],"releaseid")
# pairs <- merge(pairs,sampleinfo,by="releaseid",all.x = T, all.y=F)
# merge in diabetes duration
# dxtime <- read.csv("./Clinical data/TODAY/PAT.csv")
# dxtime <- dxtime %>% select(releaseid,dxtime)
# pairs <- merge(pairs,dxtime,by="releaseid",all.x=T,all.y=F)
# write.csv(pairs,"./ViCTER matching/matched_pairs.csv", row.names = F, na=".")
