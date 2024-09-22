library(SomaDataIO)
library(dplyr)

if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad/TODAY subaward"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = ""
}
#setwd(home_dir)

# read in data - we want to use the fully processed, normalized file ending in "anmlSMP.adat"
soma <- read_adat("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Somalogic data raw/WUS-22-001_Somalogic_normalized/WUS-22-001_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")
analytes <- getAnalyteInfo(soma)

# remove fc mouse and no protein
drop <- analytes %>% filter(Target == "Fc_MOUSE" | Target == "No Protein" | !(Organism == "Human") | !(Type == "Protein"))
apt_drop <- drop$AptName
soma <- soma %>% select(!all_of(apt_drop))
analytes <- analytes %>% filter(!Target == "Fc_MOUSE")
analytes <- analytes %>% filter(!Target == "No Protein")
analytes <- analytes %>% filter(Organism == "Human")
analytes <- analytes %>% filter(Type == "Protein")

# read in the files that will link repository ID (column A) to Somalogic ID (column C)
ids1 <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Somalogic repository link/Omics-Petter Ancillary Samples at Colorado LEAD Center - Wash U.csv")
colnames(ids1) <- c("releaseid","material_type","current_label","MASK.ID","Date.Drawn","visnum","location")
ids1$bsi_id <- NA
ids2 <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Somalogic repository link/Omics-Petter Ancillary Samples at NIDDK-TODAY - Wash U.csv")
ids2$MASK.ID <- NA
ids3 <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Somalogic repository link/Omics-Petter Ancillary Samples at NIDDK-TODAY2 - Wash U.csv")
ids3$MASK.ID <- NA
ids <- rbind(ids1, ids2, ids3)
ids$SampleDescription <- ids$current_label

# merge IDs with soma
soma <- left_join(soma,ids,by="SampleDescription")

# samples without a release ID are QC samples
soma <- soma %>% filter(!is.na(releaseid))

# fix date drawn
soma$Date.Drawn <- as.Date(soma$Date.Drawn,format="%m/%d/%Y")

# merge in bariatric surgery data
tme <- read.csv("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Clinical data/TODAY2/TME.csv")
tme <- tme %>% filter(TMETYPE==1)
tme <- tme %>% select(RELEASEID,TMETYPE,DAYSTOTME)
colnames(tme) <- c("releaseid","TMETYPE","DAYSTOTME")
soma <- left_join(soma, tme, by="releaseid")
# get rid of the 10 yr visit if DAYSTOTME<=3650
# need to be sure to keep the baseline visit for those people
# first add a variable that is sequential within ppt so we can tell baseline vs. followup
temp <- as.data.frame(soma[,"releaseid"])
names(temp) <- "releaseid"
temp <- temp %>% group_by(releaseid) 
temp <- temp %>% mutate(visit=seq_along(releaseid)) 
soma$visit <- temp$visit
soma <- soma %>% filter(!(visit==2 & !is.na(DAYSTOTME) & DAYSTOTME<3652))
# drop bariatric surgery variables
soma <- soma %>% select(-c(TMETYPE,DAYSTOTME))

# Save
save(soma,file = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Somalogic data raw/soma.Rdata")
save(analytes,file = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Somalogic data raw/analytes.Rdata")
