library(SomaDataIO)
library(stringr)
library(dplyr)
library(arsenal)

if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/"
}
setwd(home_dir)

# read in data - we want to use the fully processed, normalized file ending in "anmlSMP.adat"
soma <- read_adat("./Local cohort Somalogic data/WUS-22-002/WUS-22-002_v4.1_EDTAPlasma_hybNorm_medNormInt_plateScale_calibrate_anmlQC_qcCheck_anmlSMP.adat")
soma <- soma %>% select(-Optional2)
soma <- soma %>% filter(!is.na(SampleDescription))
analytes <- getAnalyteInfo(soma)
analytes <- analytes %>% select(AptName,SeqId,SeqIdVersion,SomaId,TargetFullName,Target,UniProt,EntrezGeneID,EntrezGeneSymbol,Organism,Units,Type)

# remove fc mouse and no protein
drop <- analytes %>% filter(Target == "Fc_MOUSE" | Target == "No Protein")
apt_drop <- drop$AptName
soma <- soma %>% select(!all_of(apt_drop))
analytes <- analytes %>% filter(!Target == "Fc_MOUSE")
analytes <- analytes %>% filter(!Target == "No Protein")

# read in 2nd dataset
soma2 <- read_adat("./Local cohort Somalogic data/WUS-23-004/WUS_23_004_2023-11-15/WUS_23_004_v4.1_EDTAPlasma.hybNorm_medNormInt_plateScale_calibrate_anmlQC_qcCheck_anmlSMP.adat")
soma2 <- soma2 %>% filter(!is.na(SampleDescription))
# delete extraneous stuff in the sample description
soma2$SampleDescription <- str_sub(soma2$SampleDescription, 12)
soma2$SampleDescription <- str_remove(soma2$SampleDescription, " ")
soma2$SampleDescription <- str_remove_all(soma2$SampleDescription, "\\(\\S+")
# there is one duplicate sample - delete the second result
soma2 <- soma2 %>% filter(!SampleDescription=="RH-14-O")
# remove FC_Mouse and no protein
soma2 <- soma2 %>% select(!all_of(apt_drop))
analytes2 <- getAnalyteInfo(soma2)
# 2nd analytes file is esssentially the same as the first except for some batch specific information we don't need
# will keep the first file

olink_plasma = read.csv("./Olink Data/Data_Clean/plasma_cleaned.csv")
olink_urine = read.csv("./Olink Data/Data_Clean/urine_cleaned.csv")

# filter out Q/C samples
soma <- rbind(soma,soma2)
# delete Pima data
soma <- soma %>% filter(!str_detect(SampleDescription,"CKDS"))

# create dataframes for each study
croc_soma <- soma %>% filter(str_detect(SampleDescription,"CRC"))
croc_olink_plasma = olink_plasma %>% filter(grepl("CRC",record_id))
croc_olink_urine = olink_urine %>% filter(grepl("CRC",record_id))

improve_soma <- soma %>% filter(str_detect(SampleDescription,"IT2D"))
improve_olink_plasma = olink_plasma %>% filter(grepl("IT",record_id))
improve_olink_urine = olink_urine %>% filter(grepl("IT",record_id))

rh_soma <- soma %>% filter(str_detect(SampleDescription,"RH"))
rh_olink_plasma = olink_plasma %>% filter(grepl("RH",record_id))
rh_olink_urine = olink_urine %>% filter(grepl("RH",record_id))

pen_soma <- soma %>% filter(str_detect(SampleDescription,"PEN"))

knight_soma <- soma %>% filter(str_detect(SampleDescription,"KGHT") | str_detect(SampleDescription,"SHB"))
knight_soma$SampleDescription <- str_trim(knight_soma$SampleDescription)

# Save CROCODILE
save(croc_soma,file = "./CROCODILE/Somalogic data/croc_soma.Rdata")
save(analytes,file = "./CROCODILE/Somalogic data/analytes.Rdata")
save(croc_olink_plasma,file = "./CROCODILE/Olink Data/croc_olink_plasma.Rdata")
save(croc_olink_urine,file = "./CROCODILE/Olink Data/croc_olink_urine.Rdata")

# Save IMPROVE
save(improve_soma,file = "./IMPROVE T2D/Somalogic data/improve_soma.Rdata")
save(analytes,file = "./IMPROVE T2D/Somalogic data/analytes.Rdata")
save(improve_olink_plasma,file = "./IMPROVE T2D/Olink Data/improve_olink_plasma.Rdata")
save(improve_olink_urine,file = "./IMPROVE T2D/Olink Data/improve_olink_urine.Rdata")

# save RH
save(rh_soma,file = "./Renal HERITAGE/Somalogic data/rh_soma.Rdata")
save(analytes,file = "./Renal HERITAGE/Somalogic data/analytes.Rdata")
save(rh_olink_plasma,file = "./Renal HERITAGE/Olink Data/rh_olink_plasma.Rdata")
save(rh_olink_urine,file = "./Renal HERITAGE/Olink Data/rh_olink_urine.Rdata")

# save PENGUIN
save(pen_soma,file = "./PENGUIN/Somalogic data/penguin_soma.Rdata")
save(analytes,file = "./PENGUIN/Somalogic data/analytes.Rdata")

# save KNIGHT
save(knight_soma,file = "./KNIGHT/Somalogic data/knight_soma.Rdata")
save(analytes,file = "./KNIGHT/Somalogic data/analytes.Rdata")
