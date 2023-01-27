library(SomaDataIO)
library(stringr)
library(readxl)

if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = "/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad"
}
setwd(home_dir)

# read in data - we want to use the fully processed, normalized file ending in "anmlSMP.adat"
soma <- read_adat("./Local cohort Somalogic data/WUS-22-002_v4.1_EDTAPlasma_hybNorm_medNormInt_plateScale_calibrate_anmlQC_qcCheck_anmlSMP.adat")
analytes <- getAnalyteInfo(soma)

olink_plasma = read_excel("./Olink Data/CROCODILE, IMPROVE, Renal HEIR_Plasma_Summary _1212_2022.xlsx")
olink_urine = read_excel("./Olink Data/CROCODILE, IMPROVE, Renal HEIR_Urine_Summary _1212_2022.xlsx")

# filter out Q/C samples
soma <- soma %>% filter(!is.na(SampleDescription))

# create dataframes for each study
croc_soma <- soma %>% filter(str_sub(SampleDescription,1,3)=="CRC")
croc_olink_plasma = olink_plasma %>% filter(STUDY=="CROCODILE")
croc_olink_urine = olink_urine %>% filter(STUDY=="CROCODILE")

improve_soma <- soma %>% filter(str_sub(SampleDescription,1,4)=="IT2D")
improve_olink_plasma = olink_plasma %>% filter(STUDY=="IMPROVE")
improve_olink_urine = olink_urine %>% filter(STUDY=="IMPROVE")

rh_soma <- soma %>% filter(str_sub(SampleDescription,1,2)=="RH")
rh_olink_plasma = olink_plasma %>% filter(STUDY=="Renal HEIR")
rh_olink_urine = olink_urine %>% filter(STUDY=="Renal HEIR")

pima_soma <- soma %>% filter(str_sub(SampleDescription,1,4)=="CKDS")

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

# save Pima
save(pima_soma,file = "./Pima/Somalogic data/pima_soma.Rdata")
save(analytes,file = "./Pima/Somalogic data/analytes.Rdata")
