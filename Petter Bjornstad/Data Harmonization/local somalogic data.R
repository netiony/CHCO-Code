library(SomaDataIO)
library(stringr)

if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad/Local cohort Somalogic data"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Local cohort Somalogic data/"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = ""
}
setwd(home_dir)

# read in data - we want to use the fully processed, normalized file ending in "anmlSMP.adat"
soma <- read_adat("./WUS-22-002_v4.1_EDTAPlasma_hybNorm_medNormInt_plateScale_calibrate_anmlQC_qcCheck_anmlSMP.adat")
analytes <- getAnalyteInfo(soma)

# filter out Q/C samples
soma <- soma %>% filter(!is.na(SampleDescription))

# create dataframes for each study
# not sure which study IDs CKDSxxx belong to
croc_soma <- soma %>% filter(str_sub(SampleDescription,1,3)=="CRC")
improve_soma <- soma %>% filter(str_sub(SampleDescription,1,4)=="IT2D")
rh_soma <- soma %>% filter(str_sub(SampleDescription,1,2)=="RH")
pima_soma <- soma %>% filter(str_sub(SampleDescription,1,4)=="CKDS")

# Save CROCODILE
if(Sys.info()["sysname"] == "Windows"){
  dir = ""
} else if (Sys.info()["sysname"] == "Linux"){
  dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/CROCODILE/Somalogic data/"
} else if (Sys.info()["sysname"] == "Darwin"){
  dir = ""
}
setwd(dir)
save(croc_soma,file = "./croc_soma.Rdata")
save(analytes,file = "./analytes.Rdata")

# Save IMPROVE
if(Sys.info()["sysname"] == "Windows"){
  dir = ""
} else if (Sys.info()["sysname"] == "Linux"){
  dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/IMPROVE T2D/Somalogic data/"
} else if (Sys.info()["sysname"] == "Darwin"){
  dir = ""
}
setwd(dir)
save(improve_soma,file = "./improve_soma.Rdata")
save(analytes,file = "./analytes.Rdata")

# save RH
if(Sys.info()["sysname"] == "Windows"){
  dir = ""
} else if (Sys.info()["sysname"] == "Linux"){
  dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Renal HERITAGE/Somalogic data/"
} else if (Sys.info()["sysname"] == "Darwin"){
  dir = ""
}
setwd(dir)
save(rh_soma,file = "./rh_soma.Rdata")
save(analytes,file = "./analytes.Rdata")

# save Pima
if(Sys.info()["sysname"] == "Windows"){
  dir = ""
} else if (Sys.info()["sysname"] == "Linux"){
  dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Pima/Somalogic data/"
} else if (Sys.info()["sysname"] == "Darwin"){
  dir = ""
}
setwd(dir)
save(rh_soma,file = "./pima_soma.Rdata")
save(analytes,file = "./analytes.Rdata")
