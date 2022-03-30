if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad/TODAY subaward"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/BDC/Projects"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = "/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward"
}

setwd(home_dir)

# COMORB dataset
comorb <- read.csv("./Clinical data/COMORB.csv")
comorb$MIC.OR.MAC <- ifelse(comorb$MAC==1 | comorb$MIC==1,1,
                            ifelse(is.na(comorb$MAC) & is.na(comorb$MIC), NA, 0))

# Save
save(comorb,file = "./Clinical data/comorb.Rdata")
