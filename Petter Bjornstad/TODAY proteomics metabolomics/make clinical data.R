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

# read in the files that will link repository ID (column A) to Somalogic ID (column C)
ids1 <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at Colorado LEAD Center - Wash U.csv")
ids2 <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at NIDDK-TODAY - Wash U.csv")
ids3 <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at NIDDK-TODAY2 - Wash U.csv")
ids <- rbind(ids1, ids2, ids3)

# Save
df = as.data.frame(comorb)
save(comorb,file = "./comorb.Rdata")
