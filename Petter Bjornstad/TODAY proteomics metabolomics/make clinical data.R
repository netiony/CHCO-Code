if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad/TODAY subaward/Clinical data"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/BDC/Projects"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = "/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Clinical data"
}

setwd(home_dir)

# COMORB dataset
comorb <- read.csv("./COMORB.csv")

# Save
df = as.data.frame(comorb)
save(comorb,file = "./comorb.Rdata")
