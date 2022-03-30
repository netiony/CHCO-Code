if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad/TODAY subaward"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = ""
}
setwd(home_dir)

# read in data - we want to use the fully processed, normalized file ending in "anmlSMP.adat"
soma <- read_adat("./Somalogic data raw/WUS-22-001_Somalogic_normalized/WUS-22-001_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")

# read in the files that will link repository ID (column A) to Somalogic ID (column C)
ids1 <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at Colorado LEAD Center - Wash U.csv")
colnames(ids1) <- c("releaseid","material_type","current_label","MASK.ID","Date.Drawn","visnum","location")
ids1$bsi_id <- NA
ids2 <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at NIDDK-TODAY - Wash U.csv")
ids2$MASK.ID <- NA
ids3 <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at NIDDK-TODAY2 - Wash U.csv")
ids3$MASK.ID <- NA
ids <- rbind(ids1, ids2, ids3)
ids$SampleDescription <- ids$current_label

# merge IDs with soma
soma <- merge(soma,ids,by="SampleDescription",all.x = T,all.y = F)

# samples without a release ID are QC samples
soma <- soma %>% filter(!is.na(releaseid))

# fix date drawn
soma$Date.Drawn <- as.Date(soma$Date.Drawn,format="%m/%d/%Y")

# Save
save(soma,file = "./Somalogic data raw/soma.Rdata")
