library(dplyr)
library(berryFunctions)
library(stringr)

if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad/TODAY subaward"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = ""
}
setwd(home_dir)

####################
# normalized urine #
####################

# read in NIH samples
nih_urine <- read.csv("./Metabolomic data/NIDDK_AA_20220427_20220620_Normalized_AF_urine.csv")

# read in LEAD samples
lead_urine <- read.csv("./Metabolomic data/Lead_AA_20220322_20220620_Normalized_AF_urine.csv")

# merge
urine <- rbind(nih_urine,lead_urine)

####################
# plasma           #
####################

# read in NIH samples
nih_plasma <- read.csv("./Metabolomic data/NIDDK_AA_20220427_20220620_Normalized_AF_plasma.csv")

# read in LEAD samples
lead_plasma <- read.csv("./Metabolomic data/Lead_AA_20220322_20220620_Normalized_AF_plasma.csv")

# merge
plasma <- rbind(nih_plasma,lead_plasma)

######################
# link IDs and merge #
######################

# read in the files that will link repository ID to sample ID
ids1 <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at Colorado LEAD Center - Wash U.csv")
colnames(ids1) <- c("releaseid","material_type","current_label","MASK.ID","Date.Drawn","visnum","location")
ids1$bsi_id <- NA
ids2 <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at NIDDK-TODAY - Wash U.csv")
ids2$MASK.ID <- NA
ids3 <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at NIDDK-TODAY2 - Wash U.csv")
ids3$MASK.ID <- NA
ids <- rbind(ids1, ids2, ids3)
ids$SampleDescription <- ids$current_label

# merge to urine
urine$current_label <- urine$Freezerworks.ID
urine <- merge(urine,ids,by="current_label",all.x = T, all.y = T)

# merge to plasma
plasma$current_label <- plasma$Freezerworks.ID
plasma <- merge(plasma,ids,by="current_label",all.x = T, all.y = T)

# remove Q/C samples

## old code below - from proteomics

# merge IDs with soma
soma <- merge(soma,ids,by="SampleDescription",all.x = T,all.y = F)

# samples without a release ID are QC samples
soma <- soma %>% filter(!is.na(releaseid))

# fix date drawn
soma$Date.Drawn <- as.Date(soma$Date.Drawn,format="%m/%d/%Y")

# Save
save(soma,file = "./Somalogic data raw/soma.Rdata")
save(analytes,file = "./Somalogic data raw/analytes.Rdata")
