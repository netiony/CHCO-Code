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
nih_urine$site <- "NIH"

# read in LEAD samples
lead_urine <- read.csv("./Metabolomic data/Lead_AA_20220322_20220620_Normalized_AF_urine.csv")
lead_urine$site <- "LEAD"

####################
# plasma           #
####################

# read in NIH samples
nih_plasma <- read.csv("./Metabolomic data/NIDDK_AA_20220427_20220620_Normalized_AF_plasma.csv")

# read in LEAD samples
lead_plasma <- read.csv("./Metabolomic data/Lead_AA_20220322_20220620_Normalized_AF_plasma.csv")


######################
# link IDs and merge #
######################

# read in the files that will link repository ID to sample ID
ids_lead <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at Colorado LEAD Center - UT Health San Antonio.csv")
ids_lead$bsi_id <- NA
ids_niddk_today <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at NIDDK-TODAY - UT Health San Antonio.csv")
ids_niddk_today$MASK.ID <- NA
ids_niddk_today2 <- read.csv("./Somalogic repository link/Omics-Petter Ancillary Samples at NIDDK-TODAY2 - UT Health San Antonio.csv")
ids_niddk_today2$MASK.ID <- NA
ids_niddk <- rbind(ids_niddk_today, ids_niddk_today2)

# merge to urine
# first NIH - merge Freezerworks ID in results file to current label in NIDDK file
nih_urine$current_label <- nih_urine$Freezerworks.ID
nih_urine <- merge(nih_urine,ids_niddk,by="current_label",all.x = T, all.y = F)
nih_urine$SAMPLE_ID <- NA
# for LEAD, need to merge MASK.ID in ID file to sample name in results file
ids_lead$Sample.Name <- str_sub(ids_lead$MASK.ID, 1, 8)
ids_lead$Sample.Name <- str_replace(ids_lead$Sample.Name, "-", "_")
ids_lead$t <- ifelse(str_trim(ids_lead$material_type)=="Urine","_U","")
ids_lead$Sample.Name <- paste0(ids_lead$Sample.Name,ids_lead$t)
ids_lead$t <- NULL
lead_urine <- merge(lead_urine,ids_lead,by="Sample.Name",all.x=T, all.y = F)
lead_urine$current_label <- NA
# combine NIH and LEAD
urine <- rbind(nih_urine,lead_urine)

# merge to plasma
# first NIH
nih_plasma$current_label <- nih_plasma$Freezerworks.ID
nih_plasma <- merge(nih_plasma,ids_niddk,by="current_label",all.x = T, all.y = F)
nih_plasma$SAMPLE_ID <- NA
# LEAD
lead_plasma <- merge(lead_plasma,ids_lead,by="Sample.Name",all.x=T, all.y = F)
lead_plasma$current_label <- NA
# combine NIH and LEAD
plasma <- rbind(nih_plasma,lead_plasma)

# Save
save(urine,file = "./Metabolomic data/urine.Rdata")
save(plasma,file = "./Metabolomic data/plasma.Rdata")
