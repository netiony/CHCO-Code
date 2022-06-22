library(dplyr)
library(berryFunctions)

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
nih_urine <- read.csv("~/Metabolomic data/NIDDK_AA_20220427_20220526_Normalized_AF_urine_norm_creatinine.csv")

# read in LEAD samples
lead_urine <- read.csv("~/Metabolomic data/Lead_AA_20220322_20220510_Normalized_AF_urine_norm_creatinine.csv")

# different number of columns
a <- as.data.frame(colnames(nih_urine))
a <- as.data.frame(a[3:nrow(a),])
colnames(a) <- "a"
a <- a %>% arrange("a")
a <- insertRows(a, 31:32 , new = NA)
b <- as.data.frame(colnames(lead_urine))
b <- as.data.frame(b[3:nrow(b),])
colnames(b) <- "b"
b <- b %>% arrange("b")
c <- cbind(a,b)

# merge
# urine <- rbind(nih,lead)

# remove Q/C samples


####################
# plasma           #
####################

# read in NIH samples
nih_plasma <- read.csv("~/Metabolomic data/NIDDK_AA_20220427_20220526_Normalized_AF_plasma.csv")

# read in LEAD samples
lead_plasma <- read.csv("~/Metabolomic data/Lead_AA_20220427_20220526_Normalized_AF_plasma.csv")

# different number of columns
a <- as.data.frame(colnames(nih_urine))
a <- as.data.frame(a[3:nrow(a),])
colnames(a) <- "a"
a <- a %>% arrange("a")
a <- insertRows(a, 31:32 , new = NA)
b <- as.data.frame(colnames(lead_urine))
b <- as.data.frame(b[3:nrow(b),])
colnames(b) <- "b"
b <- b %>% arrange("b")
c <- cbind(a,b)

# merge

# remove Q/C samples


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
save(analytes,file = "./Somalogic data raw/analytes.Rdata")
