library(SomaDataIO)
library(stringr)
library(arsenal)
library(Hmisc)
library(dplyr)
library(tidyr)
library(purrr)
library(limma)

if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/"
}
setwd(home_dir)

###################
# ANML ADAT FILES #
###################

# read in data - we want to use the fully processed, normalized file ending in "anmlSMP.adat"
soma <- read_adat("./Local cohort Somalogic data/WUS-22-002/WUS-22-002_v4.1_EDTAPlasma_hybNorm_medNormInt_plateScale_calibrate_anmlQC_qcCheck_anmlSMP.adat")
soma <- soma %>% select(-Optional2)
soma <- soma %>% filter(!is.na(SampleDescription))
analytes <- getAnalyteInfo(soma)
analytes <- analytes %>% select(AptName,SeqId,SeqIdVersion,SomaId,TargetFullName,Target,UniProt,EntrezGeneID,EntrezGeneSymbol,Organism,Units,Type)
# remove fc mouse and no protein
drop <- analytes %>% filter(Target == "Fc_MOUSE" | Target == "No Protein" | !(Organism == "Human") | !(Type == "Protein"))
apt_drop <- drop$AptName
soma <- soma %>% select(!all_of(apt_drop))
analytes <- analytes %>% filter(!Target == "Fc_MOUSE")
analytes <- analytes %>% filter(!Target == "No Protein")
analytes <- analytes %>% filter(Organism == "Human")
analytes <- analytes %>% filter(Type == "Protein")

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

# read in PANTHER data
panther <- read_adat("./Local cohort Somalogic data/PANTHER/20240126_597_Bjornstad_SOMAscan7k_WUS-24-002_data_export/WUS-24-002_2024-01-26_Somalogic_standardized_files/WUS_24_002_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")
panther <- panther %>% filter(!is.na(SampleDescription))
panther <- panther %>% select(-Optional1)
# remove fc mouse and no protein
panther <- panther %>% select(!all_of(apt_drop))
# one ID changed
panther$SampleDescription <- ifelse(panther$SampleDescription == "PAN-14-C", "PAN-14-O", panther$SampleDescription)

# filter out Q/C samples
soma <- rbind(soma,soma2,panther)
# delete Pima data
soma <- soma %>% filter(!str_detect(SampleDescription,"CKDS"))
# fix sample IDs on a few RH2 participants who changed groups
soma <- soma %>%  mutate(
  SampleDescription = case_when(
    SampleDescription == "RH2-13-T" ~ "RH2-13-O",
    SampleDescription == "RH2-06-T" ~ "RH2-06-O",
    SampleDescription == "RH2-07-T" ~ "RH2-07-O",
    SampleDescription == "RH2-38-O" ~ "RH2-38-T",
    .default = SampleDescription
  )
)

soma_combined <- soma 

# pull out KNIGHT samples only
soma_anml_keep <- soma_combined %>% filter(grepl("KGHT",SampleDescription) | grepl("SHB",SampleDescription))

# read in linkage file for Natalie's studies
link <- read.csv("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/KNIGHT/Somalogic data/BCF-23-091 Linker File 10.05.2023.csv")
link$Study.ID <- link$Sequence.Number
redcap <- read.csv("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/KNIGHT/Somalogic data/200572DataKNIGHT_DATA_LABELS_2023-12-10_1551.csv")
redcap <- redcap %>% filter(Event.Name == "Screening")
redcap <- redcap %>% select(Study.ID, Sex.assigned.at.birth)
link <- merge(link, redcap, by="Study.ID", all.x = T, all.y = F)
link$SampleDescription <- link$Barcode
link <- link %>% select(SampleDescription, Study.ID, Timepoint.Label, Sex.assigned.at.birth)
soma_anml_keep <- full_join(soma_anml_keep, link, by="SampleDescription")
soma_anml_keep$SampleDescription <- str_trim(soma_anml_keep$SampleDescription)

# create new variable for harmonized time point
soma_anml_keep$time <-  
  case_when(
    soma_anml_keep$TimePoint == "V1" ~ 1,
    soma_anml_keep$TimePoint == "V2" ~ 2,
    soma_anml_keep$Timepoint.Label == "Baseline" ~ 1,
    soma_anml_keep$Timepoint.Label == "M3" ~ 2,
    .default = NA
  )

soma_anml_keep$group <-  
  case_when(
    str_sub(soma_anml_keep$SampleDescription, 1, 2) == "MV" ~ "MTF",
    str_sub(soma_anml_keep$SampleDescription, 1, 2) == "VM" ~ "FTM",
    soma_anml_keep$Sex.assigned.at.birth == "Male" ~ "MTF",
    soma_anml_keep$Sex.assigned.at.birth == "Female" ~ "FTM",
    .default = NA
  )

# identify columns corresponding to proteins
#is_seq <- function(.x) grepl("^seq\\.[0-9]{4}", .x) # regex for analytes
is_seq <- function(.x) grepl("seq", .x)
seq <- is_seq(names(soma_anml_keep))
# log transform
soma_anml_keep <- soma_anml_keep %>% modify_if(is_seq(names(.)), log)

# Calculate change in each protein
df_diff <- soma_anml_keep %>%
  arrange(Study.ID, time) %>%
  group_by(Study.ID, group) %>%
  reframe(across(contains("seq"), ~ diff(as.numeric(.x))))

# First MTF
# Limma setup
# Y needs "genes" in rows
y <- df_diff %>%
  filter(group == "MTF") %>%
  dplyr::select(contains("seq")) %>%
  t()
# Design matrix
design_mat <- model.matrix(~1, data = df_diff[df_diff$group == "MTF",])
# Fit
fit <- lmFit(y, design_mat)
fit <- eBayes(fit)
res <- topTable(fit, coef = 1, number = dim(y)[1], sort.by = "p")
res$Target <- analytes$Target[match(rownames(res), analytes$AptName)]
res$TargetFullName <- analytes$TargetFullName[match(rownames(res), analytes$AptName)]
res$UniProt <- analytes$UniProt[match(rownames(res), analytes$AptName)]
write.csv(res, file = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/KNIGHT/Testing ANML vs median/limma_MTF_results_ANML.csv", row.names = F)

# Then FTM
# Limma setup
# Y needs "genes" in rows
y <- df_diff %>%
  filter(group == "FTM") %>%
  dplyr::select(contains("seq")) %>%
  t()
# Design matrix
design_mat <- model.matrix(~1, data = df_diff[df_diff$group == "FTM",])
# Fit
fit <- lmFit(y, design_mat)
fit <- eBayes(fit)
res <- topTable(fit, coef = 1, number = dim(y)[1], sort.by = "p")
res$Target <- analytes$Target[match(rownames(res), analytes$AptName)]
res$TargetFullName <- analytes$TargetFullName[match(rownames(res), analytes$AptName)]
res$UniProt <- analytes$UniProt[match(rownames(res), analytes$AptName)]
write.csv(res, file = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/KNIGHT/Testing ANML vs median/limma_FTM_results_ANML.csv", row.names = F)

###########################
# MEDIAN NORMALIZED FILES #
###########################

# read in data
soma_med1 <- read.delim("./Local cohort Somalogic data/WUS-23-004/WUS_23_004_GTAC_analyzed/linear.RFU_WUS_23_004_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_quantile.txt")
# read in samplesheet
samplesheet <- read.delim("./Local cohort Somalogic data/WUS-23-004/WUS_23_004_GTAC_analyzed/samplesheet_WUS_23_004.txt")
samplesheet$ArrayId <- samplesheet$X

# need to link samplesheet$X to soma_med1 column names
soma_med1 <- soma_med1 %>% select(!(X))
soma_med1 <- soma_med1 %>% select(!(EntrezGeneSymbol:Dilution))
rownames(soma_med1) <- soma_med1$AptName
soma_med1 <- soma_med1 %>% select(!(AptName))
soma_med1 <- as.data.frame(t(soma_med1))
soma_med1$ArrayId <- rownames(soma_med1)
soma_med1$ArrayId <- str_remove(soma_med1$ArrayId, "X")
soma_med <- left_join(soma_med1, samplesheet, by="ArrayId")
# pull out KNIGHT samples only
soma_med_keep <- soma_med %>% filter(grepl("KGHT",SampleDescription) | grepl("SHB",SampleDescription))
soma_med_keep$SampleDescription <- str_replace(soma_med_keep$SampleDescription , '(.*?)-(.*?)', '')

# merge with link file
soma_med_keep <- full_join(soma_med_keep, link, by="SampleDescription")
soma_med_keep$SampleDescription <- str_trim(soma_med_keep$SampleDescription)

# create new variable for harmonized time point
soma_med_keep$time <-  
  case_when(
    soma_med_keep$TimePoint == "V1" ~ 1,
    soma_med_keep$TimePoint == "V2" ~ 2,
    soma_med_keep$Timepoint.Label == "Baseline" ~ 1,
    soma_med_keep$Timepoint.Label == "M3" ~ 2,
    .default = NA
  )

soma_med_keep$group <-  
  case_when(
    str_sub(soma_med_keep$SampleDescription, 1, 2) == "MV" ~ "MTF",
    str_sub(soma_med_keep$SampleDescription, 1, 2) == "VM" ~ "FTM",
    soma_med_keep$Sex.assigned.at.birth == "Male" ~ "MTF",
    soma_med_keep$Sex.assigned.at.birth == "Female" ~ "FTM",
    .default = NA
  )

# identify columns corresponding to proteins
#is_seq <- function(.x) grepl("^seq\\.[0-9]{4}", .x) # regex for analytes
is_seq <- function(.x) grepl("seq", .x)
seq <- is_seq(names(soma_med_keep))
# log transform
soma_med_keep <- soma_med_keep %>% modify_if(is_seq(names(.)), log)

# Calculate change in each protein
df_diff <- soma_med_keep %>%
  arrange(Study.ID, time) %>%
  group_by(Study.ID, group) %>%
  reframe(across(contains("seq"), ~ diff(as.numeric(.x))))

# First MTF
# Limma setup
# Y needs "genes" in rows
y <- df_diff %>%
  filter(group == "MTF") %>%
  dplyr::select(contains("seq")) %>%
  t()
# Design matrix
design_mat <- model.matrix(~1, data = df_diff[df_diff$group == "MTF",])
# Fit
fit <- lmFit(y, design_mat)
fit <- eBayes(fit)
res <- topTable(fit, coef = 1, number = dim(y)[1], sort.by = "p")
res$Target <- analytes$Target[match(rownames(res), analytes$AptName)]
res$TargetFullName <- analytes$TargetFullName[match(rownames(res), analytes$AptName)]
res$UniProt <- analytes$UniProt[match(rownames(res), analytes$AptName)]
write.csv(res, file = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/KNIGHT/Testing ANML vs median/limma_MTF_results_median.csv", row.names = F)

# Then FTM
# Limma setup
# Y needs "genes" in rows
y <- df_diff %>%
  filter(group == "FTM") %>%
  dplyr::select(contains("seq")) %>%
  t()
# Design matrix
design_mat <- model.matrix(~1, data = df_diff[df_diff$group == "FTM",])
# Fit
fit <- lmFit(y, design_mat)
fit <- eBayes(fit)
res <- topTable(fit, coef = 1, number = dim(y)[1], sort.by = "p")
res$Target <- analytes$Target[match(rownames(res), analytes$AptName)]
res$TargetFullName <- analytes$TargetFullName[match(rownames(res), analytes$AptName)]
res$UniProt <- analytes$UniProt[match(rownames(res), analytes$AptName)]
write.csv(res, file = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/KNIGHT/Testing ANML vs median/limma_FTM_results_median.csv", row.names = F)
