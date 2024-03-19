library(SomaDataIO)
library(stringr)
library(arsenal)
library(Hmisc)
library(dplyr)
library(tidyr)
library(purrr)

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

# read in harmonized dataset to get group information
# df <- read.csv("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", 
#                na.strings = c(" ", "", "-9999",-9999))
# coenroll_id <- read.csv("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Renal HERITAGE/Data_Cleaned/coenrolled_ids.csv") %>%
#   pivot_longer(cols = 'improve_id':'crc_id',
#                values_to = "record_id") %>% 
#   dplyr::select(merged_id, record_id) %>%
#   filter(record_id != "")
# df <- df %>%
#   dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
#                    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
#                    .by = c(record_id, visit)) %>%
#   filter(participation_status!="Removed"|is.na(participation_status)) %>%
#   left_join(coenroll_id)
# df <- df %>% 
#   dplyr::select(record_id, co_enroll_id, merged_id, study, visit, date, group)
# # Merge
# soma_harmonized$record_id <- soma_harmonized$SampleDescription
# soma_harmonized <- soma_harmonized %>% mutate(
#   record_id = case_when(
#     str_detect(record_id, "IT2D") ~ str_replace(record_id, "IT2D-", "IT_"),
#     !str_detect(record_id, "IT2D") ~ record_id
#   )
# )
# soma_harmonized <- soma_harmonized %>%  mutate(
#   visit = case_when(
#     TimePoint == "BL" ~ "baseline",
#     TimePoint == "Baseline" ~ "baseline",
#     TimePoint == "3M" ~ "3_months_post_surgery",
#     TimePoint == "12M" ~ "12_months_post_surgery"
#   )
# )
# soma_harmonized <- left_join(soma_harmonized, df, by = c("record_id", "visit"))

# some extra spaces in the sample IDs
soma_combined$SampleDescription <- str_trim(soma_combined$SampleDescription)
# add labels
labels <- paste0("log(",analytes$EntrezGeneSymbol[match(colnames(soma_combined), analytes$AptName)],")")
labels[labels == "log(NA)"] <- ""
label(soma_combined) <- as.list(labels)

# write copy of the entire local SomaScan dataframe to the data harmonization folder
save(soma_combined, file = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Combined SomaScan/soma_combined_anml.RData")
save(analytes, file = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Combined SomaScan/analytes.RData")



###########################
# MEDIAN NORMALIZED FILES #
###########################

# read in data
soma_med1 <- read.delim("./Local cohort Somalogic data/WUS-22-002/linear.RFU_WUS-22-002_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_quantile.txt")
# read in samplesheet
samplesheet <- read.delim("./Local cohort Somalogic data/WUS-22-002/samplesheet_WUS-22-002.txt")
samplesheet$ArrayId <- samplesheet$X

# need to link samplesheet$X to soma_med1 column names
soma_med1 <- soma_med1 %>% select(!(X))
soma_med1 <- soma_med1 %>% select(!(EntrezGeneSymbol:Organism))
rownames(soma_med1) <- soma_med1$AptName
soma_med1 <- soma_med1 %>% select(!(AptName))
soma_med1 <- as.data.frame(t(soma_med1))
soma_med1$ArrayId <- rownames(soma_med1)
soma_med1$ArrayId <- str_remove(soma_med1$ArrayId, "X")
soma_med <- left_join(soma_med1, samplesheet, by="ArrayId")

# apply all the same processing steps as above
# remove fc mouse and no protein
soma_med <- soma_med %>% select(!all_of(apt_drop))


# will have to create an adat file?
# then figure out which study to run a comparison on median vs ANML
# probably KNIGHT since we won't have any harmonization issues