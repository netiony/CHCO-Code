library(SomaDataIO)
library(stringr)
library(dplyr)
library(arsenal)

if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/"
}
setwd(home_dir)

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

# filter out Q/C samples
soma <- rbind(soma,soma2)
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

# df named "soma" has all studies including KNIGHT, but some don't need to be merged w/ harmonized data
# make another df with only local studies to be merged w/ harmonized data
soma_harmonized <- soma %>% filter(!str_detect(SampleDescription,"KGHT") & !str_detect(SampleDescription,"SHB"))

# read in harmonized dataset to get group information
df <- read.csv("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", 
               na.strings = c(" ", "", "-9999",-9999))
coenroll_id <- read.csv("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Renal HERITAGE/Data_Cleaned/coenrolled_ids.csv") %>%
  pivot_longer(cols = 'improve_id':'crc_id',
               values_to = "record_id") %>% 
  dplyr::select(merged_id, record_id) %>%
  filter(record_id != "")
df <- df %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  filter(participation_status!="Removed"|is.na(participation_status)) %>%
  left_join(coenroll_id)
df <- df %>% 
  dplyr::select(record_id, co_enroll_id, merged_id, study, visit, date, group)
# Merge
soma_harmonized$record_id <- soma_harmonized$SampleDescription
soma_harmonized <- soma_harmonized %>% mutate(
  record_id = case_when(
    str_detect(record_id, "IT2D") ~ str_replace(record_id, "IT2D-", "IT_"),
    !str_detect(record_id, "IT2D") ~ record_id
  )
)
soma_harmonized <- soma_harmonized %>%  mutate(
  visit = case_when(
    TimePoint == "BL" ~ "baseline",
    TimePoint == "Baseline" ~ "baseline",
    TimePoint == "3M" ~ "3_months_post_surgery",
    TimePoint == "12M" ~ "12_months_post_surgery"
  )
)
soma_harmonized <- left_join(soma_harmonized, df, by = c("record_id", "visit"))
check <- soma_harmonized %>% select(record_id, merged_id, visit, group, study)
check <- check %>% filter(visit == "baseline")
check <- check %>% filter(study %in% c("IMPROVE", "RENAL-HEIR", "RENAL-HEIRitage"))
write.csv(check, "/Users/pylell/Documents/Temp/check_local_somalogic.csv", row.names = F)

# write copy of the entire local SomaScan dataframe to the data harmonization folder
save(soma_harmonized, file = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Harmonized SomaScan/soma_harmonized.RData")
save(analytes, file = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Harmonized SomaScan/analytes.RData")

# read in Olink data
olink_plasma = read.csv("./Olink Data/Data_Clean/plasma_cleaned.csv")
olink_urine = read.csv("./Olink Data/Data_Clean/urine_cleaned.csv")

# create dataframes for each study
croc_soma <- soma %>% filter(str_detect(SampleDescription,"CRC"))
croc_olink_plasma = olink_plasma %>% filter(grepl("CRC",record_id))
croc_olink_urine = olink_urine %>% filter(grepl("CRC",record_id))
improve_soma <- soma %>% filter(str_detect(SampleDescription,"IT2D"))
improve_olink_plasma = olink_plasma %>% filter(grepl("IT",record_id))
improve_olink_urine = olink_urine %>% filter(grepl("IT",record_id))
rh_soma <- soma %>% filter(str_detect(SampleDescription,"RH"))
rh_olink_plasma = olink_plasma %>% filter(grepl("RH",record_id))
rh_olink_urine = olink_urine %>% filter(grepl("RH",record_id))
pen_soma <- soma %>% filter(str_detect(SampleDescription,"PEN"))
knight_soma <- soma %>% filter(str_detect(SampleDescription,"KGHT") | str_detect(SampleDescription,"SHB"))
knight_soma$SampleDescription <- str_trim(knight_soma$SampleDescription)

# Save individual files for each study
# Save CROCODILE
save(croc_soma,file = "./CROCODILE/Somalogic data/croc_soma.Rdata")
save(analytes,file = "./CROCODILE/Somalogic data/analytes.Rdata")
save(croc_olink_plasma,file = "./CROCODILE/Olink Data/croc_olink_plasma.Rdata")
save(croc_olink_urine,file = "./CROCODILE/Olink Data/croc_olink_urine.Rdata")
# Save IMPROVE
save(improve_soma,file = "./IMPROVE T2D/Somalogic data/improve_soma.Rdata")
save(analytes,file = "./IMPROVE T2D/Somalogic data/analytes.Rdata")
save(improve_olink_plasma,file = "./IMPROVE T2D/Olink Data/improve_olink_plasma.Rdata")
save(improve_olink_urine,file = "./IMPROVE T2D/Olink Data/improve_olink_urine.Rdata")
# save RH
save(rh_soma,file = "./Renal HERITAGE/Somalogic data/rh_soma.Rdata")
save(analytes,file = "./Renal HERITAGE/Somalogic data/analytes.Rdata")
save(rh_olink_plasma,file = "./Renal HERITAGE/Olink Data/rh_olink_plasma.Rdata")
save(rh_olink_urine,file = "./Renal HERITAGE/Olink Data/rh_olink_urine.Rdata")
# save PENGUIN
save(pen_soma,file = "./PENGUIN/Somalogic data/penguin_soma.Rdata")
save(analytes,file = "./PENGUIN/Somalogic data/analytes.Rdata")
# save KNIGHT
save(knight_soma,file = "./KNIGHT/Somalogic data/knight_soma.Rdata")
save(analytes,file = "./KNIGHT/Somalogic data/analytes.Rdata")
