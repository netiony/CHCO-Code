library(dplyr)
library(REDCapR)
library(ATCapiR)
library(striprtf)

# REDCap token for RPC2
rpc2_tok <- api_tok[api_tok$Study == "RPC2",]$Token
# uri for REDCap
uri <- "https://redcap.ucdenver.edu/api/"
# read from REDCap
rpc2 <- redcap_read(redcap_uri = uri, token = rpc2_tok)$data

# read ATC codes for meds
rxnorm <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Ye Ji Choi/RXNORM.csv") 

# read meds from Phoom
meds_biopsy <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/Meds/biopsy_meds_100324.csv")%>%
  dplyr::select(c(-mrn)) %>%
  melt(id = 1:3, variable.name = "Preferred.Label", value.name = "med_yn") %>%
  filter(!is.na(med_yn)) %>%
  filter(med_yn != "") %>%
  filter(grepl("RPC-", record_id)) %>% 
  mutate(age_current=NA, 
         gender=NA, race=NA, ethnicity=NA, 
         medname=NA, medindication=NA, meddose=NA, 
         medfreq=NA, medstart=NA, medcontinues=NA, medstop=NA,
         source = "EPIC") %>%
  rename(subject_id = "record_id")

# subset to data of interest
# IDs requested by Tomas
rpc2_ids<- c("RPC-02", "RPC-03", "RPC-05", "RPC-07", "RPC-16", "RPC-17", 
             "RPC-21", "RPC-25", "RPC-26", "RPC-27", "RPC-29")

race_vars <- c("race___1", "race___2", "race___3", "race___4", "race___5", "race___6", "race___7")

rpc2_subset <- rpc2 %>%
  filter(subject_id %in% rpc2_ids) %>%
  rowwise() %>%
  mutate(race = case_when(sum(c_across(all_of(race_vars))) > 1 ~ "More than one",
                          race___1 == 1 ~ "American Indian or Alaskan Native",
                          race___2 == 1 ~ "Asian",
                          race___3 == 1 ~ "Hawaiian or Pacific Islander",
                          race___4 == 1 ~ "Black or African American",
                          race___5 == 1 ~ "White",
                          race___6 == 1 ~ "Unknown",
                          race___7 == 1 ~ "Other"),
         ethnicity = case_when(ethnicity___1 == 1 ~ "Hispanic",
                               ethnicity___2 == 1 ~ "Non-Hispanic",
                               ethnicity___3 == 1 ~ "Unknown/Not reported"),
         medname = as.character(medname)) %>% 
  left_join(rxnorm, by = c("medname" = "Class.ID")) %>%
  dplyr::select(subject_id, age_current, gender, race, ethnicity, 
                medname, Preferred.Label, medindication, meddose, 
                medfreq, medstart, medcontinues, medstop) %>%
  mutate(Preferred.Label = case_when(medname == "798928" ~ "acetaminophen 500 MG Oral Tablet [Tactinal]", 
                                     T ~ Preferred.Label),
         gender = case_when(gender == 0 ~ "Female", gender == 1 ~ "Male", gender == 0 ~ "Other"),
         kidneybx_1 = NA, kidneybx_2 = NA, med_yn = "Yes",
         source = "Self-reported") %>%
  rbind(meds_biopsy) %>%
  ungroup() %>% group_by(subject_id) %>%
  fill(age_current, gender, race, ethnicity, kidneybx_1, kidneybx_2, .direction = "updown") %>%
  filter(!is.na(Preferred.Label)) %>%
  arrange(subject_id)
  
last_columns <- tail(names(rpc2_subset), 4)

rpc2_subset <- rpc2_subset %>% 
  dplyr::select(1:5, all_of(last_columns), setdiff(names(rpc2_subset), c(1:5, last_columns)))

# Petter wants just the EPIC pulled meds and in wide format
rpc2_subset_epic <- rpc2_subset %>%
  filter(source == "EPIC") %>%
  dplyr::select(1:11, -medname, -source) %>%
  pivot_wider(names_from = "Preferred.Label", values_from = "med_yn") %>%
  left_join(subset(rpc2, select = c(subject_id, screen_egfr)))
  

write.csv(rpc2_subset_epic, "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/rpc2_biopsy_meds.csv", row.names = F, na = "")
