######################################
######################################

# DECODE KIDNEY Harmonizing Code

# List of studies:
# PANDA, RH, RH2, CASPER, CROC, IMPROVE
# PENGUIN, ATTEMPT, TODAY/TOODAY2,
# TEEN-LABS, T1-DISCO, REMODEL-T1D, RPC2

######################################
######################################

library(Hmisc)
library(dplyr)
library(labelled)

# eGFR function
# Floor age and assign qcr based on age & sex
# qcr <- floor(data[[age]])
# 
# # Assign qcr values
# qcr[qcr == 8] <- 0.46
# qcr[qcr == 9] <- 0.49
# qcr[qcr == 10] <- 0.51
# qcr[qcr == 11] <- 0.53
# qcr[qcr == 12] <- 0.57
# qcr[qcr == 13] <- 0.59
# qcr[qcr == 14] <- 0.61
# 
# qcr[data[[sex]] == "F" & qcr == 15] <- 0.64
# qcr[data[[sex]] == "F" & qcr == 16] <- 0.67
# qcr[data[[sex]] == "F" & qcr == 17] <- 0.69
# qcr[data[[sex]] == "F" & qcr == 18] <- 0.69
# qcr[data[[sex]] == "F" & qcr >= 19] <- 0.70
# 
# qcr[data[[sex]] == "M" & qcr == 15] <- 0.72
# qcr[data[[sex]] == "M" & qcr == 16] <- 0.78
# qcr[data[[sex]] == "M" & qcr == 17] <- 0.82
# qcr[data[[sex]] == "M" & qcr == 18] <- 0.85
# qcr[data[[sex]] == "M" & qcr == 19] <- 0.88
# qcr[data[[sex]] == "M" & qcr > 19] <- 0.90

egfr_fas_cr <- function(serum_creatinine, qcr) {
  return(107.3 / (serum_creatinine / qcr))
}

egfr_fas_cr_cysc <- function(serum_creatinine, qcr, cystatin_c, alpha = 0.5) {
  f1 <- serum_creatinine / qcr
  f2 <- 1 - alpha
  f3 <- cystatin_c / 0.82
  return(107.3 / ((0.5 * f1) + (f2 * f3)))
}

egfr_zappitelli <- function(height, cystatin_c, serum_creatinine) {
  return((507.76 * exp(0.003 * height)) /
           ((cystatin_c ^ 0.635) * ((serum_creatinine * 88.4) ^ 0.547)))
}

egfr_schwartz <- function(height, serum_creatinine, cystatin_c, bun, sex) {
  m <- ifelse(sex == "M", 1, 0)
  return(39.1 * ((height / serum_creatinine) ^ 0.516) * 
           ((1.8 / cystatin_c) ^ 0.294) * 
           ((30 / bun) ^ 0.169) * 
           (1.099 ^ m) * ((height / 1.4) ^ 0.188))
}

egfr_bedside_schwartz <- function(height, serum_creatinine) {
  return((41.3 * (height / 100)) / serum_creatinine)
}

egfr_ckd_epi <- function(serum_creatinine, age, sex) {
  f <- ifelse(sex == "F", 1, 0)
  a <- ifelse(sex == "M", -0.302, -0.241)
  k <- ifelse(sex == "M", 0.9, 0.7)
  
  return(142 * 
           (pmin(serum_creatinine / k, 1) ^ a) * 
           (pmax(serum_creatinine / k, 1) ^ -1.200) * 
           (0.9938 ^ age) * (1.012 * f + (1 - f)))
}

# omit ongoing studies for now: REMODEL-T1D, T1D-DISCO and RPC2

# pull existing harmonized dataset
harm_dat <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = "")

harm_dat <- harm_dat %>%
  filter(study %in% c("PANDA", "RENAL-HEIR", "RENAL-HEIRitage", "IMPROVE",
                      "CASPER", "CROCODILE", "PENGUIN")) # omit ATTEMPT from this for now since pulling from Antoine's data

### To pull from Antoine's dataset (merged_data)
load("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Data Clean/ATTEMPT_AC.RData")
attempt <- merged_data %>%
  select(-contains("date")) %>%
  mutate(eGFR_CKD_epi = egfr_ckd_epi((creatinine_s/100), age, sex),
         study = "ATTEMPT",
         group = "Type 1 Diabetes")

## TODAY/TODAY2
# load comorbidity data (comorb)
# load("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/TODAY subaward/Clinical data/comorb.Rdata")

# load baseline risk factors (baserisk)
load("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/TODAY subaward/Clinical data/TODAY/baserisk.Rdata")

today <- baserisk %>%
  select(-race, -sex) %>%
  rename(record_id = releaseid, 
         age = AGEBASE,
         race = racedesc,
         cystatin_c_s = serumcystc,
         creatinine_s = SerumCreat,
         hba1c = HbA1c,
         sex = sex_char) %>%
  mutate(sex = case_when(sex == "M" ~ "Male", 
                         sex == "F" ~ "Female"),
         eGFR_CKD_epi = egfr_ckd_epi(creatinine_s, age, sex),
         study = "TODAY/TODAY2",
         visit = "baseline",
         group = "Type 2 Diabetes")

## TEEN-LABS
# TEENLABS
load("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Teen Labs/Data_Cleaned/analysis_dataset.RData")
df <- remove_labels(df)

teenlabs <- df %>%
  rename(record_id = ID,
         sex = SEX,
         ethnicity = ETHN,
         race = RACE,
         eGFR_fas_cr = eGFR.fas_cr,
         eGFR_fas_cr_cysc = eGFR.fas_cr_cysc,
         acr_u = UACRATIO,
         cystatin_c_s = CYSC,
         creatinine_s = CREAS) %>%
  mutate(group = case_when(diab == "Yes" ~ "Type 2 Diabetes",
                           diab == "No" ~ "Obese Control"),
         acr_u = acr_u*100,
         eGFR_CKD_epi = egfr_ckd_epi(creatinine_s, age, sex),
         race = as.character(race),
         sex = as.character(sex),
         ethnicity = as.character(ethnicity),
         study = "TEEN-LABS") %>%
  select(-(starts_with("seq")))

## T1-DISCO
# omit for now (ongoing)

## REMODEL-T1D
# omit for now (ongoing)

## RPC2
# omit for now (ongoing)

######################################

# Merge

######################################

decode <- full_join(today, teenlabs) %>%
  full_join(harm_dat) %>%
  full_join(attempt)

save(decode, file = "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/DECODE/Data Clean/decode_harmonized.RData")
