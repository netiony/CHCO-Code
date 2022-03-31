library(knitr)
library(redcapAPI)
library(tidyREDCap)
library(tidyverse)
if(Sys.info()["sysname"] == "Windows"){
  home_dir = "S:/Laura/Peds Endo/Petter Bjornstad/Data Harmonization"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = "/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization"
}
setwd(home_dir)
rm(home_dir)
# API import
tokens = read.csv("api_tokens.csv")
uri = "https://redcap.ucdenver.edu/api/"
renalheir = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "Renal-HEIR"]
  )
)
penguin = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "PENGUIN"]
  )
)
crocodile = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "CROCODILE"]
  )
)
coffee = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "COFFEE"]
  )
)
casper = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "CASPER"]
  )
)
improve = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "IMPROVE"]
  )
)
###############################################################################
# Screening variables
###############################################################################
screen_vars = c("insulin","insulin_pump","screen_height","screen_weight","screen_bmi",
                "screen_bmi_percentile","waist_circumference","hip_circumference",
                "sys_bp","dys_bp","map","pulse","activity_factor","schofield",
                "hba1c","hemoglobin","screen_hematocrit","screen_serum_creatinine",
                "screen_urine_mab","screen_urine_cre","screen_urine_acr",
                "screen_pregnant","screening_labs_complete")
# RENAL HEIR
# insulin
renalheir$insulin = renalheir$diabetes_med___2
levels(renalheir$insulin) = c("No","Yes")
# Insulin pump
renalheir$insulin_pump = NA
# activity_factor
renalheir = renalheir %>% 
  unite(.,activity_factor,fem_activity_factor,male_activity_factor,na.rm = T) %>%
  mutate(schofield = schofield_female+schofield_male)
renalheir$schofield[renalheir$schofield == 0] = NA
# PENGUIN
# Rename
penguin = penguin %>% 
  rename(screen_height = phys_height,screen_weight = phys_weight,
         screen_bmi = phys_bmi,waist_circumference = phys_waistcm,
         hip_circumference = phys_hipcm,sys_bp = phys_sysbp,dys_bp = phys_diasbp,
         map = phys_map,pulse = phys_pulse,hba1c = bl_a1c,
         screen_serum_creatinine = bl_creatinine_s,screen_urine_mab = u24_mab,
         screen_urine_cre = bl_creatinine_u,screen_urine_acr = bl_uacr,
         screen_pregnant = eligibility_preg)
# Missing
penguin[,c("insulin","insulin_pump","screen_bmi_percentile",
           "activity_factor","schofield","hemoglobin","screen_hematocrit")] = NA
# CROCODILE
# insulin_pump
crocodile$insulin_pump = crocodile$diabetes_tx___1
levels(crocodile$insulin_pump) = c("No","Yes")
# insulin
crocodile$insulin = "No"
crocodile$insulin[which(crocodile$diabetes_tx___1 == "Checked"|
                          crocodile$diabetes_tx___2 == "Checked")] = "Yes"
crocodile$insulin = factor(crocodile$insulin,levels = c("No","Yes"))
# Rename
crocodile = crocodile %>% 
  rename(screen_height = phys_height,screen_weight = phys_weight,
         screen_bmi = phys_bmi,waist_circumference = phys_waistcm,
         hip_circumference = phys_hipcm,sys_bp = phys_sysbp,dys_bp = phys_diasbp,
         map = phys_map,pulse = phys_pulse,hba1c = bl_a1c,
         screen_serum_creatinine = bl_creatinine_s,screen_urine_mab = u24_mab,
         screen_urine_cre = bl_creatinine_u,screen_urine_acr = bl_uacr,
         screen_pregnant = eligibility_preg)


