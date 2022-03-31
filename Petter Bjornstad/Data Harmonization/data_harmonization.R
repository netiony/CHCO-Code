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
renalheir$Study = "RENAL-HEIR"
penguin = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "PENGUIN"]
  )
)
penguin$Study = "PENGUIN"
crocodile = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "CROCODILE"]
  )
)
crocodile$Study = "CROCODILE"
coffee = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "COFFEE"]
  )
)
coffee$Study = "COFFEE"
casper = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "CASPER"]
  )
)
casper$Study = "CASPER"
improve = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "IMPROVE"]
  )
)
improve$Study = "IMPROVE"
###############################################################################
# Demographic variables
###############################################################################



###############################################################################
# Screening variables
###############################################################################
screen_vars = c("subject_id","Study","insulin","insulin_pump","screen_height",
                "screen_weight","screen_bmi","screen_bmi_percentile",
                "waist_circumference","hip_circumference","sys_bp","dys_bp",
                "map","pulse","activity_factor","schofield","hba1c","hemoglobin",
                "screen_hematocrit","screen_serum_creatinine",
                "screen_urine_mab","screen_urine_cre","screen_urine_acr",
                "screen_pregnant","screening_labs_complete")
# RENAL HEIR
# insulin
renalheir$insulin = renalheir$diabetes_med___2
levels(renalheir$insulin) = c("No","Yes")
# Insulin pump
renalheir$insulin_pump = NA
# Pregnancy
levels(renalheir$screen_pregnant) = c("Not Pregnant","Pregnant")
# activity_factor and schofield
renalheir = renalheir %>% 
  unite(.,activity_factor,fem_activity_factor,male_activity_factor,na.rm = T) %>%
  mutate(schofield = schofield_female+schofield_male)
renalheir$schofield[renalheir$schofield == 0] = NA
# PENGUIN
# Rename
penguin = penguin %>% 
  rename(subject_id = record_id,screen_height = phys_height,screen_weight = phys_weight,
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
# insulin
crocodile$insulin = "No"
crocodile$insulin[which(crocodile$diabetes_tx___1 == "Checked"|
                          crocodile$diabetes_tx___2 == "Checked")] = "Yes"
crocodile$insulin = factor(crocodile$insulin,levels = c("No","Yes"))
# insulin_pump
crocodile$insulin_pump = crocodile$diabetes_tx___1
levels(crocodile$insulin_pump) = c("No","Yes")
# Rename
crocodile = crocodile %>% 
  rename(subject_id = record_id,screen_height = phys_height,screen_weight = phys_weight,
         screen_bmi = phys_bmi,waist_circumference = phys_waistcm,
         hip_circumference = phys_hipcm,sys_bp = phys_sysbp,dys_bp = phys_diasbp,
         map = phys_map,pulse = phys_pulse,hba1c = bl_a1c,
         screen_serum_creatinine = bl_creatinine_s,screen_urine_mab = u24_mab,
         screen_urine_cre = bl_creatinine_u,screen_urine_acr = bl_uacr,
         screen_pregnant = screen_upt)
# Missing
crocodile[,screen_vars[which(!screen_vars %in% colnames(crocodile))]] = NA
# COFFEE
# insulin
coffee$insulin = "No"
coffee$insulin[which(coffee$diabetes_med___1 == "Checked"|
                          coffee$diabetes_med___2 == "Checked")] = "Yes"
coffee$insulin = factor(coffee$insulin,levels = c("No","Yes"))
# insulin_pump
coffee$insulin_pump = coffee$diabetes_med___1
levels(coffee$insulin_pump) = c("No","Yes")
# activity_factor and schofield
coffee = coffee %>% 
  unite(.,activity_factor,fem_activity_factor,male_activity_factor,na.rm = T) %>%
  mutate(schofield = rowSums(.[c("schofield_female","schofield_male")],na.rm = T))
coffee$schofield[coffee$schofield == 0] = NA
# Other
coffee$screening_labs_complete = coffee$screening_labs_casper_complete
# CASPER
# insulin
casper$insulin = "No"
casper$insulin[which(casper$diabetes_med___1 == "Checked"|
                       casper$diabetes_med___2 == "Checked")] = "Yes"
casper$insulin = factor(casper$insulin,levels = c("No","Yes"))
# insulin_pump
casper$insulin_pump = casper$diabetes_med___1
levels(casper$insulin_pump) = c("No","Yes")
# activity_factor and schofield
casper = casper %>% 
  unite(.,activity_factor,fem_activity_factor,male_activity_factor,na.rm = T) %>%
  mutate(schofield = rowSums(.[c("schofield_female","schofield_male")],na.rm = T))
casper$schofield[casper$schofield == 0] = NA
# IMPROVE
# insulin
improve$insulin = improve$diabetes_med___2
improve$insulin = factor(improve$insulin,levels = c("No","Yes"))
# insulin_pump
improve$insulin_pump = NA
# activity_factor and schofield
improve = improve %>% 
  unite(.,activity_factor,activity_factor_female,activity_factor_male,na.rm = T) %>%
  mutate(schofield = rowSums(.[c("schofield_female","schofield_male")],na.rm = T))
improve$schofield[improve$schofield == 0] = NA
# Merge screening
screening = do.call(rbind,list(renalheir[,screen_vars],penguin[,screen_vars],
                               crocodile[,screen_vars],coffee[,screen_vars],
                               casper[,screen_vars],improve[,screen_vars]))




