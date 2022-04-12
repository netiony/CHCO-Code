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
# API import
tokens = read.csv("api_tokens.csv")
uri = "https://redcap.ucdenver.edu/api/"
renalheir = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "Renal-HEIR"]
  ),labels = F
)
renalheir$study = "RENAL-HEIR"

penguin = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "PENGUIN"]
  ),labels = F
)
penguin$study = "PENGUIN"

crocodile = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "CROCODILE"]
  ),labels = F
)
crocodile$study = "CROCODILE"

coffee = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "COFFEE"]
  ),labels = F
)
coffee$study = "COFFEE"

casper = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "CASPER"]
  ),labels = F
)
casper$study = "CASPER"

improve = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "IMPROVE"]
  ),labels = F
)
improve$study = "IMPROVE"
rm(home_dir,tokens,uri)
###############################################################################
# Demographic variables
###############################################################################
demographic_vars = c("subject_id","study","group","dob","age_current","gender",
                     "race","ethnicity","length_of_diabetes","age_at_diabetes_dx")
# RENAL HEIR
levels(renalheir$group) = c("T2D","Obese Control","Lean Control")
renalheir$age_at_diabetes_dx = renalheir$diabetes_age
# Race
races = c("American Indian or Alaskan Native","Asian",
          "Hawaiian or Pacific Islander","Black or African American",
          "White","Unknown","Other")
renalheir$race = apply(renalheir,1,function(r){
  race = r[paste0("race___",1:7)]
  w = which(race == "Checked")
  return(paste0(sort(races[w]),collapse = "/"))
})
# Ethnicity
eths = c("Hispanic","Non-Hispanic","Unknown/Not Reported")
renalheir$ethnicity = apply(renalheir,1,function(r){
  eth = r[paste0("ethnicity___",1:3)]
  w = which(eth == "Checked")
  return(paste0(eths[w],collapse = "/"))
})
# PENGUIN
penguin$subject_id = penguin$record_id
penguin$group = "PKD"
penguin$age_current = penguin$age_consent 
penguin$gender = penguin$sex
# Race
races = c("American Indian or Alaskan Native","Asian",
          "Hawaiian or Pacific Islander","Black or African American",
          "White","Other","Unknown")
penguin$race = apply(penguin,1,function(r){
  race = r[paste0("race___",1:7)]
  w = which(race == "Checked")
  return(paste0(sort(races[w]),collapse = "/"))
})
# Ethnicity
penguin$ethnicity = apply(penguin,1,function(r){
  eth = r[paste0("ethnicity___",1:3)]
  w = which(eth == "Checked")
  return(paste0(eths[w],collapse = "/"))
})
# Diabetes info
penguin[,c("length_of_diabetes","age_at_diabetes_dx")] = NA
# CROCODILE
crocodile$subject_id = crocodile$record_id
levels(crocodile$group) = c("T1D","Lean Control")
crocodile$age_current = crocodile$age_consent 
crocodile$gender = crocodile$sex
# Race
crocodile$race = apply(crocodile,1,function(r){
  race = r[paste0("race___",1:7)]
  w = which(race == "Checked")
  return(paste0(sort(races[w]),collapse = "/"))
})
# Ethnicity
crocodile$ethnicity = apply(crocodile,1,function(r){
  eth = r[paste0("ethnicity___",1:3)]
  w = which(eth == "Checked")
  return(paste0(eths[w],collapse = "/"))
})
# Diabetes info
crocodile$length_of_diabetes = crocodile$diabetes_duration
crocodile$age_at_diabetes_dx = crocodile$diabetes_dx_age
# COFFEE
coffee$group = "T1D"
# Race
races = c("American Indian or Alaskan Native","Asian",
          "Hawaiian or Pacific Islander","Black or African American",
          "White","Unknown","Other")
coffee$race = apply(coffee,1,function(r){
  race = r[paste0("race___",1:7)]
  w = which(race == "Checked")
  return(paste0(sort(races[w]),collapse = "/"))
})
# Ethnicity
coffee$ethnicity = apply(coffee,1,function(r){
  eth = r[paste0("ethnicity___",1:3)]
  w = which(eth == "Checked")
  return(paste0(eths[w],collapse = "/"))
})
# Diabetes info
coffee$age_at_diabetes_dx = coffee$diabetes_age
# CASPER
casper$group = "T1D"
# Race
casper$race = apply(casper,1,function(r){
  race = r[paste0("race___",1:7)]
  w = which(race == "Checked")
  return(paste0(sort(races[w]),collapse = "/"))
})
# Ethnicity
casper$ethnicity = apply(casper,1,function(r){
  eth = r[paste0("ethnicity___",1:3)]
  w = which(eth == "Checked")
  return(paste0(eths[w],collapse = "/"))
})
# Diabetes info
casper$age_at_diabetes_dx = casper$diabetes_age
# IMPROVE
improve$group = "T2D"
# Race
improve$race = apply(improve,1,function(r){
  race = r[paste0("race___",1:7)]
  w = which(race == "Checked")
  return(paste0(sort(races[w]),collapse = "/"))
})
# Ethnicity
improve$ethnicity = apply(improve,1,function(r){
  eth = r[paste0("ethnicity___",1:3)]
  w = which(eth == "Checked")
  return(paste0(eths[w],collapse = "/"))
})
# Diabetes info
improve$age_at_diabetes_dx = improve$diabetes_age
# Merge demographics
demographics = do.call(rbind,list(renalheir[,demographic_vars],penguin[,demographic_vars],
                                  crocodile[,demographic_vars],coffee[,demographic_vars],
                                  casper[,demographic_vars],improve[,demographic_vars]))
###############################################################################
# Screening variables
###############################################################################
screen_vars = c("subject_id","study","insulin","insulin_pump","screen_height",
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
  rename(screen_height = phys_height,screen_weight = phys_weight,
         screen_bmi = phys_bmi,waist_circumference = phys_waistcm,
         hip_circumference = phys_hipcm,sys_bp = phys_sysbp,dys_bp = phys_diasbp,
         map = phys_map,pulse = phys_pulse,hba1c = bl_a1c,
         screen_serum_creatinine = bl_creatinine_s,
         screen_urine_cre = bl_creatinine_u,screen_urine_acr = bl_uacr,
         screen_pregnant = eligibility_preg)
penguin$screen_urine_mab = NA
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
  rename(screen_height = phys_height,screen_weight = phys_weight,
         screen_bmi = phys_bmi,waist_circumference = phys_waistcm,
         hip_circumference = phys_hipcm,sys_bp = phys_sysbp,dys_bp = phys_diasbp,
         map = phys_map,pulse = phys_pulse,hba1c = bl_a1c,
         screen_serum_creatinine = bl_creatinine_s,
         screen_urine_cre = bl_creatinine_u,screen_urine_acr = bl_uacr,
         screen_pregnant = screen_upt)
crocodile$screen_urine_mab = NA
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
###############################################################################
# DXA variables
###############################################################################

dxa_vars = c("subject_id","study","dexa_date","body_fat","lean_mass",
             "trunk_mass","fat_kg","lean_kg","trunk_kg","bone_mineral_density",
             "body_composition_dxa_complete")

# RENAL-HEIR
# PENGUIN
penguin = penguin %>% 
  rename(dexa_date = dxa_date,body_fat = bodyfat_percent,lean_kg = leanmass_kg,
         fat_kg = fatmass_kg,lean_mass = leanmass_percent,trunk_mass = trunkmass_percent,
         trunk_kg = trunkmass_kg,bone_mineral_density = bmd,
         body_composition_dxa_complete = study_visit_dxa_scan_complete)
# CROCODILE
crocodile = crocodile %>%
  rename(dexa_date = dxa_date,body_fat = bodyfat_percent,lean_kg = leanmass_kg,
         fat_kg = fatmass_kg,lean_mass = leanmass_percent,trunk_mass = trunkmass_percent,
         trunk_kg = trunkmass_kg,bone_mineral_density = bmd,
         body_composition_dxa_complete = study_visit_dxa_scan_complete)
# COFFEE
coffee[,dxa_vars[3:length(dxa_vars)]] = NA
# CASPER
# IMPROVE
improve = improve %>%
  rename(dexa_date = bodcomp_date,body_fat = dxa_body_fat,lean_kg = dxa_lean_kg,
         fat_kg = dxa_fat_kg,lean_mass = dxa_lean_mass,trunk_mass = dxa_trunk_mass,
         trunk_kg = dxa_trunk_kg,bone_mineral_density = dxa_bmd,
         body_composition_dxa_complete = dxa_complete)
# Merge DXA
dxa = do.call(rbind,list(renalheir[,dxa_vars],penguin[,dxa_vars],
                               crocodile[,dxa_vars],coffee[,dxa_vars],
                               casper[,dxa_vars],improve[,dxa_vars]))
###############################################################################
# Clamp vitals
###############################################################################

clamp_vitals = c("subject_id","study","clamp_date","clamp_height",
                 "clamp_weight","clamp_sbp","clamp_dbp","clamp_map","clamp_pls")

# RENAL-HEIR - already correct
# PENGUIN
penguin = penguin %>%
  rename(clamp_height = clamp_ht,clamp_weight = clamp_wt)
# CROCODILE
crocodile = crocodile %>% 
  rename(clamp_height = clamp_ht,clamp_weight = clamp_wt)
# COFFEE
coffee$clamp_date = coffee$cf_clamp_date
# CASPER - already correct
# IMPROVE - already correct

# Merge 
# clamp_vitals = do.call(rbind,list(renalheir[,clamp_vitals],penguin[,clamp_vitals],
#                          crocodile[,clamp_vitals],coffee[,clamp_vitals],
#                          casper[,clamp_vitals],improve[,clamp_vitals]))

###############################################################################
# Clamp labs
###############################################################################

clamp_labs = c("subject_id","study","cholesterol","hdl","ldl","triglycerides",
               "total_protein","serum_sodium","cystatin_c","serum_creatinine",
               "clamp_urine_mab_baseline","clamp_urine_cre_baseline",
               "clamp_acr_baseline","clamp_urine_sodium","clamp_glucose_bl",
               "urine_glucose","clamp_urine_mab_250","clamp_urine_cre_250",
               "clamp_acr_250","clamp_urine_vol")

# RENAL-HEIR - already correct
# PENGUIN

# CROCODILE

# COFFEE
coffee[c("cholesterol","hdl","ldl","triglycerides")] = NA
coffee$clamp_glucose_bl = coffee$fbg
# CASPER - already correct
# IMPROVE- already correct

###############################################################################
# 24 hour urine labs
###############################################################################

# Only available for CROCODILE, PENGUIN, and PANTHER
urine_labs = c("subject_id","study","u24_labs","u24_na","u24_mab","u24_volume","u24_hours")

# RENAL-HEIR
renalheir[,tail(urine_labs,-2)] = NA
# PENGUIN
penguin$u24_hours = NA
penguin$u24_volume = NA
# CROCODILE
crocodile$u24_hours = NA
crocodile$u24_volume = NA
# COFFEE
coffee[,tail(urine_labs,-2)] = NA
# CASPER
casper[,tail(urine_labs,-2)] = NA
# IMPROVE
improve[,tail(urine_labs,-2)] = NA

# Merge 
urine_labs = do.call(rbind,list(renalheir[,urine_labs],penguin[,urine_labs],
                         crocodile[,urine_labs],coffee[,urine_labs],
                         casper[,urine_labs],improve[,urine_labs]))

###############################################################################
# Hyperinsulinemic-euglycemic clamp data
###############################################################################

# Only available for CROCODILE and PENGUIN							
he_clamp = c("subject_id","study","p1_raw_m","p1_raw_leanm","p1_gc_m",
             "p1_gc_leanm","p2_raw_m","p2_raw_leanm","p2_gc_m","p2_gc_leanm")

# RENAL-HEIR
renalheir[,tail(he_clamp,-2)] = NA
# COFFEE
coffee[,tail(he_clamp,-2)] = NA
# CASPER
casper[,tail(he_clamp,-2)] = NA
# IMPROVE
improve[,tail(he_clamp,-2)] = NA

# Merge 
he_clamp = do.call(rbind,list(renalheir[,he_clamp],penguin[,he_clamp],
                                crocodile[,he_clamp],coffee[,he_clamp],
                                casper[,he_clamp],improve[,he_clamp]))

###############################################################################
# FFA data from clamps
###############################################################################

# RENAL-HEIR
renalheir_ffa = renalheir[,c("subject_id","study",colnames(renalheir)[grep("ffa_.*\\d{1,}",colnames(renalheir))])]
# PENGUIN
penguin_ffa = penguin[,c("subject_id","study",colnames(penguin)[grep("ffa_.*\\d{1,}",colnames(penguin))])]
colnames(penguin_ffa) = sub("minus","minus_",colnames(penguin_ffa))
# CROCODILE
crocodile_ffa = crocodile[,c("subject_id","study",colnames(crocodile)[grep("ffa_.*\\d{1,}",colnames(crocodile))])]
colnames(crocodile_ffa) = sub("minus","minus_",colnames(crocodile_ffa))
# COFFEE - none
# CASPER - none

# Merge
ffa = full_join(renalheir_ffa,penguin_ffa)
ffa = full_join(ffa,crocodile_ffa)
# Sort columns
ffa_cols = colnames(ffa)[3:ncol(ffa)]
ffa_cols = ffa_cols[order(as.numeric(sub("ffa_","",sub("minus_","-",ffa_cols))))]
ffa = ffa[,c("subject_id","study",ffa_cols)]

## IMPROVE
improve_ffa = improve[,c("subject_id","study",colnames(improve)[grep("ffa_.*\\d{1,}",colnames(improve))])]
colnames(improve_ffa) = sub("neg","minus",colnames(improve_ffa))
colnames(improve_ffa)[grep("ffa",colnames(improve_ffa))] = 
  paste0(colnames(improve_ffa)[grep("ffa",colnames(improve_ffa))],"_mmtt")
# Add IMPROVE
ffa = full_join(ffa,improve_ffa)

###############################################################################
# C-peptide clamp data
###############################################################################

# RENAL-HEIR
renalheir_cpep = renalheir[,c("subject_id","study",colnames(renalheir)[grep("cpeptide_.*\\d{1,}",colnames(renalheir))])]
# IMPROVE
improve_cpep = improve[,c("subject_id","study",colnames(improve)[grep("cpep_.*\\d{1,}",colnames(improve))])]
colnames(improve_cpep) = sub("mmtt_","",colnames(improve_cpep))
colnames(improve_cpep)[3:ncol(improve_cpep)] = 
  paste0(colnames(improve_cpep)[3:ncol(improve_cpep)],"_mmtt")

# Merge
cpep = full_join(renalheir_cpep,improve_cpep)

###############################################################################
# Insulin clamp data
###############################################################################



# RENAL-HEIR

# PENGUIN

# CROCODILE

# COFFEE

# CASPER

# IMPROVE

###############################################################################
# BG clamp data
###############################################################################



# RENAL-HEIR

# PENGUIN

# CROCODILE

# COFFEE

# CASPER

# IMPROVE

###############################################################################
# eGFR
###############################################################################



# RENAL-HEIR

# PENGUIN

# CROCODILE

# COFFEE

# CASPER

# IMPROVE

###############################################################################
# Kidney function
###############################################################################



# RENAL-HEIR

# PENGUIN

# CROCODILE

# COFFEE

# CASPER

# IMPROVE

###############################################################################
# Intraglomerular Hemodynamics
###############################################################################



# RENAL-HEIR

# PENGUIN

# CROCODILE

# COFFEE

# CASPER

# IMPROVE

###############################################################################
# Kidney MRI
###############################################################################



# RENAL-HEIR

# PENGUIN

# CROCODILE

# COFFEE

# CASPER

# IMPROVE

###############################################################################
# PET/CT
###############################################################################



# RENAL-HEIR

# PENGUIN

# CROCODILE

# COFFEE

# CASPER

# IMPROVE

###############################################################################
# Kidney biopsy
###############################################################################



# RENAL-HEIR

# PENGUIN

# CROCODILE

# COFFEE

# CASPER

# IMPROVE

###############################################################################
# Merge everything together
###############################################################################

df = full_join(screening,demographics)
df = full_join(df,dxa)

###############################################################################
# Final formatting and calculations
###############################################################################



