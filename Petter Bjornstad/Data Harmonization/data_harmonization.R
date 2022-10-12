library(knitr)
library(redcapAPI)
library(childsds)
library(tidyverse)
library(parsedate)
source("~/GitHub/shared-resources/Data Cleaning/Calculated Variables/eGFR.R")
source("~/GitHub/shared-resources/Data Cleaning/Calculated Variables/hemodynamics.R")
setwd("~/Dropbox/Work/CHCO/Petter Bjornstad/Data Harmonization")
# API import
tokens <- read.csv("api_tokens.csv")
uri <- "https://redcap.ucdenver.edu/api/"
renalheir <- exportRecords(
  redcapConnection(
    url = uri, token = tokens$Token[tokens$Study == "Renal-HEIR"]
  ),
  labels = F
)
renalheir$study <- "RENAL-HEIR"
renalheir$visit <- "Baseline"
renalheir[renalheir == -99] <- NA
renalheir[renalheir == -999] <- NA
renalheir[renalheir == -9999] <- NA

penguin <- exportRecords(
  redcapConnection(
    url = uri, token = tokens$Token[tokens$Study == "PENGUIN"]
  ),
  labels = F
)
penguin$study <- "PENGUIN"
penguin$visit <- "Baseline"
penguin[penguin == -99] <- NA
penguin[penguin == -999] <- NA
penguin[penguin == -9999] <- NA

crocodile <- exportRecords(
  redcapConnection(
    url = uri, token = tokens$Token[tokens$Study == "CROCODILE"]
  ),
  labels = F
)
crocodile$study <- "CROCODILE"
crocodile$visit <- "Baseline"
crocodile[crocodile == -99] <- NA
crocodile[crocodile == -999] <- NA
crocodile[crocodile == -9999] <- NA

coffee <- exportRecords(
  redcapConnection(
    url = uri, token = tokens$Token[tokens$Study == "COFFEE"]
  ),
  labels = F
)
coffee$study <- "COFFEE"
coffee$visit <- "Baseline"
coffee[coffee == -99] <- NA
coffee[coffee == -999] <- NA
coffee[coffee == -9999] <- NA

casper <- exportRecords(
  redcapConnection(
    url = uri, token = tokens$Token[tokens$Study == "CASPER"]
  ),
  labels = F
)
casper$study <- "CASPER"
casper$visit <- "Baseline"
casper[casper == -99] <- NA
casper[casper == -999] <- NA
casper[casper == -9999] <- NA

improve <- exportRecords(
  redcapConnection(
    url = uri, token = tokens$Token[tokens$Study == "IMPROVE"]
  ),
  labels = F
)
improve$study <- "IMPROVE"
improve$visit=as.character(improve$study_visit)
improve$visit[is.na(improve$visit)] = "Baseline"
improve[improve == -99] <- NA
improve[improve == -999] <- NA
improve[improve == -9999] <- NA


###############################################################################
# Demographic variables
###############################################################################
demographic_vars <- c(
  "subject_id", "study", "co_enroll", "co_enroll_id", "visit",
  "group", "dob", "diagnosis_date", "gender", "race", "ethnicity",
  "age_consent"
)
# RENAL HEIR
levels(renalheir$group) <- c("T2D", "Obese Control", "Lean Control")
renalheir$age_at_diabetes_dx <- renalheir$diabetes_age
renalheir$diagnosis_date <- renalheir$diagnosis
renalheir$age_consent = renalheir$age_current
renalheir$co_enroll_id = sub("IT2D","IT",renalheir$co_enroll_id)
# Race
races <- c(
  "American Indian or Alaskan Native", "Asian",
  "Hawaiian or Pacific Islander", "Black or African American",
  "White", "Unknown", "Other"
)
renalheir$race <- apply(renalheir, 1, function(r) {
  race <- r[paste0("race___", 1:7)]
  w <- which(race == "Checked")
  return(paste0(sort(races[w]), collapse = "/"))
})
# Ethnicity
eths <- c("Hispanic", "Non-Hispanic", "Unknown/Not Reported")
renalheir$ethnicity <- apply(renalheir, 1, function(r) {
  eth <- r[paste0("ethnicity___", 1:3)]
  w <- which(eth == "Checked")
  return(paste0(eths[w], collapse = "/"))
})
# PENGUIN
penguin$subject_id <- penguin$record_id
penguin$group <- "PKD"
penguin$age_current <- penguin$age_consent
penguin$gender <- penguin$sex
# Race
races <- c(
  "American Indian or Alaskan Native", "Asian",
  "Hawaiian or Pacific Islander", "Black or African American",
  "White", "Other", "Unknown"
)
penguin$race <- apply(penguin, 1, function(r) {
  race <- r[paste0("race___", 1:7)]
  w <- which(race == "Checked")
  return(paste0(sort(races[w]), collapse = "/"))
})
# Ethnicity
penguin$ethnicity <- apply(penguin, 1, function(r) {
  eth <- r[paste0("ethnicity___", 1:3)]
  w <- which(eth == "Checked")
  return(paste0(eths[w], collapse = "/"))
})
# Missing columns
penguin[, c("diagnosis_date", "co_enroll", "co_enroll_id")] <- NA
# CROCODILE
crocodile$subject_id <- crocodile$record_id
levels(crocodile$group) <- c("T1D", "Lean Control")
crocodile$age_current <- crocodile$age_consent
crocodile$gender <- crocodile$sex
crocodile$diagnosis_date <- crocodile$diabetes_dx_date
# Race
crocodile$race <- apply(crocodile, 1, function(r) {
  race <- r[paste0("race___", 1:7)]
  w <- which(race == "Checked")
  return(paste0(sort(races[w]), collapse = "/"))
})
# Ethnicity
crocodile$ethnicity <- apply(crocodile, 1, function(r) {
  eth <- r[paste0("ethnicity___", 1:3)]
  w <- which(eth == "Checked")
  return(paste0(eths[w], collapse = "/"))
})
# Missing columns
crocodile[, c("co_enroll", "co_enroll_id")] <- NA
# COFFEE
coffee$group <- "T1D"
coffee$age_consent = coffee$age_current
# Race
races <- c(
  "American Indian or Alaskan Native", "Asian",
  "Hawaiian or Pacific Islander", "Black or African American",
  "White", "Unknown", "Other"
)
coffee$race <- apply(coffee, 1, function(r) {
  race <- r[paste0("race___", 1:7)]
  w <- which(race == "Checked")
  return(paste0(sort(races[w]), collapse = "/"))
})
# Ethnicity
coffee$ethnicity <- apply(coffee, 1, function(r) {
  eth <- r[paste0("ethnicity___", 1:3)]
  w <- which(eth == "Checked")
  return(paste0(eths[w], collapse = "/"))
})
# Diabetes info
coffee$diagnosis_date <- coffee$diagnosis
coffee[, c("co_enroll", "co_enroll_id")] <- NA
# CASPER
casper$group <- "T1D"
casper$age_consent = casper$age_current
# Race
casper$race <- apply(casper, 1, function(r) {
  race <- r[paste0("race___", 1:7)]
  w <- which(race == "Checked")
  return(paste0(sort(races[w]), collapse = "/"))
})
# Ethnicity
casper$ethnicity <- apply(casper, 1, function(r) {
  eth <- r[paste0("ethnicity___", 1:3)]
  w <- which(eth == "Checked")
  return(paste0(eths[w], collapse = "/"))
})
# Diabetes info
casper$diagnosis_date <- casper$diagnosis
casper[, c("co_enroll", "co_enroll_id")] <- NA
# IMPROVE
improve$group <- "T2D"
improve$age_consent = improve$age_current
# Race
improve$race <- apply(improve, 1, function(r) {
  race <- r[paste0("race___", 1:7)]
  w <- which(race == "Checked")
  return(paste0(sort(races[w]), collapse = "/"))
})
# Ethnicity
improve$ethnicity <- apply(improve, 1, function(r) {
  eth <- r[paste0("ethnicity___", 1:3)]
  w <- which(eth == "Checked")
  return(paste0(eths[w], collapse = "/"))
})
# Diabetes info
improve$diagnosis_date <- improve$diagnosis
# Just first visit
improve_dem <- improve[grep("screening_arm", improve$redcap_event_name), ]
# Merge demographics
demographics <- do.call(rbind, list(
  renalheir[, demographic_vars], penguin[, demographic_vars],
  crocodile[, demographic_vars], coffee[, demographic_vars],
  casper[, demographic_vars], improve_dem[, demographic_vars]
))

###############################################################################
# Screening variables
###############################################################################
screen_vars <- c(
  "subject_id", "study", "visit", "date_of_screen",
  "screen_height", "screen_weight", "screen_bmi",
  "waist_circumference", "hip_circumference", "sys_bp", "dys_bp",
  "map", "pulse", "activity_factor", "schofield", "hba1c", "hemoglobin",
  "screen_hematocrit", "screen_serum_creatinine",
  "screen_urine_mab", "screen_urine_cre", "screen_urine_acr",
  "screen_pregnant", "screening_labs_complete"
)

# RENAL HEIR
# Diabetes medications
renalheir$insulin <- factor(renalheir$diabetes_med___2,
  levels = c("Unchecked", "Checked"),
  labels = c("No", "Yes")
)
renalheir$insulin_pump <- NA
renalheir$metformin <- factor(renalheir$diabetes_med___1,
  levels = c("Unchecked", "Checked"),
  labels = c("No", "Yes")
)
renalheir$short_act_insulin <- factor(renalheir$insulin_type___1,
  levels = c("Unchecked", "Checked"),
  labels = c("No", "Yes")
)
renalheir$long_act_insulin <- factor(renalheir$insulin_type___2,
  levels = c("Unchecked", "Checked"),
  labels = c("No", "Yes")
)
renalheir$tzds <- factor(renalheir$diabetes_med_other___1,
  levels = c("Unchecked", "Checked"),
  labels = c("No", "Yes")
)
renalheir$glp1_agonists <- factor(renalheir$diabetes_med_other___2,
  levels = c("Unchecked", "Checked"),
  labels = c("No", "Yes")
)
renalheir$sglt2_inhibitors <- factor(renalheir$diabetes_med_other___3,
  levels = c("Unchecked", "Checked"),
  labels = c("No", "Yes")
)
renalheir$other_diabetes_meds <- factor(renalheir$diabetes_med_other___4,
  levels = c("Unchecked", "Checked"),
  labels = c("No", "Yes")
)
# Pregnancy
levels(renalheir$screen_pregnant) <- c("Not Pregnant", "Pregnant")
# activity_factor and schofield
renalheir <- renalheir %>%
  unite(., activity_factor, fem_activity_factor, male_activity_factor, na.rm = T) %>%
  mutate(schofield = schofield_female + schofield_male)
renalheir$schofield[renalheir$schofield == 0] <- NA
# PENGUIN
# Rename
penguin <- penguin %>%
  rename(
    screen_height = phys_height, screen_weight = phys_weight,
    screen_bmi = phys_bmi, waist_circumference = phys_waistcm,
    hip_circumference = phys_hipcm, sys_bp = phys_sysbp, dys_bp = phys_diasbp,
    map = phys_map, pulse = phys_pulse, hba1c = bl_a1c,
    screen_serum_creatinine = bl_creatinine_s,
    screen_urine_cre = bl_creatinine_u, screen_urine_acr = bl_uacr,
    screen_pregnant = eligibility_preg, date_of_screen = labs_date
  )
penguin$screen_urine_mab <- NA
# Missing
penguin[, c(
  "insulin", "insulin_pump", "metformin", "short_act_insulin",
  "long_act_insulin", "tzds", "glp1_agonists", "sglt2_inhibitors",
  "other_diabetes_meds", "activity_factor", "schofield", "hemoglobin",
  "screen_hematocrit"
)] <- NA
# CROCODILE
# Diabetes medications
crocodile$insulin <- "No"
crocodile$insulin[which(crocodile$diabetes_tx___1 == "Checked" |
  crocodile$diabetes_tx___2 == "Checked")] <- "Yes"
crocodile$insulin <- factor(crocodile$insulin, levels = c("No", "Yes"))
crocodile$insulin_pump <- factor(crocodile$diabetes_tx___1,
                              levels = c("Unchecked", "Checked"),
                              labels = c("No", "Yes")
)
crocodile$short_act_insulin <- "No"
crocodile$short_act_insulin[which(rowSums(crocodile[,paste0("injections_short_acting___",1:3)]=="Checked")>0)] <- "Yes"
crocodile$long_act_insulin <- "No"
crocodile$long_act_insulin[which(rowSums(crocodile[,paste0("injections_long_acting___",1:6)]=="Checked")>0)] <- "Yes"
crocodile$metformin <- factor(crocodile$diabetes_meds_other___1,
                              levels = c("Unchecked", "Checked"),
                              labels = c("No", "Yes")
)
crocodile$tzds <- factor(crocodile$diabetes_meds_other___2,
                         levels = c("Unchecked", "Checked"),
                         labels = c("No", "Yes")
)
crocodile$glp1_agonists <- factor(crocodile$diabetes_meds_other___3,
                                  levels = c("Unchecked", "Checked"),
                                  labels = c("No", "Yes")
)
crocodile$sglt2_inhibitors <- factor(crocodile$diabetes_meds_other___4,
                                     levels = c("Unchecked", "Checked"),
                                     labels = c("No", "Yes")
)
crocodile$other_diabetes_meds <- factor(crocodile$diabetes_meds_other___5,
                                        levels = c("Unchecked", "Checked"),
                                        labels = c("No", "Yes")
)
# Rename
crocodile <- crocodile %>%
  rename(
    screen_height = phys_height, screen_weight = phys_weight,
    screen_bmi = phys_bmi, waist_circumference = phys_waistcm,
    hip_circumference = phys_hipcm, sys_bp = phys_sysbp, dys_bp = phys_diasbp,
    map = phys_map, pulse = phys_pulse, hba1c = bl_a1c,
    screen_serum_creatinine = bl_creatinine_s,
    screen_urine_cre = bl_creatinine_u, screen_urine_acr = screen_uacr,
    screen_pregnant = screen_upt, date_of_screen = labs_date
  )
crocodile$screen_urine_mab <- NA
# Missing
crocodile[, screen_vars[which(!screen_vars %in% colnames(crocodile))]] <- NA
# COFFEE
# Diabetes medications
coffee$insulin <- "No"
coffee$insulin[which(coffee$diabetes_tx___1 == "Checked" |
                          coffee$diabetes_tx___2 == "Checked")] <- "Yes"
coffee$insulin <- factor(coffee$insulin, levels = c("No", "Yes"))
coffee$insulin_pump <- factor(coffee$diabetes_med___1,
                                 levels = c("Unchecked", "Checked"),
                                 labels = c("No", "Yes")
)
coffee$short_act_insulin <- factor(coffee$insulin_type___1,
                           levels = c("Unchecked", "Checked"),
                           labels = c("No", "Yes")
)
coffee$long_act_insulin <- factor(coffee$insulin_type___2,
                                   levels = c("Unchecked", "Checked"),
                                   labels = c("No", "Yes")
)
coffee$metformin <- factor(coffee$diabetes_med_other___1,
                              levels = c("Unchecked", "Checked"),
                              labels = c("No", "Yes")
)
coffee$tzds <- factor(coffee$diabetes_med_other___2,
                         levels = c("Unchecked", "Checked"),
                         labels = c("No", "Yes")
)
coffee$glp1_agonists <- factor(coffee$diabetes_med_other___3,
                                  levels = c("Unchecked", "Checked"),
                                  labels = c("No", "Yes")
)
coffee$sglt2_inhibitors <- factor(coffee$diabetes_med_other___4,
                                     levels = c("Unchecked", "Checked"),
                                     labels = c("No", "Yes")
)
coffee$other_diabetes_meds <- factor(coffee$diabetes_med_other___5,
                                        levels = c("Unchecked", "Checked"),
                                        labels = c("No", "Yes")
)
# activity_factor and schofield
coffee <- coffee %>%
  unite(., activity_factor, fem_activity_factor, male_activity_factor, na.rm = T) %>%
  mutate(schofield = rowSums(.[c("schofield_female", "schofield_male")], na.rm = T))
coffee$schofield[coffee$schofield == 0] <- NA
# Other
coffee$screening_labs_complete <- coffee$screening_labs_casper_complete

# CASPER
# Diabetes medications
casper$insulin <- "No"
casper$insulin[which(casper$diabetes_med___1 == "Checked" |
                       casper$diabetes_med___2 == "Checked")] <- "Yes"
casper$insulin <- factor(casper$insulin, levels = c("No", "Yes"))
casper$insulin_pump <- factor(casper$diabetes_med___1,
                              levels = c("Unchecked", "Checked"),
                              labels = c("No", "Yes")
)
casper$short_act_insulin <- factor(casper$insulin_type___1,
                                   levels = c("Unchecked", "Checked"),
                                   labels = c("No", "Yes")
)
casper$long_act_insulin <- factor(casper$insulin_type___2,
                                  levels = c("Unchecked", "Checked"),
                                  labels = c("No", "Yes")
)
casper$metformin <- factor(casper$diabetes_med_other___1,
                           levels = c("Unchecked", "Checked"),
                           labels = c("No", "Yes")
)
casper$tzds <- factor(casper$diabetes_med_other___2,
                      levels = c("Unchecked", "Checked"),
                      labels = c("No", "Yes")
)
casper$glp1_agonists <- factor(casper$diabetes_med_other___3,
                               levels = c("Unchecked", "Checked"),
                               labels = c("No", "Yes")
)
casper$sglt2_inhibitors <- factor(casper$diabetes_med_other___4,
                                  levels = c("Unchecked", "Checked"),
                                  labels = c("No", "Yes")
)
casper$other_diabetes_meds <- factor(casper$diabetes_med_other___5,
                                     levels = c("Unchecked", "Checked"),
                                     labels = c("No", "Yes")
)
# activity_factor and schofield
casper <- casper %>%
  unite(., activity_factor, fem_activity_factor, male_activity_factor, na.rm = T) %>%
  mutate(schofield = rowSums(.[c("schofield_female", "schofield_male")], na.rm = T))
casper$schofield[casper$schofield == 0] <- NA

# IMPROVE
# Diabetes medications
improve$insulin <- factor(improve$diabetes_med___2,
                            levels = c("Unchecked", "Checked"),
                            labels = c("No", "Yes")
)
improve$insulin_pump <- NA
improve$short_act_insulin <- factor(improve$insulin_type___1,
                                   levels = c("Unchecked", "Checked"),
                                   labels = c("No", "Yes")
)
improve$long_act_insulin <- factor(improve$insulin_type___2,
                                  levels = c("Unchecked", "Checked"),
                                  labels = c("No", "Yes")
)
improve$metformin <- factor(improve$diabetes_med___1,
                            levels = c("Unchecked", "Checked"),
                            labels = c("No", "Yes")
)
improve$tzds <- factor(improve$diabetes_med_other___1,
                      levels = c("Unchecked", "Checked"),
                      labels = c("No", "Yes")
)
improve$glp1_agonists <- factor(improve$diabetes_med_other___2,
                               levels = c("Unchecked", "Checked"),
                               labels = c("No", "Yes")
)
improve$sglt2_inhibitors <- factor(improve$diabetes_med_other___3,
                                  levels = c("Unchecked", "Checked"),
                                  labels = c("No", "Yes")
)
improve$other_diabetes_meds <- factor(improve$diabetes_med_other___4,
                                     levels = c("Unchecked", "Checked"),
                                     labels = c("No", "Yes")
)
# activity_factor and schofield
improve <- improve %>%
  unite(., activity_factor, activity_factor_female, activity_factor_male, na.rm = T) %>%
  mutate(schofield = rowSums(.[c("schofield_female", "schofield_male")], na.rm = T))
improve$schofield[improve$schofield == 0] <- NA
# Just first visit
improve_screen <- improve[grep("screening_arm", improve$redcap_event_name), ]
# Merge screening
screening <- do.call(rbind, list(
  renalheir[, screen_vars], penguin[, screen_vars],
  crocodile[, screen_vars], coffee[, screen_vars],
  casper[, screen_vars], improve_screen[, screen_vars]
))

###############################################################################
# DXA variables
###############################################################################

dxa_vars <- c(
  "subject_id", "study", "visit", 
  # Diabetes medications were formatted above because I (Tim) thought they were screening 
  # variables. But these variables change at each IMPROVE visit, so they are added 
  # here instead (the screening variables are limited to visit 1).
  "insulin", "insulin_pump", "metformin", "short_act_insulin", "long_act_insulin",
  "tzds", "glp1_agonists", "sglt2_inhibitors", "other_diabetes_meds",
  "dexa_date", "body_fat", "lean_mass",
  "trunk_mass", "fat_kg", "lean_kg", "trunk_kg", "bone_mineral_density",
  "body_composition_dxa_complete"
)

# RENAL-HEIR
# PENGUIN
penguin <- penguin %>%
  rename(
    dexa_date = dxa_date, body_fat = bodyfat_percent, lean_kg = leanmass_kg,
    fat_kg = fatmass_kg, lean_mass = leanmass_percent, trunk_mass = trunkmass_percent,
    trunk_kg = trunkmass_kg, bone_mineral_density = bmd,
    body_composition_dxa_complete = study_visit_dxa_scan_complete
  )
# CROCODILE
crocodile <- crocodile %>%
  rename(
    dexa_date = dxa_date, body_fat = bodyfat_percent, lean_kg = leanmass_kg,
    fat_kg = fatmass_kg, lean_mass = leanmass_percent, trunk_mass = trunkmass_percent,
    trunk_kg = trunkmass_kg, bone_mineral_density = bmd,
    body_composition_dxa_complete = study_visit_dxa_scan_complete
  )
# COFFEE
coffee[, dxa_vars[4:length(dxa_vars)]] <- NA
# CASPER
# IMPROVE
improve <- improve %>%
  rename(
    dexa_date = bodcomp_date, body_fat = dxa_body_fat, lean_kg = dxa_lean_kg,
    fat_kg = dxa_fat_kg, lean_mass = dxa_lean_mass, trunk_mass = dxa_trunk_mass,
    trunk_kg = dxa_trunk_kg, bone_mineral_density = dxa_bmd,
    body_composition_dxa_complete = dxa_complete
  )
# Merge DXA
dxa <- do.call(rbind, list(
  renalheir[, dxa_vars], penguin[, dxa_vars],
  crocodile[, dxa_vars], coffee[, dxa_vars],
  casper[, dxa_vars], improve[, dxa_vars]
))
###############################################################################
# Clamp vitals
###############################################################################

clamp_vitals <- c(
  "subject_id", "study", "visit", "clamp_date", "clamp_height",
  "clamp_weight", "clamp_sbp", "clamp_dbp", "clamp_map", "clamp_pls"
)

# RENAL-HEIR - already correct
# PENGUIN
penguin <- penguin %>%
  rename(clamp_height = clamp_ht, clamp_weight = clamp_wt)
penguin$clamp_date <- penguin$visit_date
penguin$clamp_sbp <- penguin$sys_bp
penguin$clamp_dbp <- penguin$dys_bp
penguin$clamp_map <- penguin$map
penguin$clamp_pls <- penguin$pulse
# CROCODILE
crocodile <- crocodile %>%
  rename(clamp_height = clamp_ht, clamp_weight = clamp_wt)
crocodile$clamp_date <- crocodile$visit_date
crocodile$clamp_sbp <- crocodile$sys_bp
crocodile$clamp_dbp <- crocodile$dys_bp
crocodile$clamp_map <- crocodile$map
crocodile$clamp_pls <- crocodile$pulse
# COFFEE
coffee$clamp_date <- coffee$cf_clamp_date
# CASPER - already correct
# IMPROVE - already correct
# Merge
clamp_vitals <- do.call(rbind, list(
  renalheir[, clamp_vitals], penguin[, clamp_vitals],
  crocodile[, clamp_vitals], coffee[, clamp_vitals],
  casper[, clamp_vitals], improve[, clamp_vitals]
))

###############################################################################
# Clamp labs
###############################################################################

clamp_labs <- c(
  "subject_id", "study", "visit", "cholesterol", "hdl", "ldl", "triglycerides",
  "total_protein", "serum_sodium", "cystatin_c", "serum_creatinine",
  "clamp_urine_mab_baseline", "clamp_urine_cre_baseline",
  "clamp_acr_baseline", "clamp_urine_sodium", "clamp_glucose_bl",
  "urine_glucose", "clamp_urine_mab_250", "clamp_urine_cre_250",
  "clamp_acr_250", "clamp_urine_vol", "hematocrit_minus_10",
  "hematocrit_minus_5", "hematocrit_90", "hematocrit_120"
)

# RENAL-HEIR
renalheir$hematocrit_minus_5 <- NA
# PENGUIN
penguin <- penguin %>% rename(
  cholesterol = bl_cholesterol, hdl = bl_hdl,
  ldl = bl_ldl, triglycerides = bl_triglycerides,
  total_protein = bl_tot_protein, serum_sodium = bl_na_s,
  cystatin_c = bl_cystatin_c_s,
  serum_creatinine = screen_serum_creatinine,
  clamp_urine_cre_baseline = screen_urine_cre,
  clamp_acr_baseline = screen_urine_acr, clamp_urine_sodium = bl_na_u,
  clamp_glucose_bl = bl_glucose_u
)
penguin[, c(
  "clamp_urine_mab_baseline", "urine_glucose", "clamp_urine_mab_250",
  "clamp_urine_cre_250", "clamp_acr_250", "clamp_urine_vol",
  "hematocrit_minus_10", "hematocrit_minus_5", "hematocrit_90", "hematocrit_120"
)] <- NA
# CROCODILE
crocodile <- crocodile %>% rename(
  cholesterol = bl_cholesterol, hdl = bl_hdl,
  ldl = bl_ldl, triglycerides = bl_triglycerides,
  total_protein = bl_tot_protein, serum_sodium = bl_na_s,
  serum_creatinine = screen_serum_creatinine,
  clamp_urine_cre_baseline = screen_urine_cre,
  clamp_acr_baseline = screen_urine_acr, clamp_urine_sodium = bl_na_u,
  clamp_glucose_bl = bl_glucose_u
)
crocodile[, c(
  "clamp_urine_mab_baseline", "urine_glucose", "clamp_urine_mab_250",
  "clamp_urine_cre_250", "clamp_acr_250", "clamp_urine_vol", "cystatin_c",
  "hematocrit_minus_10", "hematocrit_minus_5", "hematocrit_90", "hematocrit_120"
)] <- NA
# COFFEE
coffee[c("cholesterol", "hdl", "ldl", "triglycerides", "hematocrit_minus_10")] <- NA
coffee$clamp_glucose_bl <- coffee$fbg
# CASPER
casper$hematocrit_minus_10 <- NA
# IMPROVE
improve$hematocrit_minus_5 <- NA
# Merge
clamp_labs <- do.call(rbind, list(
  renalheir[, clamp_labs], penguin[, clamp_labs],
  crocodile[, clamp_labs], coffee[, clamp_labs],
  casper[, clamp_labs], improve[, clamp_labs]
))

###############################################################################
# 24 hour urine labs
###############################################################################

# Only available for CROCODILE, PENGUIN, and PANTHER
urine_labs <- c("subject_id", "study", "visit", "u24_labs", "u24_na", "u24_mab", "u24_volume", "u24_hours")

# RENAL-HEIR
renalheir[, tail(urine_labs, -3)] <- NA
# PENGUIN
penguin <- penguin %>% rename(u24_volume = u24_vl, u24_hours = u24_hrs)
# CROCODILE
crocodile$u24_hours <- NA
crocodile$u24_volume <- NA
# COFFEE
coffee[, tail(urine_labs, -3)] <- NA
# CASPER
casper[, tail(urine_labs, -3)] <- NA
# IMPROVE
improve[, tail(urine_labs, -3)] <- NA

# Merge
urine_labs <- do.call(rbind, list(
  renalheir[, urine_labs], penguin[, urine_labs],
  crocodile[, urine_labs], coffee[, urine_labs],
  casper[, urine_labs], improve[, urine_labs]
))

###############################################################################
# Hyperinsulinemic-euglycemic clamp data
###############################################################################

# Only available for CROCODILE and PENGUIN
he_clamp <- c(
  "subject_id", "study", "visit", "p1_raw_m", "p1_raw_leanm", "p1_gc_m",
  "p1_gc_leanm", "p2_raw_m", "p2_raw_leanm", "p2_gc_m", "p2_gc_leanm"
)

# RENAL-HEIR
renalheir[, tail(he_clamp, -3)] <- NA
# COFFEE
coffee[, tail(he_clamp, -3)] <- NA
# CASPER
casper[, tail(he_clamp, -3)] <- NA
# IMPROVE
improve[, tail(he_clamp, -3)] <- NA

# Merge
he_clamp <- do.call(rbind, list(
  renalheir[, he_clamp], penguin[, he_clamp],
  crocodile[, he_clamp], coffee[, he_clamp],
  casper[, he_clamp], improve[, he_clamp]
))

###############################################################################
# FFA data from clamps
###############################################################################

# RENAL-HEIR
renalheir_ffa <- renalheir[, c("subject_id", "study", "visit", colnames(renalheir)[grep("ffa_.*\\d{1,}", colnames(renalheir))])]
# PENGUIN
penguin_ffa <- penguin[, c("subject_id", "study", "visit", colnames(penguin)[grep("ffa_.*\\d{1,}", colnames(penguin))])]
colnames(penguin_ffa) <- sub("minus", "minus_", colnames(penguin_ffa))
# CROCODILE
crocodile_ffa <- crocodile[, c("subject_id", "study", "visit", colnames(crocodile)[grep("ffa_.*\\d{1,}", colnames(crocodile))])]
colnames(crocodile_ffa) <- sub("minus", "minus_", colnames(crocodile_ffa))
# COFFEE - none
# CASPER - none

# Merge
ffa <- full_join(renalheir_ffa, penguin_ffa)
ffa <- full_join(ffa, crocodile_ffa)
# Sort columns
ffa_cols <- colnames(ffa)[4:ncol(ffa)]
ffa_cols <- ffa_cols[order(as.numeric(sub("ffa_", "", sub("minus_", "-", ffa_cols))))]
ffa <- ffa[, c("subject_id", "study", "visit", ffa_cols)]

# IMPROVE
improve_ffa <- improve[, c("subject_id", "study", "visit", colnames(improve)[grep("ffa_.*\\d{1,}", colnames(improve))])]
colnames(improve_ffa) <- sub("neg", "minus", colnames(improve_ffa))
colnames(improve_ffa)[grep("ffa", colnames(improve_ffa))] <-
  paste0(colnames(improve_ffa)[grep("ffa", colnames(improve_ffa))], "_mmtt")
# Add IMPROVE
ffa <- full_join(ffa, improve_ffa)

###############################################################################
# C-peptide clamp data
###############################################################################

# RENAL-HEIR
renalheir_cpep <- renalheir[, c("subject_id", "study", "visit", colnames(renalheir)[grep("cpeptide_.*\\d{1,}", colnames(renalheir))])]
# IMPROVE
improve_cpep <- improve[, c("subject_id", "study", "visit", colnames(improve)[grep("cpep_.*\\d{1,}", colnames(improve))])]
colnames(improve_cpep) <- sub("mmtt_", "", colnames(improve_cpep))
colnames(improve_cpep)[4:ncol(improve_cpep)] <-
  paste0(colnames(improve_cpep)[4:ncol(improve_cpep)], "_mmtt")

# Merge
cpep <- full_join(renalheir_cpep, improve_cpep)

###############################################################################
# Insulin clamp data
###############################################################################

# RENAL-HEIR
renalheir_ins <- renalheir[, c("subject_id", "study", "visit", colnames(renalheir)[grep("insulin_.*\\d{1,}", colnames(renalheir))])]
renalheir_ins$insulin_type___1 <- NULL
renalheir_ins$insulin_type___2 <- NULL
# PENGUIN
penguin_ins <- penguin[, c("subject_id", "study", "visit", colnames(penguin)[grep("insulin_.*\\d{1,}", colnames(penguin))])]
colnames(penguin_ins) <- sub("minus", "minus_", colnames(penguin_ins))
# CROCODILE
crocodile_ins <- crocodile[, c("subject_id", "study", "visit", colnames(crocodile)[grep("insulin_.*\\d{1,}", colnames(crocodile))])]
colnames(crocodile_ins) <- sub("minus", "minus_", colnames(crocodile_ins))
# COFFEE - none
# CASPER - none

# Merge
ins <- full_join(renalheir_ins, penguin_ins)
ins <- full_join(ins, crocodile_ins)
# Sort columns
ins_cols <- colnames(ins)[4:ncol(ins)]
ins_cols <- ins_cols[order(as.numeric(sub("insulin_", "", sub("minus_", "-", ins_cols))))]
ins <- ins[, c("subject_id", "study", "visit", ins_cols)]

# IMPROVE
improve_ins <- improve[, c("subject_id", "study", "visit", colnames(improve)[grep("insulin_.*\\d{1,}", colnames(improve))])]
improve_ins$insulin_type___1 <- NULL
improve_ins$insulin_type___2 <- NULL
mmtt <- grep("mmtt", colnames(improve_ins))
colnames(improve_ins)[mmtt] <- sub("mmtt_", "", colnames(improve_ins)[mmtt])
colnames(improve_ins) <- sub("neg", "minus", colnames(improve_ins))
colnames(improve_ins)[mmtt] <-
  paste0(colnames(improve_ins)[mmtt], "_mmtt")
# Add IMPROVE
ins <- full_join(ins, improve_ins)


###############################################################################
# BG clamp data
###############################################################################

# RENAL-HEIR
renalheir_glu <- renalheir[, c("subject_id", "study", "visit", colnames(renalheir)[grep("glucose_.*\\d{1,}", colnames(renalheir))])]
# PENGUIN
penguin_glu <- penguin[, c("subject_id", "study", "visit", colnames(penguin)[grep("bg_.*\\d{1,}", colnames(penguin))])]
colnames(penguin_glu) <- sub("minus", "minus_", colnames(penguin_glu))
colnames(penguin_glu) <- sub("bg", "glucose", colnames(penguin_glu))
# CROCODILE
crocodile_glu <- crocodile[, c("subject_id", "study", "visit", colnames(crocodile)[grep("bg_.*\\d{1,}", colnames(crocodile))])]
colnames(crocodile_glu) <- sub("minus", "minus_", colnames(crocodile_glu))
colnames(crocodile_glu) <- sub("bg", "glucose", colnames(crocodile_glu))
# COFFEE
coffee_glu <- coffee[, c("subject_id", "study", "visit", colnames(coffee)[grep("glucose_.*\\d{1,}", colnames(coffee))])]
# CASPER
casper_glu <- casper[, c("subject_id", "study", "visit", colnames(casper)[grep("glucose_.*\\d{1,}", colnames(casper))])]
# IMPROVE
improve_glu <- improve[, c("subject_id", "study", "visit", colnames(improve)[grep("glucose_.*\\d{1,}", colnames(improve))])]

# Merge
glu <- full_join(renalheir_glu, penguin_glu)
glu <- full_join(glu, crocodile_glu)
glu <- full_join(glu, coffee_glu)
glu <- full_join(glu, casper_glu)
glu <- full_join(glu, improve_glu)

###############################################################################
# Kidney function
###############################################################################

kidney_vars <- c(
  "subject_id", "study", "visit", "gfr", "gfr_bsa", "pah_clear_abs",
  "pah_clear_bsa", "rpf", "rpf_bsa", "ecv", "gfr_ecv_percent", "gfr_ecv_std"
)

# RENAL-HEIR - already correct
# PENGUIN
penguin <- penguin %>% rename(
  pah_clear_abs = pah_abs, pah_clear_bsa = pah_bsa,
  gfr_ecv_percent = gfr_ecv, gfr_ecv_std = gfr_standard,
  rpf = erpf, rpf_bsa = erpf_bsa
)
# CROCODILE
crocodile <- crocodile %>% rename(
  gfr = gfr_raw, pah_clear_abs = pah_raw,
  pah_clear_bsa = pah_bsa, rpf = erpf,
  rpf_bsa = erpfbsa
)
crocodile[, c("ecv", "gfr_ecv_percent", "gfr_ecv_std")] <- NA
# COFFEE
coffee <- coffee %>% rename(
  gfr = gfr_abs, gfr_bsa = gfr_adj, pah_clear_abs = abs_pah_clear,
  rpf = rpf_abs, rpf_bsa = rpf_adj
)
# CASPER
casper <- casper %>% rename(
  pah_clear_abs = abs_pah, pah_clear_bsa = pah_bsa,
  gfr_ecv_percent = gfr_ecv, gfr_ecv_std = gfr_standard
)
# IMPROVE
improve <- improve %>% rename(
  pah_clear_abs = abs_pah, pah_clear_bsa = pah_bsa,
  rpf_bsa = erpf_bsa
)
# Merge
kidney <- do.call(rbind, list(
  renalheir[, kidney_vars], penguin[, kidney_vars],
  crocodile[, kidney_vars], coffee[, kidney_vars],
  casper[, kidney_vars], improve[, kidney_vars]
))

###############################################################################
# Kidney MRI
###############################################################################

kidney_mri <- c(
  "subject_id", "study", "visit", "o2_sats", "pcasl3d_right", "pasl2d_right",
  "adc_right", "length_right", "width_right", "depth_right",
  "volume_right", "bold_r_bl_cortex", "bold_r_bl_medulla",
  "bold_r_bl_kidney", "bold_r_pf_cortex", "bold_r_pf_medulla",
  "bold_r_pf_kidney", "pcasl3d_left", "pasl2d_left", "adc_left",
  "length_left", "width_left", "depth_left", "volume_left",
  "bold_l_bl_cortex", "bold_l_bl_medulla", "bold_l_bl_kidney",
  "bold_l_pf_cortex", "bold_l_pf_medulla", "bold_l_pf_kidney"
)

# RENAL-HEIR
renalheir <- renalheir %>% rename(pcasl3d_right = asl_right, pcasl3d_left = asl_left)
renalheir[, c(
  "adc_right", "length_right", "width_right",
  "depth_right", "volume_right", "adc_left", "length_left",
  "width_left", "depth_left", "volume_left", "pasl2d_right", "pasl2d_left"
)] <- NA
# PENGUIN
penguin[, tail(kidney_mri, -3)] <- NA
# CROCODILE
crocodile <- crocodile %>% rename(o2_sats = o2_sat)
# COFFEE
coffee <- coffee %>% rename(pcasl3d_right = asl_right, pcasl3d_left = asl_left)
coffee[, c(
  "adc_right", "length_right", "width_right",
  "depth_right", "volume_right", "adc_left", "length_left",
  "width_left", "depth_left", "volume_left", "pasl2d_right", "pasl2d_left"
)] <- NA
# CASPER
casper <- casper %>% rename(pcasl3d_right = asl_right, pcasl3d_left = asl_left)
casper[, c(
  "adc_right", "length_right", "width_right",
  "depth_right", "volume_right", "adc_left", "length_left",
  "width_left", "depth_left", "volume_left", "pasl2d_right", "pasl2d_left"
)] <- NA
# IMPROVE
improve <- improve %>% rename(pcasl3d_right = asl_right, pcasl3d_left = asl_left)
improve[, c(
  "adc_right", "length_right", "width_right",
  "depth_right", "volume_right", "adc_left", "length_left",
  "width_left", "depth_left", "volume_left", "pasl2d_right", "pasl2d_left"
)] <- NA
# Merge
kidney_mri <- do.call(rbind, list(
  renalheir[, kidney_mri], penguin[, kidney_mri],
  crocodile[, kidney_mri], coffee[, kidney_mri],
  casper[, kidney_mri], improve[, kidney_mri]
))

# Calculate FSOC = bl_bold - pf_bold
kidney_mri = kidney_mri %>%
  mutate(fsoc_r_cortex = bold_r_bl_cortex - bold_r_pf_cortex,
         fsoc_r_medulla = bold_r_bl_medulla - bold_r_pf_medulla,
         fsoc_r_kidney = bold_r_bl_kidney - bold_r_pf_kidney,
         fsoc_l_cortex = bold_l_bl_cortex - bold_l_pf_cortex,
         fsoc_l_medulla = bold_l_bl_medulla - bold_l_pf_medulla,
         fsoc_l_kidney = bold_l_bl_kidney - bold_l_pf_kidney)

###############################################################################
# PET/CT
###############################################################################

# Only in PENGUIN and CROCODILE, names are consistent in both
pet_vars <- c(
  "subject_id", "study", "visit", "pet_rc_f", "pet_rc_k2", "pet_rc_vb", "pet_rc_k1",
  "pet_rm_f", "pet_rm_k2", "pet_rm_vb", "pet_rm_k1", "pet_lc_f", "pet_lc_k2",
  "pet_lc_vb", "pet_lc_k1", "pet_lm_f", "pet_lm_k2", "pet_lm_vb", "pet_lm_k1"
)

# RENAL-HEIR
renalheir[, tail(pet_vars, -3)] <- NA
# COFFEE
coffee[, tail(pet_vars, -3)] <- NA
# CASPER
casper[, tail(pet_vars, -3)] <- NA
# IMPROVE
improve[, tail(pet_vars, -3)] <- NA
# Merge
pet <- do.call(rbind, list(
  renalheir[, pet_vars], penguin[, pet_vars],
  crocodile[, pet_vars], coffee[, pet_vars],
  casper[, pet_vars], improve[, pet_vars]
))

###############################################################################
# Kidney biopsy
###############################################################################

kidney_biopsy <- c(
  "subject_id", "study", "visit", "bx_date", "bx_kit_id","age_biopsy",
  "vitals_height", "vitals_weight", "vitals_bmi",
  "vitals_sbp", "vitals_dbp","bun",
  "gloms", paste0("glom_enlarge___", 1:7),
  "gloms_gs", "ifta", "vessels___1", "vessels___2", "vessels_other",
  "fia", "glom_tuft_area", "glom_volume_wiggins",
  "mes_matrix_area", "mes_index", "mes_volume_wiggins",
  "glom_nuc_count", "mes_nuc_count"
)
# RENAL-HEIR - correct (names taken from this study)
renalheir$age_biopsy = renalheir$age_current
renalheir$bun = renalheir$labs_bun
# PENGUIN
penguin$bun = penguin$bl_bun_s
penguin[, setdiff(kidney_biopsy,colnames(penguin))] <- NA
# CROCODILE - all correct (glom_enlarge levels match)
crocodile$age_biopsy = crocodile$age_current
crocodile$bun = crocodile$labs_bun
crocodile[, setdiff(kidney_biopsy,colnames(crocodile))] <- NA
# COFFEE
coffee[, setdiff(kidney_biopsy,colnames(coffee))] <- NA
# CASPER
casper[, setdiff(kidney_biopsy,colnames(casper))] <- NA
# IMPROVE
improve = improve %>% group_by(subject_id) %>% fill(dob)
improve$dob = parse_date(improve$dob,approx = F)
improve$bx_date = parse_date(improve$bx_date,approx = F)
improve$age_biopsy = round(as.numeric(difftime(improve$bx_date,improve$dob,units = "days"))/365.25)
improve$bun = improve$labs_bun

kidney_biopsy <- do.call(rbind, list(
  renalheir[, kidney_biopsy], penguin[, kidney_biopsy],
  crocodile[, kidney_biopsy], coffee[, kidney_biopsy],
  casper[, kidney_biopsy], improve[, kidney_biopsy]
))
kidney_biopsy <- kidney_biopsy %>%
  rename(
    bx_sbp = vitals_sbp, bx_dbp = vitals_dbp, bx_height = vitals_height,
    bx_weight = vitals_weight
  )
# Capitalize all box IDs
kidney_biopsy$bx_kit_id <- toupper(kidney_biopsy$bx_kit_id)

###############################################################################
# Merge everything together, fill non-changing variables
###############################################################################
df <- full_join(demographics, screening)
df <- full_join(df, dxa)
df <- full_join(df, clamp_vitals)
df <- full_join(df, clamp_labs)
df <- full_join(df, urine_labs)
df <- full_join(df, he_clamp)
df <- full_join(df, ffa)
df <- full_join(df, cpep)
df <- full_join(df, ins)
df <- full_join(df, glu)
df <- full_join(df, kidney)
df <- full_join(df, kidney_mri)
df <- full_join(df, kidney_biopsy)
df <- full_join(df, pet)

# Fill
fill_vars <- c("group", "dob", "gender", "race", "ethnicity", "diagnosis_date")
df <- df %>%
  group_by(subject_id) %>%
  fill(all_of(fill_vars), .direction = "downup")

###############################################################################
# Final formatting and calculated fields
###############################################################################

# Age
df$age = coalesce(df$age_consent,df$age_biopsy)

# BMI
df$bmi = coalesce(df$screen_bmi,df$vitals_bmi)

# Diabetes duration
df$disease_duration = round(as.numeric(difftime(df$diagnosis_date,df$dob,units = "days"))/365.25)

# BMI percentile
## Excluding adults
df$screen_bmi_z <- sds(
  value = df$screen_bmi,
  age = df$age,
  sex = df$gender, male = "Male", female = "Female",
  item = "bmi", type = "SDS",
  ref = cdc.ref
)
df$screen_bmi_percentile <- sds(
  value = df$screen_bmi,
  age = df$age,
  sex = df$gender, male = "Male", female = "Female",
  item = "bmi", type = "perc",
  ref = cdc.ref
)

## Including adults - over 20 treated as age == 20
df$screen_bmi_z_all <- sds(
  value = df$screen_bmi,
  age = ifelse(df$age >= 20, 20, df$age),
  sex = df$gender, male = "Male", female = "Female",
  item = "bmi", type = "SDS",
  ref = cdc.ref
)
df$screen_bmi_percentile_all <- sds(
  value = df$screen_bmi,
  age = ifelse(df$age >= 20, 20, df$age),
  sex = df$gender, male = "Male", female = "Female",
  item = "bmi", type = "perc",
  ref = cdc.ref
)

# Various eGFRs
df <- data.frame(cbind(df, egfr_calc(
  age = df$age, serum_creatinine = df$serum_creatinine,
  cystatin_c = df$cystatin_c, height = df$clamp_height, sex = df$gender
)))

# Hemodynamics
hemodynamics <- c("Pglo", "Ra", "Re", "RVR", "FF", "RBF")

hemo <- hemodynamics_calc(
  total_protein = df$total_protein,
  gfr = df$gfr, rpf = df$rpf, map = df$map,
  hematocrit_minus_10 = df$hematocrit_minus_10,
  hematocrit_minus_5 = df$hematocrit_minus_5,
  hematocrit_90 = df$hematocrit_90,
  hematocrit_120 = df$hematocrit_120,
  group = df$group
)

hemo <- hemo %>% rename(
  Pglo = Pglo_abs, Ra = Ra_abs, Re = Re_abs, RVR = RVRAbs,
  FF = FFabs, RBF = RBFabs
)

df <- data.frame(cbind(df, hemo[, hemodynamics]))

# Co-enroll "No" values
df$co_enroll[is.na(df$co_enroll)] <- "No"

# Visit names
df$visit <- factor(df$visit,
  levels = c(
    "Baseline", "Pre-Surgery",
    "3 Months Post-Surgery", "12 Months Post-Surgery"
  )
)

# Fix CROCODILE IDs
df$subject_id[df$study=="CROCODILE"] = paste0("CRC-",str_pad(df$subject_id[df$study=="CROCODILE"],2,"left","0"))

# Sort and write!
df <- df %>% arrange(study, subject_id, visit)
write.csv(df,
  file = paste0("./Data Clean/merged_dataset_", Sys.Date(), ".csv"),
  row.names = F, na = ""
)
