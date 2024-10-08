library(redcapAPI)
library(tidyverse)
library(Hmisc)
home_dir <- switch(Sys.info()["sysname"],
  "Darwin" = "/Users/timvigers/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Vigers/BDC/Janet Snell-Bergeon/CALICO",
  "Windows" = "C:/Users/timvigers/OneDrive - The University of Colorado Denver/Vigers/BDC/Janet Snell-Bergeon/CALICO",
  "Linux" = "/home/timvigers/OneDrive/Vigers/BDC/Janet Snell-Bergeon/CALICO"
)
setwd(home_dir)
# Define variable sets
# Demographics
demo_vars <- c(
  "site", "pcosdx_age", "race___1", "race___2", "race___3", "race___4",
  "race___5", "race___60", "race___unk", "race___na", "race___oth", "race___pm",
  "ethnicity", "insur_type", "pcosdx_irregular_menses", "pcosdx_anxietydx_age"
)
# Aim 1 for Melanie
aim1_vars <- c(
  "cv_bmi", "cv_bmi_percentile", "cv_bmi_z",
  "cv_hirsutism_num", "cv_hirsutism_cat", "cv_acneface",
  "cv_acneother___1", "cv_acneother___2", "cv_acneother___0",
  "cv_acneother___unk", "cv_acneother___na", "cv_acneother___oth",
  "cv_acneother___pm", "cv_ft", "cv_ft_perc",
  "cv_tt", "cv_tt_perc", "cv_dheas",
  "cv_dheas_perc", "cv_androstendione", "cv_androstendione_perc",
  "cv_lh", "cv_fsh", "cv_amh",
  "cv_a1c", "cv_tg", "cv_hdl",
  "cv_shbg", "cv_alt", "cv_fbg",
  "cv_fastinsulin", "cv_2hrglucoseogtt", "cv_acanthosisneck",
  "cv_waist", "cv_osa_sx", "cv_medications___1",
  "cv_medications___2", "cv_medications___3", "cv_medications___4",
  "cv_medications___5", "cv_medications___6", "cv_medications___7",
  "cv_medications___8", "cv_medications___9", "cv_medications___10",
  "cv_medications___11", "cv_medications___12", "cv_medications___13",
  "cv_medications___14", "cv_medications___15", "cv_medications___16",
  "cv_medications___17", "cv_medications___18", "cv_medications___19",
  "cv_medications___20", "cv_medications___21", "cv_medications___22",
  "cv_medications___23", "cv_medications___32", "cv_medications___24",
  "cv_medications___25", "cv_medications___26", "cv_medications___27",
  "cv_medications___28", "cv_medications___29", "cv_medications___30",
  "cv_medications___31", "cv_medications___60", "cv_medications___0",
  "cv_medications___unk", "cv_medications___na", "cv_medications___oth",
  "cv_medications___pm"
)
# LARC for BCH
larc_vars <- c(
  "cv_type_of_clinic", "cv_mdc_specialties___1", "cv_mdc_specialties___2",
  "cv_mdc_specialties___3", "cv_mdc_specialties___4", "cv_mdc_specialties___5",
  "cv_mdc_specialties___6", "cv_mdc_specialties___60",
  "cv_mdc_specialties___unk", "cv_mdc_specialties___na",
  "cv_mdc_specialties___oth", "cv_mdc_specialties___pm", "cv_medications___10",
  "cv_medications___11", "cv_weight", "cv_bmi", "cv_acneface",
  "cv_hirsutism_num", "cv_hirsutism_cat", "cv_menfreq", "race___1", "race___2",
  "race___3", "race___4", "race___5", "race___60", "race___unk", "race___na",
  "race___oth", "race___pm", "ethnicity", "insur_type", "cv_genderid",
  "cv_genderid_other", "cv_sexuallyactive", "cv_sexualpref",
  "cv_monthssincepcosdx", "cv_hgb", "cv_ldl", "cv_hdl", "cv_sbp", "cv_dbp"
)
# SES variables for health disparities
ses_vars <- c(
  "pcosdx_birthweight", "pcosdx_birthweight_oz", "pcosdx_birthweight_calc",
  "pcosdx_pmh___1", "pcosdx_pmh___2", "pcosdx_pmh___3", "pcosdx_pmh___4",
  "pcosdx_pmh___5", "pcosdx_pmh___6", "pcosdx_pmh___29", "pcosdx_pmh___7",
  "pcosdx_pmh___8", "pcosdx_pmh___9", "pcosdx_pmh___10", "pcosdx_pmh___11",
  "pcosdx_pmh___12", "pcosdx_pmh___22", "pcosdx_pmh___13", "pcosdx_pmh___14",
  "pcosdx_pmh___15", "pcosdx_pmh___16", "pcosdx_pmh___23", "pcosdx_pmh___24",
  "pcosdx_pmh___17", "pcosdx_pmh___18", "pcosdx_pmh___19", "pcosdx_pmh___20",
  "pcosdx_pmh___21", "pcosdx_pmh___60", "pcosdx_pmh___0", "pcosdx_pmh___unk",
  "pcosdx_pmh___na", "pcosdx_pmh___oth", "pcosdx_pmh___pm", "pmh_other",
  "pcosdx_famhx_parent___27", "pcosdx_famhx_parent___28",
  "pcosdx_famhx_parent___30", "pcosdx_famhx_parent___5",
  "pcosdx_famhx_parent___29", "pcosdx_famhx_parent___7",
  "pcosdx_famhx_parent___24", "pcosdx_famhx_parent___8",
  "pcosdx_famhx_parent___25", "pcosdx_famhx_parent___12",
  "pcosdx_famhx_parent___22", "pcosdx_famhx_parent___13",
  "pcosdx_famhx_parent___14", "pcosdx_famhx_parent___26",
  "pcosdx_famhx_parent___15", "pcosdx_famhx_parent___16",
  "pcosdx_famhx_parent___17", "pcosdx_famhx_parent___18",
  "pcosdx_famhx_parent___31", "pcosdx_famhx_parent___32",
  "pcosdx_famhx_parent___33", "pcosdx_famhx_parent___60",
  "pcosdx_famhx_parent___0", "pcosdx_famhx_parent___unk",
  "pcosdx_famhx_parent___na", "pcosdx_famhx_parent___oth",
  "pcosdx_famhx_parent___pm", "pcosdx_famhx___27", "pcosdx_famhx___28",
  "pcosdx_famhx___30", "pcosdx_famhx___5", "pcosdx_famhx___29",
  "pcosdx_famhx___7", "pcosdx_famhx___24", "pcosdx_famhx___8",
  "pcosdx_famhx___25", "pcosdx_famhx___12", "pcosdx_famhx___22",
  "pcosdx_famhx___13", "pcosdx_famhx___14", "pcosdx_famhx___26",
  "pcosdx_famhx___15", "pcosdx_famhx___16", "pcosdx_famhx___17",
  "pcosdx_famhx___18", "pcosdx_famhx___31", "pcosdx_famhx___32",
  "pcosdx_famhx___33", "pcosdx_famhx___60", "pcosdx_famhx___0",
  "pcosdx_famhx___unk", "pcosdx_famhx___na", "pcosdx_famhx___oth",
  "pcosdx_famhx___pm", "pastmeds___1", "pastmeds___2", "pastmeds___3",
  "pastmeds___4", "pastmeds___5", "pastmeds___6", "pastmeds___7",
  "pastmeds___8", "pastmeds___9", "pastmeds___10", "pastmeds___11",
  "pastmeds___12", "pastmeds___13", "pastmeds___14", "pastmeds___15",
  "pastmeds___16", "pastmeds___60", "pastmeds___0", "pastmeds___unk",
  "pastmeds___na", "pastmeds___oth", "pastmeds___pm", "pastmeds_other",
  "pcosdx_medcomments", "pcosdx_mentalhealthcounseling___1",
  "pcosdx_mentalhealthcounseling___2", "pcosdx_mentalhealthcounseling___3",
  "pcosdx_mentalhealthcounseling___0", "pcosdx_mentalhealthcounseling___unk",
  "pcosdx_mentalhealthcounseling___na", "pcosdx_mentalhealthcounseling___oth",
  "pcosdx_mentalhealthcounseling___pm", "pcosdx_age", "pcosdx_specialty",
  "pcosdx_specialty_other", "cv_visittype", "cv_type_of_clinic", "clinic_other",
  "cv_mdc_specialties___1", "cv_mdc_specialties___2", "cv_mdc_specialties___3",
  "cv_mdc_specialties___4", "cv_mdc_specialties___5", "cv_mdc_specialties___6",
  "cv_mdc_specialties___60", "cv_mdc_specialties___unk",
  "cv_mdc_specialties___na", "cv_mdc_specialties___oth",
  "cv_mdc_specialties___pm", "cv_age", "cv_monthssincepcosdx",
  "cv_newdx___3", "cv_newdx___4", "cv_newdx___5", "cv_newdx___6",
  "cv_newdx___29", "cv_newdx___7", "cv_newdx___8", "cv_newdx___9",
  "cv_newdx___10", "cv_newdx___11", "cv_newdx___12", "cv_newdx___22",
  "cv_newdx___13", "cv_newdx___14", "cv_newdx___15", "cv_newdx___16",
  "cv_newdx___27", "cv_newdx___28", "cv_newdx___17", "cv_newdx___18",
  "cv_newdx___19", "cv_newdx___20", "cv_newdx___21", "cv_newdx___23",
  "cv_newdx___24", "cv_newdx___25", "cv_newdx___26", "cv_newdx___60",
  "cv_newdx___0", "cv_newdx___unk", "cv_newdx___na", "cv_newdx___oth",
  "cv_newdx___pm", "cv_medications___1", "cv_medications___2",
  "cv_medications___3", "cv_medications___4", "cv_medications___5",
  "cv_medications___6", "cv_medications___7", "cv_medications___8",
  "cv_medications___9", "cv_medications___10", "cv_medications___11",
  "cv_medications___12", "cv_medications___13", "cv_medications___14",
  "cv_medications___15", "cv_medications___16", "cv_medications___17",
  "cv_medications___18", "cv_medications___19", "cv_medications___20",
  "cv_medications___21", "cv_medications___22", "cv_medications___23",
  "cv_medications___32", "cv_medications___24", "cv_medications___25",
  "cv_medications___26", "cv_medications___27", "cv_medications___28",
  "cv_medications___29", "cv_medications___30", "cv_medications___31",
  "cv_medications___60", "cv_medications___0", "cv_medications___unk",
  "cv_medications___na", "cv_medications___oth", "cv_medications___pm",
  "cv_medscomments", "cv_othercomment", "cv_hormonecomments", "cv_bmi",
  "cv_bmi_percentile", "cv_bmi_z", "cv_a1c", "cv_tt_perc", "cv_ft_perc",
  "cv_shbg", "cv_fbg", "cv_2hrglucoseogtt", "cv_tg_fasting", "cv_tg", "cv_hdl",
  "cv_ldl", "cv_tc", "cv_alt", "cv_ast", "cv_newmeds___1", "cv_newmeds___2",
  "cv_newmeds___3", "cv_newmeds___4", "cv_newmeds___5", "cv_newmeds___6",
  "cv_newmeds___7", "cv_newmeds___8", "cv_newmeds___9", "cv_newmeds___10",
  "cv_newmeds___11", "cv_newmeds___12", "cv_newmeds___13", "cv_newmeds___14",
  "cv_newmeds___15", "cv_newmeds___16", "cv_newmeds___17", "cv_newmeds___18",
  "cv_newmeds___19", "cv_newmeds___20", "cv_newmeds___21", "cv_newmeds___22",
  "cv_newmeds___23", "cv_newmeds___32", "cv_newmeds___24", "cv_newmeds___25",
  "cv_newmeds___26", "cv_newmeds___27", "cv_newmeds___28", "cv_newmeds___29",
  "cv_newmeds___30", "cv_newmeds___31", "cv_newmeds___60", "cv_newmeds___na",
  "cv_newmeds___unk", "cv_newmeds___oth", "cv_newmeds___pm",
  "cv_newcont_meds_othercomment", "cv_newcont_meds_hormonecomments",
  "cv_newcont_meds_discontinued", "cv_referrals___1", "cv_referrals___2",
  "cv_referrals___3", "cv_referrals___4", "cv_referrals___5",
  "cv_referrals___6", "cv_referrals___7", "cv_referrals___8",
  "cv_referrals___9", "cv_referrals___10", "cv_referrals___11",
  "cv_referrals___12", "cv_referrals___60", "cv_referrals___na",
  "cv_referrals___unk", "cv_referrals___oth", "cv_referrals___pm",
  "cv_referrals_other", "cv_referrals_reason"
)
# Open REDCap
unlockREDCap(c(rcon = "CALICO"),
  keyring = "API_KEYs",
  envir = 1,
  url = "https://redcap.ucdenver.edu/api/"
)
# Function for converting "UNK", etc. to NA in numeric fields
isMissingSpecial <- function(x, ...) {
  is.na(x) | x == "" | x == "UNK" | x == "OTH" | x == "NA" | x == "PM"
}
df <- exportRecordsTyped(rcon,
  fields = c(demo_vars, aim1_vars, larc_vars), warn_zero_coded = F,
  api_param = list(
    filterLogic = "[history_complete] = 2 || [clinical_visit_complete] = 2"
  ),
  na = list(number = isMissingSpecial, radio = isMissingSpecial)
)
# Obesity categories
df <- df %>%
  mutate(
    Obesity_perc = case_when(
      !is.na(cv_bmi_percentile) & cv_bmi_percentile >= 85 ~ "Yes",
      !is.na(cv_bmi_percentile) & cv_bmi_percentile < 85 ~ "No",
      is.na(cv_bmi_percentile) & cv_bmi >= 25 ~ "Yes",
      is.na(cv_bmi_percentile) & cv_bmi < 25 ~ "No",
      .default = NA
    ),
    Obesity_raw = ifelse(cv_bmi >= 25, "Yes", "No")
  )
df$Obesity_perc <- factor(df$Obesity_perc, levels = c("No", "Yes"))
df$Obesity_raw <- factor(df$Obesity_raw, levels = c("No", "Yes"))
# Combined race column
races <- c(
  "Caucasian", "African American", "Asian", "Pacific Islander",
  "American Indian or Alaska Native", "Other", "Unknown/Not recorded",
  "Not applicable", "Other", "Premenarchal"
)
df$combined_race <- apply(
  df[, paste0("race___", c(1:5, 60, "unk"))], 1,
  function(r) {
    w <- which(r == "Checked")
    if (length(w) == 0) {
      return(NA)
    } else if (length(w) > 1) {
      return("More than one")
    } else {
      return(races[w])
    }
  }
)
df$combined_race <- factor(df$combined_race)
# LARC
df$larc <- df$cv_medications___10 == "Checked" |
  df$cv_medications___11 == "Checked"
df$larc <- factor(df$larc, levels = c(F, T), labels = c("No", "Yes"))
# Ever on LARC vs. never
larc_users <- unique(df$record_number[df$larc == "Yes"])
df$larc_ever <- "No"
df$larc_ever[df$record_number %in% larc_users] <- "Yes"
df$larc_ever <- factor(df$larc_ever, levels = c("No", "Yes"))
# Fill down select columns
df <- df %>%
  group_by(record_number) %>%
  fill(combined_race, ethnicity, pcosdx_age, larc_ever) %>%
  ungroup()
# Age group at diagnosis
df$age_group <- cut(df$pcosdx_age, c(-Inf, 15, Inf),
  right = F, labels = c("< 15 years", ">= 15 years")
)
# Convert columns to numeric
df$cv_a1c <- suppressWarnings(as.numeric(df$cv_a1c))
df$redcap_repeat_instance <- as.numeric(df$redcap_repeat_instance)
# Fix/add labels
label(df$ethnicity) <- "Ethnicity"
label(df$pcosdx_age) <- "Age at time of PCOS diagnosis"
label(df$age_group) <- "Age Group at Diagnosis"
label(df$combined_race) <- "Race"
label(df$Obesity_perc) <- "Obesity Status (by percentile and raw value)"
label(df$Obesity_raw) <- "Obesity Status (by raw value)"
label(df$cv_a1c) <- "HbA1C"
label(df$larc) <- "On LARC"
label(df$larc_ever) <- "LARC Ever User"
label(df[, grep("___unk", colnames(df))]) <-
  as.list(sub(
    "choice=NA", "choice=Unknown/Not recorded",
    label(df[, grep("___unk", colnames(df))])
  ))
label(df[, grep("___na", colnames(df))]) <-
  as.list(sub(
    "choice=NA", "choice=Not applicable",
    label(df[, grep("___na", colnames(df))])
  ))
label(df[, grep("___oth", colnames(df))]) <-
  as.list(sub(
    "choice=NA", "choice=Other",
    label(df[, grep("___oth", colnames(df))])
  ))
label(df[, grep("___pm", colnames(df))]) <-
  as.list(sub(
    "choice=NA", "choice=Premenarchal",
    label(df[, grep("___pm", colnames(df))])
  ))
# Save
save(df, demo_vars, aim1_vars, larc_vars, ses_vars,
  file = "./Data_Clean/analysis_data.RData"
)
