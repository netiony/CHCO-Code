library(redcapAPI)
library(tidyverse)
library(Hmisc)
home_dir <- switch(Sys.info()["sysname"],
  "Darwin" = "/Users/timvigers/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Vigers/BDC/Janet Snell-Bergeon/CALICO",
  "Linux" = "/home/timvigers/OneDrive/Vigers/BDC/Janet Snell-Bergeon/CALICO"
)
setwd(home_dir)
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
# Import history and clinic visits separately for cleanliness
history <- exportReportsTyped(rcon,
  report_id = 131532, warn_zero_coded = F,
  na = list(number = isMissingSpecial, radio = isMissingSpecial)
)
clinic <- exportReportsTyped(rcon,
  report_id = 131533, warn_zero_coded = F,
  na = list(number = isMissingSpecial, radio = isMissingSpecial)
)
# Completed  entries only, drop some columns
history <- history %>%
  filter(history_complete == "Complete") %>%
  select(-redcap_repeat_instrument, -redcap_repeat_instance, -data_entryname)
clinic <- clinic %>%
  filter(clinical_visit_complete == "Complete") %>%
  select(-redcap_repeat_instrument)
# Merge
df <- full_join(history, clinic, by = join_by(record_number)) %>%
  select(record_number, redcap_repeat_instance, everything())
# Diagnostic visits are 0 months since diagnosis
df$cv_monthssincepcosdx[df$cv_visittype == "Diagnostic Visit"] <- 0
# Define variable sets
# Demographics
demo_vars <- c(
  "site", "pcosdx_age", "combined_race", "ethnicity", "insur_type"
)
# Aim 1 for Melanie
aim1_vars <- c(
  "cv_bmi", "cv_bmi_percentile", "cv_bmi_z", "cv_hirsutism_num",
  "cv_hirsutism_cat", "cv_acneface", "cv_acneother___1", "cv_acneother___2",
  "cv_acneother___0", "cv_acneother___unk", "cv_ft", "cv_ft_perc", "cv_tt",
  "cv_tt_perc", "cv_dheas", "cv_dheas_perc", "cv_androstendione",
  "cv_androstendione_perc", "cv_lh", "cv_fsh", "cv_amh", "cv_a1c", "cv_tg",
  "cv_hdl", "cv_shbg", "cv_alt", "cv_fbg", "cv_fastinsulin",
  "cv_2hrglucoseogtt", "cv_acanthosisneck", "cv_waist", "cv_osa_sx",
  "cv_medications___1", "cv_medications___2", "cv_medications___3",
  "cv_medications___4", "cv_medications___5", "cv_medications___6",
  "cv_medications___7", "cv_medications___8", "cv_medications___9",
  "cv_medications___10", "cv_medications___11", "cv_medications___12",
  "cv_medications___13", "cv_medications___14", "cv_medications___15",
  "cv_medications___16", "cv_medications___17", "cv_medications___18",
  "cv_medications___19", "cv_medications___20", "cv_medications___21",
  "cv_medications___22", "cv_medications___23", "cv_medications___32",
  "cv_medications___24", "cv_medications___25", "cv_medications___26",
  "cv_medications___27", "cv_medications___28", "cv_medications___29",
  "cv_medications___30", "cv_medications___31", "cv_medications___60",
  "cv_medications___0", "cv_medications___unk"
)
# NCH/CHOP variables
nch_chop_vars <- c(
  "pcosdx_pmh___1", "pcosdx_pmh___2", "pcosdx_pmh___3", "pcosdx_pmh___4",
  "pcosdx_pmh___5", "pcosdx_pmh___6", "pcosdx_pmh___29", "pcosdx_pmh___7",
  "pcosdx_pmh___8", "pcosdx_pmh___9", "pcosdx_pmh___10", "pcosdx_pmh___11",
  "pcosdx_pmh___12", "pcosdx_pmh___22", "pcosdx_pmh___13", "pcosdx_pmh___14",
  "pcosdx_pmh___15", "pcosdx_pmh___16", "pcosdx_pmh___23", "pcosdx_pmh___24",
  "pcosdx_pmh___17", "pcosdx_pmh___18", "pcosdx_pmh___19", "pcosdx_pmh___20",
  "pcosdx_pmh___21", "pcosdx_pmh___60", "pcosdx_pmh___0", "pcosdx_pmh___unk",
  "pcosdx_pmh___na", "pcosdx_pmh___oth", "pcosdx_pmh___pm",
  "pcosdx_obesitydx_age", "pcosdx_obesitydx_earliestage",
  "pcosdx_t1ddx_age", "pcosdx_t2ddx_age", "pcosdx_osadx_age",
  "pcosdx_anxietydx_age", "pcosdx_depressiondx_age",
  "pcosdx_bingedx_age", "pcosdx_restrictdx_age", "pcosdx_adhddx_age",
  "pcosdx_mhhosp_age", "pcosdx_siattempt_age", "pcosdx_surghx_bariatric",
  "pcosdx_surghx_pilonidal", "pcosdx_surghx_hidradenitis", "pastmeds___1",
  "pastmeds___2", "pastmeds___3", "pastmeds___4", "pastmeds___5",
  "pastmeds___6", "pastmeds___7", "pastmeds___8", "pastmeds___9",
  "pastmeds___10", "pastmeds___11", "pastmeds___12", "pastmeds___13",
  "pastmeds___14", "pastmeds___15", "pastmeds___16", "pastmeds___60",
  "pastmeds___0", "pastmeds___unk", "pastmeds___na", "pastmeds___oth",
  "pastmeds___pm", "pcosdx_mentalhealthcounseling___1",
  "pcosdx_mentalhealthcounseling___2", "pcosdx_mentalhealthcounseling___3",
  "pcosdx_mentalhealthcounseling___0", "pcosdx_mentalhealthcounseling___unk",
  "pcosdx_mentalhealthcounseling___na", "pcosdx_mentalhealthcounseling___oth",
  "pcosdx_mentalhealthcounseling___pm", "pcosdx_dietarycounseling___1",
  "pcosdx_dietarycounseling___2", "pcosdx_dietarycounseling___0",
  "pcosdx_dietarycounseling___unk", "pcosdx_dietarycounseling___na",
  "pcosdx_dietarycounseling___oth", "pcosdx_dietarycounseling___pm",
  "pcosdx_exercising_teachingplan___1", "pcosdx_exercising_teachingplan___2",
  "pcosdx_exercising_teachingplan___0", "pcosdx_exercising_teachingplan___unk",
  "pcosdx_exercising_teachingplan___na", "pcosdx_exercising_teachingplan___oth",
  "pcosdx_exercising_teachingplan___pm", "pcosdx_age",
  "pcosdx_irregular_menses", "pcosdx_hyperandrogenism", "cv_type_of_clinic",
  "cv_mdc_specialties___1", "cv_mdc_specialties___2", "cv_mdc_specialties___3",
  "cv_mdc_specialties___4", "cv_mdc_specialties___5", "cv_mdc_specialties___6",
  "cv_mdc_specialties___60", "cv_mdc_specialties___unk",
  "cv_mdc_specialties___na", "cv_mdc_specialties___oth",
  "cv_mdc_specialties___pm", "cv_monthssincepcosdx", "cv_newdx___3",
  "cv_newdx___4", "cv_newdx___5", "cv_newdx___6", "cv_newdx___29",
  "cv_newdx___7", "cv_newdx___8", "cv_newdx___9", "cv_newdx___10",
  "cv_newdx___11", "cv_newdx___12", "cv_newdx___22", "cv_newdx___13",
  "cv_newdx___14", "cv_newdx___15", "cv_newdx___16", "cv_newdx___27",
  "cv_newdx___28", "cv_newdx___17", "cv_newdx___18", "cv_newdx___19",
  "cv_newdx___20", "cv_newdx___21", "cv_newdx___23", "cv_newdx___24",
  "cv_newdx___25", "cv_newdx___26", "cv_newdx___60", "cv_newdx___0",
  "cv_newdx___unk", "cv_newdx___na", "cv_newdx___oth", "cv_newdx___pm",
  "cv_cysticacnedx_age", "cv_obesitydx_age", "cv_anxietydx_age",
  "cv_depressiondx_age", "cv_bingedx_age", "cv_restrictdx_age", "cv_adhddx_age",
  "cv_mh_hospital_age", "cv_si_attempt_age", "cv_othermedical",
  "cv_other_age_2", "cv_medications___1", "cv_medications___2",
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
  "cv_meds_diabetes2", "cv_weight", "cv_height", "cv_bmi", "cv_bmi_percentile",
  "cv_bmi_z", "cv_sbp", "cv_dbp", "cv_genderid", "cv_hirsutism_num",
  "cv_hirsutism_cat", "cv_acneface", "cv_acneother___1", "cv_acneother___2",
  "cv_acneother___0", "cv_acneother___unk", "cv_acneother___na",
  "cv_acneother___oth", "cv_acneother___pm", "cv_acanthosisneck",
  "cv_hydradinitis", "cv_alopecia", "cv_ft", "cv_ft_uln", "cv_ft_perc", "cv_tt",
  "cv_tt_assay", "cv_tt_uln", "cv_tt_perc", "cv_a1c", "cv_mood", "cv_phq2",
  "cv_phq8", "cv_phq9", "cv_cesd20", "cv_newmeds___1", "cv_newmeds___2",
  "cv_newmeds___3", "cv_newmeds___4", "cv_newmeds___5", "cv_newmeds___6",
  "cv_newmeds___7", "cv_newmeds___8", "cv_newmeds___9", "cv_newmeds___10",
  "cv_newmeds___11", "cv_newmeds___12", "cv_newmeds___13", "cv_newmeds___14",
  "cv_newmeds___15", "cv_newmeds___16", "cv_newmeds___17", "cv_newmeds___18",
  "cv_newmeds___19", "cv_newmeds___20", "cv_newmeds___21", "cv_newmeds___22",
  "cv_newmeds___23", "cv_newmeds___32", "cv_newmeds___24", "cv_newmeds___25",
  "cv_newmeds___26", "cv_newmeds___27", "cv_newmeds___28", "cv_newmeds___29",
  "cv_newmeds___30", "cv_newmeds___31", "cv_newmeds___60", "cv_newmeds___na",
  "cv_newmeds___unk", "cv_newmeds___oth", "cv_newmeds___pm",
  "cv_newmeds_diabetes", "cv_newcont_metforminformulation",
  "cv_newcont_metformindose", "cv_newcont_metformin_otherdose",
  "cv_newcont_metforminfreq", "cv_referrals___1", "cv_referrals___2",
  "cv_referrals___3", "cv_referrals___4", "cv_referrals___5",
  "cv_referrals___6", "cv_referrals___7", "cv_referrals___8",
  "cv_referrals___9", "cv_referrals___10", "cv_referrals___11",
  "cv_referrals___12", "cv_referrals___60", "cv_referrals___na",
  "cv_referrals___unk", "cv_referrals___oth", "cv_referrals___pm",
  "cv_referrals_reason", "cv_dietarycounseling", "cv_exerciseplan"
)
# LARC for BCH
larc_vars <- c(
  "cv_type_of_clinic", "cv_mdc_specialties___1", "cv_mdc_specialties___2",
  "cv_mdc_specialties___3", "cv_mdc_specialties___4", "cv_mdc_specialties___5",
  "cv_mdc_specialties___6", "cv_mdc_specialties___60",
  "cv_medications___10",
  "cv_medications___11", "cv_weight", "cv_bmi", "cv_acneface",
  "cv_hirsutism_num", "cv_hirsutism_cat", "cv_menfreq", "race___1", "race___2",
  "race___3", "race___4", "race___5", "race___60", "race___unk",
  "ethnicity", "insur_type", "cv_genderid",
  "cv_sexuallyactive", "cv_sexualpref",
  "cv_monthssincepcosdx", "cv_hgb", "cv_ldl", "cv_hdl", "cv_sbp", "cv_dbp"
)
# SES variables for health disparities
ses_vars <- c(
  "pcosdx_birthweight_calc", "cv_bmi",
  "cv_bmi_percentile", "cv_bmi_z", "cv_a1c", "cv_tt_perc", "cv_ft_perc",
  "cv_shbg", "cv_fbg", "cv_2hrglucoseogtt", "cv_tg_fasting", "cv_tg", "cv_hdl",
  "cv_ldl", "cv_tc", "cv_alt", "cv_ast", "cv_hirsutism_num", "cv_hirsutism_cat",
  "cv_acneface", "cv_acneother___1", "cv_acneother___2", "cv_acneother___0",
  "cv_acneother___unk", "pcosdx_menarche",
  "pcosdx_pmh___1", "pcosdx_pmh___2", "pcosdx_pmh___3", "pcosdx_pmh___4",
  "pcosdx_pmh___5", "pcosdx_pmh___6", "pcosdx_pmh___29", "pcosdx_pmh___7",
  "pcosdx_pmh___8", "pcosdx_pmh___9", "pcosdx_pmh___10", "pcosdx_pmh___11",
  "pcosdx_pmh___12", "pcosdx_pmh___22", "pcosdx_pmh___13", "pcosdx_pmh___14",
  "pcosdx_pmh___15", "pcosdx_pmh___16", "pcosdx_pmh___23", "pcosdx_pmh___24",
  "pcosdx_pmh___17", "pcosdx_pmh___18", "pcosdx_pmh___19", "pcosdx_pmh___20",
  "pcosdx_pmh___21", "pcosdx_pmh___60", "pcosdx_pmh___0", "pcosdx_pmh___unk",
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
  "pcosdx_famhx___27", "pcosdx_famhx___28",
  "pcosdx_famhx___30", "pcosdx_famhx___5", "pcosdx_famhx___29",
  "pcosdx_famhx___7", "pcosdx_famhx___24", "pcosdx_famhx___8",
  "pcosdx_famhx___25", "pcosdx_famhx___12", "pcosdx_famhx___22",
  "pcosdx_famhx___13", "pcosdx_famhx___14", "pcosdx_famhx___26",
  "pcosdx_famhx___15", "pcosdx_famhx___16", "pcosdx_famhx___17",
  "pcosdx_famhx___18", "pcosdx_famhx___31", "pcosdx_famhx___32",
  "pcosdx_famhx___33", "pcosdx_famhx___60", "pcosdx_famhx___0",
  "pcosdx_famhx___unk",
  "pastmeds___1", "pastmeds___2", "pastmeds___3",
  "pastmeds___4", "pastmeds___5", "pastmeds___6", "pastmeds___7",
  "pastmeds___8", "pastmeds___9", "pastmeds___10", "pastmeds___11",
  "pastmeds___12", "pastmeds___13", "pastmeds___14", "pastmeds___15",
  "pastmeds___16", "pastmeds___60", "pastmeds___0", "pastmeds___unk",
  "pcosdx_mentalhealthcounseling___1",
  "pcosdx_mentalhealthcounseling___2", "pcosdx_mentalhealthcounseling___3",
  "pcosdx_mentalhealthcounseling___0", "pcosdx_mentalhealthcounseling___unk",
  "pcosdx_age", "pcosdx_specialty",
  "cv_visittype", "cv_type_of_clinic",
  "cv_mdc_specialties___1", "cv_mdc_specialties___2", "cv_mdc_specialties___3",
  "cv_mdc_specialties___4", "cv_mdc_specialties___5", "cv_mdc_specialties___6",
  "cv_mdc_specialties___60",
  "cv_age", "cv_monthssincepcosdx",
  "cv_newdx___3", "cv_newdx___4", "cv_newdx___5", "cv_newdx___6",
  "cv_newdx___29", "cv_newdx___7", "cv_newdx___8", "cv_newdx___9",
  "cv_newdx___10", "cv_newdx___11", "cv_newdx___12", "cv_newdx___22",
  "cv_newdx___13", "cv_newdx___14", "cv_newdx___15", "cv_newdx___16",
  "cv_newdx___27", "cv_newdx___28", "cv_newdx___17", "cv_newdx___18",
  "cv_newdx___19", "cv_newdx___20", "cv_newdx___21", "cv_newdx___23",
  "cv_newdx___24", "cv_newdx___25", "cv_newdx___26", "cv_newdx___60",
  "cv_newdx___0", "cv_newdx___unk",
  "cv_medications___1", "cv_medications___2",
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
  "cv_newmeds___1", "cv_newmeds___2",
  "cv_newmeds___3", "cv_newmeds___4", "cv_newmeds___5", "cv_newmeds___6",
  "cv_newmeds___7", "cv_newmeds___8", "cv_newmeds___9", "cv_newmeds___10",
  "cv_newmeds___11", "cv_newmeds___12", "cv_newmeds___13", "cv_newmeds___14",
  "cv_newmeds___15", "cv_newmeds___16", "cv_newmeds___17", "cv_newmeds___18",
  "cv_newmeds___19", "cv_newmeds___20", "cv_newmeds___21", "cv_newmeds___22",
  "cv_newmeds___23", "cv_newmeds___32", "cv_newmeds___24", "cv_newmeds___25",
  "cv_newmeds___26", "cv_newmeds___27", "cv_newmeds___28", "cv_newmeds___29",
  "cv_newmeds___30", "cv_newmeds___31", "cv_newmeds___60",
  "cv_referrals___1", "cv_referrals___2",
  "cv_referrals___3", "cv_referrals___4", "cv_referrals___5",
  "cv_referrals___6", "cv_referrals___7", "cv_referrals___8",
  "cv_referrals___9", "cv_referrals___10", "cv_referrals___11",
  "cv_referrals___12", "cv_referrals___60",
  "cv_referrals___unk"
)
# Weight categories
df <- df %>%
  mutate(
    weight_perc = case_when(
      !is.na(cv_bmi_percentile) & cv_bmi_percentile < 85 ~ "Normal weight",
      !is.na(cv_bmi_percentile) & cv_bmi_percentile >= 85 &
        cv_bmi_percentile < 95 ~ "Overweight",
      !is.na(cv_bmi_percentile) & cv_bmi_percentile >= 95 ~ "Obese",
      !is.na(cv_bmi) & cv_bmi < 25 ~ "Normal weight",
      !is.na(cv_bmi) & cv_bmi >= 25 & cv_bmi < 30 ~ "Overweight",
      !is.na(cv_bmi) & cv_bmi >= 30 ~ "Obese",
      .default = NA
    ),
    weight_raw = cut(cv_bmi,
      breaks = c(-Inf, 25, 30, Inf), right = F,
      labels = c("Normal weight", "Overweight", "Obese")
    )
  )
df$weight_perc <- factor(df$weight_perc,
  levels = c("Normal weight", "Overweight", "Obese")
)
df$weight_raw <- factor(df$weight_raw,
  levels = c("Normal weight", "Overweight", "Obese")
)
# Combined weight category
df$weight_cat <- coalesce(df$weight_perc, df$weight_raw)
# Overweight yes no categories
df <- df %>%
  mutate(
    overweight_perc = case_when(
      !is.na(cv_bmi_percentile) & cv_bmi_percentile >= 85 ~ "Yes",
      !is.na(cv_bmi_percentile) & cv_bmi_percentile < 85 ~ "No",
      is.na(cv_bmi_percentile) & cv_bmi >= 25 ~ "Yes",
      is.na(cv_bmi_percentile) & cv_bmi < 25 ~ "No",
      .default = NA
    ),
    overweight_raw = ifelse(cv_bmi >= 25, "Yes", "No")
  )
df$overweight_perc <- factor(df$overweight_perc, levels = c("No", "Yes"))
df$overweight_raw <- factor(df$overweight_raw, levels = c("No", "Yes"))
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
# Combined family history - luckily the column numbers for family and non-family
# match up
hist_labels <- sub(
  "Non-parent family history \\(include all history, not only from the diagnostic visit\\)\\. Including siblings, first cousins, aunts, uncles, and grandparents. \\(choice\\=",
  "", as.character(label(df[, grep("pcosdx_famhx___", colnames(df))]))
)
hist_labels <- sub("\\)", "", hist_labels)
hist_labels <- sub(" \\{pcosdx_famhx_family_otherdx\\}", "", hist_labels)
hist_labels <- paste0("Any family history of ", hist_labels, "?")
hist_vars <- sub("pcosdx_famhx___", "", colnames(df[, grep(
  "pcosdx_famhx___", colnames(df)
)]))
any_hist <- data.frame(lapply(hist_vars, function(v) {
  parent <- df[, paste0("pcosdx_famhx_parent___", v)]
  nonparent <- df[, paste0("pcosdx_famhx___", v)]
  any <- (parent == "Checked") + (nonparent == "Checked")
  any <- factor(any, levels = 0:2, labels = c("No", "Yes", "Yes"))
  return(any)
}))
colnames(any_hist) <- paste0("pcosdx_any_famhx___", hist_vars)
df <- cbind(df, any_hist)
# Mental health screening
df$mental_health_screening <-
  rowSums(!is.na(df[, c("cv_phq2", "cv_phq8", "cv_phq9", "cv_cesd20")]))
df$mental_health_screening <- factor(df$mental_health_screening,
  levels = 0:4, labels = c("No", "Yes", "Yes", "Yes", "Yes")
)
# Combine two earliest obesity noticed columns
df$pcosdx_obesitydx_age_combined <-
  coalesce(df$pcosdx_obesitydx_age, df$pcosdx_obesitydx_earliestage)
label(df$pcosdx_obesitydx_age_combined) <-
  "Earliest age that overweight/obesity was noted?"
# Categorize
df$pcosdx_obesitydx_age_combined_cat <- cut(
  df$pcosdx_obesitydx_age_combined,
  c(-Inf, 5, 10, 15, Inf),
  labels = c("<= 5", ">5 and <= 10", ">10 and <= 15", ">= 15")
)
label(df$pcosdx_obesitydx_age_combined_cat) <-
  "Earliest age group that overweight/obesity was noted?"
aim1_vars <- c(
  aim1_vars, "pcosdx_obesitydx_age_combined",
  "pcosdx_obesitydx_age_combined_cat"
)
# Region
# df$Region <- df$site
# levels(df$Region) <- c(
#   "West", "Northeast", "Southeast", "Midwest", "West", "Midwest", "Northeast",
#   "Midwest", "Northeast", "Southeast", "West", "Southeast", "Southeast",
#   "Midwest", "Southeast", "Midwest", "Southwest", "Southwest"
# )
# Mental health diagnosis variables
# Convert to numeric
df$cv_phq9[df$cv_phq9 %in%
  c("performed, score not included", "performed, results not listed")] <- NA
df$cv_phq9[df$cv_phq9 == "negative"] <- 0
df$cv_cesd20[df$cv_cesd20 == "AN"] <- NA
df$cv_phq9 <- as.numeric(df$cv_phq9)
df$cv_cesd20 <- as.numeric(df$cv_cesd20)
df <- df %>%
  mutate(
    depression = cv_newdx___16 == "Checked" | pcosdx_pmh___16 == "Checked" |
      cv_phq2 >= 3 | cv_phq8 >= 10 | cv_phq9 >= 10 | cv_cesd20 >= 16,
    anxiety =
      factor(cv_newdx___15 == "Checked" | pcosdx_pmh___15 == "Checked",
        levels = c(F, T), labels = c("No", "Yes")
      ),
    bed =
      factor(cv_newdx___17 == "Checked" | pcosdx_pmh___17 == "Checked",
        levels = c(F, T), labels = c("No", "Yes")
      ),
    red =
      factor(cv_newdx___18 == "Checked" | pcosdx_pmh___18 == "Checked",
        levels = c(F, T), labels = c("No", "Yes")
      ),
    adhd =
      factor(cv_newdx___19 == "Checked" | pcosdx_pmh___19 == "Checked",
        levels = c(F, T), labels = c("No", "Yes")
      )
  )
df$depression[is.na(df$depression)] <- F
df$depression <-
  factor(df$depression, levels = c(F, T), labels = c("No", "Yes"))
# Any mental health counseling
df$pcosdx_any_mentalhealthcounseling <-
  rowSums(df[, c(paste0("pcosdx_mentalhealthcounseling___", 1:3))] == "Checked")
df$pcosdx_any_mentalhealthcounseling <-
  factor(df$pcosdx_any_mentalhealthcounseling,
    levels = 0:3,
    labels = c("No", "Yes", "Yes", "Yes")
  )
# Labels
label(df$depression) <- "Depression?"
label(df$anxiety) <- "Anxiety?"
label(df$bed) <- "Binge Eating Disorder?"
label(df$red) <- "Restrictive Eating Disorder?"
label(df$adhd) <- "ADHD?"
label(df$pcosdx_any_mentalhealthcounseling) <-
  "Any mental health counseling in the 12 months prior to diagnosis?"
# LARC
df$larc <- df$cv_medications___10 == "Checked" |
  df$cv_medications___11 == "Checked"
df$larc <- factor(df$larc, levels = c(F, T), labels = c("No", "Yes"))
# Age at LARC start, etc.
df <- df %>%
  group_by(record_number) %>%
  mutate(
    age_first_larc = first(na.omit(cv_age[larc == "Yes"])),
    time_first_larc = first(na.omit(cv_monthssincepcosdx[larc == "Yes"]))
  ) %>%
  ungroup() %>%
  mutate(
    time_pcos_to_first_larc = time_first_larc - cv_monthssincepcosdx,
    time_menarche_to_first_larc = age_first_larc - pcosdx_menarche
  )
label(df$age_first_larc) <- "Age at First LARC"
label(df$time_pcos_to_first_larc) <- "Years From PCOS Dx to First LARC"
label(df$time_menarche_to_first_larc) <- "Years From Menarche to First LARC"
# Ever on LARC vs. never
larc_users <- unique(df$record_number[df$larc == "Yes"])
df$larc_ever <- "No"
df$larc_ever[df$record_number %in% larc_users] <- "Yes"
df$larc_ever <- factor(df$larc_ever, levels = c("No", "Yes"))
# Same as LARC, but for estrogen-containing medications (EC)
df$ec <- df$cv_medications___5 == "Checked" |
  df$cv_medications___6 == "Checked" | df$cv_medications___7 == "Checked"
df$ec <- factor(df$ec, levels = c(F, T), labels = c("No", "Yes"))
df <- df %>%
  group_by(record_number) %>%
  mutate(
    age_first_ec = first(na.omit(cv_age[ec == "Yes"])),
    time_first_ec = first(na.omit(cv_monthssincepcosdx[ec == "Yes"]))
  ) %>%
  ungroup() %>%
  mutate(
    time_pcos_to_first_ec = time_first_ec - cv_monthssincepcosdx,
    time_menarche_to_first_ec = age_first_ec - pcosdx_menarche
  )
label(df$age_first_ec) <- "Age at First EC"
label(df$time_pcos_to_first_ec) <- "Years From PCOS Dx to First EC"
label(df$time_menarche_to_first_ec) <- "Years From Menarche to First EC"
ec_users <- unique(df$record_number[df$ec == "Yes"])
df$ec_ever <- "EC Never"
df$ec_ever[df$record_number %in% ec_users] <- "EC Ever"
df$ec_ever <- factor(df$ec_ever, levels = c("EC Never", "EC Ever"))
# Same as LARC and EC, but for Metformin
df$metformin <- df$cv_newmeds___1 == "Checked"
df$metformin <- factor(df$metformin, levels = c(F, T), labels = c("No", "Yes"))
df <- df %>%
  group_by(record_number) %>%
  mutate(
    age_first_metformin = first(na.omit(cv_age[metformin == "Yes"])),
    time_first_metformin = first(na.omit(cv_monthssincepcosdx[metformin == "Yes"]))
  ) %>%
  ungroup() %>%
  mutate(
    time_pcos_to_first_metformin = age_first_metformin - pcosdx_age,
    time_menarche_to_first_metformin = time_first_metformin - cv_monthssincepcosdx
  )
label(df$age_first_metformin) <- "Age at First Metformin"
label(df$time_pcos_to_first_metformin) <-
  "Years From PCOS Dx to First Metformin"
label(df$time_menarche_to_first_metformin) <-
  "Years From Menarche to First Metformin"
metformin_users <- unique(df$record_number[df$metformin == "Yes"])
df$metformin_ever <- "Metformin Never"
df$metformin_ever[df$record_number %in% metformin_users] <- "Metformin Ever"
df$metformin_ever <- factor(df$metformin_ever,
  levels = c("Metformin Ever", "Metformin Never")
)
# Lifestyle medicine means not on Metformin or EC
df$lifestyle <- df$metformin == "No" & df$ec == "No"
df$lifestyle <- factor(df$lifestyle, levels = c(F, T), labels = c("No", "Yes"))
df <- df %>%
  group_by(record_number) %>%
  mutate(
    age_first_lifestyle = first(na.omit(cv_age[lifestyle == "Yes"])),
    time_first_lifestyle = first(na.omit(cv_monthssincepcosdx[lifestyle == "Yes"]))
  ) %>%
  ungroup() %>%
  mutate(
    time_pcos_to_first_lifestyle = age_first_lifestyle - pcosdx_age,
    time_menarche_to_first_lifestyle = age_first_lifestyle - pcosdx_menarche
  )
label(df$age_first_lifestyle) <- "Age at First Lifestyle Medicine"
label(df$time_pcos_to_first_lifestyle) <-
  "Years From PCOS Dx to First Lifestyle Medicine"
label(df$time_menarche_to_first_lifestyle) <-
  "Years From Menarche to First Lifestyle Medicine"
lifestyle_users <- unique(df$record_number[df$lifestyle == "Yes"])
df$lifestyle_ever <- "Lifestyle Never"
df$lifestyle_ever[df$record_number %in% lifestyle_users] <- "Lifestyle Ever"
df$lifestyle_ever <- factor(df$lifestyle_ever,
  levels = c("Lifestyle Never", "Lifestyle Ever")
)
# Age group at diagnosis
df$age_group <- cut(df$pcosdx_age, c(-Inf, 15, Inf),
  right = F, labels = c("< 15 years", ">= 15 years")
)
# Convert columns to numeric
df$cv_a1c <- suppressWarnings(as.numeric(df$cv_a1c))
df$redcap_repeat_instance <- as.numeric(df$redcap_repeat_instance)
df$pcosdx_age <- as.numeric(df$pcosdx_age)
# Fix/add labels
label(df$insur_type) <- "Insurance Type"
# label(df$Region) <- "Region"
label(df$ethnicity) <- "Ethnicity"
label(df$pcosdx_age) <- "Age at time of PCOS diagnosis"
label(df$age_group) <- "Age Group at Diagnosis"
label(df$combined_race) <- "Race"
label(df$weight_perc) <- "Weight Category (by percentile and raw value)"
label(df$weight_raw) <- "Weight Category (by raw value)"
label(df$weight_cat) <- "Weight Category (by percentile then raw value)"
label(df$overweight_perc) <- "Overweight Status (by percentile and raw value)"
label(df$overweight_raw) <- "Overweight Status (by raw value)"
label(df$cv_a1c) <- "HbA1C"
label(df$larc) <- "On LARC?"
label(df$larc_ever) <- "LARC Ever Used?"
label(df$ec) <- "On EC?"
label(df$ec_ever) <- "EC Ever Used?"
label(df$metformin) <- "On Metformin?"
label(df$metformin_ever) <- "Metformin Ever Used?"
label(df$lifestyle) <- "On Lifestyle Medicine?"
label(df$lifestyle_ever) <- "Lifestyle Medicine Ever Used?"
label(df$pcosdx_age) <- "Age at time of PCOS diagnosis"
label(df[, paste0("pcosdx_any_famhx___", hist_vars)]) <- as.list(hist_labels)
label(df$mental_health_screening) <-
  "PHQ-2, PHQ-8, PHQ-9 or CED-S Score Available?"
label(df[, grep("___unk", colnames(df))]) <-
  as.list(sub(
    "choice=NA", "choice=Unknown/Not recorded",
    label(df[, grep("___unk", colnames(df))])
  ))
# Drop unused levels
l <- label(df)
df <- droplevels(df)
label(df) <- as.list(l)
# Save
save(df, demo_vars, ses_vars, file = "./Data_Clean/analysis_data.RData")
