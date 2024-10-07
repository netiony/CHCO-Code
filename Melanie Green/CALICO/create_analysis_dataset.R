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
  "cv_age", "cv_bmi", "cv_bmi_percentile", "cv_bmi_z",
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
  mutate(Obesity = case_when(
    !is.na(cv_bmi_percentile) & cv_bmi_percentile >= 95 ~ "Yes",
    !is.na(cv_bmi_percentile) & cv_bmi_percentile < 95 ~ "No",
    is.na(cv_bmi_percentile) & cv_bmi >= 30 ~ "Yes",
    is.na(cv_bmi_percentile) & cv_bmi < 30 ~ "Yes",
    .default = NA
  ))
df$Obesity <- factor(df$Obesity, levels = c("No", "Yes"))
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
label(df$Obesity) <- "Obesity Status"
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
save(df, demo_vars, aim1_vars, larc_vars, file = "./Data_Clean/analysis_data.RData")
