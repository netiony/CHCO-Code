---
  title: "Subset of CROCODILE for Kumar (FFA and Insulin)"
author: "Ye Ji Choi"
date: "`r format(Sys.time(), '%d %B %Y')`"
---
  
library(dplyr)
library(tidyr)
library(magrittr)
library(tableone)
library(knitr)
library(ISLR)
library(stringr)
library(Hmisc)
library(emmeans)
library(car)

home_dir = "S:/Petter Bjornstad/Data Harmonization/Data Clean"
setwd(home_dir)

# read in raw data
dat <- read.csv("harmonized_dataset.csv", na.strings = c(" ", "", "-99"))

croc <- dat %>%
  select(record_id, study, age, sex, race, ethnicity, group, diabetes_duration, 
                    baseline_ffa, ffa_minus_20, ffa_minus_10, ffa_minus_5, ffa_0, ffa_2, ffa_4, ffa_6, ffa_8, ffa_10, 
                    ffa_30, ffa_60, ffa_70, ffa_80, ffa_90, ffa_120, ffa_180, ffa_220, ffa_230, ffa_240, ffa_250, ffa_260, 
                    ffa_270, steady_state_ffa, ffa_suppression, insulin_minus_20, insulin_minus_10, insulin_minus_5, insulin_0, 
                    insulin_2, insulin_4, insulin_6, insulin_8, insulin_10, insulin_20, insulin_30, insulin_45, insulin_60, 
                    insulin_70, insulin_80, insulin_90, insulin_120, insulin_150, insulin_180, insulin_210, insulin_220, 
                    insulin_230, insulin_240, insulin_245, insulin_249, insulin_250, insulin_252, insulin_253, insulin_254, 
                    insulin_255, insulin_260, insulin_270, steady_state_insulin, insulin_amt, insulin_labs, insulin_units, 
                    long_insulin_dose, long_insulin_time, morning_insulin, p1_gc_leanm, p1_gc_m, p1_raw_leanm, p1_raw_m, 
                    p2_gc_leanm, p2_gc_m, p2_raw_leanm, p2_raw_m) %>%
  group_by(record_id)%>%
  fill(names(croc)) %>%
  summarise_all(last)

croc_clean <- croc %>%
  select_if(~ !all(is.na(.)))

croc_summary <- croc_clean %>%
  ungroup() %>%
  group_by(study) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = T))) %>% 
  filter(study != "CASPER" & study != "COFFEE") %>%
  select(study, insulin_minus_20, insulin_minus_10, insulin_minus_5, insulin_0, ffa_minus_20, ffa_minus_10, ffa_minus_5, ffa_0)

write.csv(croc_clean, "./Data Exports/data_for_kumar_crocodile_insulin_ffa_011123.csv", row.names=F)