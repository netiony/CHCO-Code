---
title: "Proteomics and DKD"
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

home_dir = "S:/Petter Bjornstad/"
setwd(home_dir)

# read in raw data
dat <- read.csv("./Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-99"))
ids <- read.csv("./Proteomics and DKD/dkd_dataset_ids.csv", na.strings = c(" ", "", "-99"))

ids %<>% 
  mutate(visit = case_when(visit == "BL" ~ "baseline",
                           visit == "3M" ~ "3_months_post_surgery",
                           visit == "12M" ~ "12_months_post_surgery")) %>%
  mutate_all(funs(str_replace(., "IT2D-", "IT_")))

clean_dat <- dat %>%
  right_join(ids, by = c("record_id", "visit")) %>%
  select(record_id, visit, age, sex, bmi, hba1c, sbp, triglycerides, m_value, p1_gc_leanm, p1_gc_m, p1_raw_leanm, p1_raw_m, 
         p2_gc_leanm, p2_gc_m, p2_raw_leanm, p2_raw_m, eGFR_Zap, eGFR_CKD_epi, eGFR_fas_cr, eGFR_fas_cr_cysc, eGFR_Schwartz, eGFR_bedside_Schwartz,
         alb_base, gfr_bsa_plasma, gfr_ecv, gfr_raw_plasma, gfr_standard, gfr_ecv_percent, gfr_ecv_std, gfr_bsa_plasma_urine, gfr_raw_plasma_urine,
         acr_u, glom_nuc_count, glom_tuft_area, glom_volume_con, glom_volume_weibel, glom_volume_wiggins, gloms_gs, gloms) %>%
  select_if(~ !all(is.na(.))) %>%
  group_by(record_id, visit) %>%
  fill(names(clean_dat)) %>%
  summarise_all(last)


write.csv(clean_dat, "./Proteomics and DKD/Data_Clean/clean_dataset.csv", row.names=F)
