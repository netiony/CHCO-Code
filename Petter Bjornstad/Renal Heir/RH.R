library(dplyr)
library(tidyr)
library(tableone)

dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))

dat <- dat %>%
  filter(study == "RENAL-HEIR") %>%
  group_by(record_id, visit) %>%
  summarise(across(where(is.character),~last(na.omit(.x))),
            across(where(is.factor),~last(na.omit(.x))),
            across(where(is.numeric),~mean(.x,na.rm = T)),.groups = "drop") %>%
  filter(participation_status=="Participated") %>%
  mutate(race_ethnicity_condensed = case_when(race == "White" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                              race == "Black or African American" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                              ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                              T ~ "Not Hispanic or Latino Other")) %>%
  rowwise() %>%
  mutate(m_i = raw_m / steady_state_insulin)
  
table_one <- CreateTableOne(vars = c("age", 
                                     "sex", 
                                     "race_ethnicity_condensed", 
                                     "weight", 
                                     "height", 
                                     "waistcm", 
                                     "hba1c", 
                                     "sbp", 
                                     "dbp", 
                                     "cholesterol", 
                                     "ldl", 
                                     "hdl", 
                                     "triglycerides",
                                     "eGFR_CKD_epi", 
                                     "acr_u", 
                                     "metformin_timepoint", 
                                     "insulin_med_timepoint", 
                                     "sglti_timepoint", 
                                     "tzd_timepoint", 
                                     "raasi_timepoint",
                                     "statin"
                                     ), 
                            strata = "group", data = dat)
print(table_one, nonnormal = c("acr_u"), noSpaces = T)

table_two <- CreateTableOne(vars = c("dexa_lean_kg", 
                                     "dexa_fat_kg", 
                                     "dexa_body_fat", 
                                     "dexa_trunk_kg", 
                                     "fasting_ffa", 
                                     "ffa_suppression", 
                                     "steady_state_insulin", 
                                     "steady_state_cpeptide", 
                                     "raw_m", 
                                     "m_i"
), 
strata = "group", data = dat)
print(table_two, noSpaces = T)
