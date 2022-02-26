# Libraries
library(tidyverse)
library(redcapAPI)
# Import api tokens

# RENAL HEIR
renal_heir = read.csv("/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Vigers/Tim and Casey/RENALHEIR_DATA_2022-02-18_1154.csv")
renal_heir = renal_heir %>% 
  select(mr_number,dob,diagnosis,gender,scr) %>%
  rename(mrn = mr_number,sex = gender,diagnosis_date = diagnosis)
renal_heir$sex = factor(renal_heir$sex,levels = 0:2,labels = c("Male","Female","Other"))

