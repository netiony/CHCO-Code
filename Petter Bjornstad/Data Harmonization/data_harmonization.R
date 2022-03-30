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
improve = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "IMPROVE"]
  )
)
penguin = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "PENGUIN"]
  )
)
crocodile = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "CROCODILE"]
  )
)
coffee = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "COFFEE"]
  )
)
casper = exportRecords(
  redcapConnection(
    url = uri,token = tokens$Token[tokens$Study == "CASPER"]
  )
)
###############################################################################
# Screening variables
###############################################################################
screen_vars = c("insulin","screen_height","screen_weight","screen_bmi",
                "screen_bmi_percentile","waist_circumference","hip_circumference",
                "sys_bp","dys_bp","map","pulse","activity_factor","schofield",
                "hba1c","hemoglobin","screen_hematocrit","screen_serum_creatinine",
                "screen_urine_mab","screen_urine_cre","screen_urine_acr",
                "screen_pregnant","screening_labs_complete")
# insulin
renalheir$insulin = renalheir$diabetes_med___2
levels(renalheir$insulin) = c("No","Yes")
# Insulin pump

# activity_factor
renalheir = renalheir %>% 
  unite(.,activity_factor,fem_activity_factor,male_activity_factor,na.rm = T) %>%
  mutate(schofield = schofield_female+schofield_male)
renalheir$schofield[renalheir$schofield == 0] = NA
