library(readxl)
# Import sample list
needed = read_excel("/Volumes/som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/CHCO_clinical_data_needed-9-15_pbedits.xlsx")
# Get harmonized dataset
source("/Users/timvigers/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization/data_harmonization.R")
df = harmonize_data()
# Combine
new_df = left_join(needed,df,by=c("Study ID"="subject_id","Visit"="visit"))
# Unique identifier for Michigan
new_df$michigan_id = new_df$`Study ID`
new_df$michigan_id[new_df$study == "IMPROVE" & new_df$Visit == "Baseline"] = 
  paste0(new_df$michigan_id[new_df$study == "IMPROVE" & new_df$Visit == "Baseline"],"_BL")
new_df$michigan_id[new_df$study == "IMPROVE" & new_df$Visit == "12 Months Post-Surgery"] = 
  paste0(new_df$michigan_id[new_df$study == "IMPROVE" & new_df$Visit == "12 Months Post-Surgery"],"_12M")
# Michigan ID first
new_df = new_df %>% select(michigan_id,everything()) %>% select(-sglt2_inhibitors)
# Write data
write.csv(new_df,na = "",row.names = F,
          file = paste0("~/Dropbox/Sample_ID_Master_Folder_CHCO_Petter/biopsy_clinical_",
                        format(Sys.time(), "%Y_%m_%d_%H%M"),".csv"))