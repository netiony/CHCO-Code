library(readxl)
# Import sample list
needed = read_excel("/Volumes/som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/CHCO_clinical_data_needed-9-15_pbedits.xlsx")
# Get harmonized dataset
source("/Users/timvigers/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization/data_harmonization.R")
df = harmonize_data()
# Combine
new_df = left_join(needed,df,by=c("Study ID"="subject_id","Visit"="visit"))
# 
check = new_df %>% select(`Study ID`,Visit,study,group,sglt2i)
