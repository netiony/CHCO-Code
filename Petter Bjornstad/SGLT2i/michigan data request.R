library(readxl)
library(dplyr)

# read file with needed IDs
ids <- readxl::read_xlsx("E:\\Petter Bjornstad\\T2D scRNA SGLT2i\\CHCO_clinical_data_needed-9-15_pbedits.xlsx")
ids <- ids %>% filter(`Has clinical info? ( as of 9-15)` != "Y")
ids <- ids %>% select(`Study ID`,`Has clinical info? ( as of 9-15)`)

# read in harmonized data
harmonized <- read.csv("E:\\Petter Bjornstad\\Data Harmonization\\Data Clean\\merged_dataset_2022-09-26.csv")
# need to fix the IDs to match the Michigan request
harmonized$subject_id<- ifelse(harmonized$study=="CROCODILE" & nchar(harmonized$subject_id)<2,
                          paste0("CRC-0",harmonized$subject_id),
                          ifelse(harmonized$study=="CROCODILE" & nchar(harmonized$subject_id)==2,
                          paste0("CRC-",harmonized$subject_id),harmonized$subject_id))

# read in old spreadsheet to get the variables we need
old <- read.csv("E:\\Petter Bjornstad\\T2D scRNA SGLT2i\\Data clean\\data_set_without_hpi_with_pgloabs.csv")

# get the harmonized data for the IDs we need
harmonized$`Study ID` <- harmonized$subject_id
harmonized_keep <- merge(harmonized,ids,by="Study ID",all.x = F, all.y = T)

