library(openxlsx)
library(dplyr)
library(tidyverse)
dir.results <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/DECODE")

#Read in Biopsy Tracker dataset
bmt <- read.xlsx("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Data Harmonization/Biopsies/Biopsy master tracker.xlsx")
colnames(bmt) <- bmt[2,]
bmt <- bmt[-c(1,2),]

bmt_sum <- bmt %>% 
  group_by(Study, `Shipped Y/N`, `scRNA status`,`Visit ID`) %>%
  summarise(Count = n()) %>%
  ungroup

write.csv(bmt_sum,fs::path(dir.results,"Biopsies and Sequencing","Biopsies_Sequencing_by_Study_Visit.csv"))
