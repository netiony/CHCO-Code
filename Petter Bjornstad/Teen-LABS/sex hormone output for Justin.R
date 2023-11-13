load("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Teen Labs/Data_Cleaned/analysis_dataset.RData")

# keep only the proteins we need
analyte_sh <- analyte_info %>% filter(UniProt %in% c("P04278","P01215|P01222","P01215|P01225","P01215|P01229","P01236"))
apt_keep <- analyte_sh$AptName
df_justin <- df %>% select(ID, visit, all_of(apt_keep))

write.csv(df_justin, "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Teen Labs/soma_sex_hormone_Justin.csv",
          row.names = F)
write.csv(analyte_sh, "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Teen Labs/analyte_info.csv",
          row.names = F)