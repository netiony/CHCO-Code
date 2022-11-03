source("/Users/timvigers/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization/data_harmonization.R")
df = harmonize_data()
# Filter and select
mri = df %>% filter(subject_id %in% paste0("CRC-",c(10,12,13,18,20,24,28,31,33,34,35))) %>%
  select(subject_id,group,age_mri,gender,mri_bmi,contains("pasl"),contains("pcasl"),
         contains("adc"),volume_left,volume_right,contains("_bl_"),
         contains("fsoc"))
# Write CSV
write.csv(mri,na = "",row.names = F,
          file = paste0("~/Desktop/bold_mri_",
                        format(Sys.time(), "%Y_%m_%d_%H%M"),".csv"))
