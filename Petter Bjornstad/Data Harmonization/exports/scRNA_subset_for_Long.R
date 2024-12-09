library(dplyr)
library(Hmisc)
library(purrr)
st_ids <- read.csv("/run/user/1023/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Spatial transcriptomics/spatial_transcriptomics_ids.csv")
so <- readRDS("/run/user/1023/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/scRNA/data_clean/seurat_data_no_computations.RDS")

st_ids <- st_ids %>%
  mutate(visit = case_when(visit == "12_months_post_surgery" ~ "12M",
                           T ~ "BL"),
         michigan_id = case_when(grepl("IT_", record_id) ~ paste0(record_id,"_",visit),
                                 T ~ record_id))
st_so_sub <- subset(so, michigan_id %in% st_ids$michigan_id)
write_rds(st_so_sub, "/run/user/1023/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/scRNA/data_clean/seurat_data_spatial_match_subset.RDS")
