library(dplyr)
library(Hmisc)
library(purrr)
library(Seurat)

so <- readRDS("/run/user/778527649/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/scRNA/data_raw/PB_90samples_Harmony_rpca_Fadhl_PhilApproved_091024.RDS")

dat <- read.csv("/run/user/778527649/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv")
dat[dat == ""] <- NA
dat_subset <- dat %>%
  dplyr::select(mrn, kit_id, date, procedure, cryostor_id, visit, record_id, group, age, sex, race, epic_glp1ra_1, epic_ever_glp1ra_1, epic_sglti2_1, epic_ever_sglt2i_1)%>%
  dplyr::mutate(kit_id = sub("KI|kl", "KL", kit_id)) %>%
  group_by(record_id) %>%
  fill(everything(), .direction = "downup") %>%
  filter(!is.na(kit_id)) %>%
  filter(visit == "baseline") %>%
  filter(procedure == "kidney_biopsy")

so_meta <- so@meta.data %>%
  dplyr::mutate(kit_id = sub("KI|kl", "KL", kit_id))
so_meta_subset <- so_meta %>%
  distinct(kit_id, .keep_all = T) %>%
  dplyr::mutate(kit_id = sub("KI|kl", "KL", kit_id))

so_meta_combined <- left_join(so_meta, dat_subset)
row.names(so_meta_combined) <- rownames(so_meta)
write.csv(so_meta_combined, "/run/user/778527649/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/PB90_clinical_metadata.csv")
save(so_meta_combined, file = "/run/user/778527649/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Exports/PB90_clinical_metadata.RData")
