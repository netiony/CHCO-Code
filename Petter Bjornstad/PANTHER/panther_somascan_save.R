library(SomaDataIO)
library(dplyr)

# read in data - we want to use the fully processed, normalized file ending in "anmlSMP.adat"
soma <- read_adat("/Volumes/Peds Endo/Petter Bjornstad/SOMAScan/PANTHER/20240126_597_Bjornstad_SOMAscan7k_WUS-24-002_data_export/WUS-24-002_2024-01-26_Somalogic_standardized_files/WUS_24_002_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")
analytes <- getAnalyteInfo(soma)

# remove fc mouse and no protein
drop <- analytes %>% filter(Target == "Fc_MOUSE" | Target == "No Protein" | !(Organism == "Human") | !(Type == "Protein"))
apt_drop <- drop$AptName
soma <- soma %>% dplyr::select(!all_of(apt_drop))
analytes <- analytes %>% filter(!Target == "Fc_MOUSE")
analytes <- analytes %>% filter(!Target == "No Protein")
analytes <- analytes %>% filter(Organism == "Human")
analytes <- analytes %>% filter(Type == "Protein")

# samples without a SampleDescription are QC samples
soma <- soma %>% filter(!is.na(SampleDescription))

# save
save(soma,file = "/Volumes/Peds Endo/Petter Bjornstad/SOMAScan/PANTHER/panther_bl_soma.Rdata")
save(analytes,file = "/Volumes/Peds Endo/Petter Bjornstad/SOMAScan/PANTHER/panther_bl_analytes.Rdata")


# non standardized soma
soma_ns <- read.delim("/Volumes/Peds Endo/Petter Bjornstad/SOMAScan/PANTHER/20240126_597_Bjornstad_SOMAscan7k_WUS-24-002_data_export/WUS-24-002_2024-01-26_GTAC_analyzed_files/customer.samples_log2RFU_WUS_24_002_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_quantiled.txt")
