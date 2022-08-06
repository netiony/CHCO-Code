library(tidyverse)
library(SomaDataIO)
df = read_adat("/Users/timvigers/Documents/Work/Petter Bjornstad/Somalogic data/WUS-22-002_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")

df = df %>% filter(SampleType == "Sample") %>% select(SampleDescription,seq.10000.28:seq.9999.1)
croc = df %>% filter(str_detect(SampleDescription, "CRC"))
ckds = df %>% filter(str_detect(SampleDescription, "CKDS"))
improve = df %>% filter(str_detect(SampleDescription, "IT2D"))
renalhair = df %>% filter(str_detect(SampleDescription, "RH"))
nrow(croc) + nrow(ckds) + nrow(improve) + nrow(renalhair) == nrow(df)

write.csv(croc,"/Users/timvigers/Documents/Work/Petter Bjornstad/croc_proteomics.csv",row.names = F,na="")
write.csv(ckds,"/Users/timvigers/Documents/Work/Petter Bjornstad/ckds_proteomics.csv",row.names = F,na="")
write.csv(improve,"/Users/timvigers/Documents/Work/Petter Bjornstad/improve_proteomics.csv",row.names = F,na="")
write.csv(renalhair,"/Users/timvigers/Documents/Work/Petter Bjornstad/rh_proteomics.csv",row.names = F,na="")

t = getAnalyteInfo(df)
