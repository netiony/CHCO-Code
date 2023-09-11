library(readxl)

data <- read_xlsx("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Melanie Green/ORANGE/Data raw/ORANGE Multiomics for LP 3.18.2023/ORANGE Clinical Data - Targeted.xlsx")

sd(data$`Matsuda (10,000/sqrt(FPG*FPI*mean G*mean I)`, na.rm = T)
sd(data$`HOMA (FPGxFPI)/405`, na.rm = T)
sd(data$`liver fat %`, na.rm = T)
sd(data$`OGTT total testerone (ng/dL)`, na.rm = T)