library(dplyr)
library(arsenal)

# my files
cdc_hgt <- read.csv("/Volumes/Shared/Shared Projects/Laura/Peds Endo/PEDSNet/CAH/Data_raw/statage.csv")
cdc_hgt$source <- "CDC"
cdc_hgt$param <- "HEIGHTCM"
cdc_hgt <- cdc_hgt %>% select(-c(L, M, S))
cdc_hgt$age <- cdc_hgt$Agemos / 12
cdc_hgt$agemo <- cdc_hgt$Agemos
cdc_hgt <- cdc_hgt %>% filter(age >=2 & age<=18)
cdc_hgt <- cdc_hgt %>% rename(p03="P3", p05="P5", p10="P10", p25="P25", p50="P50", p75="P75", p90="P90", p95="P95", p97="P97")
cdc_hgt_males <- cdc_hgt %>% filter(Sex==1)
cdc_hgt_females <- cdc_hgt %>% filter(Sex==2)
# weight
cdc_wgt <- read.csv("/Volumes/Shared/Shared Projects/Laura/Peds Endo/PEDSNet/CAH/Data_raw/wtage.csv")
cdc_wgt$source <- "CDC"
cdc_wgt$param <- "WEIGHTKG"
cdc_wgt <- cdc_wgt %>% select(-c(L, M, S))
cdc_wgt$age <- cdc_wgt$Agemos / 12
cdc_wgt$agemo <- cdc_wgt$Agemos
cdc_wgt <- cdc_wgt %>% filter(age >=2 & age<=18)
cdc_wgt <- cdc_wgt %>% rename(p03="P3", p05="P5", p10="P10", p25="P25", p50="P50", p75="P75", p90="P90", p95="P95", p97="P97")
cdc_wgt_males <- cdc_wgt %>% filter(Sex==1)
cdc_wgt_females <- cdc_wgt %>% filter(Sex==2)
# BMI
cdc_bmi <- read.csv("/Volumes/Shared/Shared Projects/Laura/Peds Endo/PEDSNet/CAH/Data_raw/bmiagerev.csv")
cdc_bmi$source <- "CDC"
cdc_bmi$param <- "BMI"
cdc_bmi$Agemos <- as.numeric(cdc_bmi$Agemos)
cdc_bmi$age <- cdc_bmi$Agemos / 12
cdc_bmi <- cdc_bmi %>% filter(age >=2 & age<=18)
cdc_bmi$P3 <- as.numeric(cdc_bmi$P3)
cdc_bmi$P5 <- as.numeric(cdc_bmi$P5)
cdc_bmi$P10 <- as.numeric(cdc_bmi$P10)
cdc_bmi$P25 <- as.numeric(cdc_bmi$P25)
cdc_bmi$P50 <- as.numeric(cdc_bmi$P50)
cdc_bmi$P75 <- as.numeric(cdc_bmi$P75)
cdc_bmi$P90 <- as.numeric(cdc_bmi$P90)
cdc_bmi$P95 <- as.numeric(cdc_bmi$P95)
cdc_bmi$P97 <- as.numeric(cdc_bmi$P97)
cdc_bmi$age <- cdc_bmi$Agemos / 12
cdc_bmi$agemo <- cdc_bmi$Agemos
cdc_bmi <- cdc_bmi %>% select(-c(L, M, S))
cdc_bmi <- cdc_bmi %>% rename(p03="P3", p05="P5", p10="P10", p25="P25", p50="P50", p75="P75", p90="P90", p95="P95", p97="P97")
cdc_bmi_males <- cdc_bmi %>% filter(Sex==1)
cdc_bmi_females <- cdc_bmi %>% filter(Sex==2)
cdc_bmi_males$age <- round(cdc_bmi_males$age,3)
cdc_bmi_males$p03 <- round(cdc_bmi_males$p03,5)
cdc_bmi_males$p05 <- round(cdc_bmi_males$p05,5)
cdc_bmi_males$p10 <- round(cdc_bmi_males$p10,5)
cdc_bmi_males$p25 <- round(cdc_bmi_males$p25,5)
cdc_bmi_males$p50 <- round(cdc_bmi_males$p50,5)
cdc_bmi_males$p75 <- round(cdc_bmi_males$p75,5)
cdc_bmi_males$p90 <- round(cdc_bmi_males$p90,5)
cdc_bmi_males$p95 <- round(cdc_bmi_males$p95,5)
cdc_bmi_males$p97 <- round(cdc_bmi_males$p97,5)
cdc_wgt_males$age <- round(cdc_wgt_males$age,3)
cdc_wgt_males$p03 <- round(cdc_wgt_males$p03,5)
cdc_wgt_males$p05 <- round(cdc_wgt_males$p05,5)
cdc_wgt_males$p10 <- round(cdc_wgt_males$p10,5)
cdc_wgt_males$p25 <- round(cdc_wgt_males$p25,5)
cdc_wgt_males$p50 <- round(cdc_wgt_males$p50,5)
cdc_wgt_males$p75 <- round(cdc_wgt_males$p75,5)
cdc_wgt_males$p90 <- round(cdc_wgt_males$p90,5)
cdc_wgt_males$p95 <- round(cdc_wgt_males$p95,5)
cdc_wgt_males$p97 <- round(cdc_wgt_males$p97,5)

# read in Taylor's file
taylor <- read.csv("/Volumes/Shared/Shared Projects/Laura/Peds Endo/PEDSNet/CAH/Data_raw/cdc_who_data_Taylor_KS_donotuse.csv")
taylor_bmi <- taylor %>% filter(param=="BMI")
taylor_bmi <- taylor_bmi %>% select(!p85) %>% select(!length)
taylor_wgt <- taylor %>% filter(param == "WEIGHTKG")
taylor_wgt <- taylor_wgt %>% select(!p85) %>% select(!length)

a <-comparedf(taylor_bmi, cdc_bmi_males)
summary(a)

b <-comparedf(taylor_wgt, cdc_wgt_males)
summary(b)

