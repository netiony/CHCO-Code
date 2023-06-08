library(fishmethods)
library(tableone)

data <- read.csv("/Volumes/Peds Endo/Lauren Shomaker/Psychotherapy alliance MS/Data/export for SAS.csv")

clus.rho(popchar = data$BO1, cluster = data$coh)
clus.rho(popchar = data$BO16, cluster = data$coh)
clus.rho(popchar = data$BO112, cluster = data$coh)
clus.rho(popchar = data$BOTX1, cluster = data$coh)
clus.rho(popchar = data$BOTX16, cluster = data$coh)
clus.rho(popchar = data$BOTX112, cluster = data$coh)

clus.rho(popchar = data$TA1, cluster = data$coh)
clus.rho(popchar = data$TATX1, cluster = data$coh)
clus.rho(popchar = data$TA16, cluster = data$coh)
clus.rho(popchar = data$TATX16, cluster = data$coh)
clus.rho(popchar = data$TA112, cluster = data$coh)
clus.rho(popchar = data$TATX112, cluster = data$coh)

clus.rho(popchar = data$bmib, cluster = data$coh)
clus.rho(popchar = data$bmi12, cluster = data$coh)
clus.rho(popchar = data$bmip1y, cluster = data$coh)
clus.rho(popchar = data$bmip3y, cluster = data$coh)

clus.rho(popchar = data$locblg, cluster = data$coh)
clus.rho(popchar = data$loc12lg, cluster = data$coh)
clus.rho(popchar = data$locp1y, cluster = data$coh)
clus.rho(popchar = data$locp3y, cluster = data$coh)

t2 <- CreateTableOne(vars=c("TA1","TATX1","TA16","TATX16","TA112","TATX112"), data=data, strata="coh")
t2 <- print(t2,varLabels=TRUE,showAllLevels=TRUE)