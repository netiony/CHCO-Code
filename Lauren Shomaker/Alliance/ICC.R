library(fishmethods)

data <- read.csv("/Volumes/Peds Endo/Lauren Shomaker/Psychotherapy alliance MS/Data/export for SAS.csv")

clus.rho(popchar = data$BO1, cluster = data$coh)
clus.rho(popchar = data$BO16, cluster = data$coh)
clus.rho(popchar = data$BO112, cluster = data$coh)
clus.rho(popchar = data$BOTX1, cluster = data$coh)
clus.rho(popchar = data$BOTX16, cluster = data$coh)
clus.rho(popchar = data$BOTX112, cluster = data$coh)