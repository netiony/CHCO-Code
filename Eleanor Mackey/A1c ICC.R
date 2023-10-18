library(TMB)
library(fishmethods)
library(readxl)

data <- read_xlsx("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Eleanor Mackey/BREATHET1Dprelimdata.xlsx")

clus.rho(popchar = data$a1c_fu1, cluster = data$cohort)
clus.rho(popchar = data$a1c_fu2, cluster = data$cohort)
clus.rho(popchar = data$a1c_post, cluster = data$cohort)
clus.rho(popchar = data$a1c_before, cluster = data$cohort)
clus.rho(popchar = data$a1c_after, cluster = data$cohort)
clus.rho(popchar = data$a1c_late, cluster = data$cohort)

# use ICC of 0.05