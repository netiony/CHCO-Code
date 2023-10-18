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
# calculate effective sample size - Rao and Scott 1992
icc <- 0.04
M <- 12
deff <- 1 + (M - 1)*icc
n_eff <- 100/deff
# 65 per group