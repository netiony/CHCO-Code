library(TMB)
library(fishmethods)
library(readxl)
library(tableone)

data <- read_xlsx("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Lauren Shomaker/Mackey BREATHE/BREATHET1Dprelimdata.xlsx")

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

# assume 15% attrition
n_eff <- 85/deff
# 59 per group

t1 <- CreateTableOne(vars=c("a1c_fu1","a1c_fu2","a1c_post","a1c_before","a1c_after","a1c_late"), data=data)
t1 <- print(t1, nonnormal=c("a1c_fu1","a1c_fu2","a1c_post","a1c_before","a1c_after","a1c_late"))
