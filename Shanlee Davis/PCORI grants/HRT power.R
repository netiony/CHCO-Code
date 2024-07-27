library(nlme)

# read in simulated data
data <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Davis/PCORI proposals/HRT/ald_sim_totN_300.csv")

# check data
table(data$sim)
table(data$cohort)
table(data$id)
table(data$period)

summary(data$y)

# code new variable for each group
# if 2 years on a drug, code 2
# if 1 year on a drug, code 1
# if not on a drug, code 0
# WHY ARE THERE ONLY 2 GROUPS IN THE SIM DATA?
data$OCP <- ifelse(data$group %in% c(""), 2,
                   ifelse(data$group %in% c(""), 1, 0))
data$transdermal <- ifelse(data$group %in% c(""), 2,
                           ifelse(data$group %in% c(""), 1, 0))
data$oral <- ifelse(data$group %in% c(""), 2,
                           ifelse(data$group %in% c(""), 1, 0))

for (i in 1:1000) {
  model <- lme(y ~ age * group, random = ~1 | id, data = data[data$sim == i,])
}
# does group need to be a time dependent covariate? how would I generate data like that?