library(nlme)
library(stringr)
library(dplyr)

# read in simulated data
data <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Davis/PCORI proposals/HRT/ald_sim_totN_300.csv")

# check data
table(data$sim, useNA = "ifany")
table(data$cohort, useNA = "ifany")
table(data$id, useNA = "ifany")
table(data$period, useNA = "ifany")
table(data$group, useNA = "ifany")

summary(data$y)

# code new variable for each group
# if 2 years on a drug, code 2
# if 1 year on a drug, code 1
# if not on a drug, code 0
data$OCP <- str_count(data$group, "OCP")
data$transdermal <- str_count(toupper(data$group), "TRANSDERMAL")
data$oral <- str_count(toupper(data$group), "ORAL")

data <- data %>% group_by(sim, id) %>% filter(row_number() == 1)

p_ocp <- NULL
p_oral <- NULL
p_transdermal <- NULL

for (i in 1:1000) {
  model_ocp <- lme(y ~ OCP , random = ~1 | id, data = data[data$sim == i,])
  p_ocp <- c(p_ocp, summary(model_ocp)$tTable[2,5])
  model_oral <- lme(y ~ oral , random = ~1 | id, data = data[data$sim == i,])
  p_oral <- c(p_oral, summary(model_oral)$tTable[2,5])
  model_transdermal <- lme(y ~ transdermal , random = ~1 | id, data = data[data$sim == i,])
  p_transdermal <- c(p_transdermal, summary(model_transdermal)$tTable[2,5])
}

# power below is expressed as percent
sum(p_ocp<0.05)/10
sum(p_oral<0.05)/10
sum(p_transdermal<0.05)/10
