library(nlme)
library(lme4)
library(ggplot2)
library(dplyr)

# read in data that Tim simulated
data <- read.csv("S:\\Shared Projects\\Laura\\BDC\\Projects\\Petter Bjornstad\\ALD R01\\
                 Power Calculation\\Data_Raw\\simulated_ald_data_0.5_1.1_half_resid.csv")
data <- read.csv("S:\\Shared Projects\\Laura\\BDC\\Projects\\Petter Bjornstad\\ALD R01\\Power Calculation\\Data_Raw\\simulated_ald_data_0.5_1.1_quarter_resid.csv")

# model checking and plots
quick <- data[data$sim==1,]
mod1 <- lmer(y ~ group*age + (1 + age|id), data=quick)
mod1 <- lme(y ~ group*age, random = ~1|id, data=quick, control=ctrl)
quick %>% ggplot(aes(x=age,y=y,group=group)) + geom_line()

# power calcs
ctrl <- lmeControl(maxIter = 1000)
p <- NULL
for (i in 1:1000) {
  quick <- data[data$sim==i,]
  try(mod1 <- lme(y ~ group*age, random = ~1|id, data=quick, control=ctrl))
  pcur <- summary(mod1)$tTable[4,5]
  p <- c(p,pcur)
}