library(haven)

data <- read_sav("H:\\Endocrinology\\Shomaker\\Psychotherapy alliance MS\\Data\\Shomaker Truncated-v4.21.20.sav")

table(data$group,data$ther)
