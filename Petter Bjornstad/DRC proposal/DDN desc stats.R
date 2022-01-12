library(dplyr)
library(tableone)

data <- read.csv("E:\\Petter Bjornstad\\Grants\\DRC data harmonization\\DDN.csv")
data$dob <- as.Date(data$dob,format = "%m/%d/%Y")
data$f8visitdate <- as.Date(data$f8visitdate,format="%m/%d/%Y") 
data$age <- as.numeric((data$f8visitdate - data$dob)/365.25)
data$dtdiabonset <- as.Date(data$dtdiabonset,format = "%m/%d/%Y")
data$duration <- as.numeric((data$dtdiabonset-data$dob)/365.25)
data$bmi <- data$f6weight/((data$height/100)^2)
data$sex_code <- as.factor(data$sex_code)
data$f5hypertensn <- as.factor(data$f5hypertensn)

baseline <- data %>% filter(redcap_event_name=="interval_0_arm_1")

vars <- c("age","duration","hba1c","sex_code","bmi","gfr","f5hypertensn")

t1 <- CreateTableOne(vars=vars, data=baseline)

# how many people have biopsies?
biopsy <- data[!data$kidneybiopsy_date=="",]
length(unique(biopsy$record_id))
       