library(childsds)

data <- read.csv("E:\\Nokoff\\151215 study\\151215TransgenderStu-DataForBMIPercentile_DATA_LABELS_2020-12-02_1525.csv")
data$Date.of.birth <- as.Date(as.character(data$Date.of.birth),format="%Y-%m-%d")
data$Date.of.study.visit <- as.Date(as.character(data$Date.of.study.visit),format="%Y-%m-%d")
data[data$Participant.ID==5,]$Date.of.birth <- as.Date("2001-06-06",format="%Y-%m-%d")
data$ageyears <- as.numeric(floor((data$Date.of.study.visit-data$Date.of.birth)/365.25))
data$bmi <- data$Weight..kg. / ((data$Height..cm./100)^2)

#data <- data[data$Participant.ID != 17,]

data$bmipct = sds(data$bmi,
              age = data$ageyears,
              sex = data$Natal.sex, male = "Male", female =  "Female",
              ref = cdc.ref,
              item = "bmi",
              type = "perc")
data$bmipct <- round(data$bmipct*100,2)

write.csv(data,"E:\\Nokoff\\151215 study\\bmi_perc.csv")