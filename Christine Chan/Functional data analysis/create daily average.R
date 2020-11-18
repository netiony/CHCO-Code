library(data.table)
library(lubridate)
library(dplyr)
library(ggplot2)

master <- read.csv("H:\\Endocrinology\\Chan\\Functional data analysis\\Cleaned CGM Data\\master.csv")

temp <- master
temp$hour <- hour(temp$timestamp)
temp$min <- minute(temp$timestamp)
temp$sec <- second(temp$timestamp)
temp$min_floor <- temp$min - (temp$min - (5*floor(temp$min/5)))
temp$min_floor <- ifelse(temp$min_floor=="0","00",temp$min_floor)
temp$min_floor <- ifelse(temp$min_floor=="5","05",temp$min_floor)
for (i in 1:nrow(temp)) {
  temp$wind[i] <- paste0(temp$hour[i],":",temp$min_floor[i])                   
}

mean_output <- aggregate(sensorglucose~subjectid + as.factor(wind), data=temp, FUN="mean")
colnames(mean_output) <- c("subjectid","wind","sensorglucose")
mean_output <- arrange(mean_output,subjectid,wind)
for (i in 1:nrow(mean_output)) {
  mean_output$t1[i] <- as.numeric(strsplit(as.character(mean_output$wind[i]),":")[[1]][1])
  mean_output$t2[i] <- as.numeric(strsplit(as.character(mean_output$wind[i]),":")[[1]][2])
}
mean_output$minutes <- (mean_output$t1*60) + mean_output$t2
mean_output <- arrange(mean_output,subjectid,minutes)

p <- ggplot(data=mean_output,aes(x=minutes,y=sensorglucose,group=subjectid))
p + geom_point()
p + geom_line()

write.csv(mean_output,"H:\\Endocrinology\\Chan\\Functional data analysis\\Cleaned CGM Data\\mean_output.csv")
