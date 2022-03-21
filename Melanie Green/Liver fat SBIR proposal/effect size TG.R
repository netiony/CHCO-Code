library(dplyr)

data <- read.csv("E:/Melanie Green/Grants/Green Wolfe liver fat SBIR/liver data TG for R.csv")
data$delta_tg <- ifelse(data$Order=="A",data$Phase.1.TG0-data$Phase.2.TG0,data$Phase.2.TG0-data$Phase.1.TG0)

# keep only those with hepatic steatosis
data_HS <- data %>% filter(ID %in% c("18-0803-01","18-0803-02","18-0803-08","18-0803-10","18-0803-12","18-0803-13",
                                  "18-0803-15","18-0803-17","18-0803-21","18-0803-22","18-0803-25","18-0803-26"))

mean(data$delta_tg,na.rm = T)
sd(data$delta_tg,na.rm = T)
mean(data_HS$delta_tg,na.rm = T)
sd(data_HS$delta_tg,na.rm = T)

r <- cor(data$delta_tg,data$Phase.1.TG0,use = "complete.obs")
r^2

data$tg_A <- ifelse(data$Order=="A",data$Phase.1.TG0,data$Phase.2.TG0)
data$tg_P <- ifelse(data$Order=="P",data$Phase.1.TG0,data$Phase.2.TG0)
