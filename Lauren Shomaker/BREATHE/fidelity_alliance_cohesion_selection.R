library(dplyr)
library(sampling)
library(readxl)

data <- read_xlsx("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Lauren Shomaker/BREATHE U01/Session ratings/List for Intervention Ratings 9.27.23.xlsx")
data$...5 <- NULL
data <- unique(data)
data <- arrange(data,Cohort)
data$Lead_Cohort <- paste(data$`Lead Facilitator`,data$Cohort)

# first select 30% of sessions for intervention fidelity, stratified by cohort 
set.seed(3654)
numvisits <- table(data$Lead_Cohort)
numvisits_select <- floor(table(data$Lead_Cohort)*0.3)
for (i in 1:8) {
  numvisits_select[i] <- ifelse(numvisits_select[i]==0, 1, numvisits_select[i])
}
temp <- strata(data=data, stratanames = "Lead_Cohort", size=numvisits_select, method="srswor")
res <- getdata(data, temp)
res <- res[,c("Modality","Cohort","Session #","Lead Facilitator")]

write.csv(res,"/Volumes/Shared/Shared Projects/Laura/Peds Endo/Lauren Shomaker/BREATHE U01/Session ratings/Intervention fidelity ratings selected visits.csv",
          row.names = F)

# then select 30% of cohorts 2-6 for alliance and cohesion, stratified by cohort and lead facilitator
data_cohort26 <- data[data$Cohort>=2 & data$Cohort<=6,]
numvisits <- table(data_cohort26$Lead_Cohort)
numvisits_select <- floor(table(data_cohort26$Lead_Cohort)*0.3)
for (i in 1:8) {
  numvisits_select[i] <- ifelse(numvisits_select[i]==0, 1, numvisits_select[i])
}
temp26 <- strata(data=data_cohort26, stratanames = "Lead_Cohort", size=numvisits_select, method="srswor")
res26 <- getdata(data_cohort26, temp26)
res26 <- res26[,c("Modality","Cohort","Session #","Lead Facilitator")]

write.csv(res26,"/Volumes/Shared/Shared Projects/Laura/Peds Endo/Lauren Shomaker/BREATHE U01/Session ratings/Intervention alliance and cohesion ratings selected visits.csv",
          row.names = F)



