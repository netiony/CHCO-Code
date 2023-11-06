library(dplyr)
library(sampling)
library(readxl)

data <- read_xlsx("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Lauren Shomaker/BREATHE U01/Session ratings/List for Intervention Ratings 9.27.23.xlsx")
data$...5 <- NULL
data <- unique(data)
data <- arrange(data,Cohort)
data$Lead_Cohort <- paste(data$`Lead Facilitator`,data$Cohort)

# first select 30% of sessions for intervention fidelity, stratified by Lead Facilitator 
set.seed(3654)
numvisits <- table(data$`Lead Facilitator`)
numvisits_select <- floor(table(data$`Lead Facilitator`)*0.3)
for (i in 1:8) {
  numvisits_select[i] <- ifelse(numvisits_select[i]==0, 1, numvisits_select[i])
}
temp <- strata(data=data, stratanames = "Lead Facilitator", size=numvisits_select, method="srswor")
res <- getdata(data, temp)
res <- res[,c("Modality","Cohort","Session #","Lead Facilitator")]
# now assign 70% of each modality to the corresponding expert
numvisits_modality <- table(res$Modality)
numvisits_modality_select <- floor(table(res$Modality)*0.7)
res <- arrange(res,Modality)
temp_modality <- strata(data=res, stratanames = "Modality", size=numvisits_modality_select, method="srswor")
res_modality <- getdata(res, temp_modality)
res_modality$Expert <- ifelse(res_modality$Modality=="HE","Lauren",
                              ifelse(res_modality$Modality=="L2B","Trish",
                                     ifelse(res_modality$Modality=="CBT","Heather",NA)))
res_modality <- res_modality[,c("Cohort","Session #","Modality","Expert")]
# merge back to res
res <- merge(res, res_modality, by=c("Cohort","Session #","Modality"), all.x = T, all.y = T)
# now randomly assign the remaining sessions
res_unassigned <- res[is.na(res$Expert),]
res_unassigned <- arrange(res_unassigned, Cohort)
nn <- floor(nrow(res_unassigned)*0.33)
res_unassigned$Expert <- c(rep("Lauren",nn),rep("Trish",nn),rep("Heather",nrow(res_unassigned)-(2*nn)))
res_unassigned <- res_unassigned[,c("Cohort","Session #","Modality","Expert")]
# merge back to res
res <- merge(res, res_unassigned, by=c("Cohort","Session #","Modality"), all.x = T, all.y = T)
res$Expert <- ifelse(is.na(res$Expert.x),res$Expert.y,res$Expert.x)
res$Expert.x <- NULL
res$Expert.y <- NULL
write.csv(res,"/Volumes/Shared/Shared Projects/Laura/Peds Endo/Lauren Shomaker/BREATHE U01/Session ratings/Intervention fidelity ratings selected visits.csv",
          row.names = F)

# then select 30% of cohorts 2-6 for alliance and cohesion, stratified by cohort and lead facilitator
data_cohort26 <- data[data$Cohort>=2 & data$Cohort<=6,]
numvisits <- table(data_cohort26$`Lead Facilitator`)
numvisits_select <- floor(table(data_cohort26$`Lead Facilitator`)*0.3)
for (i in 1:8) {
  numvisits_select[i] <- ifelse(numvisits_select[i]==0, 1, numvisits_select[i])
}
temp26 <- strata(data=data_cohort26, stratanames = "Lead Facilitator", size=numvisits_select, method="srswor")
res26 <- getdata(data_cohort26, temp26)
res26 <- res26[,c("Modality","Cohort","Session #","Lead Facilitator")]

write.csv(res26,"/Volumes/Shared/Shared Projects/Laura/Peds Endo/Lauren Shomaker/BREATHE U01/Session ratings/Intervention alliance and cohesion ratings selected visits.csv",
          row.names = F)



