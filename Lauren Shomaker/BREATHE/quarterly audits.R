library(dplyr)
library(sampling)

data <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Lauren Shomaker/BREATHE U01/Quarterly audits/Breathe Audit List AFTER 12.15.2023.csv",
                 na.strings = c(""," ","-99"))
data$Record <- paste(data$SID, data$interval)
data <- unique(data)
data_long <- data

data_long <- arrange(data_long,SITE)

# reshape wide to long
# data_long <- reshape(data, varying = c("SV1.Date","SV2.Date"), direction = "long", v.names = "Visit_date", 
#                      idvar = c("SID","SITE"))
# data_long <- data_long %>% filter(!is.na(data_long$Visit_date))
# data_long <- data_long %>% arrange(SITE)

# for now need to manually enter # of visits by strata (in order as in dataframe)
# add as 10% of the total number of visits in each stratum
# will try to fix this later
set.seed(3654)
numvisits <- table(data_long$SITE)
# these are not in the correct order - need to reorder SITEs in numvisits to match the file
numvisits_select <- floor(table(data_long$SITE)*0.1)
for (i in 1:4) {
  numvisits_select[i] <- ifelse(numvisits_select[i]==0, 1, numvisits_select[i])
}
#temp <- strata(data=data_long, stratanames = "SITE", size=c(floor(numvisits[3]*0.1),floor(numvisits[1]*0.1),floor(numvisits[4]*0.1),
#                                                            floor(numvisits[2]*0.1)), method="srswor")
temp <- strata(data=data_long, stratanames = "SITE", size=numvisits_select, method="srswor")
res <- getdata(data_long, temp)
#res <- res[,c("SID","time","Visit_date","SITE")]
res <- res[,c("Record","SITE")]


write.csv(res,"/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Lauren Shomaker/BREATHE U01/Quarterly audits/Breathe Audit List AFTER 12.15.2023 selected visits.csv",
          row.names = F)

