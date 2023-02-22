library(dplyr)
library(sampling)

data <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Lauren Shomaker/BREATHE U01/Quarterly audits/BREATHE Quaterly Audit 2.1.23.csv",
                 na.strings = c(""," ","-99"))
data <- unique(data)
data_long <- data

# reshape wide to long
# data_long <- reshape(data, varying = c("SV1.Date","SV2.Date"), direction = "long", v.names = "Visit_date", 
#                      idvar = c("SID","SITE"))
# data_long <- data_long %>% filter(!is.na(data_long$Visit_date))
# data_long <- data_long %>% arrange(SITE)

# for now need to manually enter # of visits by strata (in order as in dataframe)
# add as 10% of the total number of visits in each stratum
# will try to fix this later
set.seed(3654)
numvisits <- table(data_long$Site)
temp <- strata(data=data_long, stratanames = "Site", size=c(2,2,1,1), method="srswor")
res <- getdata(data_long, temp)
#res <- res[,c("SID","time","Visit_date","SITE")]
res <- res[,c("Record","Site")]


write.csv(res,"/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Lauren Shomaker/BREATHE U01/Quarterly audits//BREATHE Quarterly Audit 2.1.23 selected visits.csv",
          row.names = F)
