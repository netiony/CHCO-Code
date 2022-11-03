library(dplyr)
library(sampling)

data <- read.csv("E:/Lauren Shomaker/BREATHE U01/Quarterly audits/BREATHE Quaterly Audit 10.24.22.csv",
                 na.strings = c(""," ","-99"))
data <- unique(data)
data_long <- data

# reshape wide to long
# data_long <- reshape(data, varying = c("SV1.Date","SV2.Date"), direction = "long", v.names = "Visit_date", 
#                      idvar = c("SID","SITE"))
# data_long <- data_long %>% filter(!is.na(data_long$Visit_date))
# data_long <- data_long %>% arrange(SITE)

# for now need to manually enter # of visits by strata (in order as in dataframe)
# will try to fix this later
set.seed(3654)
temp <- strata(data=data_long, stratanames = "Site", size=c(1,2,1), method="srswor")
res <- getdata(data_long, temp)
res <- res[,c("SID","time","Visit_date","SITE")]
res <- res[,c("Record","Site")]


write.csv(res,"E:/Lauren Shomaker/BREATHE U01/Quarterly audits/BREATHE Quarterly Audit 10.24.22 selected visits.csv",
          row.names = F)
