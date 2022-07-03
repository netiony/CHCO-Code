library(dplyr)
library(sampling)

data <- read.csv("E:/Lauren Shomaker/BREATHE U01/Quarterly audits/BREATHE Quaterly Audit 6.30.22.csv",
                 na.strings = c(""," ","-99"))
data <- unique(data)

# reshape wide to long
data_long <- reshape(data, varying = c("SV1.Date","SV2.Date"), direction = "long", v.names = "Visit_date", 
                     idvar = c("SID","SITE"))
data_long <- data_long %>% filter(!is.na(data_long$Visit_date))
data_long <- data_long %>% arrange(SITE)

# for now need to manually enter # of visits by strata (in order as in dataframe)
# will try to fix this later
set.seed(3654)
temp <- strata(data=data_long, stratanames = "SITE", size=c(1,1,1), method="srswor")
res <- getdata(data_long, temp)
res <- res[,c("SID","time","Visit_date","SITE")]

write.csv(res,"E:/Lauren Shomaker/BREATHE U01/Quarterly audits/BREATHE Quarterly Audit 6.30.22 selected visits.csv",
          row.names = F)
