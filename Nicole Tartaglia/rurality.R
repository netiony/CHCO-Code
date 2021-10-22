library(dplyr)
library(stringr)

babydata <- read.csv("H:\\Tartaglia\\Rurality\\EXtraordinarYBabyStu-ZipCodesRurality_DATA_2021-10-20_1602.csv",na.strings = c("NA",""," "))
babydata$ses_address_zip <- str_pad(babydata$ses_address_zip,5,side="left",pad = "0")

# keep only visits with zip codes
babydata <- babydata[!is.na(babydata$ses_address_zip),]
# keep one record per person
zips <- babydata %>% group_by(ï..study_id_extraordinary) %>% filter(row_number()==1)

# check that all the IDs have a zip code
ids <- babydata$ï..study_id_extraordinary
bad <- zips$ï..study_id_extraordinary[!(zips$ï..study_id_extraordinary %in% ids)]
numids <- length(unique(babydata$ï..study_id_extraordinary))

zips$ses_address_zip <- ifelse(nchar(zips$ses_address_zip)>7,str_sub(zips$ses_address_zip,1,5),zips$ses_address_zip)

ruca <- read.csv("H:\\Tartaglia\\Rurality\\AP RUCA Codes - Pops 07.18.2019 - Copy.csv")
ruca <- ruca[,c("ZIP.CODE","RUCA2.0")]
ruca$ZIP.CODE <- as.character(ruca$ZIP.CODE)
ruca$ZIP.CODE <- str_pad(ruca$ZIP.CODE,5,side="left",pad = "0")
colnames(ruca) <- c("ses_address_zip","RUCA2.0")

# merge
zips <- merge(zips,ruca,by="ses_address_zip",all.x = T,all.y = F)

# write
write.csv(zips,"H:\\Tartaglia\\Rurality\\baby_study_zips_rurality.csv",row.names = F)