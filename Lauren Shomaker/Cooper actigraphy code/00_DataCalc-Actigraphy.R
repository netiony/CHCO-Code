################################################
# Project: 
# PI: 
# Analyst: Emily Cooper 
# Program: 00_DataCalc-Actigraphy
# Date created: 11/01/2022
# Purpose: Convert raw actigraphy data to by-day and averages
################################################


library(readxl); library(chron); library(Rmpfr)
dtafull <- read_xlsx("INPUTFILE.xlsx")

dtafull$timestamp <- paste0(dtafull$Date, " ",substr(dtafull$Time,12,19))
dtafull$DayNum <- as.POSIXlt(dtafull$Date)$wday
dtafull$Weekend2 <- ifelse(dtafull$DayNum==0 | dtafull$DayNum==6,1,0)

all <- data.frame("ID"=c(), "night"=c(),"bed_in"=c(), "bed_out"=c(),
                  "sleep_onset"=c(), "sleep_offset"=c(), "WASO"=c(), "SMidpoint"=c(), 
                  "TIB"=c(), "SDuration"=c(), "SEff"=c(), "SOL"=c())

for(j in 1:length(unique(dtafull$ID))){
  print(j)
  dta <- data.frame(subset(dtafull, dtafull$ID==unique(dtafull$ID)[j]))
  bed=0
  bed_in <- list()
  bed_out <- list()
  sleep_onset <- list()
  sleep_offset <- list()
  WASO <- list()
  day <- list()
  for(i in 1:nrow(dta)){
    new <- subset(dta, dta$timestamp > bed)
    t <- tapply(new$timestamp,list(new$Interval.Status),min)
    wake <- tryCatch({
      t[["ACTIVE"]]
    }, error=function(e){})
    bed_out <- append(bed_out, wake)
    new <- subset(dta, dta$timestamp > wake)
    t <- tapply(new$timestamp,list(new$Interval.Status),min)
    bed <- tryCatch({
      min(t[["REST-S"]], t[["REST"]])
    }, error=function(e){})
    bed_in <- append(bed_in, bed)
    sleep <- tryCatch({
      t[["REST-S"]]
    }, error=function(e){})
    sleep_onset <- append(sleep_onset, sleep)
    weekend <- new$Weekend2[new$timestamp==bed]
    day <- append(day, weekend)
  }
  bed_out <- bed_out[-1]
  if(length(bed_in) > length(bed_out)){
    bed_in <- bed_in[-length(bed_in)]
  }
  
  for(i in 1:length(bed_out)){
    new <- subset(dta, dta$timestamp >= bed_in[i] & dta$timestamp <= bed_out[i])
    t <- tapply(new$timestamp,list(new$Interval.Status),max)
    wake <- tryCatch({
      t[["REST-S"]]
    }, error=function(e){})
    sleep_offset <- append(sleep_offset, wake)
    new <- subset(dta, dta$timestamp > sleep_onset[i] & dta$timestamp < sleep_offset[i])
    awake <- length(which(new$Sleep.Wake==1))
    WASO <- append(WASO, awake)
  }
  
  if(length(sleep_offset) < length(sleep_onset)){
    sleep_onset <- sleep_onset[-length(sleep_onset)]
  }

    #### Time corrections for daylight savings ####
      #For each daylight savings start or end date, you'll need one of the following chunks
      #In each chunk, specify the record ID(s) and the date of daylight savings start or end
      #This should shift all bed time, wake time, sleep onset, and sleep offset variables by 1 hour in the appropriate direction
      #The time of the shift should always be 2am (02:00:00)
  ## For March daylight savings start
  # if(unique(dtafull$ID)[j]=="RECORD_ID"){ #NEED TO PLUG IN RECORD ID(s)
  #   bed_in <- lapply(bed_in, function(x){
  #     ifelse(x >= "YYYY-03-DD 02:00:00", #NEED TO PLUG IN MARCH DATE
  #            as.character(as.POSIXct(x)+3600), 
  #            as.character(as.POSIXct(x)))})
  #   bed_out <- lapply(bed_out, function(x){
  #     ifelse(x >= "YYYY-03-DD 02:00:00",
  #            as.character(as.POSIXct(x)+3600), 
  #            as.character(as.POSIXct(x)))})
  #   sleep_onset <- lapply(sleep_onset, function(x){
  #     ifelse(x >= "YYYY-03-DD 02:00:00",
  #            as.character(as.POSIXct(x)+3600), 
  #            as.character(as.POSIXct(x)))})
  #   sleep_offset <- lapply(sleep_offset, function(x){
  #     ifelse(x >= "YYYY-03-DD 02:00:00",
  #            as.character(as.POSIXct(x)+3600), 
  #            as.character(as.POSIXct(x)))})
  ## For November daylight savings end
  # if(unique(dtafull$ID)[j]=="RECORD_ID"){ #NEED TO PLUG IN RECORD ID(s)
  #   bed_in <- lapply(bed_in, function(x){
  #     ifelse(x >= "YYYY-11-DD 02:00:00", #NEED TO PLUG IN NOVEMBER DATE
  #            as.character(as.POSIXct(x)-3600), 
  #            as.character(as.POSIXct(x)))})
  #   bed_out <- lapply(bed_out, function(x){
  #     ifelse(x >= "YYYY-11-DD 02:00:00",
  #            as.character(as.POSIXct(x)-3600), 
  #            as.character(as.POSIXct(x)))})
  #   sleep_onset <- lapply(sleep_onset, function(x){
  #     ifelse(x >= "YYYY-11-DD 02:00:00",
  #            as.character(as.POSIXct(x)-3600), 
  #            as.character(as.POSIXct(x)))})
  #   sleep_offset <- lapply(sleep_offset, function(x){
  #     ifelse(x >= "YYYY-11-DD 02:00:00",
  #            as.character(as.POSIXct(x)-3600), 
  #            as.character(as.POSIXct(x)))})
  
  #### Calculate: TIB (time in bed) ####  
  TIB <- mapply(function(X,Y){
    difftime(Y, X, units="mins") #difftime accounts for daylight savings
  }, X=bed_in, Y=bed_out)
  
  #### Calculate: sleep duration ####
  SDuration <- mapply(function(X,Y,Z){
    difftime(Y, X, units="mins")-Z #difftime accounts for daylight savings
  }, X=sleep_onset, Y=sleep_offset, Z=WASO)
  
  #### Calculate: SOL (sleep onset latency) ####
  SOL <- mapply(function(X,Y){
    difftime(X, Y, units="mins")
  }, X=sleep_onset, Y=bed_in)
  
    #### Calculate: sleep efficiency (sleep time/time in bed) #### 
  efficiency <- SDuration/TIB 
  
  #### Calculate: sleep midpoint ([time in bed/2] + bedtime) ####
  SM <- mapply(function(X,Y){
    as.character(as.POSIXct(X)+(difftime(Y, X, units="mins")/2))
  }, X=sleep_onset, Y=sleep_offset)
  
  #### Make summary by-day data frame ####
  by_day <- data.frame("night"=unlist(as.POSIXlt(unlist(bed_out))$wday), "bed_in"=unlist(bed_in), "bed_out"=unlist(bed_out),
                       "sleep_onset"=unlist(sleep_onset), "sleep_offset"=unlist(sleep_offset), "SM"=unlist(SM), "eff"=unlist(efficiency),
                       "WASO"=unlist(WASO), "TIB"=unlist(TIB), "SDuration"=unlist(SDuration), "SOL"=unlist(SOL))
  weekday_weekend <- data.frame("night"=c("weekday","weekend"), "bed_in"=c(NA,NA), "bed_out"=c(NA,NA),
                                "sleep_onset"=c(NA,NA), "sleep_offset"=c(NA,NA), "SM"=c(NA,NA), "eff"=c(NA,NA),
                                "WASO"=c(NA,NA), "TIB"=c(NA,NA), "SDuration"=c(NA,NA), "SOL"=c(NA,NA))
  by_day$nap <- ifelse(by_day$SDuration==min(by_day$SDuration[duplicated(substr(cbind(by_day$bed_out, by_day$bed_in), 1 , 10))|
                                          duplicated(substr(cbind(by_day$bed_out, by_day$bed_in), 1 , 10), fromLast=TRUE)]),1,0)
  by_day$night <- ifelse(by_day$night==0 | by_day$night==6, "weekend", "weekday")
  #### Averages by all days, weekdays, weekend days ####
  l <- c("bed_in","bed_out","sleep_onset","sleep_offset","SM","WASO","TIB","SDuration","eff","SOL")
  wkd <- subset(by_day, by_day$night=="weekday" & by_day$nap==0); wnd <- subset(by_day, by_day$night=="weekend" & by_day$nap==0)
  if(nrow(wkd)>0){
    weekday_weekend[1,l] <- lapply(wkd[,l], function(x){ #group by weekday/weekend day
     if(class(x)=="character"){
       ifelse(chron(mean(ifelse(format(as.POSIXct(x), "%H:%M:%S") < format(as.POSIXct(wkd$sleep_offset), "%H:%M:%S"),
                                                    times(format(as.POSIXct(x), "%H:%M:%S"))+1,
                                                    times(format(as.POSIXct(x), "%H:%M:%S")))))==1, "00:00:00", 
              substr(chron(mean(ifelse(format(as.POSIXct(x), "%H:%M:%S") < format(as.POSIXct(wkd$sleep_offset), "%H:%M:%S"),
                                       times(format(as.POSIXct(x), "%H:%M:%S"))+1,
                                       times(format(as.POSIXct(x), "%H:%M:%S"))))), 11,18))
    }
    else{
      mean(x)
    }
    })}
  if(nrow(wnd)>0){
    weekday_weekend[2,l] <- lapply(wnd[,l], function(x){ #group by weekday/weekend day
     if(class(x)=="character"){
       ifelse(chron(mean(ifelse(format(as.POSIXct(x), "%H:%M:%S") < format(as.POSIXct(wnd$sleep_offset), "%H:%M:%S"),
                                times(format(as.POSIXct(x), "%H:%M:%S"))+1,
                                times(format(as.POSIXct(x), "%H:%M:%S")))))==1, "00:00:00", 
              substr(chron(mean(ifelse(format(as.POSIXct(x), "%H:%M:%S") < format(as.POSIXct(wnd$sleep_offset), "%H:%M:%S"),
                                       times(format(as.POSIXct(x), "%H:%M:%S"))+1,
                                       times(format(as.POSIXct(x), "%H:%M:%S"))))), 11,18))
    }
    else{
      mean(x)
    }
    })}
  all_days <- subset(by_day, by_day$nap==0)
  all_days <- lapply(all_days[,l], function(x){ #for all days
    if(class(x)=="character"){
      ifelse(chron(mean(ifelse(format(as.POSIXct(x), "%H:%M:%S") < format(as.POSIXct(all_days$sleep_offset), "%H:%M:%S"),
                               times(format(as.POSIXct(x), "%H:%M:%S"))+1,
                               times(format(as.POSIXct(x), "%H:%M:%S")))))==1, "00:00:00", 
             substr(chron(mean(ifelse(format(as.POSIXct(x), "%H:%M:%S") < format(as.POSIXct(all_days$sleep_offset), "%H:%M:%S"),
                                      times(format(as.POSIXct(x), "%H:%M:%S"))+1,
                                      times(format(as.POSIXct(x), "%H:%M:%S"))))), 11,18))
      }
    else{
      mean(x)
    }
  })
  
  #### Calculate: social jetlag (weekend sleep midpoint - weekday sleep midpoint) #### 
  SJL <- ifelse(is.na(weekday_weekend$SM[1])==F & is.na(weekday_weekend$SM[2])==F,
                difftime(as.POSIXct(paste0("2000-02-01 ",weekday_weekend$SM[2])), 
                         as.POSIXct(paste0("2000-02-01 ",weekday_weekend$SM[1])), units="mins"),NA)
  
  #### Write as data frame ####
  byday_final <- data.frame("ID"=rep(unique(dta$ID), length(bed_in)), "night"=unlist(as.POSIXlt(unlist(bed_out))$wday),"bed_in"=unlist(bed_in), "bed_out"=unlist(bed_out),
                            "sleep_onset"=unlist(sleep_onset), "sleep_offset"=unlist(sleep_offset), "WASO"=unlist(WASO), "SOL"=unlist(SOL),"SMidpoint"=unlist(SM), 
                            "TIB"=unlist(TIB), "SDuration"=unlist(SDuration), "SEff"=unlist(efficiency),"SJL"=rep(NA, length(bed_in)), "nap"=by_day$nap)
  summaries <- data.frame("ID"=rep(unique(dta$ID), 3), 
                          "night"=c(weekday_weekend$night[1], weekday_weekend$night[2], "all"),
                          "bed_in"=c(weekday_weekend$bed_in[1], weekday_weekend$bed_in[2], all_days$bed_in),
                          "bed_out"=c(weekday_weekend$bed_out[1], weekday_weekend$bed_out[2], all_days$bed_out),
                          "sleep_onset"=c(weekday_weekend$sleep_onset[1], weekday_weekend$sleep_onset[2], all_days$sleep_onset),
                          "sleep_offset"=c(weekday_weekend$sleep_offset[1], weekday_weekend$sleep_offset[2], all_days$sleep_offset),
                          "WASO"=c(weekday_weekend$WASO[1], weekday_weekend$WASO[2], all_days$WASO), 
                          "SOL"=c(weekday_weekend$SOL[1], weekday_weekend$SOL[2], all_days$SOL), 
                          "SMidpoint"=c(weekday_weekend$SM[1], weekday_weekend$SM[2], all_days$SM), 
                          "TIB"=c(weekday_weekend$TIB[1], weekday_weekend$TIB[2], all_days$TIB), 
                          "SDuration"=c(weekday_weekend$SDuration[1], weekday_weekend$SDuration[2], all_days$SDuration), 
                          "SEff"=c(weekday_weekend$eff[1], weekday_weekend$eff[2], all_days$eff),
                          "SJL"=c(NA, NA, SJL), "nap"=rep(NA, 3))
  summaries <- subset(summaries, is.na(summaries$night)==F)
  temp <- data.frame(rbind(byday_final, summaries))
  all <- rbind(all, temp)
}

all$night <- ifelse(all$night==0, "Sunday",
                    ifelse(all$night==1, "Monday",
                           ifelse(all$night==2, "Tuesday",
                                  ifelse(all$night==3, "Wednesday",
                                         ifelse(all$night==4, "Thursday",
                                                ifelse(all$night==5, "Friday",
                                                       ifelse(all$night==6, "Saturday", all$night)))))))

#Convert Sleep duration, TIB, and SJL to be in hours
all$SDuration <- all$SDuration/60
all$TIB <- all$TIB/60
all$SJL <- all$SJL/60

#Convert sleep efficiency to %
all$SEff <- 100*all$SEff


#Convert date-time objects to decimal time (1=am, 2=2am ... 24=midnight, >24 indicates they went to sleep after midnight)
all$bed_out_DEC <- 24*as.numeric(times(format(as.POSIXct(gsub("  "," ", paste0("2022-01-01 ",substr(all$bed_out, nchar(all$bed_out)-8, nchar(all$bed_out)))), "%Y-%m-%d %H:%M:%S", tz="UTC"), "%H:%M:%S")))                      

all$bed_in_DEC <- ifelse(24*as.numeric(times(format(as.POSIXct(gsub("  "," ", paste0("2022-01-01 ",substr(all$bed_in, nchar(all$bed_in)-8, nchar(all$bed_in)))), "%Y-%m-%d %H:%M:%S", tz="UTC"), "%H:%M:%S"))) < all$bed_out_DEC, 
                         24+24*as.numeric(times(format(as.POSIXct(gsub("  "," ", paste0("2022-01-01 ",substr(all$bed_in, nchar(all$bed_in)-8, nchar(all$bed_in)))), "%Y-%m-%d %H:%M:%S", tz="UTC"), "%H:%M:%S"))),
                         24*as.numeric(times(format(as.POSIXct(gsub("  "," ", paste0("2022-01-01 ",substr(all$bed_in, nchar(all$bed_in)-8, nchar(all$bed_in)))), "%Y-%m-%d %H:%M:%S", tz="UTC"), "%H:%M:%S"))))

all$sleep_onset_DEC <- ifelse(24*as.numeric(times(format(as.POSIXct(gsub("  "," ", paste0("2022-01-01 ",substr(all$sleep_onset, nchar(all$sleep_onset)-8, nchar(all$sleep_onset)))), "%Y-%m-%d %H:%M:%S", tz="UTC"), "%H:%M:%S"))) < all$bed_out_DEC, 
                              24+24*as.numeric(times(format(as.POSIXct(gsub("  "," ", paste0("2022-01-01 ",substr(all$sleep_onset, nchar(all$sleep_onset)-8, nchar(all$sleep_onset)))), "%Y-%m-%d %H:%M:%S", tz="UTC"), "%H:%M:%S"))),
                              24*as.numeric(times(format(as.POSIXct(gsub("  "," ", paste0("2022-01-01 ",substr(all$sleep_onset, nchar(all$sleep_onset)-8, nchar(all$sleep_onset)))), "%Y-%m-%d %H:%M:%S", tz="UTC"), "%H:%M:%S"))))

all$sleep_offset_DEC <- 24*as.numeric(times(format(as.POSIXct(gsub("  "," ", paste0("2022-01-01 ",substr(all$sleep_offset, nchar(all$sleep_offset)-8, nchar(all$sleep_offset)))), "%Y-%m-%d %H:%M:%S", tz="UTC"), "%H:%M:%S")))                      

all$SMidpoint_DEC <- ifelse(24*as.numeric(times(format(as.POSIXct(gsub("  "," ", paste0("2022-01-01 ",substr(all$SMidpoint, nchar(all$SMidpoint)-8, nchar(all$SMidpoint)))), "%Y-%m-%d %H:%M:%S", tz="UTC"), "%H:%M:%S"))) < all$bed_out_DEC, 
                            24+24*as.numeric(times(format(as.POSIXct(gsub("  "," ", paste0("2022-01-01 ",substr(all$SMidpoint, nchar(all$SMidpoint)-8, nchar(all$SMidpoint)))), "%Y-%m-%d %H:%M:%S", tz="UTC"), "%H:%M:%S"))),
                            24*as.numeric(times(format(as.POSIXct(gsub("  "," ", paste0("2022-01-01 ",substr(all$SMidpoint, nchar(all$SMidpoint)-8, nchar(all$SMidpoint)))), "%Y-%m-%d %H:%M:%S", tz="UTC"), "%H:%M:%S"))))

write.csv(all, "OUTPUTFILE.csv", row.names=FALSE)

