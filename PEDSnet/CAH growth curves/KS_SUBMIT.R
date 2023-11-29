##############################################
library(readr)
library(quantregGrowth)
library(growthcleanr)
library(ggplot2)
library(data.table)
library(lubridate)
library(tidyverse)
library(skimr)

##############################################
# upload dt
dt <- read_csv("...PATHWAY.../KS_T.Marshall_HgtWgt.csv")  # UPDATE ME :)))))
dt <- subset(dt, kline=='Case') 
dt <- rename(dt, subjid=person_id, dx_dt=first_dx_dt, date=measurement_date)
dt <- dt %>% mutate(dob = mdy(dob), dx_dt = mdy(dx_dt), date = mdy(date))
# remove measurement dates prior to dob
dt <- dt %>% mutate(exc_dt = if_else(date >= dob | is.na(date), 0, 1)) %>% subset(exc_dt==0) # (1 measurement date before dob)
# fix incorrect measurement date
dt <- dt %>% mutate(date = if_else(subjid==120573 & date=='2050-08-21', mdy(08212005), date)) %>% ungroup() # (2050 vs 2005)
# calculate ages
dt <- dt %>% mutate(agedays = (dob %--% date) / days(1), agemo = round((dob%--%date)/months(1), digits=3), age = round((dob%--%date)/years(1), digits=3))

# create subject list w/ max age
list <- dt %>% arrange(subjid) %>% group_by(subjid) %>% mutate(max_ad = max(agedays)) %>% ungroup() 
list <- list %>% distinct(subjid, .keep_all = TRUE) %>% subset(select=c(subjid, dob, max_ad)) # n=1161

# testo
dtt <- read_csv("...PATHWAY.../KS_T.Marshall_Testos.csv")  # UPDATE ME :)))))
dtt <- rename(dtt, subjid=person_id, dx_dt=first_dx_dt, t_name=concept_name, t_st=drug_exposure_start_date, t_end=drug_exposure_end_date, t=testos)
dtt <- dtt %>% subset(select=c(subjid, dob, dx_dt, kline, t_name, t_st, t_end, t)) %>% subset(kline=='Case') 
dtt <- dtt %>% mutate(dob = mdy(dob), dx_dt = mdy(dx_dt), t_st = mdy(t_st), t_end = mdy(t_end), 
                      st_ad = (dob%--%t_st)/days(1), end_ad = (dob%--%t_end)/days(1), t = if_else(t=='Y', 1, 0))
# calculate n and row number
dtt <- dtt %>% arrange(subjid, t_st) %>% group_by(subjid) %>% mutate(n = n(), row = row_number()) %>% ungroup()
# add max agedays from dt list
dtt <- left_join(dtt, list, by = c('subjid', 'dob'))
# fill missing end_ad and st_ad with 0 to facilitate merge 
dtt <- dtt %>% mutate(end_ad = if_else(is.na(end_ad), 0, end_ad), st_ad = if_else(is.na(st_ad), 0, st_ad))
# ASSUMPTION - if on t and last rx end=start (only n>1) or has no end (any n), replace end w/ max of dt agedays vs t start (i.e., assume on t through latest date)
dtt <- dtt %>% arrange(subjid, st_ad) %>% mutate(end_ad = if_else(t==1 & row==n & n>1 & end_ad==st_ad, pmax(max_ad, st_ad), end_ad),
                                                 end_ad = if_else(t==1 & row==n & end_ad==0, pmax(max_ad, st_ad), end_ad)) 
# AUUMPTION - other than final rx, if t end=start or has no end (only n>1), replace end w/ subsequent rx start 
dtt <- dtt %>% arrange(subjid, st_ad) %>% mutate(end_ad = if_else(n>1 & row<n & end_ad==0, st_ad, end_ad)) %>% 
  group_by(subjid) %>% mutate(end_ad = if_else(n>1 & row<n & end_ad==st_ad, shift(st_ad, n=-1), end_ad)) %>% ungroup() 
# merge rows w/ overlapping t start/end of any t formulation ((NOTE: t_name summarized to allow merge))
dtt <- dtt %>% arrange(subjid, st_ad) %>% group_by(subjid) %>% mutate(bn = c(0, cumsum(lead(st_ad) > cummax(end_ad))[-n()])) %>% 
  group_by(subjid, dob, dx_dt, kline, t, max_ad, bn) %>% 
  summarise(mrg=n(), t_name=first(t_name), t_st=min(t_st), t_end=max(t_end), st_ad=min(st_ad), end_ad=max(end_ad)) %>% ungroup() 
# recalculate n and row 
dtt <- dtt %>% arrange(subjid, st_ad) %>% group_by(subjid) %>% mutate(n = n(), row = row_number()) %>% ungroup() 
# ASSUMPTION - if >12.5yo (4562d) and on t can disregard gap in use (since clinically would be rare to discontinue)
dtt <- dtt %>% group_by(subjid) %>% mutate(gap_i = if_else(row>1, st_ad - shift(end_ad), 0), # amount of time off t between subsequent rx
                                           com = if_else(n>1 & end_ad>4562, n(), row_number())) %>%
  group_by(subjid, dob, dx_dt, kline, t, max_ad, com) %>% 
  summarise(mrg=sum(mrg), t_name=first(t_name), t_st=min(t_st), t_end=max(t_end), st_ad=min(st_ad), end_ad=max(end_ad), gap_i=max(gap_i)) %>% 
  ungroup() 
# re-calculate n and row 
dtt <- dtt %>% arrange(subjid, st_ad) %>% group_by(subjid) %>% mutate(n = n(), row = row_number()) %>% ungroup() 
# ASSUMPTION - if t use > 1yr can carry forward to max age (to account for lingering effect)
dtt <- dtt %>% mutate(t1y = if_else(end_ad - st_ad >= 365 & end_ad <= pmax(max_ad, end_ad), 1, 0)) %>% group_by(subjid) %>% 
  mutate(gap_o = if_else(t1y==1, max(max_ad, end_ad) - end_ad, 0), # amount of time off t between rx end and max age
         end_ad = if_else(t1y==1, max(max_ad, end_ad), end_ad), # replace end with max age/max end if t rx length > 1yr
         bn = c(0, cumsum(lead(st_ad) > cummax(end_ad))[-n()])) %>% 
  group_by(subjid, dob, dx_dt, kline, t, max_ad, bn) %>% 
  summarise(mrg=sum(mrg), t_name=first(t_name), t_st=min(t_st), t_end=max(t_end), st_ad=min(st_ad), end_ad=max(end_ad), gap_i=max(gap_i), gap_o=max(gap_o)) %>% 
  ungroup() # n=1162 (i.e., 1 patient w/ distinct t use (60692))
# re-calculate n and row 
dtt <- dtt %>% arrange(subjid, st_ad) %>% group_by(subjid) %>% mutate(n = n(), row = row_number()) %>% ungroup() 
# convert to wide format 
dtt_w <- pivot_wider(dtt, id_cols = c(subjid, dob, dx_dt, kline, t, n), names_from = row, 
                     names_sep = "_", values_from = c(st_ad, end_ad), values_fill = 0) # 1161
# combine dt and dtt_w
dt <- left_join(dt, dtt_w, by = c('subjid', 'dob', 'dx_dt', 'kline')) 
# determine t status by age
dt <- dt %>% mutate(t1 = if_else(t==1 & agedays > st_ad_1 & agedays <= end_ad_1 & st_ad_1 > 0, 1, 0), 
                    t2 = if_else(t==1 & agedays > st_ad_2 & agedays <= end_ad_2 & st_ad_2 > 0, 1, 0)) 
dt <- dt %>% mutate(t = (t1+t2))
table(dt$t) # (2402 t=Y)

# simplify dt
dt <- subset(dt, select = c(subjid, site, dob, dx_dt, kline, race, ethnicity, date, agedays, agemo, age, wgt, hgt, t))
# remove patients with 0 measures of hgt/wgt 
dt <- subset(dt, !is.na(date))

# subset data
# exclusions - age > 20.5
dt <- subset(dt, age>=0 & age<=20.5) # 15061 (all exclusions d/t age>20.5))
# exclusions - additional genetic dx (all data)
dt <- subset(dt, subjid!='47706' & subjid!='38877' & subjid!='34438') # (+/- consider id 13652)
# exclusions - premies (only 0-2yo data)
premdt <- read_csv("...PATHWAY.../KS_IDs_premi_birth.csv")  # UPDATE ME :)))))
dt <- left_join(dt, premdt, by = c('subjid')) 
dt <- dt %>% mutate(pre2 = if_else(is.na(pre2), 0, pre2)) %>% subset(age>=2 | age<2 & pre2==0) # (premature dx) ((NOTE: 144413, 84035, 148794 not on list))
dt <- subset(dt, age>=2 | age<2 & subjid!=144413 & subjid!=84035 & subjid!=188908 & subjid!=180866 & subjid!=80580 & subjid!=49876 & 
               subjid!=95375 & subjid!=148794 & subjid!=22283) # (len<40cm and/or wgt<1.5kg); (reviewed id 183711, not premie) 

# set up height 
dth <- subset(dt, select=-c(wgt, pre2)) 
dth <- dth %>% mutate(param = 'HEIGHTCM', sex = 0)
dth <- rename(dth, measurement=hgt)
# set up weight 
dtw <- subset(dt, select=-c(hgt, pre2)) 
dtw <- dtw %>% mutate(param = 'WEIGHTKG', sex = 0)
dtw <- rename(dtw, measurement=wgt)
# combine data tables 
dt <- rbindlist(list(dth,dtw)) 

# id & remove missing data
dt <- dt %>% mutate(exc_md = if_else(is.na(measurement), 1, 0))
dt <- subset(dt, exc_md == 0) 
# id & remove duplicates 
dt <- dt %>% arrange(subjid, param, agedays) %>% group_by(subjid, param, measurement) %>% mutate(dup = n(), dup_id = row_number()) %>%
  mutate(exc_dup = if_else(param=='HEIGHTCM' & dup_id>1 & min(age)<16 | # hgt - remove subsequent duplicates if age <16
                             param=='WEIGHTKG' & dup_id>1, 1, 0)) %>% ungroup() # wgt - remove all duplicates 
dt <- subset(dt, exc_dup == 0) 
# add row number
dt <- dt %>% arrange(subjid, param, agedays) %>% group_by(param) %>% mutate(row=row_number()) %>% group_by(subjid, param) %>% mutate(srow=row_number()) %>% ungroup()
# simplify dataset & confirm structure
dt <- subset(dt, select=c(subjid, site, dob, kline, dx_dt, race, ethnicity, sex, row, srow, date, age, agemo, agedays, param, measurement, t))
skim(dt) 

# prepare data as data.table and set key (sort by subjid, param, agedays)
dt <- setDT(dt)[with(dt, order(subjid, param, agedays))]
setkey(dt, subjid, param, agedays)
# create column w/ cleangrowth flag & summarize results
dt <- dt[, gcr_result := cleangrowth(subjid, param, agedays, sex, measurement, 
                                     recover.unit.error = TRUE, 
                                     include.carryforward = TRUE, 
                                     error.load.threshold = 0.75)]
dt %>% group_by(gcr_result) %>% tally(sort=TRUE) 

# create column to review gcr exclusions & note if measurement value changed (vc)
dt <- dt %>% mutate(exc_gcr=0, vc=0) 
# remove 1 - all errors for site 1453 if first entry 
dt <- dt %>% mutate(exc_gcr = if_else(site=='site_1453' & srow==1 & gcr_result!='Include', 1, exc_gcr)) 
table(dt$exc_gcr) # n(1)=59
dt_rev <- subset(dt, exc_gcr==1)
# review 2 - unit errors
dt <- dt %>% mutate(exc_gcr = if_else(exc_gcr==0 & gcr_result %in% c('Unit-Error-High', 'Unit-Error-Low'), 2, exc_gcr)) 
table(dt$exc_gcr) # n(2)=4
dt_rev2 <- subset(dt, exc_gcr==2)
# correct data when able (4)
dt_rev2 <- dt_rev2 %>% mutate(measurement = if_else(subjid==41233 & measurement==41.801, 106.17, measurement), 
                              measurement = if_else(subjid==135071 & measurement==10.80, 23.81, measurement), 
                              measurement = if_else(subjid==118175 & measurement==30.80, 13.97, measurement), 
                              measurement = if_else(subjid==130443 & measurement==102.0, 46.27, measurement))
dt_rev2 <- dt_rev2 %>% mutate(exc_gcr=0, vc=1)
# review 3 - SD cutoff 
dt <- dt %>% mutate(exc_gcr = if_else(exc_gcr==0 & gcr_result %in% c('Exclude-SD-Cutoff'), 3, exc_gcr)) 
table(dt$exc_gcr) # n(3)=5
dt_rev3 <- subset(dt, exc_gcr==3)
# correct data when able (3) # confirm remove (2) # non-physio (1) 465.099 # same as wgt (1) 14.9
dt_rev3 <- dt_rev3 %>% mutate(measurement = if_else(subjid==10373 & measurement==2737.97, 2.737, measurement), 
                              measurement = if_else(subjid==32244 & measurement==3520, 3.52, measurement), 
                              measurement = if_else(subjid==180591 & measurement==3390, 3.39, measurement))
dt_rev3 <- dt_rev3 %>% mutate(exc_gcr = if_else(measurement %in% c(2.737, 3.52, 3.39), 0, exc_gcr), 
                              vc = if_else(measurement %in% c(2.737, 3.52, 3.39), 1, vc))
# review 4 - ewma-extreme
dt <- dt %>% mutate(exc_gcr = if_else(exc_gcr==0 & gcr_result %in% c('Exclude-EWMA-Extreme', 'Exclude-EWMA-Extreme-Pair'), 4, exc_gcr)) 
table(dt$exc_gcr) # n(4)=17
dt_rev4 <- subset(dt, exc_gcr==4)
# correct data when able (4)
dt_rev4 <- dt_rev4 %>% mutate(measurement = if_else(subjid==57924 & measurement==19.05, 48.39, measurement), 
                              measurement = if_else(subjid==57924 & measurement==7.6, 3.45, measurement), 
                              measurement = if_else(subjid==141400 & measurement==66.5, 30.16, measurement),
                              measurement = if_else(subjid==204614 & measurement==81.5, 181.5, measurement))
dt_rev4 <- dt_rev4 %>% mutate(exc_gcr = if_else(measurement %in% c(48.39, 3.45, 30.16, 181.5), 0, exc_gcr), 
                              vc = if_else(measurement %in% c(48.39, 3.45, 30.16, 181.5), 1, vc)) # confirm remove (13) # non-physio (13 - inc. 2 w/ site_1453 error but 2nd entry not 1st)
# review 5 - other
dt <- dt %>% mutate(exc_gcr = if_else(exc_gcr==0 & gcr_result %in% c('Exclude-Pair-Delta-17', 'Exclude-Pair-Delta-18'), 5, exc_gcr)) 
table(dt$exc_gcr) # n(5)=3
dt_rev5 <- subset(dt, exc_gcr==5) # remove all (non-physio 3) 
# review 6 - ewma if first entry (i.e., ewma-9)
dt <- dt %>% mutate(exc_gcr = if_else(exc_gcr==0 & gcr_result %in% c('Exclude-EWMA-9'), 6, exc_gcr))
table(dt$exc_gcr) # n(6)=22
dt_rev6 <- subset(dt, exc_gcr==6)
# keep all but 2 (most only flagging d/t inc. time between next entry) (non-physio 2)
dt_rev6 <- dt_rev6 %>% mutate(exc_gcr = if_else(subjid==122248 & measurement==53.34 | subjid==177317 & measurement==13.0, exc_gcr, 0))
# review 7 - ewma remaining (ewma-11 last entry, ewma-8 middle entry)
dt <- dt %>% mutate(exc_gcr = if_else(exc_gcr==0 & gcr_result %in% c('Exclude-EWMA-8', 'Exclude-EWMA-11'), 7, exc_gcr)) 
table(dt$exc_gcr) # n(7)=56
dt_rev7 <- subset(dt, exc_gcr==7) # remove all (non-physio 56)
# review 8 - max hgt change
dt <- dt %>% mutate(exc_gcr = if_else(exc_gcr==0 & gcr_result %in% c('Exclude-Max-Height-Change'), 8, exc_gcr)) 
table(dt$exc_gcr) # n(8)=3
dt_rev8 <- subset(dt, exc_gcr==8) # remove all (non-physio 3)
# review 9 - min hgt change NOTE: remove all (difficult to remain objective in review)
dt <- dt %>% mutate(exc_gcr = if_else(exc_gcr==0 & gcr_result %in% c('Exclude-Min-Height-Change'), 9, exc_gcr)) 
table(dt$exc_gcr) # n(9)=118
dt_rev9 <- subset(dt, exc_gcr==9) # remove all (non-physio 118)

# combine exclusion reviews
dt_ex <- rbindlist(list(dt_rev, dt_rev2, dt_rev3, dt_rev4, dt_rev5, dt_rev6, dt_rev7, dt_rev8, dt_rev9), use.names=TRUE) 
# create inclusion dt subset
dt_in <- subset(dt, exc_gcr==0) 
# create new dt w/ inclusions/exclusions 
dt <- rbindlist(list(dt_in, dt_ex), use.names=TRUE) 
# remove exclusions & re-order 
dt <- dt %>% subset(exc_gcr==0) %>% subset(select = -c(exc_gcr, row, srow)) %>% 
  arrange(subjid, param, agedays) 

# calculate age at diagnosis & make 0 if prior to birth (1 patient)
dt <- dt %>% mutate(dx_age = round((dob%--%dx_dt)/years(1), digits=3), 
                    dx_age = if_else(dx_age < 0, 0, dx_age))
# add t as character
dt <- dt %>% mutate(testo = if_else(t==1, 'Y', 'N'))
# re-calculate n, row number, & confirm structure
dt <- dt %>% mutate(n=n()) %>% group_by(param) %>% mutate(row=row_number()) %>%
  group_by(param, subjid) %>% mutate(srow=row_number()) %>% ungroup()
skim(dt) # 23902


# HEIGHT
# height - subset 
dt_hgt <- subset(dt, param=='HEIGHTCM' & age>1.5 & age<20.5) # 9567 - FUCKED! (no monotone)
# height - calculate idmw
dt_hgt <- dt_hgt %>% group_by(subjid) %>% mutate(n = n()) %>% ungroup()
dt_hgt <- dt_hgt %>% mutate(idmw = 1 / (1 + abs(median(tapply(n, subjid, median)) - n)))
# re-calculate row number
dt_hgt <- dt_hgt %>% arrange(subjid, agedays) %>% mutate(row=row_number()) %>% group_by(subjid) %>% mutate(srow=row_number()) %>% ungroup()
# height - ggplot
ggplot(data=dt_hgt, aes(x=age, y=measurement, colour=t)) + ylab("height") +
  geom_point(size=0.8,alpha=0.7)+geom_smooth(se=FALSE,linewidth=0.5,colour="red")


# WEIGHT
# weight - subset 
dt_wgt <- subset(dt, param=='WEIGHTKG' & age>1.5 & age<20.5) # 10872 - DECENT!
dt_wgt <- subset(dt_wgt, subjid!=147807) # NOTE: remove rapid dec patient if no mrg
# weight - calculate idmw
dt_wgt <- dt_wgt %>% group_by(subjid) %>% mutate(n = n()) %>% ungroup()
dt_wgt <- dt_wgt %>% mutate(idmw = 1 / (1 + abs(median(tapply(n, subjid, median)) - n)))
# re-calculate row number
dt_wgt <- dt_wgt %>% arrange(subjid, agedays) %>% mutate(row=row_number()) %>% group_by(subjid) %>% mutate(srow=row_number()) %>% ungroup()
# weight - ggplot
ggplot(data=dt_wgt, aes(x=age, y=measurement, colour=t)) + ylab("weight") +
  geom_point(size=0.8,alpha=0.7)+geom_smooth(se=FALSE,linewidth=0.5,colour="red")


# BMI
# bmi - subset
dt_bmi <- subset(dt, age>1.5 & age<20.5) # - DECENT! 
# bmi - id same day measurements 
dt_bmi <- dt_bmi %>% group_by(subjid, agedays) %>% mutate(n = n()) %>% ungroup() 
dt_bmi <- subset(dt_bmi, n==2) 
# bmi - convert to wide format
dt_bmi <- pivot_wider(dt_bmi, id_cols=c(subjid, site, dob, kline, dx_dt, race, ethnicity, sex, date, age, agemo, agedays, t, testo, dx_age), 
                      names_from=param, values_from=c(measurement)) 
dt_bmi <- dt_bmi[complete.cases(dt_bmi),] 
# bmi - calculate bmi & add param
dt_bmi <- dt_bmi %>% mutate(bmi = (WEIGHTKG/HEIGHTCM/HEIGHTCM)*10000, param = 'BMI')
# bmi - calculate idmw
dt_bmi <- dt_bmi %>% group_by(subjid) %>% mutate(n = n()) %>% ungroup()
dt_bmi <- dt_bmi %>% mutate(idmw = 1 / (1 + abs(median(tapply(n, subjid, median)) - n)))
# re-calculate row number
dt_bmi <- dt_bmi %>% arrange(subjid, agedays) %>% mutate(row=row_number()) %>% group_by(subjid) %>% mutate(srow=row_number()) %>% ungroup()
# bmi - ggplot
ggplot(data=dt_bmi, aes(x=age, y=bmi, colour=t)) + ylab("bmi") +
  geom_point(size=0.8,alpha=0.7)+geom_smooth(se=FALSE,linewidth=0.5,colour="red")


# LENGTH
# length - subset 
dt_len <- subset(dt, param=='HEIGHTCM' & age<2.5) # 1828 - DECENT! (but monotone not working)
# length - calculate idmw
dt_len <- dt_len %>% group_by(subjid) %>% mutate(n = n()) %>% ungroup()
dt_len <- dt_len %>% mutate(idmw = 1 / (1 + abs(median(tapply(n, subjid, median)) - n)))
# re-calculate row number
dt_len <- dt_len %>% arrange(subjid, agedays) %>% mutate(row=row_number()) %>% group_by(subjid) %>% mutate(srow=row_number()) %>% ungroup()
# length - ggplot
ggplot(data=dt_len, aes(x=age, y=measurement, colour=t)) + ylab("length") +
  geom_point(size=0.8,alpha=0.7)+geom_smooth(se=FALSE,linewidth=0.5,colour="red")


# BIRTHWEIGHT 
# birthweight - subset 
dt_bw <- subset(dt, param=='WEIGHTKG' & age<2.5) # 2811 - DECENT!
# birthweight - calculate idmw
dt_bw <- dt_bw %>% group_by(subjid) %>% mutate(n = n()) %>% ungroup()
dt_bw <- dt_bw %>% mutate(idmw = 1 / (1 + abs(median(tapply(n, subjid, median)) - n)))
# re-calculate row number
dt_bw <- dt_bw %>% arrange(subjid, agedays) %>% mutate(row=row_number()) %>% group_by(subjid) %>% mutate(srow=row_number()) %>% ungroup()
# birthweight - ggplot
ggplot(data=dt_bw, aes(x=age, y=measurement, colour=t)) + ylab("weight") +
  geom_point(size=0.8,alpha=0.7)+geom_smooth(se=FALSE,linewidth=0.5,colour="red")


# WEIGHT FOR LENGTH
# weight for length - subset 
dt_wl <- subset(dt, age<=2) # 4034 - DECENT!
# weight for length - id same day measurements 
dt_wl <- dt_wl %>% group_by(subjid, agedays) %>% mutate(n = n()) %>% ungroup() 
dt_wl <- subset(dt_wl, n==2) 
# weight for length - convert to wide format
dt_wl <- pivot_wider(dt_wl, id_cols=c(subjid, agemo, agedays, t), 
                     names_from=param, values_from=c(measurement)) 
dt_wl <- dt_wl[complete.cases(dt_wl),] 
# weight for length - calculate idmw
dt_wl <- dt_wl %>% group_by(subjid) %>% mutate(n = n()) %>% ungroup()
dt_wl <- dt_wl %>% mutate(idmw = 1 / (1 + abs(median(tapply(n, subjid, median)) - n)))
# re-calculate row number & add param 
dt_wl <- dt_wl %>% arrange(subjid, agedays) %>% mutate(param='WFL', row=row_number()) %>% group_by(subjid) %>% mutate(srow=row_number()) %>% ungroup()
# weight for length - ggplot
ggplot(data=dt_wl, aes(x=HEIGHTCM, y=WEIGHTKG, colour=t)) + xlab("length") + ylab("weight") +
  geom_point(size=0.8,alpha=0.7)+geom_smooth(se=FALSE,linewidth=0.5,colour="red")


# CDC & WHO REFERENCE DATA
# upload and subset CDC & WHO datasets 
cdc_who <- read_csv("...PATHWAY.../cdc_who_data.csv")  # UPDATE ME :)))))
cdc_hgt <- subset(cdc_who, source=='CDC' & param=='HEIGHTCM')
cdc_wgt <- subset(cdc_who, source=='CDC' & param=='WEIGHTKG')
cdc_bmi <- subset(cdc_who, source=='CDC' & param=='BMI')
who_len <- subset(cdc_who, source=='WHO' & param=='LENGTHCM')
who_wgt <- subset(cdc_who, source=='WHO' & param=='WEIGHTKG')
who_wl <- subset(cdc_who, source=='WHO' & param=='WFL')


# HEIGHT
# height - gcrq formula 
tau2=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
mform_dt_hgt <-gcrq(data=dt_hgt, measurement ~ ps(age, monotone=1) + t + idmw, tau=tau2) # MONOTONE NOT WORKING???
mform_dt_hgt <-gcrq(data=dt_hgt, measurement ~ ps(age) + idmw + t, tau=tau2) 
# height - plot percentiles NOTE: revisit to help the uptick at the 95th percentile
plot(mform_dt_hgt, term=1, res=TRUE, ylab='height', lwd=2)

# plot - KS HEIGHT  
plot.new()
par(mar=c(3.1, 4.1, 3.1, 2.6) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(2,20), ylim=c(72,198), xaxs="i", yaxs="i")
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="2 to 20 years: Klinefelter Syndrome Males", adj=0, line=1.75, cex.main=1.8, col.main="deepskyblue4")
title(main='Stature-for-age percentiles', adj=0, line=0.45, cex.main=1.8, col.main="deepskyblue4")
title(xlab='AGE (YEARS)', line=1.95, cex.lab=1.7, font.lab=2, col.lab="deepskyblue4")
title(ylab='STATURE (CM)', line=2.85, cex.lab=1.7, font.lab=2, col.lab="deepskyblue4")
# grid marks
axis(1, at=seq(2,20,0.5), lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=72:198, lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=2:20, lwd=0.9, tck=1, cex.axis=1.6, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=seq(75,195,5), lwd=0.9, tck=1, cex.axis=1.6, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=seq(75,195,5), lwd=0.9, tck=1, cex.axis=1.6, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# curves - NOTE: re-run code 3x to darken the plot lines (idfkw..)
plot(mform_dt_hgt, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col="deepskyblue4")
plot(mform_dt_hgt, select.tau=c(0.95, 0.5, 0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
box(lwd=1, col="deepskyblue4")

# labels - need to fix (idk how to automate location w/ line xy coordinates)
x2_20 = c(19.5)
x2_20a <- c(19.4, 19.4, 19.4, 19.4, 19.4, 19.4)
x2_20b <- c(19.6, 19.6, 19.6, 19.6, 19.6, 19.6)
yhgt <- c(168.7, 174.8, 181.2, 185.9, 189.5, 193.3)
lbl <- c("10", "25", "50", "75", "90", "95")
points(x2_20a, yhgt, pch=15, cex=2.4, col="white")
points(x2_20b, yhgt, pch=15, cex=2.4, col="white")
text(x2_20, yhgt, labels=lbl, adj=0.5, cex=1.3, font=1, col="deepskyblue4")
points(x2_20, 161.2, pch=15, cex=1.8, col="white")
points(x2_20, 161.8, pch=15, cex=1.8, col="white")
text(x2_20, 161.5, labels='5', adj=0.5, cex=1.3, font=1, col="deepskyblue4")
# save as hgt_kscurve

# cdc height - gcrq formula
cdc_hgt <- subset(cdc_hgt, select=c(source, param, age, agemo, p03, p05, p10, p25, p50, p75, p90, p95, p97))
cdc_hgt.form05 <- gcrq(data=cdc_hgt, p05 ~ ps(age), tau=0.5)
cdc_hgt.form50 <- gcrq(data=cdc_hgt, p50 ~ ps(age), tau=0.5)
cdc_hgt.form95 <- gcrq(data=cdc_hgt, p95 ~ ps(age), tau=0.5)

# plot - KS HEIGHT + CDC SHADED 
plot.new()
par(mar=c(3.1, 4.1, 3.1, 2.6) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(2,20), ylim=c(72,198), xaxs="i", yaxs="i")
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="2 to 20 years: Klinefelter Syndrome Males", adj=0, line=1.75, cex.main=1.8, col.main="deepskyblue4")
title(main='Stature-for-age percentiles', adj=0, line=0.45, cex.main=1.8, col.main="deepskyblue4")
title(xlab='AGE (YEARS)', line=1.95, cex.lab=1.7, font.lab=2, col.lab="deepskyblue4")
title(ylab='STATURE (CM)', line=2.85, cex.lab=1.7, font.lab=2, col.lab="deepskyblue4")
# grid marks
axis(1, at=seq(2,20,0.5), lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=72:198, lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=2:20, lwd=0.9, tck=1, cex.axis=1.6, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=seq(75,195,5), lwd=0.9, tck=1, cex.axis=1.6, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=seq(75,195,5), lwd=0.9, tck=1, cex.axis=1.6, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# height cdc polygon
polygon(c(cdc_hgt$age, rev(cdc_hgt$age)), c(cdc_hgt$p95, rev(cdc_hgt$p05)), add=TRUE, lty=0, col=alpha("darkorange3", 0.3))
# curves (re-run the code 3x to darken the plot lines)
plot(mform_dt_hgt, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col="deepskyblue4")
plot(mform_dt_hgt, select.tau=c(0.95, 0.5, 0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
box(lwd=1, col="deepskyblue4")

# plot - KS HEIGHT 95% + CDC
plot.new()
par(mar=c(3.1, 4.1, 3.1, 2.6) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(2,20), ylim=c(72,198), xaxs="i", yaxs="i")
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="2 to 20 years: Klinefelter Syndrome Males", adj=0, line=1.75, cex.main=1.8, col.main="deepskyblue4")
title(main='Stature-for-age percentiles', adj=0, line=0.45, cex.main=1.8, col.main="deepskyblue4")
title(xlab='AGE (YEARS)', line=1.95, cex.lab=1.7, font.lab=2, col.lab="deepskyblue4")
title(ylab='STATURE (CM)', line=2.85, cex.lab=1.7, font.lab=2, col.lab="deepskyblue4")
# grid marks
axis(1, at=seq(2,20,0.5), lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=72:198, lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=2:20, lwd=0.9, tck=1, cex.axis=1.6, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=seq(75,195,5), lwd=0.9, tck=1, cex.axis=1.6, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=seq(75,195,5), lwd=0.9, tck=1, cex.axis=1.6, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# 95%ile curves 
plot(mform_dt_hgt, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_hgt, select.tau=c(0.5, 0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_hgt, select.tau=c(0.95), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_hgt, select.tau=c(0.95), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_hgt, select.tau=c(0.95), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(cdc_hgt.form95, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_hgt.form95, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_hgt.form95, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
box(lwd=1, col="deepskyblue4")
# labels
x2_20a <- c(19.4, 19.4)
x2_20b <- c(19.6, 19.6)
yhgt <- c(188.5, 193.3)
points(x2_20a, yhgt, pch=15, cex=2.4, col="white")
points(x2_20b, yhgt, pch=15, cex=2.4, col="white")
text(19.5, 193.3, labels='95', adj=0.5, cex=1.3, font=1, col="deepskyblue4")
text(19.5, 188.5, labels='95', adj=0.5, cex=1.3, font=1, col="darkorange3")

# plot - KS HEIGHT 50% + CDC
plot.new()
par(mar=c(3.1, 4.1, 3.1, 2.6) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(2,20), ylim=c(72,198), xaxs="i", yaxs="i")
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="2 to 20 years: Klinefelter Syndrome Males", adj=0, line=1.75, cex.main=1.8, col.main="deepskyblue4")
title(main='Stature-for-age percentiles', adj=0, line=0.45, cex.main=1.8, col.main="deepskyblue4")
title(xlab='AGE (YEARS)', line=1.95, cex.lab=1.7, font.lab=2, col.lab="deepskyblue4")
title(ylab='STATURE (CM)', line=2.85, cex.lab=1.7, font.lab=2, col.lab="deepskyblue4")
# grid marks
axis(1, at=seq(2,20,0.5), lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=72:198, lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=2:20, lwd=0.9, tck=1, cex.axis=1.6, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=seq(75,195,5), lwd=0.9, tck=1, cex.axis=1.6, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=seq(75,195,5), lwd=0.9, tck=1, cex.axis=1.6, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# 50%ile curves 
plot(mform_dt_hgt, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_hgt, select.tau=c(0.95, 0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_hgt, select.tau=c(0.50), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_hgt, select.tau=c(0.50), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_hgt, select.tau=c(0.50), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(cdc_hgt.form50, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_hgt.form50, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_hgt.form50, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
box(lwd=1, col="deepskyblue4")
# labels
yhgt <- c(181.2, 176.8)
points(x2_20a, yhgt, pch=15, cex=2.4, col="white")
points(x2_20b, yhgt, pch=15, cex=2.4, col="white")
text(19.5, 181.2, labels='50', adj=0.5, cex=1.3, font=1, col="deepskyblue4")
text(19.5, 176.8, labels='50', adj=0.5, cex=1.3, font=1, col="darkorange3")

# plot - KS HEIGHT 5% + CDC
plot.new()
par(mar=c(3.1, 4.1, 3.1, 2.6) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(2,20), ylim=c(72,198), xaxs="i", yaxs="i")
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="2 to 20 years: Klinefelter Syndrome Males", adj=0, line=1.75, cex.main=1.8, col.main="deepskyblue4")
title(main='Stature-for-age percentiles', adj=0, line=0.45, cex.main=1.8, col.main="deepskyblue4")
title(xlab='AGE (YEARS)', line=1.95, cex.lab=1.7, font.lab=2, col.lab="deepskyblue4")
title(ylab='STATURE (CM)', line=2.85, cex.lab=1.7, font.lab=2, col.lab="deepskyblue4")
# grid marks
axis(1, at=seq(2,20,0.5), lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=72:198, lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=2:20, lwd=0.9, tck=1, cex.axis=1.6, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=seq(75,195,5), lwd=0.9, tck=1, cex.axis=1.6, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=seq(75,195,5), lwd=0.9, tck=1, cex.axis=1.6, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# 5%ile curves 
plot(mform_dt_hgt, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_hgt, select.tau=c(0.95, 0.50), axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_hgt, select.tau=c(0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_hgt, select.tau=c(0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_hgt, select.tau=c(0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(cdc_hgt.form05, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_hgt.form05, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_hgt.form05, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
box(lwd=1, col="deepskyblue4")
# labels
x2_20 = c(19.5, 19.5)
yhgt1 <- c(161.2, 161.8)
yhgt2 <- c(164.7, 165.3)
points(x2_20, yhgt1, pch=15, cex=1.8, col="white")
points(x2_20, yhgt2, pch=15, cex=1.8, col="white")
text(19.5, 161.5, labels='5', adj=0.5, cex=1.3, font=1, col="deepskyblue4")
text(19.5, 165, labels='5', adj=0.5, cex=1.3, font=1, col="darkorange3")


# WEIGHT 
# weight - gcrq formula - NOTE: may need to correct min/max lambda 
tau2=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
mform_dt_wgt <-gcrq(data=dt_wgt, measurement ~ ps(age, lambda=5) + t + idmw, tau=tau2)
# weight - plot percentiles 
plot(mform_dt_wgt, term=1, res=TRUE, ylab='weight', lwd=2)

# plot - KS WEIGHT 
plot.new()
par(mar=c(3,4,3,2.5) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(2,20), ylim=c(2,143), xaxs="i", yaxs="i")
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="2 to 20 years: Klinefelter Syndrome Males", adj=0, line=1.55, cex.main=1.5, col.main="deepskyblue4")
title(main='Weight-for-age percentiles', adj=0, line=0.4, cex.main=1.5, col.main="deepskyblue4")
title(xlab='AGE (YEARS)', line=1.9, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
title(ylab='WEIGHT (KG)', line=2.75, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
# grid marks
axis(1, at=seq(2,20,0.5), lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=2:143, lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=2:20, lwd=0.9, tck=1, cex.axis=1.5, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=seq(0,140,5), lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=seq(0,140,5), lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# curves - (re-run the code 3x to darken the plot lines)
plot(mform_dt_wgt, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col="deepskyblue4")
plot(mform_dt_wgt, select.tau=c(0.95, 0.5, 0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
box(lwd=1, col="deepskyblue4")

# labels - need to fix (idk how to automate location w/ line xy coordinates)
x2_20 = c(19.5)
x2_20a <- c(19.4, 19.4, 19.4)
x2_20b <- c(19.6, 19.6, 19.6)
ywgt <- c(54.5, 78.8, 136.9)
lbl <- c("5", "50", "95")
points(x2_20a, ywgt, pch=15, cex=2.4, col="white")
points(x2_20b, ywgt, pch=15, cex=2.4, col="white")
text(x2_20, ywgt, labels=lbl, adj=0.5, cex=1.3, font=1, col="deepskyblue4")
# save as wgt_kscurve

# cdc weight - gcrq formula
cdc_wgt <- subset(cdc_wgt, select=c(source, param, age, agemo, p03, p05, p10, p25, p50, p75, p90, p95, p97))
cdc_wgt.form05 <- gcrq(data=cdc_wgt, p05 ~ ps(age), tau=0.5)
cdc_wgt.form50 <- gcrq(data=cdc_wgt, p50 ~ ps(age), tau=0.5)
cdc_wgt.form95 <- gcrq(data=cdc_wgt, p95 ~ ps(age), tau=0.5)

# plot - KS WEIGHT + CDC SHADED 
plot.new()
par(mar=c(3,4,3,2.5) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(2,20), ylim=c(2,143), xaxs="i", yaxs="i")
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="2 to 20 years: Klinefelter Syndrome Males", adj=0, line=1.55, cex.main=1.5, col.main="deepskyblue4")
title(main='Weight-for-age percentiles', adj=0, line=0.4, cex.main=1.5, col.main="deepskyblue4")
title(xlab='AGE (YEARS)', line=1.9, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
title(ylab='WEIGHT (KG)', line=2.75, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
# grid marks
axis(1, at=seq(2,20,0.5), lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=2:143, lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=2:20, lwd=0.9, tck=1, cex.axis=1.5, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=seq(0,140,5), lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=seq(0,140,5), lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# weight cdc polygon
polygon(c(cdc_wgt$age, rev(cdc_wgt$age)), c(cdc_wgt$p95, rev(cdc_wgt$p05)), add=TRUE, lty=0, col=alpha("darkorange3", 0.3))
# curves (re-run the code 3x to darken the plot lines)
plot(mform_dt_wgt, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col="deepskyblue4")
plot(mform_dt_wgt, select.tau=c(0.95, 0.5, 0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
box(lwd=1, col="deepskyblue4")

# plot - KS WEIGHT 95% + CDC
plot.new()
# 95%ile curves 
plot(mform_dt_wgt, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_wgt, select.tau=c(0.5, 0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_wgt, select.tau=c(0.95), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_wgt, select.tau=c(0.95), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_wgt, select.tau=c(0.95), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(cdc_wgt.form95, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_wgt.form95, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_wgt.form95, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
box(lwd=1, col="deepskyblue4")

# plot - KS WEIGHT 50% + CDC
plot.new()
# 50%ile curves 
plot(mform_dt_wgt, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_wgt, select.tau=c(0.95, 0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_wgt, select.tau=c(0.5), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_wgt, select.tau=c(0.5), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_wgt, select.tau=c(0.5), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(cdc_wgt.form50, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_wgt.form50, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_wgt.form50, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
box(lwd=1, col="deepskyblue4")

# plot - KS WEIGHT 5% + CDC
plot.new()
# 5%ile curves 
plot(mform_dt_wgt, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_wgt, select.tau=c(0.95, 0.5), axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_wgt, select.tau=c(0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_wgt, select.tau=c(0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_wgt, select.tau=c(0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(cdc_wgt.form05, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_wgt.form05, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_wgt.form05, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
box(lwd=1, col="deepskyblue4")


# BMI
# bmi - gcrq formula - NOTE: may need to correct min/max lambda
tau3=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.85, 0.9, 0.95)
mform_dt_bmi <-gcrq(data=dt_bmi, bmi ~ ps(age, lambda=6) + t + idmw, tau=tau3) # lambda to 6 for smoothing 
# bmi - plot percentiles
plot(mform_dt_bmi, term=1, res=TRUE, ylab='bmi', lwd=2)

# plot - KS BMI
plot.new()
par(mar=c(3,4,3,2.5) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(2,20), ylim=c(11.4,44.6), xaxs="i", yaxs="i")
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="2 to 20 years: Klinefelter Syndrome Males", adj=0, line=1.55, cex.main=1.5, col.main="deepskyblue4")
title(main='Body mass index-for-age percentiles', adj=0, line=0.4, cex.main=1.5, col.main="deepskyblue4")
title(xlab='AGE (YEARS)', line=1.9, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
title(ylab='BMI', line=2.75, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
# grid marks
axis(1, at=seq(2,20,0.5), lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=seq(11.4,44.6,0.2), lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=2:20, lwd=0.9, tck=1, cex.axis=1.5, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=12:44, lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=12:44, lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# curves (re-run the code 3x to darken the plot lines)
plot(mform_dt_bmi, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col="deepskyblue4")
plot(mform_dt_bmi, select.tau=c(0.95, 0.85, 0.5, 0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
box(lwd=1, col="deepskyblue4")
# save as bmi_kscurve

# cdc bmi - gcrq formula
cdc_bmi <- subset(cdc_bmi, select=c(source, param, age, agemo, p03, p05, p10, p25, p50, p75, p85, p90, p95, p97))
cdc_bmi.form05 <- gcrq(data=cdc_bmi, p05 ~ ps(age), tau=0.5)
cdc_bmi.form50 <- gcrq(data=cdc_bmi, p50 ~ ps(age), tau=0.5)
cdc_bmi.form85 <- gcrq(data=cdc_bmi, p85 ~ ps(age), tau=0.5)
cdc_bmi.form95 <- gcrq(data=cdc_bmi, p95 ~ ps(age), tau=0.5)

# plot - KS BMI + CDC SHADED 
plot.new()
par(mar=c(3,4,3,2.5) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(2,20), ylim=c(11.4,44.6), xaxs="i", yaxs="i")
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="2 to 20 years: Klinefelter Syndrome Males", adj=0, line=1.55, cex.main=1.5, col.main="deepskyblue4")
title(main='Body mass index-for-age percentiles', adj=0, line=0.4, cex.main=1.5, col.main="deepskyblue4")
title(xlab='AGE (YEARS)', line=1.9, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
title(ylab='BMI', line=2.75, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
# grid marks
axis(1, at=seq(2,20,0.5), lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=seq(11.4,44.6,0.2), lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=2:20, lwd=0.9, tck=1, cex.axis=1.5, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=12:44, lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=12:44, lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# bmi cdc polygon
polygon(c(cdc_bmi$age, rev(cdc_bmi$age)), c(cdc_bmi$p95, rev(cdc_bmi$p05)), add=TRUE, lty=0, col=alpha("darkorange3", 0.3))
# curves (re-run the code 3x to darken the plot lines)
plot(mform_dt_bmi, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col="deepskyblue4")
plot(mform_dt_bmi, select.tau=c(0.95, 0.85, 0.5, 0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
box(lwd=1, col="deepskyblue4")

# plot - KS BMI 95% + CDC
plot.new()
# 95%ile curves 
plot(mform_dt_bmi, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_bmi, select.tau=c(0.85, 0.5, 0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_bmi, select.tau=c(0.95), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_bmi, select.tau=c(0.95), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_bmi, select.tau=c(0.95), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(cdc_bmi.form95, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_bmi.form95, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_bmi.form95, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
box(lwd=1, col="deepskyblue4")

# plot - KS BMI 85% + CDC 
plot.new()
# 85%ile curves 
plot(mform_dt_bmi, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_bmi, select.tau=c(0.95, 0.5, 0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_bmi, select.tau=c(0.85), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_bmi, select.tau=c(0.85), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_bmi, select.tau=c(0.85), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(cdc_bmi.form85, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_bmi.form85, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_bmi.form85, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
box(lwd=1, col="deepskyblue4")

# plot - KS BMI 50% + CDC
plot.new()
# 50%ile curves 
plot(mform_dt_bmi, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_bmi, select.tau=c(0.95, 0.85, 0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_bmi, select.tau=c(0.5), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_bmi, select.tau=c(0.5), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_bmi, select.tau=c(0.5), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(cdc_bmi.form50, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_bmi.form50, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_bmi.form50, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
box(lwd=1, col="deepskyblue4")

# plot - KS BMI 5% + CDC
plot.new()
# 5%ile curves
plot(mform_dt_bmi, select.tau=c(0.9, 0.75, 0.25, 0.1), axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_bmi, select.tau=c(0.95, 0.85, 0.5), axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_bmi, select.tau=c(0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_bmi, select.tau=c(0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_bmi, select.tau=c(0.05), axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(cdc_bmi.form05, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_bmi.form05, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(cdc_bmi.form05, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
box(lwd=1, col="deepskyblue4")


# LENGTH 
# length - gcrq formula
tau2=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
mform_dt_len <-gcrq(data=dt_len, measurement ~ ps(agemo, monotone=1) + idmw + t, tau=tau2) # MONOTONE NOT WORKING???
mform_dt_len <-gcrq(data=dt_len, measurement ~ ps(agemo) + idmw + t, tau=tau2) 
# length - plot percentiles 
plot(mform_dt_len, term=1, res=TRUE, ylab='length', lwd=2)

# plot - KS LENGTH
plot.new()
par(mar=c(3,4,3,2.5) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(0,24), ylim=c(37,98), xaxs="i", yaxs="i")
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="Birth to 24 months: Klinefelter Syndrome Males", adj=0, line=1.55, cex.main=1.5, col.main="deepskyblue4")
title(main='Length-for-age percentiles', adj=0, line=0.4, cex.main=1.5, col.main="deepskyblue4")
title(xlab='AGE (MONTHS)', line=1.9, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
title(ylab='LENGTH (CM)', line=2.75, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
# grid marks
axis(1, at=0:24, lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=37:98, lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=0, labels='Birth', lwd=0.9, tck=1, cex.axis=1.5, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(1, at=seq(3,24,3), lwd=0.9, tck=1, cex.axis=1.5, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=seq(40,95,5), lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=seq(40,95,5), lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# curves - (re-run the code 3x to darken the plot lines)
plot(mform_dt_len, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col="deepskyblue4")
plot(mform_dt_len, select.tau=c(0.95, 0.5, 0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
box(lwd=1, col="deepskyblue4")
# save as len_kscurve

# who length - gcrq formula
who_len <- subset(who_len, select=c(source, param, age, agemo, p03, p05, p10, p25, p50, p75, p85, p90, p95, p97))
who_len.form05 <- gcrq(data=who_len, p05 ~ ps(agemo), tau=0.5)
who_len.form50 <- gcrq(data=who_len, p50 ~ ps(agemo), tau=0.5)
who_len.form95 <- gcrq(data=who_len, p95 ~ ps(agemo), tau=0.5)

# plot - KS LENGTH + WHO SHADED  
plot.new()
par(mar=c(3,4,3,2.5) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(0,24), ylim=c(37,98), xaxs="i", yaxs="i")
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="Birth to 24 months: Klinefelter Syndrome Males", adj=0, line=1.55, cex.main=1.5, col.main="deepskyblue4")
title(main='Length-for-age percentiles', adj=0, line=0.4, cex.main=1.5, col.main="deepskyblue4")
title(xlab='AGE (MONTHS)', line=1.9, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
title(ylab='LENGTH (CM)', line=2.75, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
# grid marks
axis(1, at=0:24, lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=37:98, lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=0, labels='Birth', lwd=0.9, tck=1, cex.axis=1.5, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(1, at=seq(3,24,3), lwd=0.9, tck=1, cex.axis=1.5, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=seq(40,95,5), lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=seq(40,95,5), lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# length who polygon
polygon(c(who_len$agemo, rev(who_len$agemo)), c(who_len$p95, rev(who_len$p05)), add=TRUE, lty=0, col=alpha("darkorange3", 0.3))
# curves (re-run the code 3x to darken the plot lines)
plot(mform_dt_len, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col="deepskyblue4")
plot(mform_dt_len, select.tau=c(0.95, 0.5, 0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
box(lwd=1, col="deepskyblue4")

# plot - KS LENGTH 95% + WHO
plot.new()
# 95%ile curves 
plot(mform_dt_len, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_len, select.tau=c(0.5, 0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_len, select.tau=c(0.95), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_len, select.tau=c(0.95), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_len, select.tau=c(0.95), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(who_len.form95, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(who_len.form95, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(who_len.form95, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
box(lwd=1, col="deepskyblue4")

# plot - KS LENGTH 50% + WHO
plot.new()
# 50%ile curves 
plot(mform_dt_len, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_len, select.tau=c(0.95, 0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_len, select.tau=c(0.5), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_len, select.tau=c(0.5), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_len, select.tau=c(0.5), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(who_len.form50, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(who_len.form50, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(who_len.form50, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
box(lwd=1, col="deepskyblue4")

# plot - KS LENGTH 5% + WHO
plot.new()
# 5%ile curves 
plot(mform_dt_len, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_len, select.tau=c(0.95, 0.5), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_len, select.tau=c(0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_len, select.tau=c(0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_len, select.tau=c(0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(who_len.form05, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(who_len.form05, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
plot(who_len.form05, axes=FALSE, add=TRUE, lwd=4, lty=1, col="darkorange3")
box(lwd=1, col="deepskyblue4")


# BIRTHWEIGHT
# birthweight - gcrq formula
tau2=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
mform_dt_bw <-gcrq(data=dt_bw, measurement ~ ps(agemo, lambda=6) + idmw + t, tau=tau2) 
# birthweight - plot percentiles 
plot(mform_dt_bw, term=1, res=TRUE, ylab='weight', lwd=2)

# plot - KS BIRTHWEIGHT
plot.new()
par(mar=c(3,4,3,2.5) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(0,24), ylim=c(0.4,16.6), xaxs="i", yaxs="i")
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="Birth to 24 months: Klinefelter Syndrome Males", adj=0, line=1.55, cex.main=1.5, col.main="deepskyblue4")
title(main='Weight-for-age percentiles', adj=0, line=0.4, cex.main=1.5, col.main="deepskyblue4")
title(xlab='AGE (MONTHS)', line=1.9, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
title(ylab='WEIGHT (KG)', line=2.75, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
# grid marks
axis(1, at=0:24, lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=seq(0.4,16.6,0.2), lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=0, labels='Birth', lwd=0.9, tck=1, cex.axis=1.5, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(1, at=seq(3,24,3), lwd=0.9, tck=1, cex.axis=1.5, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=1:16, lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=1:16, lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# curves - (re-run the code 3x to darken the plot lines)
plot(mform_dt_bw, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col="deepskyblue4")
plot(mform_dt_bw, select.tau=c(0.95, 0.5, 0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
box(lwd=1, col="deepskyblue4")

# who birthweight - gcrq formulas
who_wgt <- subset(who_wgt, select=c(source, param, age, agemo, p03, p05, p10, p25, p50, p75, p85, p90, p95, p97))
who_wgt.form05 <- gcrq(data=who_wgt, p05 ~ ps(agemo), tau=0.5)
who_wgt.form50 <- gcrq(data=who_wgt, p50 ~ ps(agemo), tau=0.5)
who_wgt.form95 <- gcrq(data=who_wgt, p95 ~ ps(agemo), tau=0.5)

# plot - KS BIRTHWEIGHT + WHO SHADED  
plot.new()
par(mar=c(3,4,3,2.5) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(0,24), ylim=c(0.4,16.6), xaxs="i", yaxs="i")
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="Birth to 24 months: Klinefelter Syndrome Males", adj=0, line=1.55, cex.main=1.5, col.main="deepskyblue4")
title(main='Weight-for-age percentiles', adj=0, line=0.4, cex.main=1.5, col.main="deepskyblue4")
title(xlab='AGE (MONTHS)', line=1.9, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
title(ylab='WEIGHT (KG)', line=2.75, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
# grid marks
axis(1, at=0:24, lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=seq(0.4,16.6,0.2), lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=0, labels='Birth', lwd=0.9, tck=1, cex.axis=1.5, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(1, at=seq(3,24,3), lwd=0.9, tck=1, cex.axis=1.5, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=1:16, lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=1:16, lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# birthweight who polygon
polygon(c(who_wgt$agemo, rev(who_wgt$agemo)), c(who_wgt$p95, rev(who_wgt$p05)), add=TRUE, lty=0, col=alpha("darkorange3", 0.3))
# curves (re-run the code 3x to darken the plot lines)
plot(mform_dt_bw, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col="deepskyblue4")
plot(mform_dt_bw, select.tau=c(0.95, 0.5, 0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
box(lwd=1, col="deepskyblue4")

# plot - KS BIRTHWEIGHT 95% + WHO
plot.new()
# 95%ile curves 
plot(mform_dt_bw, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_bw, select.tau=c(0.5, 0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_bw, select.tau=c(0.95), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_bw, select.tau=c(0.95), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_bw, select.tau=c(0.95), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(who_wgt.form95, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
plot(who_wgt.form95, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
plot(who_wgt.form95, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
box(lwd=1, col="deepskyblue4")

# plot - KS BIRTHWEIGHT 50% + WHO
plot.new()
# 50%ile curves 
plot(mform_dt_bw, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_bw, select.tau=c(0.95, 0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_bw, select.tau=c(0.5), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_bw, select.tau=c(0.5), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_bw, select.tau=c(0.5), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(who_wgt.form50, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
plot(who_wgt.form50, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
plot(who_wgt.form50, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
box(lwd=1, col="deepskyblue4")

# plot - KS BIRTHWEIGHT 5% + WHO
plot.new()
# 5%ile curves 
plot(mform_dt_bw, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_bw, select.tau=c(0.95, 0.5), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_bw, select.tau=c(0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_bw, select.tau=c(0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_bw, select.tau=c(0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(who_wgt.form05, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
plot(who_wgt.form05, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
plot(who_wgt.form05, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
box(lwd=1, col="deepskyblue4")


# WEIGHT-FOR-LENGTH
# weight for length - gcrq formula
tau2=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
mform_dt_wl <-gcrq(data=dt_wl, WEIGHTKG ~ ps(HEIGHTCM) + idmw + t, tau=tau2)
# weight for length - plot percentiles 
plot(mform_dt_wl, term=1, res=TRUE, ylab='weight', xlab='length', lwd=2)

# plot - KS WEIGHT-FOR-LENGTH
plot.new()
par(mar=c(3,4,3,2.5) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(45,90), ylim=c(0.4,18.6), xaxs="i", yaxs="i") # KS range 42.6-93.9; WHO range 45-*** 
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="Birth to 24 months: Klinefelter Syndrome Males", adj=0, line=1.55, cex.main=1.5, col.main="deepskyblue4")
title(main='Weight-for-length percentiles', adj=0, line=0.4, cex.main=1.5, col.main="deepskyblue4")
title(xlab='LENGTH (CM)', line=1.9, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
title(ylab='WEIGHT (KG)', line=2.75, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
# grid marks - NOTE: need to limit x-axis to 40 otherwise percentiles fall off grid (and consider extending age range to extend length greater)
axis(1, at=45:90, lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=seq(0.4,18.6,0.2), lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=seq(45,90,5), lwd=0.9, tck=1, cex.axis=1.5, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=1:18, lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=1:18, lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# curves (re-run the code 3x to darken the plot lines)
plot(mform_dt_wl, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col="deepskyblue4")
plot(mform_dt_wl, select.tau=c(0.95, 0.5, 0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
box(lwd=1, col="deepskyblue4")

# weight for length who - gcrq formulas 
who_wl <- subset(who_wl, select=c(source, param, length, p03, p05, p10, p25, p50, p75, p85, p90, p95, p97))
who_wl.form05 <- gcrq(data=who_wl, p05 ~ ps(length), tau=0.5)
who_wl.form50 <- gcrq(data=who_wl, p50 ~ ps(length), tau=0.5)
who_wl.form95 <- gcrq(data=who_wl, p95 ~ ps(length), tau=0.5)

# plot - KS WEIGHT-FOR-LENGTH + WHO SHADED
plot.new()
par(mar=c(3,4,3,2.5) + 0.1) #for use with a title 
#par(mar=c(3,4,1,2.5) + 0.1) #for no title
plot.window(xlim=c(45,90), ylim=c(0.4,18.6), xaxs="i", yaxs="i") # KS range 42.6-93.9; WHO range 45-*** 
box(lwd=1, col="deepskyblue4")
# axis titles
title(main="Birth to 24 months: Klinefelter Syndrome Males", adj=0, line=1.55, cex.main=1.5, col.main="deepskyblue4")
title(main='Weight-for-length percentiles', adj=0, line=0.4, cex.main=1.5, col.main="deepskyblue4")
title(xlab='LENGTH (CM)', line=1.9, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
title(ylab='WEIGHT (KG)', line=2.75, cex.lab=1.6, font.lab=2, col.lab="deepskyblue4")
# grid marks - NOTE: need to limit x-axis to 40 otherwise percentiles fall off grid (and consider extending age range to extend length greater)
axis(1, at=45:90, lwd=0.3, tck=1, cex.axis=0.05, col.axis="white", col.ticks="deepskyblue4")
axis(2, at=seq(0.4,18.6,0.2), lwd=0.3, tck=1, cex.axis=0.01, las=1, col.axis="white",col.ticks="deepskyblue4")
# axis labels
axis(1, at=seq(45,90,5), lwd=0.9, tck=1, cex.axis=1.5, padj=-0.6, gap.axis = 0.5, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(2, at=1:18, lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.65, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")
axis(4, at=1:18, lwd=0.9, tck=1, cex.axis=1.5, las=1, hadj=0.35, col.axis="deepskyblue4", col.ticks="deepskyblue4", fg="deepskyblue4")

# weight-for-length who polygon
polygon(c(who_wl$length, rev(who_wl$length)), c(who_wl$p95, rev(who_wl$p05)), add=TRUE, lty=0, col=alpha("darkorange3", 0.3))
# curves (re-run the code 3x to darken the plot lines)
plot(mform_dt_wl, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col="deepskyblue4")
plot(mform_dt_wl, select.tau=c(0.95, 0.5, 0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
box(lwd=1, col="deepskyblue4")

# plot - KS WEIGHT-FOR-LENGTH 95% + WHO
plot.new()
# 95%ile curves 
plot(mform_dt_wl, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_wl, select.tau=c(0.5, 0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_wl, select.tau=c(0.95), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_wl, select.tau=c(0.95), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_wl, select.tau=c(0.95), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(who_wl.form95, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
plot(who_wl.form95, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
plot(who_wl.form95, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
box(lwd=1, col="deepskyblue4")

# plot - KS WEIGHT-FOR-LENGTH 50% + WHO
plot.new()
# 50%ile curves 
plot(mform_dt_wl, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_wl, select.tau=c(0.95, 0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_wl, select.tau=c(0.5), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_wl, select.tau=c(0.5), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_wl, select.tau=c(0.5), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(who_wl.form50, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
plot(who_wl.form50, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
plot(who_wl.form50, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
box(lwd=1, col="deepskyblue4")

# plot - KS WEIGHT-FOR-LENGTH 5% + WHO
plot.new()
# 5%ile curves 
plot(mform_dt_wl, select.tau=c(0.9, 0.75, 0.25, 0.1), term=1, axes=FALSE, add=TRUE, lwd=2.5, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_wl, select.tau=c(0.95, 0.5), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col=alpha("deepskyblue4", 0.75))
plot(mform_dt_wl, select.tau=c(0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_wl, select.tau=c(0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(mform_dt_wl, select.tau=c(0.05), term=1, axes=FALSE, add=TRUE, lwd=4, lty=1, col="deepskyblue4")
plot(who_wl.form05, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
plot(who_wl.form05, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
plot(who_wl.form05, axes=FALSE,add=TRUE,lwd=4, col='darkorange3')
box(lwd=1, col="deepskyblue4")

