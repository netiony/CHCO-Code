data <- read.csv("E:\\Lauren Shomaker\\MIND\\Clark F31\\EC MIND Baseline Condensed R 12.3.20.csv")

mod <- lm(wbisi ~ actiwatch_ave_sleep_eff*teen_maas_score, data=data)
summary(mod)

mod <- lm(wbisi ~ actiwatch_ave_sleep_eff*teen_self_compassion_score, data=data)
summary(mod)

mod <- lm(wbisi ~ actiwatch_ave_time_asleep*teen_maas_score, data=data)
summary(mod)

mod <- lm(wbisi ~ actiwatch_ave_time_asleep*teen_self_compassion_score, data=data)
summary(mod)

mod <- lm(wbisi ~ actiwatch_ave_timebed*teen_maas_score, data=data)
summary(mod)

mod <- lm(wbisi ~ actiwatch_ave_timebed*teen_self_compassion_score, data=data)
summary(mod)

mod <- lm(wbisi ~ actiwatch_ave_onsetlatency*teen_maas_score, data=data)
summary(mod)

mod <- lm(wbisi ~ actiwatch_ave_onsetlatency*teen_self_compassion_score, data=data)
summary(mod)

mod <- lm(wbisi ~ actiwatch_ave_waso*teen_maas_score, data=data)
summary(mod)

mod <- lm(wbisi ~ actiwatch_ave_waso*teen_self_compassion_score, data=data)
summary(mod)