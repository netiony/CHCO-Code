# Goal is to recruit a total of 40 per group (dapagliflozin and placebo)
# Stratified by sex and moderate or severe albuminuria
# To account for possible extra randomizations, will create a sequence of 80 randomizations for each of the 4 strata
# with permuted random blocks, don't get exactly the same N per stratum but each stratum is balanced by treatment

library(blockrand)

set.seed(3654)
male_moderate_biopsy <- blockrand(n=80, id.prefix='MMB', block.prefix='MMB',stratum='Male_moderate_albuminuria_biopsy', block.sizes = 1:4)
male_moderate_no_biopsy <- blockrand(n=80, id.prefix='MMN', block.prefix='MMN',stratum='Male_moderate_albuminuria_no_biopsy', block.sizes = 1:4)

male_severe_biopsy <- blockrand(n=80, id.prefix='MSB', block.prefix='MSB',stratum='Male_severe_albuminuria_biopsy', block.sizes = 1:4)
male_severe_no_biopsy <- blockrand(n=80, id.prefix='MSN', block.prefix='MSN',stratum='Male_severe_albuminuria_no_biopsy', block.sizes = 1:4)

female_moderate_biopsy <- blockrand(n=80, id.prefix='FMB', block.prefix='FMB',stratum='Female_moderate_albuminuria_biopsy', block.sizes = 1:4)
female_moderate_no_biopsy <- blockrand(n=80, id.prefix='FMN', block.prefix='FMN',stratum='Female_moderate_albuminuria_no_biopsy', block.sizes = 1:4)

female_severe_biopsy <- blockrand(n=80, id.prefix='FSB', block.prefix='FSB',stratum='Female_severe_albuminuria_biopsy', block.sizes = 1:4)
female_severe_no_biopsy <- blockrand(n=80, id.prefix='FSN', block.prefix='FSN',stratum='Female_severe_albuminuria_no_biopsy', block.sizes = 1:4)

rand <- rbind(male_severe_biopsy, male_severe_no_biopsy, male_moderate_biopsy, male_moderate_no_biopsy, 
              female_severe_biopsy, female_severe_no_biopsy, female_moderate_biopsy, female_moderate_no_biopsy)
rand$treatment_char <- ifelse(rand$treatment == "A", "Placebo", "Dapagliflozin")

write.csv(rand, "/Volumes/pylell/Kendrick Bjornstad KTR trial/KTR trial 23-1360 randomization with biopsy.csv", row.names = F)