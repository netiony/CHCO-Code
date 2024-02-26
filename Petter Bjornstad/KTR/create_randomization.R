# Goal is to recruit a total of 40 per group (dapagliflozin and placebo)
# Stratified by sex and moderate or severe albuminuria
# To account for possible extra randomizations, will create a sequence of 80 randomizations for each of the 4 strata
# with permuted random blocks, don't get exactly the same N per stratum but each stratum is balanced by treatment

library(blockrand)

set.seed(3654)
male_moderate <- blockrand(n=80, id.prefix='MM', block.prefix='MM',stratum='Male_moderate_albuminuria', block.sizes = 1:4)
male_severe <- blockrand(n=80, id.prefix='MS', block.prefix='MS',stratum='Male_severe_albuminuria', block.sizes = 1:4)
female_moderate <- blockrand(n=80, id.prefix='FM', block.prefix='FM',stratum='Female_moderate_albuminuria', block.sizes = 1:4)
female_severe <- blockrand(n=80, id.prefix='FS', block.prefix='FS',stratum='Female_severe_albuminuria', block.sizes = 1:4)

rand <- rbind(male_severe, male_moderate, female_severe, female_moderate)
rand$treatment_char <- ifelse(rand$treatment == "A", "Placebo", "Dapagliflozin")

write.csv(rand, "/Volumes/pylell/Kendrick Bjornstad KTR trial/KTR trial 23-1360 randomization.csv", row.names = F)