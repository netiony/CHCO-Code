library(readxl)
library(blockrand)

set.seed(1234586444)

# goal for the entire study is 30 participants per group
# adding a new stratum of BMI >= 35 - 45
# there are two other strata (<25, 25 - <35), this randomization sequence was done previously
# will generate a sequence of 80 assignments just in case

rand <- blockrand(n=80, levels = c(1:2), block.sizes = 2)

# need to update below - going to email Sam and ask him if he cares if 1 = placebo and 2 = sema 
# i.e., does he want it to match what he is already using
rand$treatment_char <- ifelse(rand$treatment == "1", "Placebo", "Semaglutide")

# output
write.csv(rand, 
          "/Volumes/pylell/T1-DISCO/randomization_sequence_BMI_35-45.csv",
          row.names = F)



