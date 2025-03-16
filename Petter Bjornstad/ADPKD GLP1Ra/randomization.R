library(readxl)
library(blockrand)
library(dplyr)

set.seed(1234586444)

# goal is to randomize 63 subjects per group
# will generate a sequence of 80 subjects per strata just in case
# strata are sex, age (<45 and >=45 yrs), and baseline BMI (>=30 and <30)
# use random permuted block sizes

# males, age <45, BMI <30
rand_males_lt45_bmilt30 <- blockrand(n=80, levels = c(1:2), block.sizes = 2)
rand_males_lt45_bmilt30$treatment_char <- ifelse(rand_males_lt45_bmilt30$treatment == "1", "Placebo", "Tirzepatide")
rand_males_lt45_bmilt30$strata <- "Male, <45 years, BMI <30"

# males, age >=45, BMI <30
rand_males_ge45_bmilt30 <- blockrand(n=80, levels = c(1:2), block.sizes = 2)
rand_males_ge45_bmilt30$treatment_char <- ifelse(rand_males_ge45_bmilt30$treatment == "1", "Placebo", "Tirzepatide")
rand_males_ge45_bmilt30$strata <- "Male, >=45 years, BMI <30"

# males, age <45, BMI >=30
rand_males_lt45_bmige30 <- blockrand(n=80, levels = c(1:2), block.sizes = 2)
rand_males_lt45_bmige30$treatment_char <- ifelse(rand_males_lt45_bmige30$treatment == "1", "Placebo", "Tirzepatide")
rand_males_lt45_bmige30$strata <- "Male, <45 years, BMI >=30"

# males, age >=45, BMI >=30
rand_males_ge45_bmige30 <- blockrand(n=80, levels = c(1:2), block.sizes = 2)
rand_males_ge45_bmige30$treatment_char <- ifelse(rand_males_ge45_bmige30$treatment == "1", "Placebo", "Tirzepatide")
rand_males_ge45_bmige30$strata <- "Male, >=45 years, BMI >=30"

# females, age <45, BMI <30
rand_females_lt45_bmilt30 <- blockrand(n=80, levels = c(1:2), block.sizes = 2)
rand_females_lt45_bmilt30$treatment_char <- ifelse(rand_females_lt45_bmilt30$treatment == "1", "Placebo", "Tirzepatide")
rand_females_lt45_bmilt30$strata <- "Female, <45 years, BMI <30"

# females, age >=45, BMI <30
rand_females_ge45_bmilt30 <- blockrand(n=80, levels = c(1:2), block.sizes = 2)
rand_females_ge45_bmilt30$treatment_char <- ifelse(rand_females_ge45_bmilt30$treatment == "1", "Placebo", "Tirzepatide")
rand_females_ge45_bmilt30$strata <- "Female, >=45 years, BMI <30"

# females, age <45, BMI >=30
rand_females_lt45_bmige30 <- blockrand(n=80, levels = c(1:2), block.sizes = 2)
rand_females_lt45_bmige30$treatment_char <- ifelse(rand_females_lt45_bmige30$treatment == "1", "Placebo", "Tirzepatide")
rand_females_lt45_bmige30$strata <- "Female, <45 years, BMI >=30"

# females, age >=45, BMI >=30
rand_females_ge45_bmige30 <- blockrand(n=80, levels = c(1:2), block.sizes = 2)
rand_females_ge45_bmige30$treatment_char <- ifelse(rand_females_ge45_bmige30$treatment == "1", "Placebo", "Tirzepatide")
rand_females_ge45_bmige30$strata <- "Female, >=45 years, BMI >=30"

rand <- rbind(rand_males_lt45_bmilt30, rand_males_ge45_bmilt30, rand_males_lt45_bmige30, rand_males_ge45_bmige30,
              rand_females_lt45_bmilt30, rand_females_ge45_bmilt30, rand_females_lt45_bmige30, rand_females_ge45_bmige30)
rand <- rand %>% select(strata, treatment_char)
colnames(rand) <- c("Strata", "Arm (char)")

# output
write.csv(rand, 
          "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Nowak/ADPKD GLP1-Ra/Randomization/rand_sequence_ADPKD_GLP1Ra.csv",
          row.names = F)



