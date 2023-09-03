data <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Melanie Green/TEAL/TEAL_randomization.csv")
table(data$random_group,data$bmi_35)

# subset the data to the number actually randomized
# I think 34
data_rand <- data[1:34,]
table(data_rand$random_group,data_rand$bmi_35)

# read in the actual randomization data from REDCap
actual_rand <- read.csv("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Melanie Green/TEAL/TEALRandomization_DATA_2023-09-02_2157.csv")
# remove first 17 who were assiged to drug
actual_rand_post17 <- actual_rand[18:nrow(actual_rand),]
table(actual_rand_post17$random_group,actual_rand_post17$bmi_35)
