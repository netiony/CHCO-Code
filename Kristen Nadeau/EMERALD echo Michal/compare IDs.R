library(stringr)

data_michal <- read.csv("E:\\Nadeau\\Nadeau EMERALD\\Echo\\Michal\\Raw data\\EMERALG_Echo_Electro_final.csv")

# correct randomizations
data_michal$randomization_m <- data_michal$randomization
data_michal$randomization_m <- ifelse(data_michal$subject_id %in% c("6377-02","6377-24"),1,data_michal$randomization_m)

# code randomization group
data_michal$randgroup_m <- ifelse(data_michal$randomization_m==1,"Metformin","Placebo")

data_michal <- data_michal[,c("subject_id","randgroup_m","randomization_m")]
data_michal$subject_id <- str_trim(data_michal$subject_id)

# read in Alex data
data_alex <- read.csv("E:\\Nadeau\\Nadeau EMERALD\\Echo\\Michal\\Clean data\\comparing_ids.csv")
data_alex$subject_id <- data_alex$Study.ID
data_alex$Study.ID <- NULL
data_alex$randomization_a <- data_alex$randomization
data_alex$randomization <- NULL
data_alex$subject_id <- str_trim(data_alex$subject_id)

# merge
check <- merge(data_michal,data_alex,by="subject_id",all.x = T, all.y = T)
check <- unique(check)
