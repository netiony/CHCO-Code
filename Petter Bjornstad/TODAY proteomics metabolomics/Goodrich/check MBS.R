library(dplyr)
library(stringr)
library(readxl)

g_ids <- read_xlsx('/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/Goodrich/Data requests/unique_release_id_in_TODAY.xlsx')

if(Sys.info()["sysname"] == "Windows"){
  home_dir = "E:/Petter Bjornstad/TODAY subaward"
} else if (Sys.info()["sysname"] == "Linux"){
  home_dir = "~/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward"
} else if (Sys.info()["sysname"] == "Darwin"){
  home_dir = "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/TODAY subaward"
}

setwd(home_dir)

# read in TME dataset to find those who had bariatric surgery
tme <- read.csv("./Clinical data/TODAY2/TME.csv")
tme <- tme %>% filter(TMETYPE==1)
tme <- tme %>% select(RELEASEID,TMETYPE,DAYSTOTME)
colnames(tme) <- c("releaseid","TMETYPE","DAYSTOTME")

g_ids <- left_join(g_ids, tme, by = "releaseid")

write.csv(g_ids, '/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/Goodrich/Data requests/MBS_PFAS.csv', row.names = F)
