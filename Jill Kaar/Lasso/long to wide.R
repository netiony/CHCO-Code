# read in long data
long <- read.csv("H:\\Endocrinology\\Kaar\\LASSO\\Raw data\\long.csv")

# separate out demographics data
demovars <- c("")

# testing long to wide on PA
PA <- long[long$redcap_repeat_instrument=="PA",]
