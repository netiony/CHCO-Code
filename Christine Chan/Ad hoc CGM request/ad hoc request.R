library(cgmanalysis)


cleandata(inputdirectory = "H:\\Endocrinology\\Chan\\Ad hoc CGM request\\Files for processing Raw data for Laura",
          outputdirectory = "H:\\Endocrinology\\Chan\\Ad hoc CGM request\\Cleaned files",id_filename = T)

cgmvariables(inputdirectory = "H:\\Endocrinology\\Chan\\Ad hoc CGM request\\Cleaned files",
             outputdirectory = "H:\\Endocrinology\\Chan\\Ad hoc CGM request\\Output",
             customintervals = list(c(70,140)))


# read in original dataset
#orig <- read.csv("H:\\Endocrinology\\Chan\\Chan CF CGM\\Repeat primary analysis\\Data\\export_tim.csv")
#group <- orig[,c("")]


