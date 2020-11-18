source('./cleandata.r')

inputdirectory <- "H:\\Endocrinology\\Chan\\Functional data analysis\\Raw CGM Data"
outputdirectory <- "H:\\Endocrinology\\Chan\\Functional data analysis\\Cleaned CGM Data"

cleandata(inputdirectory,outputdirectory,removegaps = TRUE,
          gapfill = TRUE,maximumgap = 20,id_filename = T,verbose = T)



