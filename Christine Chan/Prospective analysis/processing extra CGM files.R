library(cgmanalysis)

# need to get the new csv files in the right format - ID, date/time, glucose
# read in the data and create date time variable


cleandata(inputdirectory = "E:/Christine Chan/Prospective analysis/Data raw/Corrected outcome data/Extra CGM raw",
          outputdirectory = "E:/Christine Chan/Prospective analysis/Data raw/Corrected outcome data/Extra CGM processed")

cgmvariables(inputdirectory = "E:/Christine Chan/Prospective analysis/Data raw/Corrected outcome data/Extra CGM processed",
             outputdirectory = "E:/Christine Chan/Prospective analysis/Data raw/Corrected outcome data/Extra CGM upload")