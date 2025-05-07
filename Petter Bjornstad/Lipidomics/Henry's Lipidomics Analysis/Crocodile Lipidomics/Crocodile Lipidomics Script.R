##Install lipidr package before if necessary
##To use lipidr for your analysis using numerical matrix as input, you need 2 files:
 #Numerical table where lipids are rows and samples are columns. Lipid names 
  #should be in the first column, and sample names are in the first row. 
 #A table with the sample annotation / groups, where the sample names are in 
  #first column. Note the sample names must be identical in the two files. 

# Load Packages -----------------------------------------------------------
BiocManager::install("lipidr")  # If from Bioconductor
  library(lipidr)
install.packages("BiocManager")
BiocManager::install("BiocParallel")
install.packages("matrixStats")
library(readr)
install.packages("dplyr")
  library(dplyr)


# Set Working Directory ---------------------------------------------------
setwd("C:/Users/Henry Mangalapalli/OneDrive - UW/Documents/M3D/Bjornstad Lab Rotation/Crocodile Lipidomics")
getwd()
list.files()  # Lists all files in the current directory
              # Settings gear -> "Go to WD" if you can't see files!

# Importing ---------------------------------------------------------------
##Import table with data, where lipids are rows & samples are columns
croc_lipid_name <- read_csv("crocodile_lipidomics_name.csv")
View(croc_lipid_name)
##Import table with sample annotations
croc_lipid <- read_csv("crocodile_lipidomics.csv")
View(croc_lipid)

# Re-labeling - OPTIONAL -------------------------------------------------------------
##Use file name in working directory, NOT in 'Importing' tab above!(I didn't do
  #this here because I changed the name twice LOL)
file.rename("crocodile_lipidomics_name", "croc_lipid_name")
file.rename("crocodile_lipidomics", "croc_lipid")
##

# View & Explore ----------------------------------------------------------
##Check and Explore the Data
 #View Summary
summary(croc_lipid)
 

# Normalizing -------------------------------------------------------------
##Normalize the Data
##Something something normalizing by ignoring NA values
##Ask Hailey for the code if I can't figure it out
# Define the Min-Max normalization function
normalize_min_max <- function(x) {(x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))}

# Apply Min-Max normalization to all numeric columns of the dataset
###THIS IS THE ONE!!!###
croc_lipid_normalized <- croc_lipid %>%
  mutate(across(where(is.numeric), normalize_min_max))

  # Apply to a column or entire dataset (Maybe try another time?)
croc_lipid <- normalize_min_max(croc_lipid)

#ChatGPT says this step is necessary but this keeps killing my R sessions
install.packages("devtools")  # if not already installed
devtools::install_github("ahmohamed/lipidr")
library(lipidr)
##Clone code from GitHub

# Plot PCA
plot_pca(croc_lipid_normalized, color_by = "group")


