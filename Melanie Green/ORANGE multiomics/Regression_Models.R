rm(list = ls())

#to derive combinations
library(gtools)
library(tidyverse)
#July 2020 for Ms
dataIn = "/Volumes/RI Biostatistics Core/Shared/Shared Projects/Bothwell/Peds ENDO/Siebert - Multiomics/Data/merged multiomics data.csv"
dataOut = "/Volumes/RI Biostatistics Core/Shared/Shared Projects/Bothwell/Peds ENDO/Siebert - Multiomics/Data/omes_models_9_22_2023.csv"
dataOut = "~/Downloads/omes_models_9_29_2023.csv"


##########Begin processing
dd = read.csv(dataIn)

# Filter to treatment arm 0
dd <- dd[dd$Treatment.Arm.PLA...0..EAA.1 == 0, ]

# Data cleaning - remove rows missing > 80% data
dd <- dd[which(rowMeans(!is.na(dd)) > 0.8),]

# require liver fat data 
dd <- dd[!is.na(dd$liver.fat..),]


# Set lower threshold of GLP as 0.5
dd$GLP.1.T.1000...180..Start.of.OSTT <- ifelse(dd$GLP.1.T.1000...180..Start.of.OSTT == "<1", "0.5", dd$GLP.1.T.1000...180..Start.of.OSTT)
dd$GLP.1.T.1000...180..Start.of.OSTT <- as.numeric(dd$GLP.1.T.1000...180..Start.of.OSTT)

# Recode NA for REM
dd$REM.AHI <- ifelse(dd$REM.AHI == "N/A", NA, dd$REM.AHI)
dd$REM.AHI <- as.numeric(dd$REM.AHI)

# remove height, weight, and Fat Mass, and Lean Mass (we have fat and lean percent so we don't need both)
dd <- dd %>% dplyr::select(-c(Fat_Mass, Lean_Mass, `OGTT.Height..cm.`, `OGTT.Wt..kg.`, `OGTT.BMI.Zscore`))
# rename(Matsuda_Index = Matsuda..10.000.sqrt.FPG.FPI.mean.G.mean.I.)

##### NEED to calculate phyla relative abundance
### Calculate F:B ratio
# Calculate phyla relative abundance
smb_sum = rowSums(dd[,669:684])
dd[,669:684] = dd[,669:684]/smb_sum
# f:b ratio
dd$f_b_ratio <- dd$stool_microbiome_phyla_Bacteria.Firmicutes/dd$stool_microbiome_phyla_Bacteria.Bacteroidetes

# remove columns we don't need (extra levels of stool microbiome)
dd <- dd[,-c(307:398,619:684)]

# remove .y suffixes 
dd <- dd[, -grep("ine.y", colnames(dd))]


## convert stool_microbiome_genus_ variables into relative abundance 
smb_sum = rowSums(dd[,296:515])
dd[,296:515] = dd[,296:515]/smb_sum


# Fix prefix of names 
names(dd) <- gsub(x = names(dd), pattern = "plasma_carnitine_", replacement = "pa_")
names(dd) <- gsub(x = names(dd), pattern = "plasma_aa_L.", replacement = "pa_")
names(dd) <- gsub(x = names(dd), pattern = "plasma_lipids_abs_", replacement = "pa_")
names(dd) <- gsub(x = names(dd), pattern = "plasma_lipids_global_", replacement = "pg_")
names(dd) <- gsub(x = names(dd), pattern = "stool_lipids_norm_", replacement = "sg_")
names(dd) <- gsub(x = names(dd), pattern = "stool_microbiome_alpha_", replacement = "smb_")
names(dd) <- gsub(x = names(dd), pattern = "stool_microbiome_genus_", replacement = "smb_")

## change prefix of vars 
colnames(dd)[c(4:9, 12:36, 38:50, 52:54)] <- paste0("co_", colnames(dd)[c(4:9, 12:36, 38:50, 52:54)])
colnames(dd)[c(10, 11, 37, 51, 55, 56, 516)] <- paste0("po_", colnames(dd)[c(10, 11, 37, 51, 55, 56, 516)])


## Fix microbiome suffixes
colnames(dd)[c(298:324, 326:350, 352:464, 466:514)] <- paste0("smb_", 
      sapply(strsplit(colnames(dd)[c(298:324, 326:350, 352:464, 466:514)], split = "[.]"), tail, 2)[1,],
      "_",
      sapply(strsplit(colnames(dd)[c(298:324, 326:350, 352:464, 466:514)], split = "[.]"), tail, 2)[2,])


names(dd)[51] <- "po_Matsuda.Index"


# Remove columns with > 80% missing data or more than 33.3% of values are 0
dd <- dd[, which(colMeans(!is.na(dd)) > 0.8)]
dd <- dd[, which(colSums(dd != 0) > length(dd$X)/3)] # if a column has more than 1/3 0, then we are removing it



write.csv(dd, "/Volumes/RI Biostatistics Core/Shared/Shared Projects/Bothwell/Peds ENDO/Siebert - Multiomics/Data/merged multiomics data clean.csv")




myOmes = c("po_", "co_", "pa_", "pg_", "sg_", "smb_")

# list assay pairs
cc = combinations(n=length(myOmes),r=2)
omeCols = lapply(myOmes, function(x) grep(paste("^", x, sep=""), colnames(dd)))
#Inspect list of omes
cbind(myOmes[cc[,1]], myOmes[cc[,2]])
#This is the "within" ome comparisons, e.g. maMb at one time point compared to another


ccShort = cc

#create function regress
regress = function(oA, oB, assays){
  oA = omeCols[[oA]]
  oB = omeCols[[oB]]
  
  #TODO:  Check for adequate numbers
  
  for (h in oA)
  {	
    #print(h)
    #h = 38
    ana1 = colnames(dd)[h]
    print(ana1)
    
    for (i in oB)
    {    
      #i=60
      print(paste(h,i))
      ana2 = colnames(dd)[i]
      #print(ana2)
      #subset for complete cases.
      myCols = c(h,i)
      #dcc = dd[complete.cases(dd[, myCols]),myCols]
      dcc = dd[complete.cases(dd[, myCols]),]
      #print(nrow(dcc))
      #check for credible counts by arm
      #if((min(table(dcc$Arm))) > 10 )
      
      #{
      
      #lm handles null values without problems
      #create full model:  Mb ~ IC + Cohort + IC*Cohort
      #subset for complete cases.  TODO.  Do once instead of for every loop
      #dcc = dd[complete.cases(dd[, c(h,i,coCol)]),]
      m = lm(dcc[,h] ~ as.numeric(as.character(dcc[,i])) )
      dffitsMax = max(abs(dffits(m)))
      myCor = cor.test(dcc[,h] , as.numeric(as.character(dcc[,i])) )
      #mRed2 = lm(dcc[,h] ~ dcc[,i])
      
      #summary(m)$coefficients[4,4]
      corr = myCor$estimate
      pCorr = myCor$p.value
      pSlope = NA
      pSlope = try(summary(m)$coefficients[2,4])
      #append to output file
      write(paste(assays, ana1, ana2, corr, pCorr, dffitsMax, sep=","), dataOut, append=TRUE)
      
      #}
    }
    
    
  }
}

header = paste("assays", "ana1", "ana2", "Corr", "pCorr", "maxInfluence", sep=",")
write(header, dataOut, append=FALSE)

#for (k in 16:28)
for (k in 1:nrow(ccShort)){
	
	#k = 1
	print(paste("processing row ", k)) 
	oA = ccShort[k,1]
	oB = ccShort[k,2]
	#omeCols[[oA]]
	#omeCols[[oB]]
	try(regress(oA, oB, gsub("_", "", paste(myOmes[[oA]], myOmes[[oB]], sep="."))))
}



#TODO:  "Source" the regression function




library(table1)
table1(~ po_liver.fat.. + co_OGTT.BMI + co_OGTT.BMI.ile + po_ALT..U.L..CHC.only + co_Free.Androgen.index, data = dd)


