
library("jsonlite")
library("RColorBrewer")
cohortColors <- brewer.pal(4,"Dark2")

#For ms, July 2020
dataIn = "/Volumes/RI Biostatistics Core/Shared/Shared Projects/Bothwell/Peds ENDO/Siebert - Multiomics/Data/merged multiomics data clean.csv"
models = "/Volumes/RI Biostatistics Core/Shared/Shared Projects/Bothwell/Peds ENDO/Siebert - Multiomics/Data/omes_models_filtered_9_22_2023.csv"
jsonOut = "/Volumes/RI Biostatistics Core/Shared/Shared Projects/Bothwell/Peds ENDO/Siebert - Multiomics/Data/omes_json_models_9_22_2023.json"


##########

df =  read.csv(dataIn)
df$Cohort = "All Patients"
cohort = levels(as.factor(df$Cohort))

#retrieve assays
n = nrow(df)
#nCo override

nCo = 1

#Load regression results 

bosTop = read.csv(models)
bosTop$assays = as.factor(bosTop$assays)
nlevels(bosTop$assays)



#derive assays and colors
assayPairs = levels(bosTop$assays)
assayNames = sort(paste(unique(unlist(lapply(assayPairs, function(x) strsplit(x, ".", fixed=TRUE)))), "_", sep=""))
assays = list()
nAssays = length(assayNames)
myColors <- brewer.pal((nAssays+2),"Set3")[3:(nAssays+2)]
assays[assayNames] = myColors[1:nAssays]

 
mySummary = NULL
#data for microPlot glyph
mPlot = NULL

#sample-level data to support detailed plots
myDetail = df[, c("ID", "Cohort")]

for (myRow in 1:nrow(bosTop))
{
	
	# Label axes 
	yAna = as.character(bosTop$ana2[myRow])
	xAna = as.character(bosTop$ana1[myRow])
	myYlab = yAna
	myXlab = xAna
		
	yInd = which(colnames(df)==yAna)
	xInd = which(colnames(df)==xAna)
		
	#Determine range of "x" analyte 
	#(by Cohort)
	ra = range(df[,xInd], na.rm=TRUE)

	#If these analytes have not been accessed before, add them to myDetail
	xDetail = df[xInd]
	colnames(xDetail) = myXlab
	yDetail = df[yInd]
	colnames(yDetail) = myYlab
	if (! (myYlab %in% colnames(myDetail))) 
		myDetail = cbind(myDetail,yDetail)
	if (! (myXlab %in% colnames(myDetail))) 
		myDetail = cbind(myDetail,xDetail)

	#build regression model	
	m = lm(df[,yInd] ~ df[,xInd]) 

	#effect in terms of 1 unit increase in X 	
	estEff = summary(m)$coefficients[2,1]
	#pVal reference slope
	pCorr = summary(m)$coefficients[2,4]
	#influence data
	infl = round(abs(dffits(m)),2)
    inflLabel = paste("infl_",myRow, sep="")

	myDetail$infl=""
    myDetail$infl[as.numeric(names(infl))] = infl
    colnames(myDetail)[ncol(myDetail)] = inflLabel
	
	#extract data from model
	coefs = coef(m)
	intercept = NULL
	slope = NULL
	i1 = coefs[1] 
	s1 = coefs[2] 
	intercept = c(intercept, i1)
	slope = c(slope, s1)
	#Coefficients are ordered as follows:  intercept, slope (for reference); intercepts for remaining levels, slopes for remaining
	#Thus, for 1st non-reference: intercept, slope are at 3, 3 + (nlevels - 1); where 3 = 2 + (levelPos - 1)
	#Thus, for 2nd non-reference: intercept, slope are at 4, 4 + (nlevels - 1); where 4 = 2 + (levelPos - 2)
	#In this case, nlevels = nCo
	
	#Let c represent the cohort (levelPosition - 1)
	if (nCo > 1)
	{
		for (c in 1:(nCo - 1))
		{
			offset = c + 2
			iC = i1 + coefs[offset]
			sC = s1 + coefs[offset + (nCo -1)]
			intercept = c(intercept, iC)
			slope = c(slope, sC)
			#print(offset)
			#print(coefs[offset])
			#print(coefs[offset + (nCo - 1)])
		}
	}

	
	#Walk the line(s)--traverse the intercept/slope lists and derive the point sets (xMin, y1), (xMax, y2)
	#Group/cohort-level xMin and xMax values allow us to draw lines that are limited to the 
	#observed range of each group
	#Points are recorded as sequence of 4, with (xMin, y1), (xMax, y2) semantics assumed
	#JavaScript needs the endpoints of the lines, as opposed to slope/intercept
	
	myPoints = NULL
	for (li in 1:length(intercept)){
		#TODO--generalize
		#xMin = ra[[li]][1]
		#xMax = ra[[li]][2]
		xMin = ra[1]
		xMax = ra[2]

		y1 = unname(slope[li] * xMin + intercept[li])
		y2 = unname(slope[li] * xMax + intercept[li])
		myPoints = c(myPoints, xMin, y1, xMax, y2)
	}
	
	#derive "global" min and max and append to myPoints
	#These let us specify range of axes in the plots
	xMin = min(df[xInd], na.rm=TRUE)
	xMax = max(df[xInd], na.rm=TRUE)
	yMin = min(df[yInd], na.rm=TRUE)
	yMax = max(df[yInd], na.rm=TRUE)
	#semantics xMin, yMin, xMax, yMax
	myPoints = c(myPoints, xMin, yMin, xMax, yMax)
	
	thisRow = data.frame(matrix(myPoints, nrow=1))
	
	thisSummary = data.frame(key=myRow, Assays=bosTop[myRow, "assays"], Analyte1=myYlab, Analyte2=myXlab, pAdjFDR = bosTop[myRow, "pAdjFDR"], Corr=bosTop[myRow, "Corr"], estEff=estEff, maxInfluence=bosTop[myRow, "maxInfluence"], row.names=NULL)
		
	mySummary = rbind(mySummary, thisSummary)	
	mPlot = rbind(mPlot, thisRow)
}



#Combine the data frames, knowing they have same number of rows in same order
mySummary$mPlot = mPlot


# Filter out estimated effects that are close to 0 - Find a better way to do this
mySummary <- mySummary[!(abs(mySummary$estEff)) < 0.0000001,]

myStatement2 = sprintf("The study includes approximately %i participants.  estEff represents the estimated amount by which the outcome (Y-axis) changes for a unit increase in X. Bonferroni adjustments are made to co.pg and co.pa pairs.", n)


# myStatement = paste(myStatement, myStatement2, sep="")



myConfig = list(assayAware="true", assays=assays, influenceAware="true", cohort=cohort, colors=cohortColors, statement=myStatement2)
 
#toJSON defaults to digits=4
myObject = list(configuration=myConfig, statement=myStatement2, config=cohort, colors=cohortColors, summary=mySummary, detail=myDetail )
oJ = toJSON(myObject, pretty=TRUE, dataframe="columns", digits=3)
write(oJ, jsonOut)


