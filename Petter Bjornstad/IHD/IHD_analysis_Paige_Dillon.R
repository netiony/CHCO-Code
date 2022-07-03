library(Hmisc)
library(knitr)
library(Table1)
library(redcapAPI)
library(dplyr)
library(ISLR)
library(stringr)

#setwd("/Users/paigedillon/Desktop")
data1 <- read.csv("/Volumes/som/PEDS/RI Biostatistics Core/Shared/Shared Material/Shared Code/R/Paige_project1/Final_Dataset_01.csv", na.strings = c(" ", "", "-99"))

data2 <- read.csv("/Volumes/som/PEDS/RI Biostatistics Core/Shared/Shared Material/Shared Code/R/Paige_project1/Table01data.csv", na.strings = c(" ", "", "-99"))

data3 <- read.csv("/Volumes/som/PEDS/RI Biostatistics Core/Shared/Shared Material/Shared Code/R/Paige_project1/Final_Dataset_02.csv", na.strings = c(" ", "", "-99"))
data4 <- read.csv("/Volumes/som/PEDS/RI Biostatistics Core/Shared/Shared Material/Shared Code/R/Paige_project1/Final_Dataset_03.csv", na.strings = c(" ", "", "-99"))

#data2$Subject_ID <- data2$RENALHEIR.Table1DataCarson_DATA_LABELS_2022.03.16_2020
finaldata <- merge(data1, data2, by=("Subject_ID"), all.x = TRUE, all.y = TRUE)

finaldata2 <- merge(data3, data4, by="Subject_ID" , all.x = TRUE , all.y = TRUE)

#finaldata$X.FFA_Suppression..x. <- "X.FFA.Suppression"
#finaldata2$X.FFA.Suppression <- "X.FFA.Suppression"


merged_data <- merge(finaldata, finaldata2, by = intersect(names(finaldata), names(finaldata2))  ,  all.x = TRUE , all.y = TRUE)
merged_data$X_FFA_Supression <- paste(finaldata$X.FFA_Suppression..x. , finaldata2$X.FFA.Suppression)
merged_data$X.FFA.Suppression <-NULL
merged_data$X.FFA_Suppression..x. <-NULL
#merged_data$X_FFA_Supression <- str_c(finaldata$X.FFA_Suppression..x., '' , finaldata2$X.FFA.Suppression)

#creating table
#install.packages('Hmisc')
library(Hmisc)
setwd("/Users/paigedillon/Desktop")
source('/Volumes/som/PEDS/RI Biostatistics Core/Shared/Shared Material/Shared Code/R/Paige_project1/table1.R')
#source('/Volumes/som/PEDS/RI Biostatistics Core/Shared/Shared Material/Shared Code/R/Paige_project1/IHD_analysis.R')

#dist_check(merged_data)
merged_data$group <- ifelse(str_detect(merged_data$Subject_ID, "T"), "T2D", 
                            ifelse(str_detect(merged_data$Subject_ID, "L"), "Lean", "Obese"))
merged_data$group <- as.factor(merged_data$group)

# labels for table 01
label(merged_data$Sex)= "Sex"
label(merged_data$Age.at.time.of.Consent)= "Age at time of consent"
label(merged_data$Race..choice.American.Indian.or.Alaskan.Native.)= "Is your race American Indian or Alaskan Native?"
label(merged_data$Race..choice.Asian.)= "Is your race Asian?"
label(merged_data$Race..choice.Hawaiian.or.Pacific.Islander.)= "Is your race Hawaiian or Pacific Islander?"
label(merged_data$Race..choice.White.)= "Is your race white?"
label(merged_data$Race..choice.Unknown.)= "Race Unknown"
label(merged_data$Race..choice.Other.)= "Race Other"
label(merged_data$Race..if.Other)= "Is your race different than those previously asked?"
label(merged_data$Race..choice.Black.or.African.American.)= "Is your race Black or African American?"
label(merged_data$Ethnicity..choice.Hispanic.)= "Are you hispanic?"
label(merged_data$Ethnicity..choice.Non.Hispanic.)= "Are you non-Hispanic?"
label(merged_data$Ethnicity..choice.Unknown.Not.Reported.)= "Ethnicity Unknown or Not Reported"
label(merged_data$Screening.A1c.....)= "Screening A1C"
label(merged_data$Length.of.Diabetes.Diagnosis..months.)= "Length of diabetes diagnosis in months"
label(merged_data$Height..cm.)= "Height in Cm"
label(merged_data$Weight_kg)= "Weight in kg"
label(merged_data$BMI)= "BMI"
label(merged_data$Body.Fat....)= "Body Fat"
label(merged_data$X.Avg..Waist.Circumference..cm.)= "Average Waist Circumference"
label(merged_data$Baseline.Cholesterol..mg.dl.)= "Baseline Cholesterol (mg/dl)"
label(merged_data$Baseline.HDL..mg.dl.)= "Baseline HDL (mg/dl)"
label(merged_data$Baseline.LDL..mg.dl.)= "Baseline LDL (mg/dl)"
label(merged_data$Baseline.Triglycerides..mg.dl.)= "Baseline Triglycerides (mg/dl)"
label(merged_data$Systolic.Blood.Pressure..SBP.)= "Systolic Blood Pressure (SBP)"
label(merged_data$Diastolic.Blood.Pressure..DBP.)= "Diastolic Blood Pressure (DBP)"
label(merged_data$Mean.Arterial.Pressure..MAP.)= "Mean Arterial Pressure (MAP)"
label(merged_data$Serum.Creatinine...mg.dl.)= "Serum Creatinine (mg/dl)"
label(merged_data$Baseline.Serum.Cystatin.C..clamp.)= "Baseline Serum Cystatin"
label(merged_data$BSA.normalized.GFR)= "BSA normalized GFR"
label(merged_data$Absolute.GFR)= "Absolute GFR"
label(merged_data$Absolute.ERPF)= "Absolute ERPF"
label(merged_data$BSA.normalized.ERPF)= "BSA normalized ERPF"
label(merged_data$O2_Consumption._cortex)= "O2 consumption- cortex"
label(merged_data$O2_Consumption_medulla)= "O2 consumption medulla"
label(merged_data$O2_Consumption_Kidney)= "O2 consumption kidney"

#Labeling the other variables
label(merged_data$Subject_ID)= "Subject ID"
label(merged_data$ACR_baseline)= "ACR Baseline"
label(merged_data$GFR_BSAnormalized)= "GFR BSA Normalized"
label(merged_data$ERPF_BSAnormalized)= "ERPF BSA Normalized"
label(merged_data$Gender)= "Gender"
label(merged_data$Date.of.Diagnosis.with.Diabetes)= "Date of Diagnosis with Diabetes"
label(merged_data$Weight..kg.)= "Weight (kg)"
label(merged_data$Fat.Mass..kg.)= "Fat Mass (kg)"
label(merged_data$Notes)= "Notes"
label(merged_data$AIRg)= "AIRg"
label(merged_data$SS.insulin)= "SS. insulin"
label(merged_data$Disposition.index..M.I.xAIRg)= "Disposition Index (M.I.x AIRg)"
label(merged_data$Baseline.Serum.creatinine..clamp.)= "Baseline Serum Creatinine"
label(merged_data$Baseline.ACR)= "Baseline ACR"
label(merged_data$X_FFA_Supression)= "X.FFA. Suppression"
label(merged_data$Date.of.Birth)= "Date of Birth"
label(merged_data$Baseline.Cystatin.C..mg.L.)= "Baseline Cystatin (mg/L)"
label(merged_data$ACPRg..nmol.L.)= "ACPRg (nmol/L)"
label(merged_data$SS.CP)= " SS CP"
label(merged_data$GIR)= "GIR"
label(merged_data$group)= "Group"

#Table 01
#tab.1_bygrp<-final_table(merged_data,c('Sex','Age.at.time.of.Consent','Race..choice.American.Indian.or.Alaskan.Native.','Race..choice.Asian.', 'Race..choice.Black.or.African.American.','Race..choice.White.', 'Race..choice.Unknown.', 'Race..choice.Other.', 'Race..if.Other','Ethnicity..choice.Hispanic.', 'Ethnicity..choice.Non.Hispanic.','Screening.A1c.....', 'Length.of.Diabetes.Diagnosis..months.' ,'Height..cm.' , 'Weight_kg', 'BMI', 'Body.Fat....', 'X.Avg..Waist.Circumference..cm.','Baseline.Cholesterol..mg.dl.', 'Baseline.HDL..mg.dl.', 'Baseline.LDL..mg.dl.', 'Baseline.Triglycerides..mg.dl.', 'Systolic.Blood.Pressure..SBP.', 'Diastolic.Blood.Pressure..DBP.','Mean.Arterial.Pressure..MAP.', 'Serum.Creatinine...mg.dl.', 'Baseline.Serum.Cystatin.C..clamp.','BSA.normalized.GFR', 'Absolute.GFR', 'Absolute.ERPF','BSA.normalized.ERPF', 'O2_Consumption._cortex', 'O2_Consumption_medulla','O2_Consumption_Kidney' ),
                         #merged_data$group,margin=2,single=F,ron=1,col.names=T, summary.stat='mean')
#tab.1_bygrp

#length of diabetes in months is what was messing up the group, maybe we'll have to manually add this to the table later

tab.1_bygrp<-final_table(merged_data,c('Sex','Age.at.time.of.Consent','Race..choice.American.Indian.or.Alaskan.Native.','Race..choice.Asian.', 'Race..choice.Black.or.African.American.','Race..choice.White.', 'Race..choice.Unknown.', 'Race..choice.Other.', 'Race..if.Other','Ethnicity..choice.Hispanic.', 'Ethnicity..choice.Non.Hispanic.','Screening.A1c.....','Height..cm.', 'Weight_kg', 'BMI', 'Body.Fat....', 'X.Avg..Waist.Circumference..cm.','Baseline.Cholesterol..mg.dl.', 'Baseline.HDL..mg.dl.', 'Baseline.LDL..mg.dl.', 'Baseline.Triglycerides..mg.dl.', 'Systolic.Blood.Pressure..SBP.', 'Diastolic.Blood.Pressure..DBP.','Mean.Arterial.Pressure..MAP.', 'Serum.Creatinine...mg.dl.', 'Baseline.Serum.Cystatin.C..clamp.','BSA.normalized.GFR', 'Absolute.GFR', 'Absolute.ERPF','BSA.normalized.ERPF', 'O2_Consumption._cortex', 'O2_Consumption_medulla','O2_Consumption_Kidney'),
                         merged_data$group,margin=2,single=F,ron=1,col.names=T, summary.stat='mean')
tab.1_bygrp
#what I took out bc 0- Race..choice.Hawaiian.or.Pacific.Islander, Ethnicity..choice.Unknown.Not.Reported.
#Weight_kg and Weight...kg are different variables or at least have different # in them so I did not merge them




