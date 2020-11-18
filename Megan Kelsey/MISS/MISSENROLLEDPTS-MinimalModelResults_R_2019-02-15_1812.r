graphics.off()
#Load Hmisc library
library(Hmisc)
#Read Data
data=read.csv('H:/Endocrinology/Kelsey/Kelsey microbiome/Raw data/MISSENROLLEDPTS-MinimalModelResults_DATA_2019-02-15_1812.csv')
#Setting Labels

label(data$miss_id)="MISS ID"
label(data$date_of_study_visit_mm)="Date of Study Visit"
label(data$insulin_sensitivity)="Insulin Sensitivity (Si)"
label(data$insulin_sensitivity_fsd)="Insulin Sensitivity (Si) FSD"
label(data$insulin_secretion_mm)="Insulin Secretion (AIRg)"
label(data$disposition_index)="Disposition Index (DI)"
label(data$disposition_index_fsd)="Disposition Index (DI) FSD"
label(data$sg_mm)="Glucose Effectiveness (Sg)"
label(data$gezi_mm)="Glucose Effectiveness @ 0 Insulin (GEZI)"
label(data$kg_mm)="Glucose Tolerance (Kg)"
label(data$minimal_model_complete)="Complete?"
#Setting Units


#Setting Factors(will create new variable for factors)
data$minimal_model_complete.factor = factor(data$minimal_model_complete,levels=c("0","1","2"))

levels(data$minimal_model_complete.factor)=c("Incomplete","Unverified","Complete")
