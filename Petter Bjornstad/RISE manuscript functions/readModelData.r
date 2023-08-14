#' Read IVGTT Modelling Data
#' 
#' `readModelData()` reads IVGTT modelling data from the RISE SAS Share library
#' 
#' @description This function will read the IVGTT modelling data stored on the
#'              RISE SAS share drive and adds 'label' and 'units' attributes to 
#'              each of the variables.  These attributes are used by statTable
#'              to automatically create row labels with (units).  See example 
#'              below.
#'          
#' @return
#'
#' @examples
#' library(StatAnalysis)
#' test.data = readModelData()
#' statTable(data=test.data,
#'           columns = {stats(by=visit)},
#'           rows = {
#'           'N': .n.
#'           'Sex': sex
#'           'Height':height
#'           basal.glu
#'           'Mean glucose': mean.glu
#'           basal.ins
#'           'Mean Insulin': mean.ins
#'           isrb.gref3
#'           pot.0
#'           sec.0
#'           dos.4
#'           dos.10.6
#'           })
#'
#' @export
readModelData = function(){
  library(RODBC)
  RISEshare <-odbcConnect("RISE Share",believeNRows = FALSE)
  # View(sqlColumns(DPPshare, "risework.ogttmodel"))
  # tbllist<-sqlTables(RISEshare) # list all tables
  Modeldata <- sqlQuery(RISEshare, 'select * from risework.ogttmodel', as.is=TRUE)
  odbcCloseAll()
  rm(RISEshare)
  
  # Convert all variable names to lowercase and replace _ with .
  names(Modeldata) = tolower(names(Modeldata))
  names(Modeldata) = gsub('_','.',names(Modeldata))
  
  # Change names beginning with .
  names(Modeldata)[names(Modeldata) == '.2.h.glu'] = 'glu.2h'
  names(Modeldata)[names(Modeldata) == '.2.h.ogis'] = 'ogis.2h'
  names(Modeldata)[names(Modeldata) == '.3.h.ogis'] = 'ogis.3h'
  
  # Add labels and units attribute to variables.  These are automatically used
  # by statTable to make row labels 
  attr(Modeldata$sex, 'label') = 'Sex (M/F)'
  
  attr(Modeldata$age, 'label') = 'Age'
  attr(Modeldata$age, 'units') = 'yrs'
  
  attr(Modeldata$height, 'label') = 'Height'
  attr(Modeldata$height, 'units') = 'm'
  
  attr(Modeldata$weight, 'label') = 'Weight'
  attr(Modeldata$weight, 'units') = 'kg'
  
  attr(Modeldata$bmi, 'label') = 'BMI'
  attr(Modeldata$bmi, 'units') = 'kg/m2'
  
  attr(Modeldata$bsa, 'label') = 'Body surface area'
  attr(Modeldata$bsa, 'units') = 'm2'
  
  attr(Modeldata$basal.glu, 'label') = 'Basal glucose'
  attr(Modeldata$basal.glu, 'units') = 'mmol/L'
  
  attr(Modeldata$mean.glu, 'label') = 'Mean glucose during the test, calculated from the glucose AUC'
  attr(Modeldata$mean.glu, 'units') = 'mmol/L'
  
  attr(Modeldata$basal.ins, 'label') = 'Basal insulin'
  attr(Modeldata$basal.ins, 'units') = 'pmol/L'
  
  attr(Modeldata$mean.ins, 'label') = 'Mean insulin during the test, calculated from the insulin AUC'
  attr(Modeldata$mean.ins, 'units') = 'pmol/L'
  
  attr(Modeldata$glu.2h, 'label') = '2-h glucose'
  attr(Modeldata$glu.2h, 'units') = 'mmol/L'
  
  attr(Modeldata$ogis.2h, 'label') = 'OGIS (2-h equation)'
  attr(Modeldata$ogis.2h, 'units') = 'ml min-1m-2'
  
  attr(Modeldata$ogis.3h, 'label') = 'OGIS (3-h equation)'
  attr(Modeldata$ogis.3h, 'units') = 'ml min-1m-2'
  
  attr(Modeldata$basal.isr, 'label') ='Basal insulin secretion'
  attr(Modeldata$basal.isr, 'units') ='pmol min-1m-2'
  
  attr(Modeldata$isr.gref1, 'label') ='Insulin secretion at 5 mmol/L glucose from the dose-response'
  attr(Modeldata$isr.gref1, 'units') ='pmol min-1m-2'
  
  attr(Modeldata$isr.gref2, 'label') ='Insulin secretion at 6.5 mmol/L glucose from the dose-response	(pmol min-1m-2)'
  attr(Modeldata$isr.gref2, 'units') ='pmol min-1m-2'

  attr(Modeldata$isr.gref3, 'label')	='Insulin secretion at 7 mmol/L glucose from the dose-response'
  attr(Modeldata$isr.gref3, 'units')	='pmol min-1m-2'
  
  attr(Modeldata$isrb.gref1, 'label')='Insulin secretion at 5 mmol/L glucose from the dose-response, adjusted for basal potentiation'
  attr(Modeldata$isrb.gref1, 'units')='pmol min-1m-2'
  
  attr(Modeldata$isrb.gref2, 'label') ='Insulin secretion at 6.5 mmol/L glucose from the dose-response, adjusted for basal potentiation'
  attr(Modeldata$isrb.gref2, 'units') ='pmol min-1m-2'
  
  attr(Modeldata$isrb.gref3, 'label')	='Insulin secretion at 7 mmol/L glucose from the dose-response, adjusted for basal potentiation'
  attr(Modeldata$isrb.gref3, 'units')	='pmol min-1m-2'
  
  attr(Modeldata$glu.sens, 'label')='Glucose sensitivity: mean slope of the dose-response in the observed glucose range'
  attr(Modeldata$glu.sens, 'units')='pmol min-1m-2mM-1'
  
  attr(Modeldata$rate.sens, 'label')='Rate sensitivity: parameter of the derivative component'
  attr(Modeldata$rate.sens, 'units')='pmol m-2mM-1'
  
  attr(Modeldata$pfr1, 'label')	= 'Potentiation factor ratio: ratio between the mean value in the intervals [100-120] and [0-20] min'
  attr(Modeldata$pfr1, 'units')	= ''
  
  attr(Modeldata$pfr2, 'label')	= 'Potentiation factor ratio: ratio between the mean value in the intervals [160-180] and [0-20] min'
  attr(Modeldata$pfr2, 'units')	= ''
  
  attr(Modeldata$total.isr, 'label')=	'Integral of total insulin secretion'
  attr(Modeldata$total.isr, 'units')=	'nmol m-2'
  
  attr(Modeldata$period, 'label')	='Time period for the calculation of the integral of insulin secretion'
  attr(Modeldata$period, 'units')	='min'
  
  attr(Modeldata$glu.tol, 'label')	='Glucose tolerance (ADA1997)'

  attr(Modeldata$stumvoll, 'label')=	'Stumvolls insulin sensitivity index'
  attr(Modeldata$stumvoll, 'units')=	'mL min-1kg-1'
  
  attr(Modeldata$matsuda, 'label')	='Matsudas insulin sensitivity index'
  attr(Modeldata$matsuda, 'units')	=''
  
  attr(Modeldata$clinsb, 'label')	='basal insulin clearance (basal insulin secretion/basal insulin)'
  attr(Modeldata$clinsb, 'units')	='L min-1m-2'
  
  attr(Modeldata$clins, 'label')	='OGTT insulin clearance (insulin secretion AUC/insulin AUC)'
  attr(Modeldata$clins, 'units')	='L min-1m-2'
  
  
  attr(Modeldata$dos.4, 'label')	=	'Insulin secretion at 4 mmol glucose concentration  (pmol min-1m-2)'
  attr(Modeldata$dos.5, 'label')	=	'Insulin secretion at 5 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.6, 'label')	=	'Insulin secretion at 6 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.7, 'label')	=	'Insulin secretion at 7 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.8, 'label')	=	'Insulin secretion at 8 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.9, 'label')	=	'Insulin secretion at 9 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.10, 'label')	=	'Insulin secretion at 10 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.11, 'label')	=	'Insulin secretion at 11 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.12, 'label')	=	'Insulin secretion at 12 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.13, 'label')	=	'Insulin secretion at 13 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.14, 'label')	=	'Insulin secretion at 14 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.15, 'label')	=	'Insulin secretion at 15 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.10.2, 'label')	=	'Insulin secretion at 10.2 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.10.4, 'label')	=	'Insulin secretion at 10.4 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.10.6, 'label')	=	'Insulin secretion at 10.6 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.10.8, 'label')	=	'Insulin secretion at 10.8 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.11.2, 'label')	=	'Insulin secretion at 11.2 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.11.4, 'label')	=	'Insulin secretion at 11.4 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.11.6, 'label')	=	'Insulin secretion at 11.6 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.11.8, 'label')	=	'Insulin secretion at 11.8 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.12.2, 'label')	=	'Insulin secretion at 12.2 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.12.4, 'label')	=	'Insulin secretion at 12.4 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.12.6, 'label')	=	'Insulin secretion at 12.6 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.12.8, 'label')	=	'Insulin secretion at 12.8 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.13.2, 'label')	=	'Insulin secretion at 13.2 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.13.4, 'label')	=	'Insulin secretion at 13.4 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.13.6, 'label')	=	'Insulin secretion at 13.6 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.13.8, 'label')	=	'Insulin secretion at 13.8 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.14.2, 'label')	=	'Insulin secretion at 14.2 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.14.4, 'label')	=	'Insulin secretion at 14.4 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.14.6, 'label')	=	'Insulin secretion at 14.6 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.14.8, 'label')	=	'Insulin secretion at 14.8 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.15.2, 'label')	=	'Insulin secretion at 15.2 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.15.4, 'label')	=	'Insulin secretion at 15.4 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.4.2, 'label')	=	'Insulin secretion at 4.2 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.4.4, 'label')	=	'Insulin secretion at 4.4 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.4.6, 'label')	=	'Insulin secretion at 4.6 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.4.8, 'label')	=	'Insulin secretion at 4.8 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.5.2, 'label')	=	'Insulin secretion at 5.2 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.5.4, 'label')	=	'Insulin secretion at 5.4 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.5.6, 'label')	=	'Insulin secretion at 5.6 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.5.8, 'label')	=	'Insulin secretion at 5.8 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.6.2, 'label')	=	'Insulin secretion at 6.2 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.6.4, 'label')	=	'Insulin secretion at 6.4 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.6.6, 'label')	=	'Insulin secretion at 6.6 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.6.8, 'label')	=	'Insulin secretion at 6.8 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.7.2, 'label')	=	'Insulin secretion at 7.2 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.7.4, 'label')	=	'Insulin secretion at 7.4 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.7.6, 'label')	=	'Insulin secretion at 7.6 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.7.8, 'label')	=	'Insulin secretion at 7.8 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.8.2, 'label')	=	'Insulin secretion at 8.2 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.8.4, 'label')	=	'Insulin secretion at 8.4 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.8.6, 'label')	=	'Insulin secretion at 8.6 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.8.8, 'label')	=	'Insulin secretion at 8.8 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.9.2, 'label')	=	'Insulin secretion at 9.2 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.9.4, 'label')	=	'Insulin secretion at 9.4 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.9.6, 'label')	=	'Insulin secretion at 9.6 mmol glucose concentration (pmol min-1m-2)'
  attr(Modeldata$dos.9.8, 'label')	=	'Insulin secretion at 9.8 mmol glucose concentration (pmol min-1m-2)'
  
  attr(Modeldata$dos.4, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.5, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.6, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.7, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.8, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.9, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.10, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.11, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.12, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.13, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.14, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.15, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.10.2, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.10.4, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.10.6, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.10.8, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.11.2, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.11.4, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.11.6, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.11.8, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.12.2, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.12.4, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.12.6, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.12.8, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.13.2, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.13.4, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.13.6, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.13.8, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.14.2, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.14.4, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.14.6, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.14.8, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.15.2, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.15.4, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.4.2, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.4.4, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.4.6, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.4.8, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.5.2, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.5.4, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.5.6, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.5.8, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.6.2, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.6.4, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.6.6, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.6.8, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.7.2, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.7.4, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.7.6, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.7.8, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.8.2, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.8.4, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.8.6, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.8.8, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.9.2, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.9.4, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.9.6, 'units')	=	'pmol min-1m-2'
  attr(Modeldata$dos.9.8, 'units')	=	'pmol min-1m-2'
  
  attr(Modeldata$sec.0, 'label')		='Secretion value at time 0 min'
  attr(Modeldata$sec.5, 'label')		='Secretion value at time 5 min'
  attr(Modeldata$sec.10, 'label')		='Secretion value at time 10 min'
  attr(Modeldata$sec.15, 'label')		='Secretion value at time 15 min'
  attr(Modeldata$sec.20, 'label')		='Secretion value at time 20 min'
  attr(Modeldata$sec.25, 'label')		='Secretion value at time 25 min'
  attr(Modeldata$sec.30, 'label')		='Secretion value at time 30 min'
  attr(Modeldata$sec.35, 'label')		='Secretion value at time 35 min'
  attr(Modeldata$sec.40, 'label')		='Secretion value at time 40 min'
  attr(Modeldata$sec.45, 'label')		='Secretion value at time 45 min'
  attr(Modeldata$sec.50, 'label')		='Secretion value at time 50 min'
  attr(Modeldata$sec.55, 'label')		='Secretion value at time 55 min'
  attr(Modeldata$sec.60, 'label')		='Secretion value at time 60 min'
  attr(Modeldata$sec.65, 'label')		='Secretion value at time 65 min'
  attr(Modeldata$sec.70, 'label')		='Secretion value at time 70 min'
  attr(Modeldata$sec.75, 'label')		='Secretion value at time 75 min'
  attr(Modeldata$sec.80, 'label')		='Secretion value at time 80 min'
  attr(Modeldata$sec.85, 'label')		='Secretion value at time 85 min'
  attr(Modeldata$sec.90, 'label')		='Secretion value at time 90 min'
  attr(Modeldata$sec.95, 'label')		='Secretion value at time 95 min'
  attr(Modeldata$sec.100, 'label')		='Secretion value at time 100 min'
  attr(Modeldata$sec.105, 'label')		='Secretion value at time 105 min'
  attr(Modeldata$sec.110, 'label')		='Secretion value at time 110 min'
  attr(Modeldata$sec.115, 'label')		='Secretion value at time 115 min'
  attr(Modeldata$sec.120, 'label')		='Secretion value at time 120 min'
  attr(Modeldata$sec.125, 'label')		='Secretion value at time 125 min'
  attr(Modeldata$sec.130, 'label')		='Secretion value at time 130 min'
  attr(Modeldata$sec.135, 'label')		='Secretion value at time 135 min'
  attr(Modeldata$sec.140, 'label')		='Secretion value at time 140 min'
  attr(Modeldata$sec.145, 'label')		='Secretion value at time 145 min'
  attr(Modeldata$sec.150, 'label')		='Secretion value at time 150 min'
  attr(Modeldata$sec.155, 'label')		='Secretion value at time 155 min'
  attr(Modeldata$sec.160, 'label')		='Secretion value at time 160 min'
  attr(Modeldata$sec.165, 'label')		='Secretion value at time 165 min'
  attr(Modeldata$sec.170, 'label')		='Secretion value at time 170 min'
  attr(Modeldata$sec.175, 'label')		='Secretion value at time 175 min'
  attr(Modeldata$sec.180, 'label')		='Secretion value at time 180 min'
  
  attr(Modeldata$sec.0, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.5, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.10, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.15, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.20, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.25, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.30, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.35, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.40, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.45, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.50, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.55, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.60, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.65, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.70, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.75, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.80, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.85, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.90, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.95, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.100, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.105, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.110, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.115, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.120, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.125, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.130, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.135, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.140, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.145, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.150, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.155, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.160, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.165, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.170, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.175, 'units')		='pmol/min/m^2'
  attr(Modeldata$sec.180, 'units')		='pmol/min/m^2'
  

  attr(Modeldata$pot.0, 'label')		='Potentiation value at time 0 min'
  attr(Modeldata$pot.5, 'label')		='Potentiation value at time 5 min'
  attr(Modeldata$pot.10, 'label')		='Potentiation value at time 10 min'
  attr(Modeldata$pot.15, 'label')		='Potentiation value at time 15 min'
  attr(Modeldata$pot.20, 'label')		='Potentiation value at time 20 min'
  attr(Modeldata$pot.25, 'label')		='Potentiation value at time 25 min'
  attr(Modeldata$pot.30, 'label')		='Potentiation value at time 30 min'
  attr(Modeldata$pot.35, 'label')		='Potentiation value at time 35 min'
  attr(Modeldata$pot.40, 'label')		='Potentiation value at time 40 min'
  attr(Modeldata$pot.45, 'label')		='Potentiation value at time 45 min'
  attr(Modeldata$pot.50, 'label')		='Potentiation value at time 50 min'
  attr(Modeldata$pot.55, 'label')		='Potentiation value at time 55 min'
  attr(Modeldata$pot.60, 'label')		='Potentiation value at time 60 min'
  attr(Modeldata$pot.65, 'label')		='Potentiation value at time 65 min'
  attr(Modeldata$pot.70, 'label')		='Potentiation value at time 70 min'
  attr(Modeldata$pot.75, 'label')		='Potentiation value at time 75 min'
  attr(Modeldata$pot.80, 'label')		='Potentiation value at time 80 min'
  attr(Modeldata$pot.85, 'label')		='Potentiation value at time 85 min'
  attr(Modeldata$pot.90, 'label')		='Potentiation value at time 90 min'
  attr(Modeldata$pot.95, 'label')		='Potentiation value at time 95 min'
  attr(Modeldata$pot.100, 'label')		='Potentiation value at time 100 min'
  attr(Modeldata$pot.105, 'label')		='Potentiation value at time 105 min'
  attr(Modeldata$pot.110, 'label')		='Potentiation value at time 110 min'
  attr(Modeldata$pot.115, 'label')		='Potentiation value at time 115 min'
  attr(Modeldata$pot.120, 'label')		='Potentiation value at time 120 min'
  attr(Modeldata$pot.125, 'label')		='Potentiation value at time 125 min'
  attr(Modeldata$pot.130, 'label')		='Potentiation value at time 130 min'
  attr(Modeldata$pot.135, 'label')		='Potentiation value at time 135 min'
  attr(Modeldata$pot.140, 'label')		='Potentiation value at time 140 min'
  attr(Modeldata$pot.145, 'label')		='Potentiation value at time 145 min'
  attr(Modeldata$pot.150, 'label')		='Potentiation value at time 150 min'
  attr(Modeldata$pot.155, 'label')		='Potentiation value at time 155 min'
  attr(Modeldata$pot.160, 'label')		='Potentiation value at time 160 min'
  attr(Modeldata$pot.165, 'label')		='Potentiation value at time 165 min'
  attr(Modeldata$pot.170, 'label')		='Potentiation value at time 170 min'
  attr(Modeldata$pot.175, 'label')		='Potentiation value at time 175 min'
  attr(Modeldata$pot.180, 'label')		='Potentiation value at time 180 min'

  attr(Modeldata$pot.0, 'units')		=''
  attr(Modeldata$pot.5, 'units')		=''
  attr(Modeldata$pot.10, 'units')		=''
  attr(Modeldata$pot.15, 'units')		=''
  attr(Modeldata$pot.20, 'units')		=''
  attr(Modeldata$pot.25, 'units')		=''
  attr(Modeldata$pot.30, 'units')		=''
  attr(Modeldata$pot.35, 'units')		=''
  attr(Modeldata$pot.40, 'units')		=''
  attr(Modeldata$pot.45, 'units')		=''
  attr(Modeldata$pot.50, 'units')		=''
  attr(Modeldata$pot.55, 'units')		=''
  attr(Modeldata$pot.60, 'units')		=''
  attr(Modeldata$pot.65, 'units')		=''
  attr(Modeldata$pot.70, 'units')		=''
  attr(Modeldata$pot.75, 'units')		=''
  attr(Modeldata$pot.80, 'units')		=''
  attr(Modeldata$pot.85, 'units')		=''
  attr(Modeldata$pot.90, 'units')		=''
  attr(Modeldata$pot.95, 'units')		=''
  attr(Modeldata$pot.100, 'units')		=''
  attr(Modeldata$pot.105, 'units')		=''
  attr(Modeldata$pot.110, 'units')		=''
  attr(Modeldata$pot.115, 'units')		=''
  attr(Modeldata$pot.120, 'units')		=''
  attr(Modeldata$pot.125, 'units')		=''
  attr(Modeldata$pot.130, 'units')		=''
  attr(Modeldata$pot.135, 'units')		=''
  attr(Modeldata$pot.140, 'units')		=''
  attr(Modeldata$pot.145, 'units')		=''
  attr(Modeldata$pot.150, 'units')		=''
  attr(Modeldata$pot.155, 'units')		=''
  attr(Modeldata$pot.160, 'units')		=''
  attr(Modeldata$pot.165, 'units')		=''
  attr(Modeldata$pot.170, 'units')		=''
  attr(Modeldata$pot.175, 'units')		=''
  attr(Modeldata$pot.180, 'units')		=''
  
  return(Modeldata)
  
}