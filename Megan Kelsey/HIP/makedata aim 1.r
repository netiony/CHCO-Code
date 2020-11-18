#setwd("H:\\Endocrinology\\Kelsey\\Kesley HIP study ISS puberty\\Ergui poster\\Data for models")
setwd("H:\\Endocrinology\\Kelsey\\Kesley HIP study ISS puberty\\Final dataset cleaning\\Aim 1")

library(tidyr)
library(reshape2)

#Clear existing data and graphics
rm(list=ls())
graphics.off()
#Load Hmisc library
library(Hmisc)
#Read Data
data=read.csv('ivgtt data.csv',na.strings = c("-99"))
#Setting Labels

label(data$hip_id_screen)="HIP ID"
label(data$redcap_event_name)="Event Name"
label(data$date_of_study_visit)="Date of Study Visit "
label(data$study_visit_number_svl)="Study Visit Number"
label(data$tcholes_lv)="Total Cholesterol"
label(data$tg_lv)="Triglycerides"
label(data$hdl_lv)="HDL"
label(data$ldl_lv)="LDL"
label(data$ldl_meth_lv)="LDL Method"
label(data$dhea_s)="DHEA-S"
label(data$estradiol)="Estradiol"
label(data$total_testosterone)="Total Testosterone"
label(data$shbg)="SHBG"
label(data$adiponect_lv)="Adiponectin"
label(data$lept)="Leptin"
label(data$gluc_neg15_iv)="-15 minute glucose"
label(data$gluc_0_iv)="0 minute glucose"
label(data$gluc_2_iv)="2 minute glucose"
label(data$gluc_3_iv)="3 minute glucose"
label(data$gluc_4_iv)="4 minute glucose"
label(data$gluc_5_iv)="5minute glucose"
label(data$gluc_6_iv)="6 minute glucose"
label(data$gluc_8_iv)="8 minute glucose"
label(data$gluc_10_iv)="10 minute glucose"
label(data$gluc_12_iv)="12 minute glucose"
label(data$gluc_14_iv)="14 minute glucose"
label(data$gluc_16_iv)="16 minute glucose"
label(data$gluc_19_iv)="19 minute glucose"
label(data$gluc_22_iv)="22 minute glucose"
label(data$gluc_25_iv)="25 minute glucose"
label(data$gluc_30_iv)="30 minute glucose"
label(data$gluc_35_iv)="35 minute glucose"
label(data$gluc_40_iv)="40 minute glucose"
label(data$gluc_50_iv)="50 minute glucose"
label(data$gluc_60_iv)="60 minute glucose"
label(data$gluc_70_iv)="70 minute glucose"
label(data$gluc_80_iv)="80 minute glucose"
label(data$gluc_90_iv)="90 minute glucose"
label(data$gluc_100_iv)="100 minute glucose"
label(data$gluc_120_iv)="120 minute glucose"
label(data$gluc_140_iv)="140 minute glucose"
label(data$gluc_160_iv)="160 minute glucose"
label(data$gluc_180_iv)="180 minute glucose"
label(data$ins_neg15_iv)="-15 minute insulin"
label(data$ins_0_iv)="0 minute insulin"
label(data$ins_2_iv)="2 minute insulin"
label(data$ins_3_iv)="3 minute insulin"
label(data$ins_4_iv)="4 minute insulin"
label(data$ins_5_iv)="5 minute insulin"
label(data$ins_6_iv)="6 minute insulin"
label(data$ins_8_iv)="8 minute insulin"
label(data$ins_10_iv)="10 minute insulin"
label(data$ins_12_iv)="12 minute insulin"
label(data$ins_14_iv)="14 minute insulin"
label(data$ins_16_iv)="16 minute insulin"
label(data$ins_19_iv)="19 minute insulin"
label(data$ins_22_iv)="22 minute insulin"
label(data$ins_25_iv)="25 minute insulin"
label(data$ins_30_iv)="30 minute insulin"
label(data$ins_35_iv)="35 minute insulin"
label(data$ins_40_iv)="40 minute insulin"
label(data$ins_50_iv)="50 minute insulin"
label(data$ins_60_iv)="60 minute insulin"
label(data$ins_70_iv)="70 minute insulin"
label(data$ins_80_iv)="80 minute insulin"
label(data$ins_90_iv)="90 minute insulin"
label(data$ins_100_iv)="100 minute insulin"
label(data$ins_120_iv)="120 minute insulin"
label(data$ins_140_iv)="140 minute insulin"
label(data$ins_160_iv)="160 minute insulin"
label(data$ins_180_iv)="180 minute insulin"
label(data$study_visit_adult_core_labs_complete)="Complete?"
label(data$dt_of_study_visit)="Date of Study Visit "
label(data$sv_num)="Study Visit Number"
label(data$ast)="AST"
label(data$alt)="ALT"
label(data$crp_lv)="CRP"
label(data$igf_lv)="IGF-1"
label(data$hba)="HbA1c"
label(data$study_visit_other_labs_complete)="Complete?"
label(data$dov_lv)="Date of Study Visit "
label(data$study_visit_number)="Study Visit Number"
label(data$urine_sv)="Was a urine sample collected at this visit? "
label(data$height_lv)="Height (cm)"
label(data$weight_lv)="Weight (kilograms)"
label(data$bmi_screen_v1)="BMI"
label(data$bmi_pct_screen_v1)="BMI Percentile"
label(data$zscore_lv)="z-score"
label(data$wc_1_lv)="Waist Circumference 1"
label(data$wc_2_lv)="Waist Circumference 2"
label(data$wc_3_lv)="Waist Circumference 3"
label(data$wc_avg_lv)="Waist Circumference Average"
label(data$study_visit_form_complete)="Complete?"
label(data$date_of_study_visit_mm)="Date of Study Visit"
label(data$study_visit_number_mm)="Study Visit Number"
label(data$insulin_sensitivity)="Insulin Sensitivity (Si)"
label(data$insulin_sensitivity_fsd)="Insulin Sensitivity (Si) FSD"
label(data$insulin_secretion_mm)="Insulin Secretion (AIRg)"
label(data$disposition_index)="Disposition Index (DI)"
label(data$disposition_index_fsd)="Disposition Index (DI) FSD"
label(data$sg_mm)="Glucose Effectiveness (Sg)"
label(data$gezi_mm)="Glucose Effectiveness @ 0 Insulin (GEZI)"
label(data$kg_mm)="Glucose Tolerance (Kg)"
label(data$minimal_model_complete)="Complete?"
label(data$dov_dexa)="Date of DEXA"
label(data$study_visit_number_dexa)="Study Visit Number"
label(data$lean_mass_dexa)="Lean Mass"
label(data$fat_mass_dexa)="Fat Mass"
label(data$bone_density_dexa)="Bone Mineral Density"
label(data$bone_content_dexa)="Bone Mineral Content"
label(data$fat_percentage_dexa)="% Fat"
label(data$lean_percentage_dexa)="% Lean"
label(data$dexa_complete)="Complete?"
label(data$date_pdpar)="Date of Study Visit"
label(data$sv_num_pdpar)="Study Visit Number"
label(data$mets)="METS"
label(data$pdpar_complete)="Complete?"
label(data$date_urine)="Date of Study Visit "
label(data$type_visit_urine)="Type of Visit"
label(data$visit_number)="IVGTT/Exam Visit Number"
label(data$ur_cr)="Urine creatinine"
label(data$e1c_cr)="E1c/Cr"
label(data$pdg_cr)="Pdg/Cr"
label(data$lh_cr)="LH/Cr"
label(data$fsh_cr)="FSH/Cr"
label(data$urine_labs_complete)="Complete?"
#Setting Units

units(data$tcholes_lv)="mg/dL"
units(data$tg_lv)="mg/dL"
units(data$hdl_lv)="mg/dL"
units(data$ldl_lv)="mg/dL"
units(data$adiponect_lv)="ug/mL"
units(data$crp_lv)="mg/dL"
units(data$igf_lv)="ng/mL"
units(data$height_lv)="cm"
units(data$weight_lv)="kilograms"
units(data$bmi_screen_v1)="kilograms/m2"
units(data$bmi_pct_screen_v1)="percent"
units(data$wc_1_lv)="cm"
units(data$wc_2_lv)="cm"
units(data$wc_3_lv)="cm"
units(data$wc_avg_lv)="cm"
units(data$lean_mass_dexa)="g"
units(data$fat_mass_dexa)="g"
units(data$bone_density_dexa)="g/cm3"
units(data$bone_content_dexa)="g"
units(data$fat_percentage_dexa)="%"
units(data$lean_percentage_dexa)="%"

#Setting Factors(will create new variable for factors)
data$redcap_event_name.factor = factor(data$redcap_event_name,levels=c("ivgtt_visit_1_arm_1","exam_visit_1_arm_1","exam_visit_2_arm_1","ivgtt_visit_2_arm_1","exam_visit_3_arm_1","exam_visit_4_arm_1","ivgtt_visit_3_arm_1","exam_visit_5_arm_1","exam_visit_6_arm_1","exam_visit_7_arm_1","exam_visit_8_arm_1","exam_visit_9_arm_1","exam_visit_10_arm_1","exam_visit_11_arm_1","exam_visit_12_arm_1","ivgtt_visit_1_arm_2","med_visit_1_arm_2","exam_visit_2_arm_2","med_visit_3_arm_2","exam_visit_4_arm_2","ivgtt_visit_2_arm_2","med_visit_5_arm_2","exam_visit_6_arm_2","med_visit_7_arm_2","exam_visit_8_arm_2","ivgtt_visit_3_arm_2","med_visit_9_arm_2","exam_visit_10_arm_2","med_visit_11_arm_2","exam_visit_12_arm_2","ivgtt_visit_4_arm_2","med_visit_13_arm_2","exam_visit_14_arm_2","med_visit_15_arm_2","exam_visit_16_arm_2","med_visit_17_arm_2","exam_visit_18_arm_2","med_visit_19_arm_2","exam_visit_20_arm_2"))
data$study_visit_number_svl.factor = factor(data$study_visit_number_svl,levels=c("1","2","3","4","5"))
data$ldl_meth_lv.factor = factor(data$ldl_meth_lv,levels=c("O","1"))
data$study_visit_adult_core_labs_complete.factor = factor(data$study_visit_adult_core_labs_complete,levels=c("0","1","2"))
data$sv_num.factor = factor(data$sv_num,levels=c("1","2","3","4","5"))
data$study_visit_other_labs_complete.factor = factor(data$study_visit_other_labs_complete,levels=c("0","1","2"))
data$study_visit_number.factor = factor(data$study_visit_number,levels=c("1","2","3","4","5"))
data$urine_sv.factor = factor(data$urine_sv,levels=c("1","0"))
data$study_visit_form_complete.factor = factor(data$study_visit_form_complete,levels=c("0","1","2"))
data$study_visit_number_mm.factor = factor(data$study_visit_number_mm,levels=c("1","2","3","4","5"))
data$minimal_model_complete.factor = factor(data$minimal_model_complete,levels=c("0","1","2"))
data$study_visit_number_dexa.factor = factor(data$study_visit_number_dexa,levels=c("1","2","3","4","5"))
data$dexa_complete.factor = factor(data$dexa_complete,levels=c("0","1","2"))
data$sv_num_pdpar.factor = factor(data$sv_num_pdpar,levels=c("1","2","3","4"))
data$pdpar_complete.factor = factor(data$pdpar_complete,levels=c("0","1","2"))
data$type_visit_urine.factor = factor(data$type_visit_urine,levels=c("1","2"))
data$urine_labs_complete.factor = factor(data$urine_labs_complete,levels=c("0","1","2"))

levels(data$redcap_event_name.factor)=c("IVGTT Visit 1 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 1 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 2 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 2 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 3 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 4 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 3 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 5 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 6 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 7 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 8 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 9 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 10 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 11 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 12 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 1 (Arm 2: TREATMENT GROUP)","Med. Visit 1 (Arm 2: TREATMENT GROUP)","Exam Visit 2 (Arm 2: TREATMENT GROUP)","Med Visit 3 (Arm 2: TREATMENT GROUP)","Exam Visit 4 (Arm 2: TREATMENT GROUP)","IVGTT Visit 2 (Arm 2: TREATMENT GROUP)","Med Visit 5 (Arm 2: TREATMENT GROUP)","Exam Visit 6 (Arm 2: TREATMENT GROUP)","Med Visit 7 (Arm 2: TREATMENT GROUP)","Exam Visit 8 (Arm 2: TREATMENT GROUP)","IVGTT Visit 3 (Arm 2: TREATMENT GROUP)","Med. Visit 9 (Arm 2: TREATMENT GROUP)","Exam Visit 10 (Arm 2: TREATMENT GROUP)","Med. Visit 11 (Arm 2: TREATMENT GROUP)","Exam Visit 12 (Arm 2: TREATMENT GROUP)","IVGTT Visit 4 (Arm 2: TREATMENT GROUP)","Med. Visit 13 (Arm 2: TREATMENT GROUP)","Exam Visit 14 (Arm 2: TREATMENT GROUP)","Med. Visit 15 (Arm 2: TREATMENT GROUP)","Exam Visit 16 (Arm 2: TREATMENT GROUP)","Med. Visit 17 (Arm 2: TREATMENT GROUP)","Exam Visit 18 (Arm 2: TREATMENT GROUP)","Med. Visit 19 (Arm 2: TREATMENT GROUP)","Exam Visit 20 (Arm 2: TREATMENT GROUP)")
levels(data$study_visit_number_svl.factor)=c("1","2","3","4","5")
levels(data$ldl_meth_lv.factor)=c("Measured","Calculated")
levels(data$study_visit_adult_core_labs_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$sv_num.factor)=c("1","2","3","4","5")
levels(data$study_visit_other_labs_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$study_visit_number.factor)=c("1","2","3","4","5")
levels(data$urine_sv.factor)=c("Yes","No")
levels(data$study_visit_form_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$study_visit_number_mm.factor)=c("1","2","3","4","5")
levels(data$minimal_model_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$study_visit_number_dexa.factor)=c("1","2","3","4","5")
levels(data$dexa_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$sv_num_pdpar.factor)=c("1","2","3","4")
levels(data$pdpar_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$type_visit_urine.factor)=c("IVGTT","Exam")
levels(data$urine_labs_complete.factor)=c("Incomplete","Unverified","Complete")

# add randomization dat
rand <- read.csv("H:\\Endocrinology\\Kelsey\\Kesley HIP study ISS puberty\\Final dataset cleaning\\Randomization\\HIP Study Randomization Schema.csv")
names(rand)[1] <- "hip_id_screen"
final <- merge(data,rand,by=c("hip_id_screen"),all.x=TRUE, all.y=TRUE)
check <- c("hip_id_screen","Randomization.Group","redcap_event_name")
final[check]
dim(final)

#add leading zeroes to some IDs
final$hip_id_screen <- as.character(final$hip_id_screen)
final$hip_id_screen[nchar(final$hip_id_screen)==2] <- paste0("0", final$hip_id_screen[nchar(final$hip_id_screen)==2])
final$hip_id_screen[nchar(final$hip_id_screen)==1] <- paste0("00", final$hip_id_screen[nchar(final$hip_id_screen)==1])

# delete two participants who did not get randomized
final <- final[!(final$hip_id_screen=="078-T" | final$hip_id_screen=="140-T" ),]
dim(final)

# remove those on metformin
nomet <- final
nomet <- nomet[!(final$Randomization.Group=="1") | is.na(final$Randomization.Group),]

# delete visits with missing insulin sensitivity, secretion, and DI
nomet <- nomet[!is.na(nomet$insulin_secretion_mm) & !is.na(nomet$insulin_sensitivity) & !is.na(nomet$disposition_index),]

# need to delete those with only one visit
foranalysis <- nomet[nomet$hip_id_screen %in% names(table(nomet$hip_id_screen))[table(nomet$hip_id_screen) > 1],]

# add in baseline Tanner, race, sex
#Read Data
data=read.csv('baseline vars.csv')
#Setting Labels

# need to add leading zeroes to IDs in the baseline file
#add leading zeroes to some IDs
data$hip_id_screen <- as.character(data$hip_id_screen)
data$hip_id_screen[nchar(data$hip_id_screen)==2] <- paste0("0", data$hip_id_screen[nchar(data$hip_id_screen)==2])
data$hip_id_screen[nchar(data$hip_id_screen)==1] <- paste0("00", data$hip_id_screen[nchar(data$hip_id_screen)==1])


label(data$hip_id_screen)="HIP ID"
label(data$redcap_event_name)="Event Name"
label(data$sex)="Gender"
label(data$ethnicity)="Ethnicity"
label(data$race___1)="Race (choice=1 Native American)"
label(data$race___2)="Race (choice=2 Asian)"
label(data$race___3)="Race (choice=3 Black)"
label(data$race___4)="Race (choice=4 White)"
label(data$race___5)="Race (choice=5 Hispanic)"
label(data$race___6)="Race (choice=6 Other)"
label(data$if_other_race_ic)="If other race, please specify"
label(data$tanner_breast_exam)="Tanner Stage, Breast"
label(data$tanner_hair_exam)="Tanner Stage, Pubic Hair"
label(data$tanner_genetalia_exam)="Tanner Stage, Genetalia"
label(data$tanner_testis_vol_exam)="Tanner Stage, Testicular Volume    T1   < = 3ml    T2      >3 to 8 ml    T3      >=8 to 12 ml    T4     >=12 to < =15 ml    T5   >15 ml"
label(data$bmi_cat_final)="Final BMI Category  "
#Setting Units


#Setting Factors(will create new variable for factors)
data$redcap_event_name.factor = factor(data$redcap_event_name,levels=c("ivgtt_visit_1_arm_1","exam_visit_1_arm_1","exam_visit_2_arm_1","ivgtt_visit_2_arm_1","exam_visit_3_arm_1","exam_visit_4_arm_1","ivgtt_visit_3_arm_1","exam_visit_5_arm_1","exam_visit_6_arm_1","exam_visit_7_arm_1","exam_visit_8_arm_1","exam_visit_9_arm_1","exam_visit_10_arm_1","exam_visit_11_arm_1","exam_visit_12_arm_1","ivgtt_visit_1_arm_2","med_visit_1_arm_2","exam_visit_2_arm_2","med_visit_3_arm_2","exam_visit_4_arm_2","ivgtt_visit_2_arm_2","med_visit_5_arm_2","exam_visit_6_arm_2","med_visit_7_arm_2","exam_visit_8_arm_2","ivgtt_visit_3_arm_2","med_visit_9_arm_2","exam_visit_10_arm_2","med_visit_11_arm_2","exam_visit_12_arm_2","ivgtt_visit_4_arm_2","med_visit_13_arm_2","exam_visit_14_arm_2","med_visit_15_arm_2","exam_visit_16_arm_2","med_visit_17_arm_2","exam_visit_18_arm_2","med_visit_19_arm_2","exam_visit_20_arm_2"))
data$sex.factor = factor(data$sex,levels=c("0","1"))
data$ethnicity.factor = factor(data$ethnicity,levels=c("0","1","2"))
data$race___1.factor = factor(data$race___1,levels=c("0","1"))
data$race___2.factor = factor(data$race___2,levels=c("0","1"))
data$race___3.factor = factor(data$race___3,levels=c("0","1"))
data$race___4.factor = factor(data$race___4,levels=c("0","1"))
data$race___5.factor = factor(data$race___5,levels=c("0","1"))
data$race___6.factor = factor(data$race___6,levels=c("0","1"))
data$tanner_breast_exam.factor = factor(data$tanner_breast_exam,levels=c("1","2","3","4","5"))
data$tanner_hair_exam.factor = factor(data$tanner_hair_exam,levels=c("1","2","3","4","5"))
data$tanner_genetalia_exam.factor = factor(data$tanner_genetalia_exam,levels=c("1","2","3","4","5"))
data$tanner_testis_vol_exam.factor = factor(data$tanner_testis_vol_exam,levels=c("1","2","3","4","5"))
data$bmi_cat_final.factor = factor(data$bmi_cat_final,levels=c("1","2"))

levels(data$redcap_event_name.factor)=c("IVGTT Visit 1 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 1 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 2 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 2 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 3 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 4 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 3 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 5 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 6 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 7 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 8 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 9 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 10 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 11 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 12 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 1 (Arm 2: TREATMENT GROUP)","Med. Visit 1 (Arm 2: TREATMENT GROUP)","Exam Visit 2 (Arm 2: TREATMENT GROUP)","Med Visit 3 (Arm 2: TREATMENT GROUP)","Exam Visit 4 (Arm 2: TREATMENT GROUP)","IVGTT Visit 2 (Arm 2: TREATMENT GROUP)","Med Visit 5 (Arm 2: TREATMENT GROUP)","Exam Visit 6 (Arm 2: TREATMENT GROUP)","Med Visit 7 (Arm 2: TREATMENT GROUP)","Exam Visit 8 (Arm 2: TREATMENT GROUP)","IVGTT Visit 3 (Arm 2: TREATMENT GROUP)","Med. Visit 9 (Arm 2: TREATMENT GROUP)","Exam Visit 10 (Arm 2: TREATMENT GROUP)","Med. Visit 11 (Arm 2: TREATMENT GROUP)","Exam Visit 12 (Arm 2: TREATMENT GROUP)","IVGTT Visit 4 (Arm 2: TREATMENT GROUP)","Med. Visit 13 (Arm 2: TREATMENT GROUP)","Exam Visit 14 (Arm 2: TREATMENT GROUP)","Med. Visit 15 (Arm 2: TREATMENT GROUP)","Exam Visit 16 (Arm 2: TREATMENT GROUP)","Med. Visit 17 (Arm 2: TREATMENT GROUP)","Exam Visit 18 (Arm 2: TREATMENT GROUP)","Med. Visit 19 (Arm 2: TREATMENT GROUP)","Exam Visit 20 (Arm 2: TREATMENT GROUP)")
levels(data$sex.factor)=c("0 Female","1 Male")
levels(data$ethnicity.factor)=c("Hispanic or Latino","NOT Hispanic or Latino","Unknown / Not Reported")
levels(data$race___1.factor)=c("Unchecked","Checked")
levels(data$race___2.factor)=c("Unchecked","Checked")
levels(data$race___3.factor)=c("Unchecked","Checked")
levels(data$race___4.factor)=c("Unchecked","Checked")
levels(data$race___5.factor)=c("Unchecked","Checked")
levels(data$race___6.factor)=c("Unchecked","Checked")
levels(data$tanner_breast_exam.factor)=c("Tanner 1","Tanner 2","Tanner 3","Tanner 4","Tanner 5")
levels(data$tanner_hair_exam.factor)=c("Tanner 1","Tanner 2","Tanner 3","Tanner 4","Tanner 5")
levels(data$tanner_genetalia_exam.factor)=c("Tanner 1","Tanner 2","Tanner 3","Tanner 4","Tanner 5")
levels(data$tanner_testis_vol_exam.factor)=c("Tanner 1","Tanner 2","Tanner 3","Tanner 4","Tanner 5")
levels(data$bmi_cat_final.factor)=c("Normal Weight","Obese")

data <- data[,-2]

dim(foranalysis)
foranalysis <- merge(foranalysis,data,by=c("hip_id_screen"),all.x=TRUE, all.y=FALSE)
dim(foranalysis)

# combine race/ethnicity
foranalysis$sumrace <- foranalysis$race___1  + foranalysis$race___2 + foranalysis$race___3 + foranalysis$race___4 + foranalysis$race___5
foranalysis$race_eth[foranalysis$ethnicity==0] <- "Hispanic"
foranalysis$race_eth[(foranalysis$ethnicity==1 | foranalysis$ethnicity==2) & foranalysis$race___1==1] <- "Native American"
foranalysis$race_eth[(foranalysis$ethnicity==1 | foranalysis$ethnicity==2) & foranalysis$race___2==1] <- "Asian"
foranalysis$race_eth[(foranalysis$ethnicity==1 | foranalysis$ethnicity==2) & foranalysis$race___3==1] <- "Black"
foranalysis$race_eth[(foranalysis$ethnicity==1 | foranalysis$ethnicity==2) & foranalysis$race___4==1] <- "NHW"
foranalysis$race_eth[(foranalysis$ethnicity==1 | foranalysis$ethnicity==2) & foranalysis$race___5==1] <- "Hispanic"
foranalysis$race_eth[is.na(foranalysis$race_eth)] <- "Other"
racevars <- c("sumrace","ethnicity.factor","race___1","race___2","race___3","race___4","race___5","race_eth")

# combine Tanner
foranalysis$tanner[foranalysis$sex==0] <- foranalysis$tanner_breast_exam.factor[foranalysis$sex==0]
foranalysis$tanner[foranalysis$sex==1] <- foranalysis$tanner_testis_vol_exam.factor[foranalysis$sex==1]
tannervars <- c("sex.factor","tanner_breast_exam.factor","tanner_testis_vol_exam.factor","tanner")

# add famhx
data=read.csv('famhx.csv')
#Setting Labels

label(data$hip_id_screen)="HIP ID"
label(data$redcap_event_name)="Event Name"
label(data$mother_history_screen)="Is mothers history available?"
label(data$type_2_db_mother_screen)="Type 2 Diabetes"
label(data$father_history_screen)="Is fathers history available?"
label(data$type_2_db_father_screen)="Type 2 Diabetes"
label(data$sib_no)="How many siblings does participant have?"
label(data$sib_1_history_sv)="Is sibling 1 history available?"
label(data$sb1_full_or_half_sv)="Is sibling 1 half or full sibling? "
label(data$type2_db_sb1_screen)="Type 2 Diabetes"
label(data$sb2_history_screen)="Is sibling 2 history available?"
label(data$sb2_full_or_half_screen)="Is sibling 2 half or full sibling?"
label(data$type2_db_sb2_screen)="Type 2 Diabetes"
label(data$sb3_history_screen)="Is the sibling 3 history available?"
label(data$full_sb3_half_or_full_sv)="Is sibling 3 half or full sibling?"
label(data$type2_db_sb3_screen)="Type 2 Diabetes"
label(data$sb4_history_screen)="Is sibling 4 history available?"
label(data$sb4_half_or_full_screen)="Is sibling 4 half or full sibling?"
label(data$type2_db_sb4_screen)="Type 2 Diabetes"
label(data$sb5_history_screen)="Is sibling 5 history available?"
label(data$sb5_half_or_full_screening)="Is sibling 5 half or full sibling?"
label(data$type2_db_sb5_screen)="Type 2 Diabetes"
#Setting Units


#Setting Factors(will create new variable for factors)
data$redcap_event_name.factor = factor(data$redcap_event_name,levels=c("ivgtt_visit_1_arm_1","exam_visit_1_arm_1","exam_visit_2_arm_1","ivgtt_visit_2_arm_1","exam_visit_3_arm_1","exam_visit_4_arm_1","ivgtt_visit_3_arm_1","exam_visit_5_arm_1","exam_visit_6_arm_1","exam_visit_7_arm_1","exam_visit_8_arm_1","exam_visit_9_arm_1","exam_visit_10_arm_1","exam_visit_11_arm_1","exam_visit_12_arm_1","ivgtt_visit_1_arm_2","med_visit_1_arm_2","exam_visit_2_arm_2","med_visit_3_arm_2","exam_visit_4_arm_2","ivgtt_visit_2_arm_2","med_visit_5_arm_2","exam_visit_6_arm_2","med_visit_7_arm_2","exam_visit_8_arm_2","ivgtt_visit_3_arm_2","med_visit_9_arm_2","exam_visit_10_arm_2","med_visit_11_arm_2","exam_visit_12_arm_2","ivgtt_visit_4_arm_2","med_visit_13_arm_2","exam_visit_14_arm_2","med_visit_15_arm_2","exam_visit_16_arm_2","med_visit_17_arm_2","exam_visit_18_arm_2","med_visit_19_arm_2","exam_visit_20_arm_2"))
data$mother_history_screen.factor = factor(data$mother_history_screen,levels=c("1","0"))
data$type_2_db_mother_screen.factor = factor(data$type_2_db_mother_screen,levels=c("1","0"))
data$father_history_screen.factor = factor(data$father_history_screen,levels=c("1","0"))
data$type_2_db_father_screen.factor = factor(data$type_2_db_father_screen,levels=c("1","0"))
data$sib_1_history_sv.factor = factor(data$sib_1_history_sv,levels=c("1","0"))
data$sb1_full_or_half_sv.factor = factor(data$sb1_full_or_half_sv,levels=c("1","2"))
data$type2_db_sb1_screen.factor = factor(data$type2_db_sb1_screen,levels=c("1","0"))
data$sb2_history_screen.factor = factor(data$sb2_history_screen,levels=c("1","0"))
data$sb2_full_or_half_screen.factor = factor(data$sb2_full_or_half_screen,levels=c("1","2"))
data$type2_db_sb2_screen.factor = factor(data$type2_db_sb2_screen,levels=c("1","0"))
data$sb3_history_screen.factor = factor(data$sb3_history_screen,levels=c("1","0"))
data$full_sb3_half_or_full_sv.factor = factor(data$full_sb3_half_or_full_sv,levels=c("1","2"))
data$type2_db_sb3_screen.factor = factor(data$type2_db_sb3_screen,levels=c("1","0"))
data$sb4_history_screen.factor = factor(data$sb4_history_screen,levels=c("1","0"))
data$sb4_half_or_full_screen.factor = factor(data$sb4_half_or_full_screen,levels=c("1","2"))
data$type2_db_sb4_screen.factor = factor(data$type2_db_sb4_screen,levels=c("1","0"))
data$sb5_history_screen.factor = factor(data$sb5_history_screen,levels=c("1","0"))
data$sb5_half_or_full_screening.factor = factor(data$sb5_half_or_full_screening,levels=c("1","2"))
data$type2_db_sb5_screen.factor = factor(data$type2_db_sb5_screen,levels=c("1","0"))

levels(data$redcap_event_name.factor)=c("IVGTT Visit 1 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 1 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 2 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 2 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 3 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 4 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 3 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 5 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 6 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 7 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 8 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 9 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 10 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 11 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 12 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 1 (Arm 2: TREATMENT GROUP)","Med. Visit 1 (Arm 2: TREATMENT GROUP)","Exam Visit 2 (Arm 2: TREATMENT GROUP)","Med Visit 3 (Arm 2: TREATMENT GROUP)","Exam Visit 4 (Arm 2: TREATMENT GROUP)","IVGTT Visit 2 (Arm 2: TREATMENT GROUP)","Med Visit 5 (Arm 2: TREATMENT GROUP)","Exam Visit 6 (Arm 2: TREATMENT GROUP)","Med Visit 7 (Arm 2: TREATMENT GROUP)","Exam Visit 8 (Arm 2: TREATMENT GROUP)","IVGTT Visit 3 (Arm 2: TREATMENT GROUP)","Med. Visit 9 (Arm 2: TREATMENT GROUP)","Exam Visit 10 (Arm 2: TREATMENT GROUP)","Med. Visit 11 (Arm 2: TREATMENT GROUP)","Exam Visit 12 (Arm 2: TREATMENT GROUP)","IVGTT Visit 4 (Arm 2: TREATMENT GROUP)","Med. Visit 13 (Arm 2: TREATMENT GROUP)","Exam Visit 14 (Arm 2: TREATMENT GROUP)","Med. Visit 15 (Arm 2: TREATMENT GROUP)","Exam Visit 16 (Arm 2: TREATMENT GROUP)","Med. Visit 17 (Arm 2: TREATMENT GROUP)","Exam Visit 18 (Arm 2: TREATMENT GROUP)","Med. Visit 19 (Arm 2: TREATMENT GROUP)","Exam Visit 20 (Arm 2: TREATMENT GROUP)")
levels(data$mother_history_screen.factor)=c("Yes","No")
levels(data$type_2_db_mother_screen.factor)=c("Yes","No")
levels(data$father_history_screen.factor)=c("Yes","No")
levels(data$type_2_db_father_screen.factor)=c("Yes","No")
levels(data$sib_1_history_sv.factor)=c("Yes","No")
levels(data$sb1_full_or_half_sv.factor)=c("0 Half","1 Full")
levels(data$type2_db_sb1_screen.factor)=c("Yes","No")
levels(data$sb2_history_screen.factor)=c("Yes","No")
levels(data$sb2_full_or_half_screen.factor)=c("0 Half","1 Full")
levels(data$type2_db_sb2_screen.factor)=c("Yes","No")
levels(data$sb3_history_screen.factor)=c("Yes","No")
levels(data$full_sb3_half_or_full_sv.factor)=c("0 Half","1 Full")
levels(data$type2_db_sb3_screen.factor)=c("Yes","No")
levels(data$sb4_history_screen.factor)=c("Yes","No")
levels(data$sb4_half_or_full_screen.factor)=c("0 Half","1 Full")
levels(data$type2_db_sb4_screen.factor)=c("Yes","No")
levels(data$sb5_history_screen.factor)=c("Yes","No")
levels(data$sb5_half_or_full_screening.factor)=c("0 Half","1 Full")
levels(data$type2_db_sb5_screen.factor)=c("Yes","No")

# create variable for family history of T2D
# per Megan, if one parent hx available, code per parent 
# first parent history
data$parhx[data$mother_history_screen==0 & data$father_history_screen==0] <- NA
data$parhx[(data$mother_history_screen==1 | data$father_history_screen==1) &
                    (data$type_2_db_mother_screen==1 | data$type_2_db_father_screen==1)] <- 1
data$parhx[(data$mother_history_screen==1 & data$father_history_screen==0) &
                    (data$type_2_db_mother_screen==1) ] <- 1
data$parhx[(data$mother_history_screen==0 & data$father_history_screen==1) &
                    (data$type_2_db_father_screen==1) ] <- 1
data$parhx[(data$mother_history_screen==1 | data$father_history_screen==1) &
                    (data$type_2_db_mother_screen==0 & data$type_2_db_father_screen==0)] <- 0
data$parhx[(data$mother_history_screen==1 & data$father_history_screen==0) &
                    (data$type_2_db_mother_screen==0) ] <- 0
data$parhx[(data$mother_history_screen==0 & data$father_history_screen==1) &
                    (data$type_2_db_father_screen==0) ] <- 0
# then sib history
# per Megan, if one sib hx available, code per sib
# no siblings
data$sibhx[data$sib_no==0] <- 0
# at least one sibling and all missing
data$sibhx[data$sib_no>=1 & data$sib_1_history_sv==0 & 
                    (data$sb2_history_screen==0 | is.na(data$sb2_history_screen)) &
                    (data$sb3_history_screen==0 | is.na(data$sb3_history_screen)) &
                    (data$sb4_history_screen==0 | is.na(data$sb4_history_screen)) &
                    (data$sb5_history_screen==0 | is.na(data$sb5_history_screen)) 
                    ] <- NA
# at least one sibling and all negative
data$sibhx[data$sib_no>=1 & 
                    (data$type2_db_sb1_screen==0 | is.na(data$type2_db_sb1_screen)) &
                  (data$type2_db_sb2_screen==0 | is.na(data$type2_db_sb2_screen)) &
                     (data$type2_db_sb3_screen==0 | is.na(data$type2_db_sb3_screen)) &
                        (data$type2_db_sb4_screen==0 | is.na(data$type2_db_sb4_screen)) &
                           (data$type2_db_sb5_screen==0 | is.na(data$type2_db_sb5_screen)) 
                  ] <- 0
# at least one sibling and at least one positive
data$sibhx[data$sib_no>=1 & 
                    ( data$type2_db_sb1_screen==1 |
                    data$type2_db_sb2_screen==1 |
                    data$type2_db_sb3_screen==1 |
                    data$type2_db_sb4_screen==1 |
                    data$type2_db_sb5_screen==1) 
                  ] <- 1
# now combine parent and sib hx
data$famhxt2d[is.na(data$parhx) & is.na(data$sibhx)] <- NA
data$famhxt2d[data$parhx==1 | data$sibhx==1] <- 1
data$famhxt2d[data$parhx==0 & data$sibhx==0] <- 0
data$famhxt2d[data$parhx==0 & is.na(data$sibhx)] <- 0
data$famhxt2d[is.na(data$parhx) & data$sibhx==0] <- 0
famhx <- data[data$redcap_event_name.factor=="IVGTT Visit 1 (Arm 1: NOT TREATMENT GROUP)" |
                data$redcap_event_name.factor=="IVGTT Visit 1 (Arm 2: TREATMENT GROUP)",
              c("hip_id_screen","famhxt2d")]

# write final long dataset
write.csv(foranalysis,file="for_analysis.csv")

# make a wide dataset
keep <- c("hip_id_screen","date_of_study_visit","sv_num","insulin_sensitivity","insulin_secretion_mm",
          "disposition_index","fat_percentage_dexa","bmi_cat_final.factor","sex.factor","race_eth","tanner")
temp <- foranalysis[keep]
wide <- reshape(temp, timevar="sv_num", idvar= c("hip_id_screen","bmi_cat_final.factor","sex.factor","race_eth","tanner"),direction="wide")

# calculate deltas
# use visit 3 as last visit
wide$delta_si[nchar(wide$hip_id_screen)==5] <- wide$insulin_sensitivity.3[nchar(wide$hip_id_screen)==5] - wide$insulin_sensitivity.1[nchar(wide$hip_id_screen)==5]
wide$delta_si[nchar(wide$hip_id_screen)!=5] <- wide$insulin_sensitivity.3[nchar(wide$hip_id_screen)!=5] - wide$insulin_sensitivity.1[nchar(wide$hip_id_screen)!=5]
wide$delta_inssec[nchar(wide$hip_id_screen)==5] <- wide$insulin_secretion_mm.3[nchar(wide$hip_id_screen)==5] - wide$insulin_secretion_mm.1[nchar(wide$hip_id_screen)==5]
wide$delta_inssec[nchar(wide$hip_id_screen)!=5] <- wide$insulin_secretion_mm.3[nchar(wide$hip_id_screen)!=5] - wide$insulin_secretion_mm.1[nchar(wide$hip_id_screen)!=5]
wide$delta_di[nchar(wide$hip_id_screen)==5] <- wide$disposition_index.3[nchar(wide$hip_id_screen)==5] - wide$disposition_index.1[nchar(wide$hip_id_screen)==5]
wide$delta_di[nchar(wide$hip_id_screen)!=5] <- wide$disposition_index.3[nchar(wide$hip_id_screen)!=5] - wide$disposition_index.1[nchar(wide$hip_id_screen)!=5]
wide$delta_bodyfat[nchar(wide$hip_id_screen)==5] <- wide$fat_percentage_dexa.3[nchar(wide$hip_id_screen)==5] - wide$fat_percentage_dexa.1[nchar(wide$hip_id_screen)==5]
wide$delta_bodyfat[nchar(wide$hip_id_screen)!=5] <- wide$fat_percentage_dexa.3[nchar(wide$hip_id_screen)!=5] - wide$fat_percentage_dexa.1[nchar(wide$hip_id_screen)!=5]

# delete people with missing delta Si, insulin secretion or DI
dim(wide)
wide <- wide[!is.na(wide$delta_di) & !is.na(wide$delta_inssec) & !is.na(wide$delta_di),]
dim(wide)

# set labels
label(wide$delta_di)="Change in DI"
label(wide$delta_inssec)="Change in insulin secretion"
label(wide$delta_si)="Change in insulin sensitivity"
label(wide$delta_bodyfat)="Change in percent body fat"

# add baseline predictors to wide dataset
# need to make sure only visit 1
base <- c("hip_id_screen","redcap_event_name","lept","adiponect_lv","fat_percentage_dexa",
          "crp_lv","tg_lv","hdl_lv", "mets","dhea_s","igf_lv","gluc_0_iv","hba")
base <- foranalysis[base]
base <- base[base$redcap_event_name=="ivgtt_visit_1_arm_1" | base$redcap_event_name=="ivgtt_visit_1_arm_2",]
#  now get rid of redcap event name
base <- base[,-2]
# calculate tg:HDL
base$tg_hdl <- base$tg_lv/base$hdl_lv   

# merge family hx with base
dim(base)
base <- merge(base,famhx,by=c("hip_id_screen"),all.x=TRUE, all.y=FALSE)
dim(base)
colnames(base) <- c("hip_id_screen","lept0","adipo0","fat_perc0","crp0","tg0","hdl0","mets0","dheas0","igf0","gluc00","hba1c0","tg_hdl0","famhxt2d0")

# merge base with wide
dim(wide)
wide <- merge(wide,base,by="hip_id_screen",all.x=TRUE, all.y=FALSE)
dim(wide)

wide$lept0 <- as.numeric(wide$lept0)
wide$crp0 <- as.numeric(wide$crp0)
label(wide$lept0)="Baseline leptin"
label(wide$adipo0)="Baseline adiponectin"
label(wide$fat_perc0)="Baseline % fat"
label(wide$crp0)="Baseline CRP"
label(wide$tg0)="Baseline TG"
label(wide$hdl0)="Baseline HDL"
label(wide$mets0)="Baseline mets"
label(wide$dheas0)="Baseline DHEAS"
label(wide$igf0)="Baseline IGF"
label(wide$gluc00)="Baseline fasting glucose"
label(wide$hba1c0)="Baseline HbA1c"
label(wide$tg_hdl0)="Baseline HDL"
label(wide$famhxt2d0)="Family history T2D"

