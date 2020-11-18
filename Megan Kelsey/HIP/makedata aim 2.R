#setwd("H:\\Endocrinology\\Kelsey\\Kesley HIP study ISS puberty\\Ergui poster\\Data for models")
#setwd("C:\\Temp\\Nashville trip\\Kesley HIP study ISS puberty\\Ergui poster\\Data for models")
#setwd("C:\\Temp\\Nashville trip\\Kesley HIP study ISS puberty\\Final dataset cleaning")
setwd("H:\\Endocrinology\\Kelsey\\Kesley HIP study ISS puberty\\Final dataset cleaning")

library(tidyr)
library(reshape2)
library(dplyr)
library(pracma)

#Clear existing data and graphics
rm(list=ls())
graphics.off()
#Load Hmisc library
library(Hmisc)
#Read Data
data=read.csv('IVGTT visits clean 082818.csv',na.strings = c("-99"," ",""))
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
#label(data$study_visit_adult_core_labs_complete)="Complete?"
label(data$dt_of_study_visit)="Date of Study Visit "
label(data$sv_num)="Study Visit Number"
label(data$ast)="AST"
label(data$alt)="ALT"
label(data$crp_lv)="CRP"
label(data$igf_lv)="IGF-1"
label(data$hba)="HbA1c"
#label(data$study_visit_other_labs_complete)="Complete?"
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
# data$redcap_event_name.factor = factor(data$redcap_event_name,levels=c("ivgtt_visit_1_arm_1","exam_visit_1_arm_1","exam_visit_2_arm_1","ivgtt_visit_2_arm_1","exam_visit_3_arm_1","exam_visit_4_arm_1","ivgtt_visit_3_arm_1","exam_visit_5_arm_1","exam_visit_6_arm_1","exam_visit_7_arm_1","exam_visit_8_arm_1","exam_visit_9_arm_1","exam_visit_10_arm_1","exam_visit_11_arm_1","exam_visit_12_arm_1","ivgtt_visit_1_arm_2","med_visit_1_arm_2","exam_visit_2_arm_2","med_visit_3_arm_2","exam_visit_4_arm_2","ivgtt_visit_2_arm_2","med_visit_5_arm_2","exam_visit_6_arm_2","med_visit_7_arm_2","exam_visit_8_arm_2","ivgtt_visit_3_arm_2","med_visit_9_arm_2","exam_visit_10_arm_2","med_visit_11_arm_2","exam_visit_12_arm_2","ivgtt_visit_4_arm_2","med_visit_13_arm_2","exam_visit_14_arm_2","med_visit_15_arm_2","exam_visit_16_arm_2","med_visit_17_arm_2","exam_visit_18_arm_2","med_visit_19_arm_2","exam_visit_20_arm_2"))
# data$study_visit_number_svl.factor = factor(data$study_visit_number_svl,levels=c("1","2","3","4","5"))
# data$ldl_meth_lv.factor = factor(data$ldl_meth_lv,levels=c("O","1"))
# data$study_visit_adult_core_labs_complete.factor = factor(data$study_visit_adult_core_labs_complete,levels=c("0","1","2"))
# data$sv_num.factor = factor(data$sv_num,levels=c("1","2","3","4","5"))
# data$study_visit_other_labs_complete.factor = factor(data$study_visit_other_labs_complete,levels=c("0","1","2"))
# data$study_visit_number.factor = factor(data$study_visit_number,levels=c("1","2","3","4","5"))
# data$urine_sv.factor = factor(data$urine_sv,levels=c("1","0"))
# data$study_visit_form_complete.factor = factor(data$study_visit_form_complete,levels=c("0","1","2"))
# data$study_visit_number_mm.factor = factor(data$study_visit_number_mm,levels=c("1","2","3","4","5"))
# data$minimal_model_complete.factor = factor(data$minimal_model_complete,levels=c("0","1","2"))
# data$study_visit_number_dexa.factor = factor(data$study_visit_number_dexa,levels=c("1","2","3","4","5"))
# data$dexa_complete.factor = factor(data$dexa_complete,levels=c("0","1","2"))
# data$sv_num_pdpar.factor = factor(data$sv_num_pdpar,levels=c("1","2","3","4"))
# data$pdpar_complete.factor = factor(data$pdpar_complete,levels=c("0","1","2"))
# data$type_visit_urine.factor = factor(data$type_visit_urine,levels=c("1","2"))
# data$urine_labs_complete.factor = factor(data$urine_labs_complete,levels=c("0","1","2"))
# 
# levels(data$redcap_event_name.factor)=c("IVGTT Visit 1 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 1 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 2 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 2 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 3 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 4 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 3 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 5 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 6 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 7 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 8 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 9 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 10 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 11 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 12 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 1 (Arm 2: TREATMENT GROUP)","Med. Visit 1 (Arm 2: TREATMENT GROUP)","Exam Visit 2 (Arm 2: TREATMENT GROUP)","Med Visit 3 (Arm 2: TREATMENT GROUP)","Exam Visit 4 (Arm 2: TREATMENT GROUP)","IVGTT Visit 2 (Arm 2: TREATMENT GROUP)","Med Visit 5 (Arm 2: TREATMENT GROUP)","Exam Visit 6 (Arm 2: TREATMENT GROUP)","Med Visit 7 (Arm 2: TREATMENT GROUP)","Exam Visit 8 (Arm 2: TREATMENT GROUP)","IVGTT Visit 3 (Arm 2: TREATMENT GROUP)","Med. Visit 9 (Arm 2: TREATMENT GROUP)","Exam Visit 10 (Arm 2: TREATMENT GROUP)","Med. Visit 11 (Arm 2: TREATMENT GROUP)","Exam Visit 12 (Arm 2: TREATMENT GROUP)","IVGTT Visit 4 (Arm 2: TREATMENT GROUP)","Med. Visit 13 (Arm 2: TREATMENT GROUP)","Exam Visit 14 (Arm 2: TREATMENT GROUP)","Med. Visit 15 (Arm 2: TREATMENT GROUP)","Exam Visit 16 (Arm 2: TREATMENT GROUP)","Med. Visit 17 (Arm 2: TREATMENT GROUP)","Exam Visit 18 (Arm 2: TREATMENT GROUP)","Med. Visit 19 (Arm 2: TREATMENT GROUP)","Exam Visit 20 (Arm 2: TREATMENT GROUP)")
# levels(data$study_visit_number_svl.factor)=c("1","2","3","4","5")
# levels(data$ldl_meth_lv.factor)=c("Measured","Calculated")
# levels(data$study_visit_adult_core_labs_complete.factor)=c("Incomplete","Unverified","Complete")
# levels(data$sv_num.factor)=c("1","2","3","4","5")
# levels(data$study_visit_other_labs_complete.factor)=c("Incomplete","Unverified","Complete")
# levels(data$study_visit_number.factor)=c("1","2","3","4","5")
# levels(data$urine_sv.factor)=c("Yes","No")
# levels(data$study_visit_form_complete.factor)=c("Incomplete","Unverified","Complete")
# levels(data$study_visit_number_mm.factor)=c("1","2","3","4","5")
# levels(data$minimal_model_complete.factor)=c("Incomplete","Unverified","Complete")
# levels(data$study_visit_number_dexa.factor)=c("1","2","3","4","5")
# levels(data$dexa_complete.factor)=c("Incomplete","Unverified","Complete")
# levels(data$sv_num_pdpar.factor)=c("1","2","3","4")
# levels(data$pdpar_complete.factor)=c("Incomplete","Unverified","Complete")
# levels(data$type_visit_urine.factor)=c("IVGTT","Exam")
# levels(data$urine_labs_complete.factor)=c("Incomplete","Unverified","Complete")
# 
# add randomization dat
#rand <- read.csv("H:\\Endocrinology\\Kelsey\\Kesley HIP study ISS puberty\\Ergui poster\\Randomization\\HIP Study Randomization Schema.csv")
#rand <- read.csv("C:\\Temp\\Nashville trip\\Kesley HIP study ISS puberty\\Ergui poster\\Randomization\\HIP Study Randomization Schema.csv")
rand <- read.csv("H:\\Endocrinology\\Kelsey\\Kesley HIP study ISS puberty\\Final dataset cleaning\\Randomization\\HIP Study Randomization Schema.csv")

names(rand)[1] <- "hip_id_screen"
final <- merge(data,rand,by=c("hip_id_screen"),all.x=TRUE, all.y=TRUE)
check <- c("hip_id_screen","Randomization.Group","redcap_event_name")


#add leading zeroes to some IDs
final$hip_id_screen <- as.character(final$hip_id_screen)
final$hip_id_screen[nchar(final$hip_id_screen)==2] <- paste0("0", final$hip_id_screen[nchar(final$hip_id_screen)==2])
final$hip_id_screen[nchar(final$hip_id_screen)==1] <- paste0("00", final$hip_id_screen[nchar(final$hip_id_screen)==1])

# delete two participants who did not get randomized
final <- final[!(final$hip_id_screen=="078-T" | final$hip_id_screen=="140-T" ),]
dim(final)

# remove those on metformin
#nomet <- final
#nomet <- nomet[!(final$Randomization.Group=="1") | is.na(final$Randomization.Group),]

# delete visits with missing insulin sensitivity, secretion, and DI
#nomet <- nomet[!is.na(nomet$insulin_secretion_mm) & !is.na(nomet$insulin_sensitivity) & !is.na(nomet$disposition_index),]

# need to delete those with only one visit
foranalysis <- final[final$hip_id_screen %in% names(table(final$hip_id_screen))[table(final$hip_id_screen) > 1],]

# add in baseline Tanner, race, sex
#Read Data
#data=read.csv('baseline vars.csv')
#Setting Labels

# need to add leading zeroes to IDs in the baseline file
#add leading zeroes to some IDs
#data$hip_id_screen <- as.character(data$hip_id_screen)
#data$hip_id_screen[nchar(data$hip_id_screen)==2] <- paste0("0", data$hip_id_screen[nchar(data$hip_id_screen)==2])
#data$hip_id_screen[nchar(data$hip_id_screen)==1] <- paste0("00", data$hip_id_screen[nchar(data$hip_id_screen)==1])


label(foranalysis$hip_id_screen)="HIP ID"
label(foranalysis$redcap_event_name)="Event Name"
label(foranalysis$sex)="Gender"
label(foranalysis$ethnicity)="Ethnicity"
label(foranalysis$race___1)="Race (choice=1 Native American)"
label(foranalysis$race___2)="Race (choice=2 Asian)"
label(foranalysis$race___3)="Race (choice=3 Black)"
label(foranalysis$race___4)="Race (choice=4 White)"
label(foranalysis$race___5)="Race (choice=5 Hispanic)"
label(foranalysis$race___6)="Race (choice=6 Other)"
label(foranalysis$if_other_race_ic)="If other race, please specify"
label(foranalysis$tanner_breast_exam)="Tanner Stage, Breast"
label(foranalysis$tanner_hair_exam)="Tanner Stage, Pubic Hair"
label(foranalysis$tanner_genetalia_exam)="Tanner Stage, Genetalia"
label(foranalysis$tanner_testis_vol_exam)="Tanner Stage, Testicular Volume    T1   < = 3ml    T2      >3 to 8 ml    T3      >=8 to 12 ml    T4     >=12 to < =15 ml    T5   >15 ml"
label(foranalysis$bmi_cat_final)="Final BMI Category  "
#Setting Units


#Setting Factors(will create new variable for factors)
# foranalysis$redcap_event_name.factor = factor(foranalysis$redcap_event_name,levels=c("ivgtt_visit_1_arm_1","exam_visit_1_arm_1","exam_visit_2_arm_1","ivgtt_visit_2_arm_1","exam_visit_3_arm_1","exam_visit_4_arm_1","ivgtt_visit_3_arm_1","exam_visit_5_arm_1","exam_visit_6_arm_1","exam_visit_7_arm_1","exam_visit_8_arm_1","exam_visit_9_arm_1","exam_visit_10_arm_1","exam_visit_11_arm_1","exam_visit_12_arm_1","ivgtt_visit_1_arm_2","med_visit_1_arm_2","exam_visit_2_arm_2","med_visit_3_arm_2","exam_visit_4_arm_2","ivgtt_visit_2_arm_2","med_visit_5_arm_2","exam_visit_6_arm_2","med_visit_7_arm_2","exam_visit_8_arm_2","ivgtt_visit_3_arm_2","med_visit_9_arm_2","exam_visit_10_arm_2","med_visit_11_arm_2","exam_visit_12_arm_2","ivgtt_visit_4_arm_2","med_visit_13_arm_2","exam_visit_14_arm_2","med_visit_15_arm_2","exam_visit_16_arm_2","med_visit_17_arm_2","exam_visit_18_arm_2","med_visit_19_arm_2","exam_visit_20_arm_2"))
# foranalysis$sex.factor = factor(foranalysis$sex,levels=c("0 Female","1 Male"))
# foranalysis$ethnicity.factor = factor(foranalysis$ethnicity,levels=c("0","1","2"))
# foranalysis$race___1.factor = factor(foranalysis$race___1,levels=c("0","1"))
# foranalysis$race___2.factor = factor(foranalysis$race___2,levels=c("0","1"))
# foranalysis$race___3.factor = factor(foranalysis$race___3,levels=c("0","1"))
# foranalysis$race___4.factor = factor(foranalysis$race___4,levels=c("0","1"))
# foranalysis$race___5.factor = factor(foranalysis$race___5,levels=c("0","1"))
# foranalysis$race___6.factor = factor(foranalysis$race___6,levels=c("0","1"))
# foranalysis$tanner_breast_exam.factor = factor(foranalysis$tanner_breast_exam,levels=c("1","2","3","4","5"))
# foranalysis$tanner_hair_exam.factor = factor(foranalysis$tanner_hair_exam,levels=c("1","2","3","4","5"))
# foranalysis$tanner_genetalia_exam.factor = factor(foranalysis$tanner_genetalia_exam,levels=c("1","2","3","4","5"))
# foranalysis$tanner_testis_vol_exam.factor = factor(foranalysis$tanner_testis_vol_exam,levels=c("1","2","3","4","5"))
# foranalysis$bmi_cat_final.factor = factor(foranalysis$bmi_cat_final,levels=c("1","2"))
# 
# levels(foranalysis$redcap_event_name.factor)=c("IVGTT Visit 1 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 1 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 2 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 2 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 3 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 4 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 3 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 5 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 6 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 7 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 8 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 9 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 10 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 11 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 12 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 1 (Arm 2: TREATMENT GROUP)","Med. Visit 1 (Arm 2: TREATMENT GROUP)","Exam Visit 2 (Arm 2: TREATMENT GROUP)","Med Visit 3 (Arm 2: TREATMENT GROUP)","Exam Visit 4 (Arm 2: TREATMENT GROUP)","IVGTT Visit 2 (Arm 2: TREATMENT GROUP)","Med Visit 5 (Arm 2: TREATMENT GROUP)","Exam Visit 6 (Arm 2: TREATMENT GROUP)","Med Visit 7 (Arm 2: TREATMENT GROUP)","Exam Visit 8 (Arm 2: TREATMENT GROUP)","IVGTT Visit 3 (Arm 2: TREATMENT GROUP)","Med. Visit 9 (Arm 2: TREATMENT GROUP)","Exam Visit 10 (Arm 2: TREATMENT GROUP)","Med. Visit 11 (Arm 2: TREATMENT GROUP)","Exam Visit 12 (Arm 2: TREATMENT GROUP)","IVGTT Visit 4 (Arm 2: TREATMENT GROUP)","Med. Visit 13 (Arm 2: TREATMENT GROUP)","Exam Visit 14 (Arm 2: TREATMENT GROUP)","Med. Visit 15 (Arm 2: TREATMENT GROUP)","Exam Visit 16 (Arm 2: TREATMENT GROUP)","Med. Visit 17 (Arm 2: TREATMENT GROUP)","Exam Visit 18 (Arm 2: TREATMENT GROUP)","Med. Visit 19 (Arm 2: TREATMENT GROUP)","Exam Visit 20 (Arm 2: TREATMENT GROUP)")
# levels(foranalysis$sex.factor)=c("0 Female","1 Male")
# levels(foranalysis$ethnicity.factor)=c("Hispanic or Latino","NOT Hispanic or Latino","Unknown / Not Reported")
# levels(foranalysis$race___1.factor)=c("Unchecked","Checked")
# levels(foranalysis$race___2.factor)=c("Unchecked","Checked")
# levels(foranalysis$race___3.factor)=c("Unchecked","Checked")
# levels(foranalysis$race___4.factor)=c("Unchecked","Checked")
# levels(foranalysis$race___5.factor)=c("Unchecked","Checked")
# levels(foranalysis$race___6.factor)=c("Unchecked","Checked")
# levels(foranalysis$tanner_breast_exam.factor)=c("Tanner 1","Tanner 2","Tanner 3","Tanner 4","Tanner 5")
# levels(foranalysis$tanner_hair_exam.factor)=c("Tanner 1","Tanner 2","Tanner 3","Tanner 4","Tanner 5")
# levels(foranalysis$tanner_genetalia_exam.factor)=c("Tanner 1","Tanner 2","Tanner 3","Tanner 4","Tanner 5")
# levels(foranalysis$tanner_testis_vol_exam.factor)=c("Tanner 1","Tanner 2","Tanner 3","Tanner 4","Tanner 5")
# levels(foranalysis$bmi_cat_final.factor)=c("Normal Weight","Obese")

#data <- data[,-2]

#dim(foranalysis)
#foranalysis <- merge(foranalysis,data,by=c("hip_id_screen"),all.x=TRUE, all.y=FALSE)
#dim(foranalysis)

# combine race/ethnicity
#foranalysis$sumrace <- as.numeric(foranalysis$race___1)  + as.numeric(foranalysis$race___2) + 
#  as.numeric(foranalysis$race___3) + as.numeric(foranalysis$race___4) + as.numeric(foranalysis$race___5)
foranalysis$race_eth[foranalysis$ethnicity==0] <- "Hispanic"
foranalysis$race_eth[(foranalysis$ethnicity==1 | foranalysis$ethnicity==2) & foranalysis$race___1==1] <- "Native American"
foranalysis$race_eth[(foranalysis$ethnicity==1 | foranalysis$ethnicity==2) & foranalysis$race___2==1] <- "Asian"
foranalysis$race_eth[(foranalysis$ethnicity==1 | foranalysis$ethnicity==2) & foranalysis$race___3==1] <- "Black"
foranalysis$race_eth[(foranalysis$ethnicity==1 | foranalysis$ethnicity==2) & foranalysis$race___4==1] <- "NHW"
foranalysis$race_eth[(foranalysis$ethnicity==1 | foranalysis$ethnicity==2) & foranalysis$race___5==1] <- "Hispanic"
foranalysis$race_eth[is.na(foranalysis$race_eth)] <- "Other"
racevars <- c("sumrace","ethnicity.factor","race___1","race___2","race___3","race___4","race___5","race_eth")

# combine Tanner
foranalysis$tanner[foranalysis$sex=="0 Female"] <- levels(foranalysis$tanner_breast_exam)[foranalysis$tanner_breast_exam[foranalysis$sex=="0 Female"]]
foranalysis$tanner[foranalysis$sex=="1 Male"] <- levels(foranalysis$tanner_testis_vol_exam)[foranalysis$tanner_testis_vol_exam[foranalysis$sex=="1 Male"]]
tannervars <- c("sex","tanner_breast_exam","tanner_testis_vol_exam","tanner")

# add famhx
#data=read.csv('famhx.csv')
#Setting Labels

#label(data$hip_id_screen)="HIP ID"
#label(data$redcap_event_name)="Event Name"
label(foranalysis$mother_history_screen)="Is mothers history available?"
label(foranalysis$type_2_db_mother_screen)="Type 2 Diabetes"
label(foranalysis$father_history_screen)="Is fathers history available?"
label(foranalysis$type_2_db_father_screen)="Type 2 Diabetes"
label(foranalysis$sib_no)="How many siblings does participant have?"
label(foranalysis$sib_1_history_sv)="Is sibling 1 history available?"
label(foranalysis$sb1_full_or_half_sv)="Is sibling 1 half or full sibling? "
label(foranalysis$type2_db_sb1_screen)="Type 2 Diabetes"
label(foranalysis$sb2_history_screen)="Is sibling 2 history available?"
label(foranalysis$sb2_full_or_half_screen)="Is sibling 2 half or full sibling?"
label(foranalysis$type2_db_sb2_screen)="Type 2 Diabetes"
label(foranalysis$sb3_history_screen)="Is the sibling 3 history available?"
label(foranalysis$full_sb3_half_or_full_sv)="Is sibling 3 half or full sibling?"
label(foranalysis$type2_db_sb3_screen)="Type 2 Diabetes"
label(foranalysis$sb4_history_screen)="Is sibling 4 history available?"
label(foranalysis$sb4_half_or_full_screen)="Is sibling 4 half or full sibling?"
label(foranalysis$type2_db_sb4_screen)="Type 2 Diabetes"
label(foranalysis$sb5_history_screen)="Is sibling 5 history available?"
label(foranalysis$sb5_half_or_full_screening)="Is sibling 5 half or full sibling?"
label(foranalysis$type2_db_sb5_screen)="Type 2 Diabetes"
#Setting Units


#Setting Factors(will create new variable for factors)
# foranalysis$redcap_event_name.factor = factor(foranalysis$redcap_event_name,levels=c("ivgtt_visit_1_arm_1","exam_visit_1_arm_1","exam_visit_2_arm_1","ivgtt_visit_2_arm_1","exam_visit_3_arm_1","exam_visit_4_arm_1","ivgtt_visit_3_arm_1","exam_visit_5_arm_1","exam_visit_6_arm_1","exam_visit_7_arm_1","exam_visit_8_arm_1","exam_visit_9_arm_1","exam_visit_10_arm_1","exam_visit_11_arm_1","exam_visit_12_arm_1","ivgtt_visit_1_arm_2","med_visit_1_arm_2","exam_visit_2_arm_2","med_visit_3_arm_2","exam_visit_4_arm_2","ivgtt_visit_2_arm_2","med_visit_5_arm_2","exam_visit_6_arm_2","med_visit_7_arm_2","exam_visit_8_arm_2","ivgtt_visit_3_arm_2","med_visit_9_arm_2","exam_visit_10_arm_2","med_visit_11_arm_2","exam_visit_12_arm_2","ivgtt_visit_4_arm_2","med_visit_13_arm_2","exam_visit_14_arm_2","med_visit_15_arm_2","exam_visit_16_arm_2","med_visit_17_arm_2","exam_visit_18_arm_2","med_visit_19_arm_2","exam_visit_20_arm_2"))
# foranalysis$mother_history_screen.factor = factor(foranalysis$mother_history_screen,levels=c("1","0"))
# foranalysis$type_2_db_mother_screen.factor = factor(foranalysis$type_2_db_mother_screen,levels=c("1","0"))
# foranalysis$father_history_screen.factor = factor(foranalysis$father_history_screen,levels=c("1","0"))
# foranalysis$type_2_db_father_screen.factor = factor(foranalysis$type_2_db_father_screen,levels=c("1","0"))
# foranalysis$sib_1_history_sv.factor = factor(foranalysis$sib_1_history_sv,levels=c("1","0"))
# foranalysis$sb1_full_or_half_sv.factor = factor(foranalysis$sb1_full_or_half_sv,levels=c("1","2"))
# foranalysis$type2_db_sb1_screen.factor = factor(foranalysis$type2_db_sb1_screen,levels=c("1","0"))
# foranalysis$sb2_history_screen.factor = factor(foranalysis$sb2_history_screen,levels=c("1","0"))
# foranalysis$sb2_full_or_half_screen.factor = factor(foranalysis$sb2_full_or_half_screen,levels=c("1","2"))
# foranalysis$type2_db_sb2_screen.factor = factor(foranalysis$type2_db_sb2_screen,levels=c("1","0"))
# foranalysis$sb3_history_screen.factor = factor(foranalysis$sb3_history_screen,levels=c("1","0"))
# foranalysis$full_sb3_half_or_full_sv.factor = factor(foranalysis$full_sb3_half_or_full_sv,levels=c("1","2"))
# foranalysis$type2_db_sb3_screen.factor = factor(foranalysis$type2_db_sb3_screen,levels=c("1","0"))
# foranalysis$sb4_history_screen.factor = factor(foranalysis$sb4_history_screen,levels=c("1","0"))
# foranalysis$sb4_half_or_full_screen.factor = factor(foranalysis$sb4_half_or_full_screen,levels=c("1","2"))
# foranalysis$type2_db_sb4_screen.factor = factor(foranalysis$type2_db_sb4_screen,levels=c("1","0"))
# foranalysis$sb5_history_screen.factor = factor(foranalysis$sb5_history_screen,levels=c("1","0"))
# foranalysis$sb5_half_or_full_screening.factor = factor(foranalysis$sb5_half_or_full_screening,levels=c("1","2"))
# foranalysis$type2_db_sb5_screen.factor = factor(foranalysis$type2_db_sb5_screen,levels=c("1","0"))
# 
# levels(foranalysis$redcap_event_name.factor)=c("IVGTT Visit 1 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 1 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 2 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 2 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 3 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 4 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 3 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 5 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 6 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 7 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 8 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 9 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 10 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 11 (Arm 1: NOT TREATMENT GROUP)","Exam Visit 12 (Arm 1: NOT TREATMENT GROUP)","IVGTT Visit 1 (Arm 2: TREATMENT GROUP)","Med. Visit 1 (Arm 2: TREATMENT GROUP)","Exam Visit 2 (Arm 2: TREATMENT GROUP)","Med Visit 3 (Arm 2: TREATMENT GROUP)","Exam Visit 4 (Arm 2: TREATMENT GROUP)","IVGTT Visit 2 (Arm 2: TREATMENT GROUP)","Med Visit 5 (Arm 2: TREATMENT GROUP)","Exam Visit 6 (Arm 2: TREATMENT GROUP)","Med Visit 7 (Arm 2: TREATMENT GROUP)","Exam Visit 8 (Arm 2: TREATMENT GROUP)","IVGTT Visit 3 (Arm 2: TREATMENT GROUP)","Med. Visit 9 (Arm 2: TREATMENT GROUP)","Exam Visit 10 (Arm 2: TREATMENT GROUP)","Med. Visit 11 (Arm 2: TREATMENT GROUP)","Exam Visit 12 (Arm 2: TREATMENT GROUP)","IVGTT Visit 4 (Arm 2: TREATMENT GROUP)","Med. Visit 13 (Arm 2: TREATMENT GROUP)","Exam Visit 14 (Arm 2: TREATMENT GROUP)","Med. Visit 15 (Arm 2: TREATMENT GROUP)","Exam Visit 16 (Arm 2: TREATMENT GROUP)","Med. Visit 17 (Arm 2: TREATMENT GROUP)","Exam Visit 18 (Arm 2: TREATMENT GROUP)","Med. Visit 19 (Arm 2: TREATMENT GROUP)","Exam Visit 20 (Arm 2: TREATMENT GROUP)")
# levels(foranalysis$mother_history_screen.factor)=c("Yes","No")
# levels(foranalysis$type_2_db_mother_screen.factor)=c("Yes","No")
# levels(foranalysis$father_history_screen.factor)=c("Yes","No")
# levels(foranalysis$type_2_db_father_screen.factor)=c("Yes","No")
# levels(foranalysis$sib_1_history_sv.factor)=c("Yes","No")
# levels(foranalysis$sb1_full_or_half_sv.factor)=c("0 Half","1 Full")
# levels(foranalysis$type2_db_sb1_screen.factor)=c("Yes","No")
# levels(foranalysis$sb2_history_screen.factor)=c("Yes","No")
# levels(foranalysis$sb2_full_or_half_screen.factor)=c("0 Half","1 Full")
# levels(foranalysis$type2_db_sb2_screen.factor)=c("Yes","No")
# levels(foranalysis$sb3_history_screen.factor)=c("Yes","No")
# levels(foranalysis$full_sb3_half_or_full_sv.factor)=c("0 Half","1 Full")
# levels(foranalysis$type2_db_sb3_screen.factor)=c("Yes","No")
# levels(foranalysis$sb4_history_screen.factor)=c("Yes","No")
# levels(foranalysis$sb4_half_or_full_screen.factor)=c("0 Half","1 Full")
# levels(foranalysis$type2_db_sb4_screen.factor)=c("Yes","No")
# levels(foranalysis$sb5_history_screen.factor)=c("Yes","No")
# levels(foranalysis$sb5_half_or_full_screening.factor)=c("0 Half","1 Full")
# levels(foranalysis$type2_db_sb5_screen.factor)=c("Yes","No")

# create variable for family history of T2D
# per Megan, if one parent hx available, code per parent 
# first parent history
foranalysis$parhx[foranalysis$mother_history_screen=="No" & foranalysis$father_history_screen=="No"] <- NA
foranalysis$parhx[(foranalysis$mother_history_screen=="Yes" | foranalysis$father_history_screen=="Yes") &
                    (foranalysis$type_2_db_mother_screen=="Yes" | foranalysis$type_2_db_father_screen=="Yes")] <- "Yes"
foranalysis$parhx[(foranalysis$mother_history_screen=="Yes" & foranalysis$father_history_screen=="No") &
                    (foranalysis$type_2_db_mother_screen=="Yes") ] <- "Yes"
foranalysis$parhx[(foranalysis$mother_history_screen=="No" & foranalysis$father_history_screen=="Yes") &
                    (foranalysis$type_2_db_father_screen=="Yes") ] <- "Yes"
foranalysis$parhx[(foranalysis$mother_history_screen=="Yes" | foranalysis$father_history_screen=="Yes") &
                    (foranalysis$type_2_db_mother_screen=="No" & foranalysis$type_2_db_father_screen=="No")] <- "No"
foranalysis$parhx[(foranalysis$mother_history_screen=="Yes" & foranalysis$father_history_screen=="No") &
                    (foranalysis$type_2_db_mother_screen=="No") ] <- "No"
foranalysis$parhx[(foranalysis$mother_history_screen=="No" & foranalysis$father_history_screen=="Yes") &
                    (foranalysis$type_2_db_father_screen=="No") ] <- "No"
# then sib history
# per Megan, if one sib hx available, code per sib
# no siblings
foranalysis$sibhx[foranalysis$sib_no==0] <- "No"
# at least one sibling and all missing
foranalysis$sibhx[foranalysis$sib_no >= 1 & foranalysis$sib_1_history_sv=="No" & 
                    (foranalysis$sb2_history_screen=="No" | is.na(foranalysis$sb2_history_screen)) &
                    (foranalysis$sb3_history_screen=="No" | is.na(foranalysis$sb3_history_screen)) &
                    (foranalysis$sb4_history_screen=="No" | is.na(foranalysis$sb4_history_screen)) &
                    (foranalysis$sb5_history_screen=="No" | is.na(foranalysis$sb5_history_screen)) 
                    ] <- NA
#View(foranalysis[c("sib_no","sib_1_history_sv","sb2_history_screen","type2_db_sb1_screen",
                #   "type2_db_sb2_screen","sibhx","parhx","famhxt2d")])
# at least one sibling and all negative
foranalysis$sibhx[foranalysis$sib_no>=1 & 
                    (foranalysis$type2_db_sb1_screen=="No" | is.na(foranalysis$type2_db_sb1_screen)) &
                    (foranalysis$type2_db_sb2_screen=="No" | is.na(foranalysis$type2_db_sb2_screen)) &
                    (foranalysis$type2_db_sb3_screen=="No" | is.na(foranalysis$type2_db_sb3_screen)) &
                    (foranalysis$type2_db_sb4_screen=="No" | is.na(foranalysis$type2_db_sb4_screen)) &
                    (foranalysis$type2_db_sb5_screen=="No" | is.na(foranalysis$type2_db_sb5_screen)) 
                  ] <- "No"
# at least one sibling and at least one positive
foranalysis$sibhx[foranalysis$sib_no>=1 & 
                    ( foranalysis$type2_db_sb1_screen=="Yes" |
                    foranalysis$type2_db_sb2_screen=="Yes" |
                    foranalysis$type2_db_sb3_screen=="Yes" |
                    foranalysis$type2_db_sb4_screen=="Yes" |
                    foranalysis$type2_db_sb5_screen=="Yes") 
                  ] <- "Yes"
# now combine parent and sib hx
foranalysis$famhxt2d[is.na(foranalysis$parhx) & is.na(foranalysis$sibhx)] <- NA
foranalysis$famhxt2d[foranalysis$parhx=="Yes" | foranalysis$sibhx=="Yes"] <- "Yes"
foranalysis$famhxt2d[foranalysis$parhx=="No" & foranalysis$sibhx=="No"] <- "No"
foranalysis$famhxt2d[foranalysis$parhx=="No" & is.na(foranalysis$sibhx)] <- "No"
foranalysis$famhxt2d[is.na(foranalysis$parhx) & foranalysis$sibhx=="No"] <- "No"
famhx <- foranalysis[foranalysis$redcap_event_name.factor=="IVGTT Visit 1 (Arm 1: NOT TREATMENT GROUP)" |
                foranalysis$redcap_event_name.factor=="IVGTT Visit 1 (Arm 2: TREATMENT GROUP)",
              c("hip_id_screen","famhxt2d")]

# write final long dataset
write.csv(foranalysis,file="for_analysis.csv", na="")

# make a wide dataset
keep <- c("hip_id_screen","date_of_study_visit","sv_num","insulin_sensitivity","insulin_secretion_mm",
          "disposition_index","fat_percentage_dexa","bmi_cat_final","sex","race_eth","tanner",
          "Randomization.Group","lept","adiponect_lv","tcholes_lv","hdl_lv","ldl_meth_lv","ldl_lv","tg_lv",
          "crp_lv","fat_percentage_dexa","zscore_lv","alt","ast","dhea_s","wc_avg_lv","igf_lv",
          "bp_systol_exam","bp_diastol_exam","famhxt2d","BMI","date_of_study_visit","dob")
temp <- foranalysis[keep]

wide <- reshape(temp, timevar="sv_num", idvar= c("hip_id_screen"),direction="wide")
dim(wide)
# rename baseline variables
wide$bmi_cat_final <- wide$bmi_cat_final.1
wide$sex <- wide$sex.1
wide$race_eth <- wide$race_eth.1
wide$tanner <- wide$tanner.1
wide$Randomization.Group <- wide$Randomization.Group.1
wide$famhxt2d <- wide$famhxt2d.1
wide$base_bmi <- wide$BMI.1
wide$baseline_date <- wide$date_of_study_visit.1
wide$dob <- wide$dob.1

wide <- select(wide,-c("bmi_cat_final.1","bmi_cat_final.2","bmi_cat_final.3","bmi_cat_final.4","sex.1",
                       "sex.2","sex.3","sex.4","race_eth.1","race_eth.2","race_eth.3","race_eth.4",
                       "tanner.1","tanner.2","tanner.3","tanner.4","Randomization.Group.1",
                       "Randomization.Group.2","Randomization.Group.3","Randomization.Group.4",
                       "famhxt2d.1","famhxt2d.2","famhxt2d.3","famhxt2d.4","BMI.1","BMI.2","BMI.3","BMI.4",
                       "date_of_study_visit.1","date_of_study_visit.2","date_of_study_visit.3","date_of_study_visit.4",
                       "date_of_study_visit.1.1","date_of_study_visit.1.2","date_of_study_visit.1.3","date_of_study_visit.1.4",
                       "dob.1","dob.2","dob.3","dob.4"))

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

wide <- select(wide,-c("fat_percentage_dexa.1.1","fat_percentage_dexa.1.2","fat_percentage_dexa.1.3","fat_percentage_dexa.1.4"))

# delete people with missing delta Si, insulin secretion or DI
dim(wide)
wide <- wide[!is.na(wide$delta_di) & !is.na(wide$delta_inssec) & !is.na(wide$delta_di),]
dim(wide)

# delete those who were not randomized
#dim(wide)
#wide <- wide[!is.na(wide$Randomization.Group),]
#dim(wide)

#delete normal weight kids
wide <- wide[wide$bmi_cat_final=="Obese",]

# set labels
label(wide$delta_di)="Change in DI"
label(wide$delta_inssec)="Change in insulin secretion"
label(wide$delta_si)="Change in insulin sensitivity"
label(wide$delta_bodyfat)="Change in percent body fat"

# calculate age
wide$dob <- as_date(wide$dob)
wide$temp<- as.Date(wide$baseline_date,format = "%m/%d/%Y")
wide$age_base <- as.numeric(floor(difftime(wide$temp,wide$dob,units="weeks")/52))


