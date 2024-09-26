# Read Data
data <- read.csv("/Users/timvigers/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Vigers/BDC/Janet Snell-Bergeon/CALICO/Data_Raw/CALICOMulticenterPCO-NCHCHOPProject_DATA_2024-09-25_0932.csv",na.strings = "")
# Setting Labels
label(data$record_number) <- "Record number"
label(data$redcap_repeat_instrument) <- "Repeat Instrument"
label(data$redcap_repeat_instance) <- "Repeat Instance"
label(data$pcosdx_pmh___1) <- "Past medical history (do not include diagnoses made this visit) (choice=Premature adrenarche)"
label(data$pcosdx_pmh___2) <- "Past medical history (do not include diagnoses made this visit) (choice=Central precocious puberty)"
label(data$pcosdx_pmh___3) <- "Past medical history (do not include diagnoses made this visit) (choice=Pregnancy)"
label(data$pcosdx_pmh___4) <- "Past medical history (do not include diagnoses made this visit) (choice=Cystic acne)"
label(data$pcosdx_pmh___5) <- "Past medical history (do not include diagnoses made this visit) (choice=Overweight/Obesity)"
label(data$pcosdx_pmh___6) <- "Past medical history (do not include diagnoses made this visit) (choice=Pre-diabetes)"
label(data$pcosdx_pmh___29) <- "Past medical history (do not include diagnoses made this visit) (choice=Type 1 diabetes)"
label(data$pcosdx_pmh___7) <- "Past medical history (do not include diagnoses made this visit) (choice=Type 2 diabetes)"
label(data$pcosdx_pmh___8) <- "Past medical history (do not include diagnoses made this visit) (choice=OSA)"
label(data$pcosdx_pmh___9) <- "Past medical history (do not include diagnoses made this visit) (choice=Use of CPAP)"
label(data$pcosdx_pmh___10) <- "Past medical history (do not include diagnoses made this visit) (choice=Pseudotumor)"
label(data$pcosdx_pmh___11) <- "Past medical history (do not include diagnoses made this visit) (choice=GERD)"
label(data$pcosdx_pmh___12) <- "Past medical history (do not include diagnoses made this visit) (choice=Nonalcoholic fatty liver disease)"
label(data$pcosdx_pmh___22) <- "Past medical history (do not include diagnoses made this visit) (choice=Hypothyroidism)"
label(data$pcosdx_pmh___13) <- "Past medical history (do not include diagnoses made this visit) (choice=Hypertension)"
label(data$pcosdx_pmh___14) <- "Past medical history (do not include diagnoses made this visit) (choice=Hyperlipidemia)"
label(data$pcosdx_pmh___15) <- "Past medical history (do not include diagnoses made this visit) (choice=Anxiety)"
label(data$pcosdx_pmh___16) <- "Past medical history (do not include diagnoses made this visit) (choice=Depression)"
label(data$pcosdx_pmh___23) <- "Past medical history (do not include diagnoses made this visit) (choice=Suicide attempt)"
label(data$pcosdx_pmh___24) <- "Past medical history (do not include diagnoses made this visit) (choice=Mental health hospitalization)"
label(data$pcosdx_pmh___17) <- "Past medical history (do not include diagnoses made this visit) (choice=Binge eating disorder (including bulemia))"
label(data$pcosdx_pmh___18) <- "Past medical history (do not include diagnoses made this visit) (choice=Restricting eating disorder (including anorexia))"
label(data$pcosdx_pmh___19) <- "Past medical history (do not include diagnoses made this visit) (choice=ADHD)"
label(data$pcosdx_pmh___20) <- "Past medical history (do not include diagnoses made this visit) (choice=Asthma)"
label(data$pcosdx_pmh___21) <- "Past medical history (do not include diagnoses made this visit) (choice=Migraines)"
label(data$pcosdx_pmh___60) <- "Past medical history (do not include diagnoses made this visit) (choice=Other {pmh_other})"
label(data$pcosdx_pmh___0) <- "Past medical history (do not include diagnoses made this visit) (choice=None)"
label(data$pcosdx_pmh___unk) <- "Past medical history (do not include diagnoses made this visit) (choice=Unknown/Not recorded)"
label(data$pcosdx_pmh___na) <- "Past medical history (do not include diagnoses made this visit) (choice=Not applicable)"
label(data$pcosdx_pmh___oth) <- "Past medical history (do not include diagnoses made this visit) (choice=Other)"
label(data$pcosdx_pmh___pm) <- "Past medical history (do not include diagnoses made this visit) (choice=Premenarchal)"
label(data$pmh_other) <- "other"
label(data$site) <- "Study Site"
label(data$race___1) <- "Race (choice=Caucasian)"
label(data$race___2) <- "Race (choice=African American)"
label(data$race___3) <- "Race (choice=Asian)"
label(data$race___4) <- "Race (choice=Pacific Islander)"
label(data$race___5) <- "Race (choice=American Indian or Alaska Native)"
label(data$race___60) <- "Race (choice=Other {race_other})"
label(data$race___unk) <- "Race (choice=Unknown/Not recorded)"
label(data$race___na) <- "Race (choice=Not applicable)"
label(data$race___oth) <- "Race (choice=Other)"
label(data$race___pm) <- "Race (choice=Premenarchal)"
label(data$race_other) <- "Other"
label(data$ethnicity) <- "Ethnicity "
label(data$insur_type) <- "Insurance Type"
label(data$insur_other) <- "Other"
label(data$pcosdx_obesitydx_age) <- "Age when overweight/obesity diagnosed (years)"
label(data$pcosdx_obesitydx_earliestage) <- "If answered UNK for the previous question, when was the earliest age that overweight/obesity was noted? "
label(data$pcosdx_t1ddx_age) <- "Age when type 1 diabetes diagnosed (years)"
label(data$pcosdx_t2ddx_age) <- "Age when type 2 diabetes diagnosed (years)"
label(data$pcosdx_osadx_age) <- "Age OSA diagnosed (years)"
label(data$pcosdx_anxietydx_age) <- "Age when anxiety diagnosed (years)"
label(data$pcosdx_depressiondx_age) <- "Age when depression diagnosed (years)"
label(data$pcosdx_bingedx_age) <- "Age when binge eating noted (years)"
label(data$pcosdx_restrictdx_age) <- "Age when restrictive eating noted (years)"
label(data$pcosdx_adhddx_age) <- "Age when ADHD diagnosed (years)"
label(data$pcosdx_mhhosp_age) <- "Age when mental health hospitalization occurred (years)"
label(data$pcosdx_siattempt_age) <- "Age when suicide attempt occurred (years)"
label(data$pcosdx_surghx_bariatric) <- "Age at bariatric surgery (years)"
label(data$pcosdx_surghx_pilonidal) <- "Age at pilonidal cyst surgery (years)"
label(data$pcosdx_surghx_hidradenitis) <- "Age at incision and drainage for Hidradenitis (years)"
label(data$pastmeds___1) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Accutane)"
label(data$pastmeds___2) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Estrogen-containing pill)"
label(data$pastmeds___3) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Estrogen-containing patch)"
label(data$pastmeds___4) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Estrogen-containing ring)"
label(data$pastmeds___5) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Cyclic progesterone)"
label(data$pastmeds___6) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Daily progesterone)"
label(data$pastmeds___7) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Progesterone implant)"
label(data$pastmeds___8) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Progesterone IUD)"
label(data$pastmeds___9) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Topiramate for weight loss)"
label(data$pastmeds___10) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Phentermine for weight loss)"
label(data$pastmeds___11) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Topirimate/phentermine for weight loss)"
label(data$pastmeds___12) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Liraglutiade (saxcenda) for weight loss)"
label(data$pastmeds___13) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Semaglutide (Wygovy/Ozempic) for weight loss)"
label(data$pastmeds___14) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Metformin)"
label(data$pastmeds___15) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Vitamin D)"
label(data$pastmeds___16) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Multivitamin)"
label(data$pastmeds___60) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Other {pastmeds_other})"
label(data$pastmeds___0) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=None)"
label(data$pastmeds___unk) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Unknown/Not recorded)"
label(data$pastmeds___na) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Not applicable)"
label(data$pastmeds___oth) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Other)"
label(data$pastmeds___pm) <- "Medications taken during 12 months prior to PCOS diagnostic visit (may or may not be currently taking) (choice=Premenarchal)"
label(data$pastmeds_other) <- "Other"
label(data$pcosdx_mentalhealthcounseling___1) <- "Any mental health counseling in the 12 months prior to this visit (choice=Group therapy)"
label(data$pcosdx_mentalhealthcounseling___2) <- "Any mental health counseling in the 12 months prior to this visit (choice=Individual therapy)"
label(data$pcosdx_mentalhealthcounseling___3) <- "Any mental health counseling in the 12 months prior to this visit (choice=Inpatient therapy)"
label(data$pcosdx_mentalhealthcounseling___0) <- "Any mental health counseling in the 12 months prior to this visit (choice=Documented none)"
label(data$pcosdx_mentalhealthcounseling___unk) <- "Any mental health counseling in the 12 months prior to this visit (choice=Unknown/Not recorded)"
label(data$pcosdx_mentalhealthcounseling___na) <- "Any mental health counseling in the 12 months prior to this visit (choice=Not applicable)"
label(data$pcosdx_mentalhealthcounseling___oth) <- "Any mental health counseling in the 12 months prior to this visit (choice=Other)"
label(data$pcosdx_mentalhealthcounseling___pm) <- "Any mental health counseling in the 12 months prior to this visit (choice=Premenarchal)"
label(data$pcosdx_dietarycounseling___1) <- "Participated in dietary counseling sessions or met with a dietician prior to this visit (choice=Group therapy)"
label(data$pcosdx_dietarycounseling___2) <- "Participated in dietary counseling sessions or met with a dietician prior to this visit (choice=Individual therapy)"
label(data$pcosdx_dietarycounseling___0) <- "Participated in dietary counseling sessions or met with a dietician prior to this visit (choice=Documented none)"
label(data$pcosdx_dietarycounseling___unk) <- "Participated in dietary counseling sessions or met with a dietician prior to this visit (choice=Unknown/Not recorded)"
label(data$pcosdx_dietarycounseling___na) <- "Participated in dietary counseling sessions or met with a dietician prior to this visit (choice=Not applicable)"
label(data$pcosdx_dietarycounseling___oth) <- "Participated in dietary counseling sessions or met with a dietician prior to this visit (choice=Other)"
label(data$pcosdx_dietarycounseling___pm) <- "Participated in dietary counseling sessions or met with a dietician prior to this visit (choice=Premenarchal)"
label(data$pcosdx_exercising_teachingplan___1) <- "Received an exercise plan or met with an exercise specialist prior to this visit (choice=Group therapy)"
label(data$pcosdx_exercising_teachingplan___2) <- "Received an exercise plan or met with an exercise specialist prior to this visit (choice=Individual therapy)"
label(data$pcosdx_exercising_teachingplan___0) <- "Received an exercise plan or met with an exercise specialist prior to this visit (choice=Documented none)"
label(data$pcosdx_exercising_teachingplan___unk) <- "Received an exercise plan or met with an exercise specialist prior to this visit (choice=Unknown/Not recorded)"
label(data$pcosdx_exercising_teachingplan___na) <- "Received an exercise plan or met with an exercise specialist prior to this visit (choice=Not applicable)"
label(data$pcosdx_exercising_teachingplan___oth) <- "Received an exercise plan or met with an exercise specialist prior to this visit (choice=Other)"
label(data$pcosdx_exercising_teachingplan___pm) <- "Received an exercise plan or met with an exercise specialist prior to this visit (choice=Premenarchal)"
label(data$pcosdx_age) <- "Age at time of PCOS diagnosis"
label(data$pcosdx_irregular_menses) <- "Irregular menses defined"
label(data$pcosdx_hyperandrogenism) <- "Hyperandrogenism diagnostic method"
label(data$cv_type_of_clinic) <- "Type of clinic visit"
label(data$cv_mdc_specialties___1) <- "Specialties included in multidisciplinary clinic  (choice=Endocrine)"
label(data$cv_mdc_specialties___2) <- "Specialties included in multidisciplinary clinic  (choice=Gynecology)"
label(data$cv_mdc_specialties___3) <- "Specialties included in multidisciplinary clinic  (choice=Nutrition)"
label(data$cv_mdc_specialties___4) <- "Specialties included in multidisciplinary clinic  (choice=Psychiatry)"
label(data$cv_mdc_specialties___5) <- "Specialties included in multidisciplinary clinic  (choice=Dermatology)"
label(data$cv_mdc_specialties___6) <- "Specialties included in multidisciplinary clinic  (choice=Adolescent Medicine)"
label(data$cv_mdc_specialties___60) <- "Specialties included in multidisciplinary clinic  (choice=Other {cv_mdc_other})"
label(data$cv_mdc_specialties___unk) <- "Specialties included in multidisciplinary clinic  (choice=Unknown)"
label(data$cv_mdc_specialties___na) <- "Specialties included in multidisciplinary clinic  (choice=Not applicable)"
label(data$cv_mdc_specialties___oth) <- "Specialties included in multidisciplinary clinic  (choice=Other)"
label(data$cv_mdc_specialties___pm) <- "Specialties included in multidisciplinary clinic  (choice=Premenarchal)"
label(data$cv_monthssincepcosdx) <- "Months since PCOS diagnosis"
label(data$cv_cysticacnedx_age) <- "Months since PCOS dx when cystic acne diagnosed"
label(data$cv_obesitydx_age) <- "Months since PCOS dx when obesity/overweight diagnosed"
label(data$cv_anxietydx_age) <- "Months since PCOS dx when anxiety diagnosed"
label(data$cv_depressiondx_age) <- "Months since PCOS dx when depression diagnosed"
label(data$cv_bingedx_age) <- "Months since PCOS dx when binge eating noted "
label(data$cv_restrictdx_age) <- "Months since PCOS dx when restrictive eating noted"
label(data$cv_adhddx_age) <- "Months since PCOS dx when ADHD diagnosed "
label(data$cv_mh_hospital_age) <- "Months since PCOS dx when mental health hospitalization occurred"
label(data$cv_si_attempt_age) <- "Months since PCOS dx when suicide attempt occurred"
label(data$cv_othermedical) <- "Other diagnosis"
label(data$cv_other_age_2) <- "Months since PCOS dx when [cv_othermedical] diagnosed "
label(data$cv_medications___1) <- "Current treatment (not prescribed this visit) (choice=Metformin)"
label(data$cv_medications___2) <- "Current treatment (not prescribed this visit) (choice=Insulin)"
label(data$cv_medications___3) <- "Current treatment (not prescribed this visit) (choice=Secondary diabetes medication {cv_meds_diabetes2})"
label(data$cv_medications___4) <- "Current treatment (not prescribed this visit) (choice=Reflux medications {cv_reflux_meds})"
label(data$cv_medications___5) <- "Current treatment (not prescribed this visit) (choice=Estrogen-containing pill {cv_meds_pill})"
label(data$cv_medications___6) <- "Current treatment (not prescribed this visit) (choice=Estrogen-containing patch)"
label(data$cv_medications___7) <- "Current treatment (not prescribed this visit) (choice=Estrogen-containing ring)"
label(data$cv_medications___8) <- "Current treatment (not prescribed this visit) (choice=Cyclic progesterone)"
label(data$cv_medications___9) <- "Current treatment (not prescribed this visit) (choice=Daily progesterone)"
label(data$cv_medications___10) <- "Current treatment (not prescribed this visit) (choice=Progesterone implant)"
label(data$cv_medications___11) <- "Current treatment (not prescribed this visit) (choice=Progesterone IUD)"
label(data$cv_medications___12) <- "Current treatment (not prescribed this visit) (choice=Progesterone injection)"
label(data$cv_medications___13) <- "Current treatment (not prescribed this visit) (choice=Topical acne medication)"
label(data$cv_medications___14) <- "Current treatment (not prescribed this visit) (choice=Oral acne antibiotic)"
label(data$cv_medications___15) <- "Current treatment (not prescribed this visit) (choice=Accutane)"
label(data$cv_medications___16) <- "Current treatment (not prescribed this visit) (choice=Spironolactone)"
label(data$cv_medications___17) <- "Current treatment (not prescribed this visit) (choice=Atypical antipsychotic)"
label(data$cv_medications___18) <- "Current treatment (not prescribed this visit) (choice=Antidepressant or anti-anxiety)"
label(data$cv_medications___19) <- "Current treatment (not prescribed this visit) (choice=ADHD medication)"
label(data$cv_medications___20) <- "Current treatment (not prescribed this visit) (choice=Sleep aids {cv_sleep_meds})"
label(data$cv_medications___21) <- "Current treatment (not prescribed this visit) (choice=Omega-3 fatty acids)"
label(data$cv_medications___22) <- "Current treatment (not prescribed this visit) (choice=Anti-hypertensive)"
label(data$cv_medications___23) <- "Current treatment (not prescribed this visit) (choice=Vitamin D)"
label(data$cv_medications___32) <- "Current treatment (not prescribed this visit) (choice=Multivitamin)"
label(data$cv_medications___24) <- "Current treatment (not prescribed this visit) (choice=Lipid-lowering)"
label(data$cv_medications___25) <- "Current treatment (not prescribed this visit) (choice=Topiramate for weight loss)"
label(data$cv_medications___26) <- "Current treatment (not prescribed this visit) (choice=Phentermine for weight)"
label(data$cv_medications___27) <- "Current treatment (not prescribed this visit) (choice=Topirimate/phentermine for weight loss)"
label(data$cv_medications___28) <- "Current treatment (not prescribed this visit) (choice=Liraglutiade (saxcenda) for weight loss)"
label(data$cv_medications___29) <- "Current treatment (not prescribed this visit) (choice=Semaglutide (Wygovy/Ozempic) for weight loss)"
label(data$cv_medications___30) <- "Current treatment (not prescribed this visit) (choice=Levothyroxine)"
label(data$cv_medications___31) <- "Current treatment (not prescribed this visit) (choice=Systemic steroids (ex: prednisone, dexamethasone, solumedrol))"
label(data$cv_medications___60) <- "Current treatment (not prescribed this visit) (choice=Other {cv_othercomment})"
label(data$cv_medications___0) <- "Current treatment (not prescribed this visit) (choice=None)"
label(data$cv_medications___unk) <- "Current treatment (not prescribed this visit) (choice=Unknown/Not recorded)"
label(data$cv_medications___na) <- "Current treatment (not prescribed this visit) (choice=Not applicable)"
label(data$cv_medications___oth) <- "Current treatment (not prescribed this visit) (choice=Other)"
label(data$cv_medications___pm) <- "Current treatment (not prescribed this visit) (choice=Premenarchal)"
label(data$cv_meds_diabetes2) <- "Secondary diabetes meds"
label(data$cv_weight) <- "Weight (kg)"
label(data$cv_height) <- "Height (cm)"
label(data$cv_bmi) <- "BMI"
label(data$cv_bmi_percentile) <- "BMI percentile"
label(data$cv_bmi_z) <- "BMI z-score"
label(data$cv_sbp) <- "Systolic BP (mmHg) "
label(data$cv_dbp) <- "Diastolic BP (mmHg) "
label(data$cv_genderid) <- "Patients gender identity"
label(data$cv_hirsutism_num) <- "Numerical FGS score "
label(data$cv_hirsutism_cat) <- "FGS score category"
label(data$cv_acneface) <- "Acne severity on the face"
label(data$cv_acneother___1) <- "Acne elsewhere on the body (choice=Back)"
label(data$cv_acneother___2) <- "Acne elsewhere on the body (choice=Chest)"
label(data$cv_acneother___0) <- "Acne elsewhere on the body (choice=None)"
label(data$cv_acneother___unk) <- "Acne elsewhere on the body (choice=Unknown/Not recorded)"
label(data$cv_acneother___na) <- "Acne elsewhere on the body (choice=Not applicable)"
label(data$cv_acneother___oth) <- "Acne elsewhere on the body (choice=Other)"
label(data$cv_acneother___pm) <- "Acne elsewhere on the body (choice=Premenarchal)"
label(data$cv_acanthosisneck) <- "Acanthosis at the neck"
label(data$cv_hydradinitis) <- "Any hydradinitis (may be described as abscess or red lesions in axilla, pannus or groin)"
label(data$cv_alopecia) <- "Androgenic alopecia (hair loss at front of scalp)"
label(data$cv_ft) <- "Free testosterone (pg/mL) "
label(data$cv_ft_uln) <- "Free Testosterone (pg/mL) value for upper limit of normal for sex, age/tanner stage  for assay"
label(data$cv_ft_perc) <- "Free testosterone (pg/mL) percentage of upper limit of normal"
label(data$cv_tt) <- "Total testosterone (ng/dL) "
label(data$cv_tt_assay) <- "Assay for total testosterone"
label(data$cv_tt_uln) <- "Total testosterone (ng/dL) value for upper limit of normal for assay"
label(data$cv_tt_perc) <- "Total testosterone percentage of upper limit of normal"
label(data$cv_a1c) <- "HbA1C (%)"
label(data$cv_mood) <- "Psychology screening or mood assessment at this visit"
label(data$cv_phq2) <- "PHQ-2 result, if performed at visit"
label(data$cv_phq8) <- "PHQ-8 result, if performed at visit"
label(data$cv_phq9) <- "PHQ-9 result, if performed at visit"
label(data$cv_cesd20) <- "CESD-20 result, if performed at visit"
label(data$cv_newmeds___1) <- "New or continued medications/supplements prescribed at this visit  (choice=Metformin)"
label(data$cv_newmeds___2) <- "New or continued medications/supplements prescribed at this visit  (choice=Insulin)"
label(data$cv_newmeds___3) <- "New or continued medications/supplements prescribed at this visit  (choice=Secondary diabetes medication {cv_newmeds_diabetes})"
label(data$cv_newmeds___4) <- "New or continued medications/supplements prescribed at this visit  (choice=Reflux medications {cv_reflux_newmeds})"
label(data$cv_newmeds___5) <- "New or continued medications/supplements prescribed at this visit  (choice=Estrogen-containing pill {cv_newmeds_pill})"
label(data$cv_newmeds___6) <- "New or continued medications/supplements prescribed at this visit  (choice=Estrogen-containing patch)"
label(data$cv_newmeds___7) <- "New or continued medications/supplements prescribed at this visit  (choice=Estrogen-containing ring)"
label(data$cv_newmeds___8) <- "New or continued medications/supplements prescribed at this visit  (choice=Cyclic progesterone)"
label(data$cv_newmeds___9) <- "New or continued medications/supplements prescribed at this visit  (choice=Daily progesterone)"
label(data$cv_newmeds___10) <- "New or continued medications/supplements prescribed at this visit  (choice=Progesterone implant)"
label(data$cv_newmeds___11) <- "New or continued medications/supplements prescribed at this visit  (choice=Progesterone IUD)"
label(data$cv_newmeds___12) <- "New or continued medications/supplements prescribed at this visit  (choice=Progesterone injection)"
label(data$cv_newmeds___13) <- "New or continued medications/supplements prescribed at this visit  (choice=Topical acne medication)"
label(data$cv_newmeds___14) <- "New or continued medications/supplements prescribed at this visit  (choice=Oral acne antibiotic)"
label(data$cv_newmeds___15) <- "New or continued medications/supplements prescribed at this visit  (choice=Accutane)"
label(data$cv_newmeds___16) <- "New or continued medications/supplements prescribed at this visit  (choice=Spironolactone)"
label(data$cv_newmeds___17) <- "New or continued medications/supplements prescribed at this visit  (choice=Atypical antipsychotic)"
label(data$cv_newmeds___18) <- "New or continued medications/supplements prescribed at this visit  (choice=Antidepressant or anti-anxiety)"
label(data$cv_newmeds___19) <- "New or continued medications/supplements prescribed at this visit  (choice=ADHD medication)"
label(data$cv_newmeds___20) <- "New or continued medications/supplements prescribed at this visit  (choice=Sleep aids {cv_sleep_newmeds})"
label(data$cv_newmeds___21) <- "New or continued medications/supplements prescribed at this visit  (choice=Omega-3 fatty acids)"
label(data$cv_newmeds___22) <- "New or continued medications/supplements prescribed at this visit  (choice=Anti-hypertensive)"
label(data$cv_newmeds___23) <- "New or continued medications/supplements prescribed at this visit  (choice=Vitamin D)"
label(data$cv_newmeds___32) <- "New or continued medications/supplements prescribed at this visit  (choice=Multivitamin)"
label(data$cv_newmeds___24) <- "New or continued medications/supplements prescribed at this visit  (choice=Lipid-lowering)"
label(data$cv_newmeds___25) <- "New or continued medications/supplements prescribed at this visit  (choice=Topiramate for weight loss)"
label(data$cv_newmeds___26) <- "New or continued medications/supplements prescribed at this visit  (choice=Phentermine for weight)"
label(data$cv_newmeds___27) <- "New or continued medications/supplements prescribed at this visit  (choice=Topirimate/phentermine for weight loss)"
label(data$cv_newmeds___28) <- "New or continued medications/supplements prescribed at this visit  (choice=Liraglutiade (saxcenda) for weight loss)"
label(data$cv_newmeds___29) <- "New or continued medications/supplements prescribed at this visit  (choice=Semaglutide (Wygovy/Ozempic) for weight loss)"
label(data$cv_newmeds___30) <- "New or continued medications/supplements prescribed at this visit  (choice=Levothyroxine)"
label(data$cv_newmeds___31) <- "New or continued medications/supplements prescribed at this visit  (choice=Systemic steroids (ex: prednisone, dexamethasone, solumedrol))"
label(data$cv_newmeds___60) <- "New or continued medications/supplements prescribed at this visit  (choice=Other)"
label(data$cv_newmeds___na) <- "New or continued medications/supplements prescribed at this visit  (choice=None)"
label(data$cv_newmeds___unk) <- "New or continued medications/supplements prescribed at this visit  (choice=Unknown)"
label(data$cv_newmeds___oth) <- "New or continued medications/supplements prescribed at this visit  (choice=Other)"
label(data$cv_newmeds___pm) <- "New or continued medications/supplements prescribed at this visit  (choice=Premenarchal)"
label(data$cv_newmeds_diabetes) <- "Secondary diabetes meds"
label(data$cv_newcont_metforminformulation) <- "Metformin formulation"
label(data$cv_newcont_metformindose) <- "Metformin total daily dose"
label(data$cv_newcont_metformin_otherdose) <- "Metformin other dose"
label(data$cv_newcont_metforminfreq) <- "Metformin frequency"
label(data$cv_referrals___1) <- "Specialty referral from this appointment (choice=Gynecology)"
label(data$cv_referrals___2) <- "Specialty referral from this appointment (choice=Endocrinology)"
label(data$cv_referrals___3) <- "Specialty referral from this appointment (choice=Dermatology)"
label(data$cv_referrals___4) <- "Specialty referral from this appointment (choice=Cardiology)"
label(data$cv_referrals___5) <- "Specialty referral from this appointment (choice=Gastroenterology)"
label(data$cv_referrals___6) <- "Specialty referral from this appointment (choice=Sleep or pulmonary medicine)"
label(data$cv_referrals___7) <- "Specialty referral from this appointment (choice=Sleep study)"
label(data$cv_referrals___8) <- "Specialty referral from this appointment (choice=Bariatric surgery)"
label(data$cv_referrals___9) <- "Specialty referral from this appointment (choice=Obesity management program)"
label(data$cv_referrals___10) <- "Specialty referral from this appointment (choice=Sleep study)"
label(data$cv_referrals___11) <- "Specialty referral from this appointment (choice=Psychological evaluation)"
label(data$cv_referrals___12) <- "Specialty referral from this appointment (choice=Adolescent Medicine)"
label(data$cv_referrals___60) <- "Specialty referral from this appointment (choice=Other {cv_referrals_other})"
label(data$cv_referrals___na) <- "Specialty referral from this appointment (choice=None)"
label(data$cv_referrals___unk) <- "Specialty referral from this appointment (choice=Unknown/not recorded)"
label(data$cv_referrals___oth) <- "Specialty referral from this appointment (choice=Other)"
label(data$cv_referrals___pm) <- "Specialty referral from this appointment (choice=Premenarchal)"
label(data$cv_referrals_other) <- "Other specialty referral type"
label(data$cv_referrals_reason) <- "Reason for specialty referral"
label(data$cv_dietarycounseling) <- "Received dietary counseling/goals at visit"
label(data$cv_exerciseplan) <- "Received an exercise plan at visit"
label(data$history_complete) <- "Complete?"
label(data$clinical_visit_complete) <- "Complete?"
# Setting Factors(will create new variable for factors)
data$redcap_repeat_instrument.factor <- factor(data$redcap_repeat_instrument, levels = c("clinical_visit"))
data$pcosdx_pmh___1.factor <- factor(data$pcosdx_pmh___1, levels = c("0", "1"))
data$pcosdx_pmh___2.factor <- factor(data$pcosdx_pmh___2, levels = c("0", "1"))
data$pcosdx_pmh___3.factor <- factor(data$pcosdx_pmh___3, levels = c("0", "1"))
data$pcosdx_pmh___4.factor <- factor(data$pcosdx_pmh___4, levels = c("0", "1"))
data$pcosdx_pmh___5.factor <- factor(data$pcosdx_pmh___5, levels = c("0", "1"))
data$pcosdx_pmh___6.factor <- factor(data$pcosdx_pmh___6, levels = c("0", "1"))
data$pcosdx_pmh___29.factor <- factor(data$pcosdx_pmh___29, levels = c("0", "1"))
data$pcosdx_pmh___7.factor <- factor(data$pcosdx_pmh___7, levels = c("0", "1"))
data$pcosdx_pmh___8.factor <- factor(data$pcosdx_pmh___8, levels = c("0", "1"))
data$pcosdx_pmh___9.factor <- factor(data$pcosdx_pmh___9, levels = c("0", "1"))
data$pcosdx_pmh___10.factor <- factor(data$pcosdx_pmh___10, levels = c("0", "1"))
data$pcosdx_pmh___11.factor <- factor(data$pcosdx_pmh___11, levels = c("0", "1"))
data$pcosdx_pmh___12.factor <- factor(data$pcosdx_pmh___12, levels = c("0", "1"))
data$pcosdx_pmh___22.factor <- factor(data$pcosdx_pmh___22, levels = c("0", "1"))
data$pcosdx_pmh___13.factor <- factor(data$pcosdx_pmh___13, levels = c("0", "1"))
data$pcosdx_pmh___14.factor <- factor(data$pcosdx_pmh___14, levels = c("0", "1"))
data$pcosdx_pmh___15.factor <- factor(data$pcosdx_pmh___15, levels = c("0", "1"))
data$pcosdx_pmh___16.factor <- factor(data$pcosdx_pmh___16, levels = c("0", "1"))
data$pcosdx_pmh___23.factor <- factor(data$pcosdx_pmh___23, levels = c("0", "1"))
data$pcosdx_pmh___24.factor <- factor(data$pcosdx_pmh___24, levels = c("0", "1"))
data$pcosdx_pmh___17.factor <- factor(data$pcosdx_pmh___17, levels = c("0", "1"))
data$pcosdx_pmh___18.factor <- factor(data$pcosdx_pmh___18, levels = c("0", "1"))
data$pcosdx_pmh___19.factor <- factor(data$pcosdx_pmh___19, levels = c("0", "1"))
data$pcosdx_pmh___20.factor <- factor(data$pcosdx_pmh___20, levels = c("0", "1"))
data$pcosdx_pmh___21.factor <- factor(data$pcosdx_pmh___21, levels = c("0", "1"))
data$pcosdx_pmh___60.factor <- factor(data$pcosdx_pmh___60, levels = c("0", "1"))
data$pcosdx_pmh___0.factor <- factor(data$pcosdx_pmh___0, levels = c("0", "1"))
data$pcosdx_pmh___unk.factor <- factor(data$pcosdx_pmh___unk, levels = c("0", "1"))
data$pcosdx_pmh___na.factor <- factor(data$pcosdx_pmh___na, levels = c("0", "1"))
data$pcosdx_pmh___oth.factor <- factor(data$pcosdx_pmh___oth, levels = c("0", "1"))
data$pcosdx_pmh___pm.factor <- factor(data$pcosdx_pmh___pm, levels = c("0", "1"))
data$site.factor <- factor(data$site, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"))
data$race___1.factor <- factor(data$race___1, levels = c("0", "1"))
data$race___2.factor <- factor(data$race___2, levels = c("0", "1"))
data$race___3.factor <- factor(data$race___3, levels = c("0", "1"))
data$race___4.factor <- factor(data$race___4, levels = c("0", "1"))
data$race___5.factor <- factor(data$race___5, levels = c("0", "1"))
data$race___60.factor <- factor(data$race___60, levels = c("0", "1"))
data$race___unk.factor <- factor(data$race___unk, levels = c("0", "1"))
data$race___na.factor <- factor(data$race___na, levels = c("0", "1"))
data$race___oth.factor <- factor(data$race___oth, levels = c("0", "1"))
data$race___pm.factor <- factor(data$race___pm, levels = c("0", "1"))
data$ethnicity.factor <- factor(data$ethnicity, levels = c("0", "1", "UNK"))
data$insur_type.factor <- factor(data$insur_type, levels = c("1", "2", "3", "0", "60", "UNK"))
data$pastmeds___1.factor <- factor(data$pastmeds___1, levels = c("0", "1"))
data$pastmeds___2.factor <- factor(data$pastmeds___2, levels = c("0", "1"))
data$pastmeds___3.factor <- factor(data$pastmeds___3, levels = c("0", "1"))
data$pastmeds___4.factor <- factor(data$pastmeds___4, levels = c("0", "1"))
data$pastmeds___5.factor <- factor(data$pastmeds___5, levels = c("0", "1"))
data$pastmeds___6.factor <- factor(data$pastmeds___6, levels = c("0", "1"))
data$pastmeds___7.factor <- factor(data$pastmeds___7, levels = c("0", "1"))
data$pastmeds___8.factor <- factor(data$pastmeds___8, levels = c("0", "1"))
data$pastmeds___9.factor <- factor(data$pastmeds___9, levels = c("0", "1"))
data$pastmeds___10.factor <- factor(data$pastmeds___10, levels = c("0", "1"))
data$pastmeds___11.factor <- factor(data$pastmeds___11, levels = c("0", "1"))
data$pastmeds___12.factor <- factor(data$pastmeds___12, levels = c("0", "1"))
data$pastmeds___13.factor <- factor(data$pastmeds___13, levels = c("0", "1"))
data$pastmeds___14.factor <- factor(data$pastmeds___14, levels = c("0", "1"))
data$pastmeds___15.factor <- factor(data$pastmeds___15, levels = c("0", "1"))
data$pastmeds___16.factor <- factor(data$pastmeds___16, levels = c("0", "1"))
data$pastmeds___60.factor <- factor(data$pastmeds___60, levels = c("0", "1"))
data$pastmeds___0.factor <- factor(data$pastmeds___0, levels = c("0", "1"))
data$pastmeds___unk.factor <- factor(data$pastmeds___unk, levels = c("0", "1"))
data$pastmeds___na.factor <- factor(data$pastmeds___na, levels = c("0", "1"))
data$pastmeds___oth.factor <- factor(data$pastmeds___oth, levels = c("0", "1"))
data$pastmeds___pm.factor <- factor(data$pastmeds___pm, levels = c("0", "1"))
data$pcosdx_mentalhealthcounseling___1.factor <- factor(data$pcosdx_mentalhealthcounseling___1, levels = c("0", "1"))
data$pcosdx_mentalhealthcounseling___2.factor <- factor(data$pcosdx_mentalhealthcounseling___2, levels = c("0", "1"))
data$pcosdx_mentalhealthcounseling___3.factor <- factor(data$pcosdx_mentalhealthcounseling___3, levels = c("0", "1"))
data$pcosdx_mentalhealthcounseling___0.factor <- factor(data$pcosdx_mentalhealthcounseling___0, levels = c("0", "1"))
data$pcosdx_mentalhealthcounseling___unk.factor <- factor(data$pcosdx_mentalhealthcounseling___unk, levels = c("0", "1"))
data$pcosdx_mentalhealthcounseling___na.factor <- factor(data$pcosdx_mentalhealthcounseling___na, levels = c("0", "1"))
data$pcosdx_mentalhealthcounseling___oth.factor <- factor(data$pcosdx_mentalhealthcounseling___oth, levels = c("0", "1"))
data$pcosdx_mentalhealthcounseling___pm.factor <- factor(data$pcosdx_mentalhealthcounseling___pm, levels = c("0", "1"))
data$pcosdx_dietarycounseling___1.factor <- factor(data$pcosdx_dietarycounseling___1, levels = c("0", "1"))
data$pcosdx_dietarycounseling___2.factor <- factor(data$pcosdx_dietarycounseling___2, levels = c("0", "1"))
data$pcosdx_dietarycounseling___0.factor <- factor(data$pcosdx_dietarycounseling___0, levels = c("0", "1"))
data$pcosdx_dietarycounseling___unk.factor <- factor(data$pcosdx_dietarycounseling___unk, levels = c("0", "1"))
data$pcosdx_dietarycounseling___na.factor <- factor(data$pcosdx_dietarycounseling___na, levels = c("0", "1"))
data$pcosdx_dietarycounseling___oth.factor <- factor(data$pcosdx_dietarycounseling___oth, levels = c("0", "1"))
data$pcosdx_dietarycounseling___pm.factor <- factor(data$pcosdx_dietarycounseling___pm, levels = c("0", "1"))
data$pcosdx_exercising_teachingplan___1.factor <- factor(data$pcosdx_exercising_teachingplan___1, levels = c("0", "1"))
data$pcosdx_exercising_teachingplan___2.factor <- factor(data$pcosdx_exercising_teachingplan___2, levels = c("0", "1"))
data$pcosdx_exercising_teachingplan___0.factor <- factor(data$pcosdx_exercising_teachingplan___0, levels = c("0", "1"))
data$pcosdx_exercising_teachingplan___unk.factor <- factor(data$pcosdx_exercising_teachingplan___unk, levels = c("0", "1"))
data$pcosdx_exercising_teachingplan___na.factor <- factor(data$pcosdx_exercising_teachingplan___na, levels = c("0", "1"))
data$pcosdx_exercising_teachingplan___oth.factor <- factor(data$pcosdx_exercising_teachingplan___oth, levels = c("0", "1"))
data$pcosdx_exercising_teachingplan___pm.factor <- factor(data$pcosdx_exercising_teachingplan___pm, levels = c("0", "1"))
data$pcosdx_irregular_menses.factor <- factor(data$pcosdx_irregular_menses, levels = c("1", "2", "4", "5", "6", "UNK"))
data$pcosdx_hyperandrogenism.factor <- factor(data$pcosdx_hyperandrogenism, levels = c("1", "2", "3", "UNK"))
data$cv_type_of_clinic.factor <- factor(data$cv_type_of_clinic, levels = c("1", "2", "4", "5", "6", "3"))
data$cv_mdc_specialties___1.factor <- factor(data$cv_mdc_specialties___1, levels = c("0", "1"))
data$cv_mdc_specialties___2.factor <- factor(data$cv_mdc_specialties___2, levels = c("0", "1"))
data$cv_mdc_specialties___3.factor <- factor(data$cv_mdc_specialties___3, levels = c("0", "1"))
data$cv_mdc_specialties___4.factor <- factor(data$cv_mdc_specialties___4, levels = c("0", "1"))
data$cv_mdc_specialties___5.factor <- factor(data$cv_mdc_specialties___5, levels = c("0", "1"))
data$cv_mdc_specialties___6.factor <- factor(data$cv_mdc_specialties___6, levels = c("0", "1"))
data$cv_mdc_specialties___60.factor <- factor(data$cv_mdc_specialties___60, levels = c("0", "1"))
data$cv_mdc_specialties___unk.factor <- factor(data$cv_mdc_specialties___unk, levels = c("0", "1"))
data$cv_mdc_specialties___na.factor <- factor(data$cv_mdc_specialties___na, levels = c("0", "1"))
data$cv_mdc_specialties___oth.factor <- factor(data$cv_mdc_specialties___oth, levels = c("0", "1"))
data$cv_mdc_specialties___pm.factor <- factor(data$cv_mdc_specialties___pm, levels = c("0", "1"))
data$cv_medications___1.factor <- factor(data$cv_medications___1, levels = c("0", "1"))
data$cv_medications___2.factor <- factor(data$cv_medications___2, levels = c("0", "1"))
data$cv_medications___3.factor <- factor(data$cv_medications___3, levels = c("0", "1"))
data$cv_medications___4.factor <- factor(data$cv_medications___4, levels = c("0", "1"))
data$cv_medications___5.factor <- factor(data$cv_medications___5, levels = c("0", "1"))
data$cv_medications___6.factor <- factor(data$cv_medications___6, levels = c("0", "1"))
data$cv_medications___7.factor <- factor(data$cv_medications___7, levels = c("0", "1"))
data$cv_medications___8.factor <- factor(data$cv_medications___8, levels = c("0", "1"))
data$cv_medications___9.factor <- factor(data$cv_medications___9, levels = c("0", "1"))
data$cv_medications___10.factor <- factor(data$cv_medications___10, levels = c("0", "1"))
data$cv_medications___11.factor <- factor(data$cv_medications___11, levels = c("0", "1"))
data$cv_medications___12.factor <- factor(data$cv_medications___12, levels = c("0", "1"))
data$cv_medications___13.factor <- factor(data$cv_medications___13, levels = c("0", "1"))
data$cv_medications___14.factor <- factor(data$cv_medications___14, levels = c("0", "1"))
data$cv_medications___15.factor <- factor(data$cv_medications___15, levels = c("0", "1"))
data$cv_medications___16.factor <- factor(data$cv_medications___16, levels = c("0", "1"))
data$cv_medications___17.factor <- factor(data$cv_medications___17, levels = c("0", "1"))
data$cv_medications___18.factor <- factor(data$cv_medications___18, levels = c("0", "1"))
data$cv_medications___19.factor <- factor(data$cv_medications___19, levels = c("0", "1"))
data$cv_medications___20.factor <- factor(data$cv_medications___20, levels = c("0", "1"))
data$cv_medications___21.factor <- factor(data$cv_medications___21, levels = c("0", "1"))
data$cv_medications___22.factor <- factor(data$cv_medications___22, levels = c("0", "1"))
data$cv_medications___23.factor <- factor(data$cv_medications___23, levels = c("0", "1"))
data$cv_medications___32.factor <- factor(data$cv_medications___32, levels = c("0", "1"))
data$cv_medications___24.factor <- factor(data$cv_medications___24, levels = c("0", "1"))
data$cv_medications___25.factor <- factor(data$cv_medications___25, levels = c("0", "1"))
data$cv_medications___26.factor <- factor(data$cv_medications___26, levels = c("0", "1"))
data$cv_medications___27.factor <- factor(data$cv_medications___27, levels = c("0", "1"))
data$cv_medications___28.factor <- factor(data$cv_medications___28, levels = c("0", "1"))
data$cv_medications___29.factor <- factor(data$cv_medications___29, levels = c("0", "1"))
data$cv_medications___30.factor <- factor(data$cv_medications___30, levels = c("0", "1"))
data$cv_medications___31.factor <- factor(data$cv_medications___31, levels = c("0", "1"))
data$cv_medications___60.factor <- factor(data$cv_medications___60, levels = c("0", "1"))
data$cv_medications___0.factor <- factor(data$cv_medications___0, levels = c("0", "1"))
data$cv_medications___unk.factor <- factor(data$cv_medications___unk, levels = c("0", "1"))
data$cv_medications___na.factor <- factor(data$cv_medications___na, levels = c("0", "1"))
data$cv_medications___oth.factor <- factor(data$cv_medications___oth, levels = c("0", "1"))
data$cv_medications___pm.factor <- factor(data$cv_medications___pm, levels = c("0", "1"))
data$cv_meds_diabetes2.factor <- factor(data$cv_meds_diabetes2, levels = c("1", "2", "3"))
data$cv_genderid.factor <- factor(data$cv_genderid, levels = c("1", "2", "3", "60", "UNK"))
data$cv_hirsutism_cat.factor <- factor(data$cv_hirsutism_cat, levels = c("0", "1", "2", "3", "4", "NA", "UNK"))
data$cv_acneface.factor <- factor(data$cv_acneface, levels = c("1", "2", "3", "0", "4", "UNK"))
data$cv_acneother___1.factor <- factor(data$cv_acneother___1, levels = c("0", "1"))
data$cv_acneother___2.factor <- factor(data$cv_acneother___2, levels = c("0", "1"))
data$cv_acneother___0.factor <- factor(data$cv_acneother___0, levels = c("0", "1"))
data$cv_acneother___unk.factor <- factor(data$cv_acneother___unk, levels = c("0", "1"))
data$cv_acneother___na.factor <- factor(data$cv_acneother___na, levels = c("0", "1"))
data$cv_acneother___oth.factor <- factor(data$cv_acneother___oth, levels = c("0", "1"))
data$cv_acneother___pm.factor <- factor(data$cv_acneother___pm, levels = c("0", "1"))
data$cv_acanthosisneck.factor <- factor(data$cv_acanthosisneck, levels = c("1", "2", "3", "4", "0", "5", "UNK"))
data$cv_hydradinitis.factor <- factor(data$cv_hydradinitis, levels = c("1", "0", "UNK"))
data$cv_alopecia.factor <- factor(data$cv_alopecia, levels = c("1", "0", "UNK"))
data$cv_tt_assay.factor <- factor(data$cv_tt_assay, levels = c("1", "OTH", "UNK"))
data$cv_mood.factor <- factor(data$cv_mood, levels = c("1", "3", "0", "2"))
data$cv_newmeds___1.factor <- factor(data$cv_newmeds___1, levels = c("0", "1"))
data$cv_newmeds___2.factor <- factor(data$cv_newmeds___2, levels = c("0", "1"))
data$cv_newmeds___3.factor <- factor(data$cv_newmeds___3, levels = c("0", "1"))
data$cv_newmeds___4.factor <- factor(data$cv_newmeds___4, levels = c("0", "1"))
data$cv_newmeds___5.factor <- factor(data$cv_newmeds___5, levels = c("0", "1"))
data$cv_newmeds___6.factor <- factor(data$cv_newmeds___6, levels = c("0", "1"))
data$cv_newmeds___7.factor <- factor(data$cv_newmeds___7, levels = c("0", "1"))
data$cv_newmeds___8.factor <- factor(data$cv_newmeds___8, levels = c("0", "1"))
data$cv_newmeds___9.factor <- factor(data$cv_newmeds___9, levels = c("0", "1"))
data$cv_newmeds___10.factor <- factor(data$cv_newmeds___10, levels = c("0", "1"))
data$cv_newmeds___11.factor <- factor(data$cv_newmeds___11, levels = c("0", "1"))
data$cv_newmeds___12.factor <- factor(data$cv_newmeds___12, levels = c("0", "1"))
data$cv_newmeds___13.factor <- factor(data$cv_newmeds___13, levels = c("0", "1"))
data$cv_newmeds___14.factor <- factor(data$cv_newmeds___14, levels = c("0", "1"))
data$cv_newmeds___15.factor <- factor(data$cv_newmeds___15, levels = c("0", "1"))
data$cv_newmeds___16.factor <- factor(data$cv_newmeds___16, levels = c("0", "1"))
data$cv_newmeds___17.factor <- factor(data$cv_newmeds___17, levels = c("0", "1"))
data$cv_newmeds___18.factor <- factor(data$cv_newmeds___18, levels = c("0", "1"))
data$cv_newmeds___19.factor <- factor(data$cv_newmeds___19, levels = c("0", "1"))
data$cv_newmeds___20.factor <- factor(data$cv_newmeds___20, levels = c("0", "1"))
data$cv_newmeds___21.factor <- factor(data$cv_newmeds___21, levels = c("0", "1"))
data$cv_newmeds___22.factor <- factor(data$cv_newmeds___22, levels = c("0", "1"))
data$cv_newmeds___23.factor <- factor(data$cv_newmeds___23, levels = c("0", "1"))
data$cv_newmeds___32.factor <- factor(data$cv_newmeds___32, levels = c("0", "1"))
data$cv_newmeds___24.factor <- factor(data$cv_newmeds___24, levels = c("0", "1"))
data$cv_newmeds___25.factor <- factor(data$cv_newmeds___25, levels = c("0", "1"))
data$cv_newmeds___26.factor <- factor(data$cv_newmeds___26, levels = c("0", "1"))
data$cv_newmeds___27.factor <- factor(data$cv_newmeds___27, levels = c("0", "1"))
data$cv_newmeds___28.factor <- factor(data$cv_newmeds___28, levels = c("0", "1"))
data$cv_newmeds___29.factor <- factor(data$cv_newmeds___29, levels = c("0", "1"))
data$cv_newmeds___30.factor <- factor(data$cv_newmeds___30, levels = c("0", "1"))
data$cv_newmeds___31.factor <- factor(data$cv_newmeds___31, levels = c("0", "1"))
data$cv_newmeds___60.factor <- factor(data$cv_newmeds___60, levels = c("0", "1"))
data$cv_newmeds___na.factor <- factor(data$cv_newmeds___na, levels = c("0", "1"))
data$cv_newmeds___unk.factor <- factor(data$cv_newmeds___unk, levels = c("0", "1"))
data$cv_newmeds___oth.factor <- factor(data$cv_newmeds___oth, levels = c("0", "1"))
data$cv_newmeds___pm.factor <- factor(data$cv_newmeds___pm, levels = c("0", "1"))
data$cv_newmeds_diabetes.factor <- factor(data$cv_newmeds_diabetes, levels = c("1", "2", "3"))
data$cv_newcont_metforminformulation.factor <- factor(data$cv_newcont_metforminformulation, levels = c("1", "2", "3", "60"))
data$cv_newcont_metformindose.factor <- factor(data$cv_newcont_metformindose, levels = c("500", "750", "850", "1000", "1500", "1700", "2000", "60"))
data$cv_newcont_metforminfreq.factor <- factor(data$cv_newcont_metforminfreq, levels = c("1", "2", "3"))
data$cv_referrals___1.factor <- factor(data$cv_referrals___1, levels = c("0", "1"))
data$cv_referrals___2.factor <- factor(data$cv_referrals___2, levels = c("0", "1"))
data$cv_referrals___3.factor <- factor(data$cv_referrals___3, levels = c("0", "1"))
data$cv_referrals___4.factor <- factor(data$cv_referrals___4, levels = c("0", "1"))
data$cv_referrals___5.factor <- factor(data$cv_referrals___5, levels = c("0", "1"))
data$cv_referrals___6.factor <- factor(data$cv_referrals___6, levels = c("0", "1"))
data$cv_referrals___7.factor <- factor(data$cv_referrals___7, levels = c("0", "1"))
data$cv_referrals___8.factor <- factor(data$cv_referrals___8, levels = c("0", "1"))
data$cv_referrals___9.factor <- factor(data$cv_referrals___9, levels = c("0", "1"))
data$cv_referrals___10.factor <- factor(data$cv_referrals___10, levels = c("0", "1"))
data$cv_referrals___11.factor <- factor(data$cv_referrals___11, levels = c("0", "1"))
data$cv_referrals___12.factor <- factor(data$cv_referrals___12, levels = c("0", "1"))
data$cv_referrals___60.factor <- factor(data$cv_referrals___60, levels = c("0", "1"))
data$cv_referrals___na.factor <- factor(data$cv_referrals___na, levels = c("0", "1"))
data$cv_referrals___unk.factor <- factor(data$cv_referrals___unk, levels = c("0", "1"))
data$cv_referrals___oth.factor <- factor(data$cv_referrals___oth, levels = c("0", "1"))
data$cv_referrals___pm.factor <- factor(data$cv_referrals___pm, levels = c("0", "1"))
data$cv_dietarycounseling.factor <- factor(data$cv_dietarycounseling, levels = c("1", "0", "2", "UNK"))
data$cv_exerciseplan.factor <- factor(data$cv_exerciseplan, levels = c("1", "0", "2", "UNK"))
data$history_complete.factor <- factor(data$history_complete, levels = c("0", "1", "2"))
data$clinical_visit_complete.factor <- factor(data$clinical_visit_complete, levels = c("0", "1", "2"))

levels(data$redcap_repeat_instrument.factor) <- c("Clinical Visit")
levels(data$pcosdx_pmh___1.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___2.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___3.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___4.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___5.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___6.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___29.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___7.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___8.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___9.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___10.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___11.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___12.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___22.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___13.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___14.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___15.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___16.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___23.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___24.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___17.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___18.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___19.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___20.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___21.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___60.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___0.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___unk.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___na.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___oth.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_pmh___pm.factor) <- c("Unchecked", "Checked")
levels(data$site.factor) <- c("Childrens Hospital Colorado, Denver (DEN)", "NIH/Childrens National (DC)", "Childrens Hospital of Philadelphia (CHOP)", "PENN", "Childrens Hospital of Pittsburgh (PIT)", "University of Florida, Gainesville (UFPCOS)", "Mercy Childrens Hospital-Kansas City (KC)", "Childrens Hospital Los Angeles (L)", "Cook County Hospital, Chicago (CHI)", "New York/Columbia (COL)", "Northwestern/Lurie Childrens (LC)", "Boston Childrens (BCH)", "University of Alabama (UAB)", "Los Angeles /Cedars Sinai (CSMC)", "Virginia Commonwealth University (VCU)", "University of Virginia (UVA)", "Nationwide Childrens Hospital, Ohio (NCH)", "Vanderbilt (VAN)", "UChicago Medicine (UCM)", "Texas A&M (TAMU)", "Driscoll Childrens Hospital (DCH)")
levels(data$race___1.factor) <- c("Unchecked", "Checked")
levels(data$race___2.factor) <- c("Unchecked", "Checked")
levels(data$race___3.factor) <- c("Unchecked", "Checked")
levels(data$race___4.factor) <- c("Unchecked", "Checked")
levels(data$race___5.factor) <- c("Unchecked", "Checked")
levels(data$race___60.factor) <- c("Unchecked", "Checked")
levels(data$race___unk.factor) <- c("Unchecked", "Checked")
levels(data$race___na.factor) <- c("Unchecked", "Checked")
levels(data$race___oth.factor) <- c("Unchecked", "Checked")
levels(data$race___pm.factor) <- c("Unchecked", "Checked")
levels(data$ethnicity.factor) <- c("Non-Hispanic", "Hispanic", "Unknown")
levels(data$insur_type.factor) <- c("Public", "Private", "Military", "None", "Other {insur_other}", "Unknown/Not recorded")
levels(data$pastmeds___1.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___2.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___3.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___4.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___5.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___6.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___7.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___8.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___9.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___10.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___11.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___12.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___13.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___14.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___15.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___16.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___60.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___0.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___unk.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___na.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___oth.factor) <- c("Unchecked", "Checked")
levels(data$pastmeds___pm.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_mentalhealthcounseling___1.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_mentalhealthcounseling___2.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_mentalhealthcounseling___3.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_mentalhealthcounseling___0.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_mentalhealthcounseling___unk.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_mentalhealthcounseling___na.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_mentalhealthcounseling___oth.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_mentalhealthcounseling___pm.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_dietarycounseling___1.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_dietarycounseling___2.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_dietarycounseling___0.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_dietarycounseling___unk.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_dietarycounseling___na.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_dietarycounseling___oth.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_dietarycounseling___pm.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_exercising_teachingplan___1.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_exercising_teachingplan___2.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_exercising_teachingplan___0.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_exercising_teachingplan___unk.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_exercising_teachingplan___na.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_exercising_teachingplan___oth.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_exercising_teachingplan___pm.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_irregular_menses.factor) <- c(">1 year post-menarche: > 90 days for any one cycle", ">1 year post-menarche: > or = 6 week cycles (< or = 8 cycles a year)", "Primary amenorrhea by age 15", "Primary amenorrhea > 3 years post-thelarche", "Excess non-ovulatory menses", "Unknown/Not recorded")
levels(data$pcosdx_hyperandrogenism.factor) <- c("Only clinical", "Only biochemical", "Both clinical and biochemical", "Unknown/Not recorded")
levels(data$cv_type_of_clinic.factor) <- c("Endocrine-only appointment", "Multidisciplinary clinic", "Adolescent Medicine only", "Gynecology only", "Pediatrics/Family Medicine", "Other {clinic_other}")
levels(data$cv_mdc_specialties___1.factor) <- c("Unchecked", "Checked")
levels(data$cv_mdc_specialties___2.factor) <- c("Unchecked", "Checked")
levels(data$cv_mdc_specialties___3.factor) <- c("Unchecked", "Checked")
levels(data$cv_mdc_specialties___4.factor) <- c("Unchecked", "Checked")
levels(data$cv_mdc_specialties___5.factor) <- c("Unchecked", "Checked")
levels(data$cv_mdc_specialties___6.factor) <- c("Unchecked", "Checked")
levels(data$cv_mdc_specialties___60.factor) <- c("Unchecked", "Checked")
levels(data$cv_mdc_specialties___unk.factor) <- c("Unchecked", "Checked")
levels(data$cv_mdc_specialties___na.factor) <- c("Unchecked", "Checked")
levels(data$cv_mdc_specialties___oth.factor) <- c("Unchecked", "Checked")
levels(data$cv_mdc_specialties___pm.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___1.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___2.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___3.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___4.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___5.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___6.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___7.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___8.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___9.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___10.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___11.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___12.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___13.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___14.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___15.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___16.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___17.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___18.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___19.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___20.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___21.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___22.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___23.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___32.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___24.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___25.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___26.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___27.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___28.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___29.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___30.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___31.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___60.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___0.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___unk.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___na.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___oth.factor) <- c("Unchecked", "Checked")
levels(data$cv_medications___pm.factor) <- c("Unchecked", "Checked")
levels(data$cv_meds_diabetes2.factor) <- c("GLP1", "SGLT2", "TZD")
levels(data$cv_genderid.factor) <- c("Female (cisgender)", "Male (transgender)", "Nonbinary", "Other {cv_genderid_other}", "Unknown/Not recorded")
levels(data$cv_hirsutism_cat.factor) <- c("None (0-6)", "Mild (7-9)", "Moderate (10-15)", "Severe (>15)", "Noted by clinician, severity unknown.", "Not applicable, patient had permanent hair removal (i.e. laser hair removal)", "Unknown/Not recorded")
levels(data$cv_acneface.factor) <- c("Mild", "Moderate (pustular)", "Severe (nodular, cystic, scarring)", "None", "Acne noted by clinician, severity unknown.", "Unknown/Not recorded")
levels(data$cv_acneother___1.factor) <- c("Unchecked", "Checked")
levels(data$cv_acneother___2.factor) <- c("Unchecked", "Checked")
levels(data$cv_acneother___0.factor) <- c("Unchecked", "Checked")
levels(data$cv_acneother___unk.factor) <- c("Unchecked", "Checked")
levels(data$cv_acneother___na.factor) <- c("Unchecked", "Checked")
levels(data$cv_acneother___oth.factor) <- c("Unchecked", "Checked")
levels(data$cv_acneother___pm.factor) <- c("Unchecked", "Checked")
levels(data$cv_acanthosisneck.factor) <- c("Back of neck, barely visible (minimal)", "Back of neck, obvious (mild)", "Lateral sides of neck (moderate)", "Circumferential (severe)", "None", "Noted by clinician, severity unknown.", "Unknown/Not recorded")
levels(data$cv_hydradinitis.factor) <- c("Yes (abscess or red lesions in axilla, pannus or groin)", "No", "Unknown/Not recorded")
levels(data$cv_alopecia.factor) <- c("Yes", "No", "Unknown/Not recorded")
levels(data$cv_tt_assay.factor) <- c("LCMSMS", "Other {cv_tt_assay_oth}", "Not recorded/unknown")
levels(data$cv_mood.factor) <- c("Yes", "Yes, asked by clinician, but no scale used.", "No", "No, but referred to a therapist")
levels(data$cv_newmeds___1.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___2.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___3.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___4.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___5.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___6.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___7.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___8.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___9.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___10.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___11.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___12.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___13.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___14.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___15.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___16.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___17.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___18.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___19.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___20.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___21.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___22.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___23.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___32.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___24.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___25.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___26.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___27.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___28.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___29.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___30.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___31.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___60.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___na.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___unk.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___oth.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds___pm.factor) <- c("Unchecked", "Checked")
levels(data$cv_newmeds_diabetes.factor) <- c("GLP1", "SGLT2", "TZD")
levels(data$cv_newcont_metforminformulation.factor) <- c("Standard release", "Extended release", "Riomet (liquid)", "Other {cv_newcont_metformin_other}")
levels(data$cv_newcont_metformindose.factor) <- c("500", "750", "850", "1000", "1500", "1700", "2000", "Other dose {cv_newcont_metformin_otherdose}")
levels(data$cv_newcont_metforminfreq.factor) <- c("1", "2", "3")
levels(data$cv_referrals___1.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___2.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___3.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___4.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___5.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___6.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___7.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___8.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___9.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___10.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___11.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___12.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___60.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___na.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___unk.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___oth.factor) <- c("Unchecked", "Checked")
levels(data$cv_referrals___pm.factor) <- c("Unchecked", "Checked")
levels(data$cv_dietarycounseling.factor) <- c("Yes", "No", "No, but referred to a dietitian", "Unknown/not recorded")
levels(data$cv_exerciseplan.factor) <- c("Yes", "No", "No, but referred to a exercise specialist", "Unknown/not recorded")
levels(data$history_complete.factor) <- c("Incomplete", "Unverified", "Complete")
levels(data$clinical_visit_complete.factor) <- c("Incomplete", "Unverified", "Complete")
