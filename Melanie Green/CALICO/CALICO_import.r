# Read Data
data <- read.csv("/Users/timvigers/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Vigers/BDC/Janet Snell-Bergeon/CALICO/Data_Raw/CALICOMulticenterPCO_DATA_2024-10-01_0840.csv", na.strings = "")
# Setting Labels
label(data$record_number) <- "Record number"
label(data$redcap_repeat_instrument) <- "Repeat Instrument"
label(data$redcap_repeat_instance) <- "Repeat Instance"
label(data$redcap_data_access_group) <- "Data Access Group"
label(data$site) <- "Study Site"
label(data$data_entryname) <- "Person entering data"
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
label(data$pcosdx_famhx_parent___27) <- "Parent history (include all history, not only from the diagnostic visit) (choice=PCOS)"
label(data$pcosdx_famhx_parent___28) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Reproductive cancer)"
label(data$pcosdx_famhx_parent___30) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Endometriosis)"
label(data$pcosdx_famhx_parent___5) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Overweight/Obesity)"
label(data$pcosdx_famhx_parent___29) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Type 1 diabetes)"
label(data$pcosdx_famhx_parent___7) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Type 2 diabetes)"
label(data$pcosdx_famhx_parent___24) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Gestational Diabetes)"
label(data$pcosdx_famhx_parent___8) <- "Parent history (include all history, not only from the diagnostic visit) (choice=OSA)"
label(data$pcosdx_famhx_parent___25) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Bariatric surgery)"
label(data$pcosdx_famhx_parent___12) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Nonalcoholic fatty liver disease)"
label(data$pcosdx_famhx_parent___22) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Hypothyroidism)"
label(data$pcosdx_famhx_parent___13) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Hypertension)"
label(data$pcosdx_famhx_parent___14) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Hyperlipidemia)"
label(data$pcosdx_famhx_parent___26) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Myocardial infarction prior to age 50)"
label(data$pcosdx_famhx_parent___15) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Anxiety)"
label(data$pcosdx_famhx_parent___16) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Depression)"
label(data$pcosdx_famhx_parent___17) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Binge eating disorder (including bulemia))"
label(data$pcosdx_famhx_parent___18) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Restricting eating disorder (including anorexia))"
label(data$pcosdx_famhx_parent___31) <- "Parent history (include all history, not only from the diagnostic visit) (choice=IBS)"
label(data$pcosdx_famhx_parent___32) <- "Parent history (include all history, not only from the diagnostic visit) (choice=POTS)"
label(data$pcosdx_famhx_parent___33) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Chronic fatigue)"
label(data$pcosdx_famhx_parent___60) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Other {pcosdx_famhx_parent_otherdx})"
label(data$pcosdx_famhx_parent___0) <- "Parent history (include all history, not only from the diagnostic visit) (choice=None)"
label(data$pcosdx_famhx_parent___unk) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Unknown/Not recorded)"
label(data$pcosdx_famhx_parent___na) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Not applicable)"
label(data$pcosdx_famhx_parent___oth) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Other)"
label(data$pcosdx_famhx_parent___pm) <- "Parent history (include all history, not only from the diagnostic visit) (choice=Premenarchal)"
label(data$pcosdx_famhx_parent_otherdx) <- "Other parent diagnosis"
label(data$pcosdx_famhx___27) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=PCOS)"
label(data$pcosdx_famhx___28) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Reproductive cancer)"
label(data$pcosdx_famhx___30) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Endometriosis)"
label(data$pcosdx_famhx___5) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Overweight/Obesity)"
label(data$pcosdx_famhx___29) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Type 1 diabetes)"
label(data$pcosdx_famhx___7) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Type 2 diabetes)"
label(data$pcosdx_famhx___24) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Gestational Diabetes)"
label(data$pcosdx_famhx___8) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=OSA)"
label(data$pcosdx_famhx___25) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Bariatric surgery)"
label(data$pcosdx_famhx___12) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Nonalcoholic fatty liver disease)"
label(data$pcosdx_famhx___22) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Hypothyroidism)"
label(data$pcosdx_famhx___13) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Hypertension)"
label(data$pcosdx_famhx___14) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Hyperlipidemia)"
label(data$pcosdx_famhx___26) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Myocardial infarction prior to age 50)"
label(data$pcosdx_famhx___15) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Anxiety)"
label(data$pcosdx_famhx___16) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Depression)"
label(data$pcosdx_famhx___17) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Binge eating disorder (including bulemia))"
label(data$pcosdx_famhx___18) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Restricting eating disorder (including anorexia))"
label(data$pcosdx_famhx___31) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=IBS)"
label(data$pcosdx_famhx___32) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=POTS)"
label(data$pcosdx_famhx___33) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Chronic fatigue)"
label(data$pcosdx_famhx___60) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Other {pcosdx_famhx_family_otherdx})"
label(data$pcosdx_famhx___0) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=None)"
label(data$pcosdx_famhx___unk) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Unknown/Not recorded)"
label(data$pcosdx_famhx___na) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Not applicable)"
label(data$pcosdx_famhx___oth) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Other)"
label(data$pcosdx_famhx___pm) <- "Non-parent family history (include all history, not only from the diagnostic visit). Including siblings, first cousins, aunts, uncles, and grandparents. (choice=Premenarchal)"
label(data$pcosdx_famhx_pcos_specify___1) <- "Please indicate which non-parent family member has PCOS (choice=Sister)"
label(data$pcosdx_famhx_pcos_specify___2) <- "Please indicate which non-parent family member has PCOS (choice=Maternal other family member)"
label(data$pcosdx_famhx_pcos_specify___3) <- "Please indicate which non-parent family member has PCOS (choice=Paternal other family member)"
label(data$pcosdx_famhx_pcos_specify___unk) <- "Please indicate which non-parent family member has PCOS (choice=Unknown/Not recorded)"
label(data$pcosdx_famhx_pcos_specify___na) <- "Please indicate which non-parent family member has PCOS (choice=Not applicable)"
label(data$pcosdx_famhx_pcos_specify___oth) <- "Please indicate which non-parent family member has PCOS (choice=Other)"
label(data$pcosdx_famhx_pcos_specify___pm) <- "Please indicate which non-parent family member has PCOS (choice=Premenarchal)"
label(data$pcosdx_famhx_family_otherdx) <- "Other non-parent family diagnosis"
label(data$pcosdx_birthweight) <- "Birthweight (lb)"
label(data$pcosdx_birthweight_oz) <- "Birthweight (oz)"
label(data$pcosdx_birthweight_calc) <- "Birthweight"
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
label(data$pcosdx_padx_age_noted) <- "Age premature adrenarche first noted by family"
label(data$pcosdx_padx_age) <- "Age when premature adrenarche diagnosed (years)"
label(data$pcosdx_padx_age_appt) <- "Age at diagnostic appointment for premature adrenarche"
label(data$pcosdx_padx_body_odor) <- "Body odor present at diagnostic appointment for premature adrenarche"
label(data$pcosdx_padx_axillary_hair) <- "Axillary hair present at diagnostic appointment for premature adrenarche"
label(data$pcosdx_padx_tanner_pubic) <- "Pubic hair Tanner stage at diagnostic appointment for premature adrenarche"
label(data$pcosdx_padx_tanner_breast) <- "Breast Tanner stage at diagnostic appointment for premature adrenarche"
label(data$pcosdx_padx_acne) <- "Acne present at diagnostic appointment for premature adrenarche"
label(data$pcosdx_padx_acne_severity) <- "Acne severity at diagnostic appointment for premature adrenarche"
label(data$pcosdx_padx_dheas) <- "DHEAS (µg/dL) at diagnosis of premature adrenarche"
label(data$pcosdx_padx_dheas_uln) <- "DHEAS (µg/dL) value for upper limit of normal for sex, age/tanner stage for assay "
label(data$pcosdx_padx_tt) <- "Total Testosterone (ng/dL) at diagnosis of premature adrenarche"
label(data$pcosdx_padx_tt_uln) <- "Total Testosterone (ng/dL) value for upper limit of normal for sex, age/tanner stage for assay"
label(data$pcosdx_padx_ft) <- "Free Testosterone (pg/mL) at diagnosis of premature adrenarche"
label(data$pcosdx_padx_ft_uln) <- "Free Testosterone (pg/mL) value for upper limit of normal for sex, age/tanner stage for assay"
label(data$pcosdx_padx_bmi) <- "BMI at diagnosis of premature adrenarche"
label(data$pcosdx_padx_bmi_perc) <- "BMI percentile at diagnosis of premature adrenarche"
label(data$pcosdx_padx_bone_age) <- "Bone age "
label(data$pcosdx_padx_age_bone_age) <- "Age at bone age"
label(data$pcosdx_cppdx_age) <- "Age when central precocious puberty diagnosed (years)"
label(data$pcosdx_pregdx_age) <- "Age when pregnant (years)"
label(data$pcosdx_cysticacne) <- "Age when cystic acne diagnosed (years)"
label(data$pcosdx_obesitydx_age) <- "Age when overweight/obesity diagnosed (years)"
label(data$pcosdx_obesitydx_earliestage) <- "If answered UNK for the previous question, when was the earliest age that overweight/obesity was noted? "
label(data$pcosdx_prediabetesdx_age) <- "Age when pre-diabetes diagnosed (years)"
label(data$pcosdx_t1ddx_age) <- "Age when type 1 diabetes diagnosed (years)"
label(data$pcosdx_t2ddx_age) <- "Age when type 2 diabetes diagnosed (years)"
label(data$pcosdx_osadx_age) <- "Age OSA diagnosed (years)"
label(data$pcosdx_cpapdx_age) <- "Age when CPAP started being used (years)"
label(data$pcosdx_pseudodx_age) <- "Age when pseudotumor diagnosed"
label(data$pcosdx_gerddx_age) <- "Age when GERD diagnosed (years)"
label(data$pcosdx_flddx_age) <- "Age when nonalcoholic fatty liver disease diagnosed (years)"
label(data$pcosdx_htndx_age) <- "Age when hypertension diagnosed (years)"
label(data$pcosdx_hlddx_age) <- "Age when hyperlipidemia diagnosed (years)"
label(data$pcosdx_anxietydx_age) <- "Age when anxiety diagnosed (years)"
label(data$pcosdx_depressiondx_age) <- "Age when depression diagnosed (years)"
label(data$pcosdx_siattempt_age) <- "Age when suicide attempt occurred (years)"
label(data$pcosdx_mhhosp_age) <- "Age when mental health hospitalization occurred (years)"
label(data$pcosdx_bingedx_age) <- "Age when binge eating noted (years)"
label(data$pcosdx_restrictdx_age) <- "Age when restrictive eating noted (years)"
label(data$pcosdx_adhddx_age) <- "Age when ADHD diagnosed (years)"
label(data$pcosdx_asthmadx_age) <- "Age when asthma diagnosed (years)"
label(data$pcosdx_migrainesdx_age) <- "Age when migraines diagnosed (years)"
label(data$pcosdx_hypothydx_age) <- "Age when hypothyroidism diagnosed (years)"
label(data$pcosdx_other) <- "Other diagnosis, please specify age of diagnosis (years) if available"
label(data$pcosdx_surghx___1) <- "Surgical history (at time of PCOS diagnosis) (choice=Tonsillectomy)"
label(data$pcosdx_surghx___2) <- "Surgical history (at time of PCOS diagnosis) (choice=Adenoidectomy)"
label(data$pcosdx_surghx___3) <- "Surgical history (at time of PCOS diagnosis) (choice=Bariatric surgery)"
label(data$pcosdx_surghx___4) <- "Surgical history (at time of PCOS diagnosis) (choice=Pilonidal cyst)"
label(data$pcosdx_surghx___5) <- "Surgical history (at time of PCOS diagnosis) (choice=Incision and drainage for Hidradenitis)"
label(data$pcosdx_surghx___6) <- "Surgical history (at time of PCOS diagnosis) (choice=Cholycystectomy)"
label(data$pcosdx_surghx___7) <- "Surgical history (at time of PCOS diagnosis) (choice=Liver biopsy)"
label(data$pcosdx_surghx___8) <- "Surgical history (at time of PCOS diagnosis) (choice=Ovarian cyst removal)"
label(data$pcosdx_surghx___9) <- "Surgical history (at time of PCOS diagnosis) (choice=Ovarian torsion)"
label(data$pcosdx_surghx___10) <- "Surgical history (at time of PCOS diagnosis) (choice=Pelvic exploratory laparotemy)"
label(data$pcosdx_surghx___60) <- "Surgical history (at time of PCOS diagnosis) (choice=Other {surghx_other})"
label(data$pcosdx_surghx___0) <- "Surgical history (at time of PCOS diagnosis) (choice=None)"
label(data$pcosdx_surghx___unk) <- "Surgical history (at time of PCOS diagnosis) (choice=Unknown/Not recorded)"
label(data$pcosdx_surghx___na) <- "Surgical history (at time of PCOS diagnosis) (choice=Not applicable)"
label(data$pcosdx_surghx___oth) <- "Surgical history (at time of PCOS diagnosis) (choice=Other)"
label(data$pcosdx_surghx___pm) <- "Surgical history (at time of PCOS diagnosis) (choice=Premenarchal)"
label(data$pcosdx_surghx_tonsil) <- "Age at tonsillectomy (years)"
label(data$pcosdx_surghx_adenoid) <- "Age at adenoidectomy (years)"
label(data$pcosdx_surghx_bariatric) <- "Age at bariatric surgery (years)"
label(data$pcosdx_surghx_pilonidal) <- "Age at pilonidal cyst surgery (years)"
label(data$pcosdx_surghx_hidradenitis) <- "Age at incision and drainage for Hidradenitis (years)"
label(data$pcosdx_surghx_cholycyst) <- "Age at cholycystectomy (years)"
label(data$pcosdx_surghx_liverbiopsy) <- "Age at liver biopsy (years)"
label(data$pcosdx_surghx_ovarycyst) <- "Age at ovarian cyst removal (years)"
label(data$pcosdx_surghx_ovarytorsion) <- "Age at ovarian torsion surgery (years)"
label(data$pcosdx_surghx_pelvicexplor) <- "Age at pelvic exploratory laparotemy (years)"
label(data$pcosdx_surghx_otherage) <- "Other surgery, please specify age of surgery (years) if available"
label(data$surghx_other) <- "Other surgery"
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
label(data$pcosdx_medcomments) <- "Medication Comments"
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
label(data$pcosdx_schooltype) <- "Type of school at time of PCOS diagnosis"
label(data$schooltype_other) <- "other"
label(data$pcosdx_schoolability) <- "School performance at time of PCOS diagnosis"
label(data$pcosdx_learning) <- "Learning disability or individual educational plan (IEP) at time of PCOS diagnosis"
label(data$pcosdx_age) <- "Age at time of PCOS diagnosis"
label(data$pcosdx_menarche) <- "Menarche age (i.e. 10 years and 6 months = 10.5)"
label(data$pcosdx_menarchetime) <- "Gynecologic age at PCOS diagnosis (months since menarche)"
label(data$pcosdx_specialty) <- "What medical specialty made the diagnosis of PCOS"
label(data$pcosdx_specialty_other) <- "Other medical specialty that made the diagnosis of PCOS"
label(data$pcosdx_irregular_menses) <- "Irregular menses defined"
label(data$pcosdx_hyperandrogenism) <- "Hyperandrogenism diagnostic method"
label(data$pcosdx_lmp_monthsago) <- "Last menstrual period"
label(data$pcosdx_lmp_fromdx) <- "Months between LMP and PCOS diagnosis"
label(data$pcosdx_menfreq) <- "Approximate menses in the last 12 months"
label(data$pcosdx_us) <- "Pelvic ultrasound"
label(data$pcosdx_us_age) <- "Age during pelvic ultrasound (yrs)"
label(data$pcosdx_us_rt_follicle) <- "Right ovary follicle count"
label(data$pcosdx_us_rt_size) <- "Right ovary volume (mL)"
label(data$pcosdx_us_lt_follicle) <- "Left ovary follicle count"
label(data$pcosdx_us_lt_size) <- "Left ovary volume (mL)"
label(data$pcosdx_us_uterussize) <- "Uterus volume (mL)"
label(data$pcosdx_us_uterussize3d) <- "Uterus size (cm x cm x cm)"
label(data$pcosdx_us_stripe) <- "Thickness of endometrial stripe (mm)"
label(data$pcosdx_us_comments) <- "Ultrasound comments"
label(data$history_complete) <- "Complete?"
label(data$cv_data_entry_personnel) <- "Data entered by: "
label(data$cv_visittype) <- "Type of Clinical Visit"
label(data$cv_type_of_clinic) <- "Type of clinic visit"
label(data$clinic_other) <- "Clinic Other"
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
label(data$cv_mdc_other) <- "other"
label(data$cv_monthssincepcosdx) <- "Months since PCOS diagnosis"
label(data$cv_age) <- "Age at follow-up"
label(data$cv_newdx___3) <- "New diagnosis since last visit and including any made this visit (choice=Pregnancy)"
label(data$cv_newdx___4) <- "New diagnosis since last visit and including any made this visit (choice=Cystic acne)"
label(data$cv_newdx___5) <- "New diagnosis since last visit and including any made this visit (choice=Overweight/Obesity)"
label(data$cv_newdx___6) <- "New diagnosis since last visit and including any made this visit (choice=Pre-diabetes)"
label(data$cv_newdx___29) <- "New diagnosis since last visit and including any made this visit (choice=Type 1 diabetes)"
label(data$cv_newdx___7) <- "New diagnosis since last visit and including any made this visit (choice=Type 2 diabetes)"
label(data$cv_newdx___8) <- "New diagnosis since last visit and including any made this visit (choice=Obstructive sleep apnea)"
label(data$cv_newdx___9) <- "New diagnosis since last visit and including any made this visit (choice=Use of continuous positive airway pressure (CPAP))"
label(data$cv_newdx___10) <- "New diagnosis since last visit and including any made this visit (choice=Pseudotumor)"
label(data$cv_newdx___11) <- "New diagnosis since last visit and including any made this visit (choice=Gastroesophageal Reflux Disease)"
label(data$cv_newdx___12) <- "New diagnosis since last visit and including any made this visit (choice=Nonalcoholic fatty liver disease)"
label(data$cv_newdx___22) <- "New diagnosis since last visit and including any made this visit (choice=Hypothyroidism)"
label(data$cv_newdx___13) <- "New diagnosis since last visit and including any made this visit (choice=Hypertension)"
label(data$cv_newdx___14) <- "New diagnosis since last visit and including any made this visit (choice=Hyperlipidemia)"
label(data$cv_newdx___15) <- "New diagnosis since last visit and including any made this visit (choice=Anxiety)"
label(data$cv_newdx___16) <- "New diagnosis since last visit and including any made this visit (choice=Depression)"
label(data$cv_newdx___27) <- "New diagnosis since last visit and including any made this visit (choice=Suicide attempt)"
label(data$cv_newdx___28) <- "New diagnosis since last visit and including any made this visit (choice=Mental health hospitalization)"
label(data$cv_newdx___17) <- "New diagnosis since last visit and including any made this visit (choice=Binge eating disorder (including bulimia))"
label(data$cv_newdx___18) <- "New diagnosis since last visit and including any made this visit (choice=Restricting eating disorder (including anorexia))"
label(data$cv_newdx___19) <- "New diagnosis since last visit and including any made this visit (choice=Attention deficit hyperactivity disorder (ADHD or ADD))"
label(data$cv_newdx___20) <- "New diagnosis since last visit and including any made this visit (choice=Asthma)"
label(data$cv_newdx___21) <- "New diagnosis since last visit and including any made this visit (choice=Migraines)"
label(data$cv_newdx___23) <- "New diagnosis since last visit and including any made this visit (choice=Endometriosis)"
label(data$cv_newdx___24) <- "New diagnosis since last visit and including any made this visit (choice=Irritable bowel syndrome)"
label(data$cv_newdx___25) <- "New diagnosis since last visit and including any made this visit (choice=Postural and orthostatic tachycardic syndrome)"
label(data$cv_newdx___26) <- "New diagnosis since last visit and including any made this visit (choice=Chronic fatigue)"
label(data$cv_newdx___60) <- "New diagnosis since last visit and including any made this visit (choice=Other {cv_othermedical})"
label(data$cv_newdx___0) <- "New diagnosis since last visit and including any made this visit (choice=None)"
label(data$cv_newdx___unk) <- "New diagnosis since last visit and including any made this visit (choice=Unknown/Not recorded)"
label(data$cv_newdx___na) <- "New diagnosis since last visit and including any made this visit (choice=Not applicable)"
label(data$cv_newdx___oth) <- "New diagnosis since last visit and including any made this visit (choice=Other)"
label(data$cv_newdx___pm) <- "New diagnosis since last visit and including any made this visit (choice=Premenarchal)"
label(data$cv_pregdx_age) <- "Months since PCOS dx when pregnant"
label(data$cv_cysticacnedx_age) <- "Months since PCOS dx when cystic acne diagnosed"
label(data$cv_obesitydx_age) <- "Months since PCOS dx when obesity/overweight diagnosed"
label(data$cv_prediabetesdx_age) <- "Months since PCOS dx when pre-diabetes diagnosed"
label(data$cv_t1ddx_date) <- "Date type 1 diabetes diagnosed"
label(data$cv_t1ddx_age) <- "Months since PCOS dx when type 1 diabetes diagnosed"
label(data$cv_t2ddx_age) <- "Months since PCOS dx when T2D diagnosed"
label(data$cv_osadx_age) <- "Months since PCOS dx when OSA diagnosed "
label(data$cv_gerddx_age) <- "Months since PCOS dx when GERD diagnosed"
label(data$cv_cpapdx_age) <- "Months since PCOS dx when CPAP started being used"
label(data$cv_pseudodx_age) <- "Months since PCOS dx when pseudotumor diagnosed (years)"
label(data$cv_flddx_age) <- "Months since PCOS dx when nonalcoholic fatty liver disease diagnosed"
label(data$cv_htndx_age) <- "Months since PCOS dx when hypertension diagnosed"
label(data$cv_hlddx_age) <- "Months since PCOS dx when hyperlipidemia diagnosed"
label(data$cv_anxietydx_age) <- "Months since PCOS dx when anxiety diagnosed"
label(data$cv_depressiondx_age) <- "Months since PCOS dx when depression diagnosed"
label(data$cv_si_attempt_age) <- "Months since PCOS dx when suicide attempt occurred"
label(data$cv_mh_hospital_age) <- "Months since PCOS dx when mental health hospitalization occurred"
label(data$cv_bingedx_age) <- "Months since PCOS dx when binge eating noted "
label(data$cv_restrictdx_age) <- "Months since PCOS dx when restrictive eating noted"
label(data$cv_adhddx_age) <- "Months since PCOS dx when ADHD diagnosed "
label(data$cv_asthmadx_age) <- "Months since PCOS dx when asthma diagnosed"
label(data$cv_migrainesdx_age) <- "Months since PCOS dx when migraines diagnosed "
label(data$cv_hypothydx_age) <- "Months since PCOS dx when hypothyroidism diagnosed "
label(data$cv_endometriosisdx_date) <- "Date confirmed endometriosis"
label(data$cv_endometriosisdx_age) <- "Months since PCOS dx when endometriosis diagnosed"
label(data$cv_ibsdx_date) <- "Date confirmed IBS"
label(data$cv_ibsdx_age) <- "Months since PCOS dx when IBS diagnosed"
label(data$cv_potsdx_date) <- "Date confirmed POTS"
label(data$cv_potsdx_age) <- "Months since PCOS dx when POTS diagnosed"
label(data$cv_fatiguedx_date) <- "Date confirmed chronic fatigue diagnosis"
label(data$cv_fatiguedx_age) <- "Months since PCOS dx when chronic fatigue diagnosed"
label(data$cv_othermedical) <- "Other diagnosis"
label(data$cv_other_age_2) <- "Months since PCOS dx when [cv_othermedical] diagnosed "
label(data$cv_surghx___1) <- "Surgery since PCOS diagnosis (choice=Tonsillectomy)"
label(data$cv_surghx___2) <- "Surgery since PCOS diagnosis (choice=Adenoidectomy)"
label(data$cv_surghx___3) <- "Surgery since PCOS diagnosis (choice=Bariatric surgery)"
label(data$cv_surghx___4) <- "Surgery since PCOS diagnosis (choice=Pilonidal cyst)"
label(data$cv_surghx___5) <- "Surgery since PCOS diagnosis (choice=Incision and drainage for Hidradenitis)"
label(data$cv_surghx___6) <- "Surgery since PCOS diagnosis (choice=Cholycystectomy)"
label(data$cv_surghx___7) <- "Surgery since PCOS diagnosis (choice=Liver biopsy)"
label(data$cv_surghx___8) <- "Surgery since PCOS diagnosis (choice=Ovarian cyst removal)"
label(data$cv_surghx___9) <- "Surgery since PCOS diagnosis (choice=Ovarian torsion)"
label(data$cv_surghx___10) <- "Surgery since PCOS diagnosis (choice=Pelvic exploratory laparotemy)"
label(data$cv_surghx___60) <- "Surgery since PCOS diagnosis (choice=Other {cv_surghx_other})"
label(data$cv_surghx___0) <- "Surgery since PCOS diagnosis (choice=None)"
label(data$cv_surghx___unk) <- "Surgery since PCOS diagnosis (choice=Unknown/Not recorded)"
label(data$cv_surghx___na) <- "Surgery since PCOS diagnosis (choice=Not applicable)"
label(data$cv_surghx___oth) <- "Surgery since PCOS diagnosis (choice=Other)"
label(data$cv_surghx___pm) <- "Surgery since PCOS diagnosis (choice=Premenarchal)"
label(data$cv_surghx_tonsil) <- "Age at tonsillectomy (years)"
label(data$cv_surghx_adenoid) <- "Age at adenoidectomy (years)"
label(data$cv_surghx_bariatric) <- "Age at bariatric surgery (years)"
label(data$cv_surghx_pilonidal) <- "Age at pilonidal cyst surgery (years)"
label(data$cv_surghx_hidradenitis) <- "Age at incision and drainage for Hidradenitis (years)"
label(data$cv_surghx_cholycyst) <- "Age at cholycystectomy (years)"
label(data$cv_surghx_liverbiopsy) <- "Age at liver biopsy (years)"
label(data$cv_surghx_ovarycyst) <- "Age at ovarian cyst removal (years)"
label(data$cv_surghx_ovarytorsion) <- "Age at ovarian torsion surgery (years)"
label(data$cv_surghx_pelvicexplor) <- "Age at pelvic exploratory laparotomy (years)"
label(data$cv_surghx_other) <- "Other surgery since PCOS diagnosis (include age if known)"
label(data$cv_surghx_other_age) <- "Age at [cv_surghx_other] (years)"
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
label(data$cv_meds_pill) <- "Brand of pill"
label(data$cv_meds_diabetes2) <- "Secondary diabetes meds"
label(data$cv_reflux_meds) <- "Reflux meds"
label(data$cv_sleep_meds) <- "Sleep aids"
label(data$cv_medscomments) <- "Medication comments"
label(data$cv_othercomment) <- "Other medications or supplements (list the other medications prescribed and provide dose, if available)"
label(data$cv_hormonecomments) <- "Hormone therapy comments"
label(data$cv_metforminformulation) <- "Metformin formulation"
label(data$metformin_form_other) <- "Other"
label(data$cv_metformintdd) <- "Metformin total daily dose "
label(data$cv_meds_metformin_otherdose) <- "Other metformin dose"
label(data$cv_metforminfreq) <- "Metformin frequency (per day)"
label(data$cv_misseddoses) <- "How many doses of metformin per week does the patient report missing?"
label(data$cv_gisymptoms) <- "Reports GI side effects with metformin?"
label(data$cv_sideeffectsreported___1) <- "What GI side effects from metformin do they report? (choice=Nausea)"
label(data$cv_sideeffectsreported___2) <- "What GI side effects from metformin do they report? (choice=Diarrhea)"
label(data$cv_sideeffectsreported___3) <- "What GI side effects from metformin do they report? (choice=Cramps)"
label(data$cv_sideeffectsreported___4) <- "What GI side effects from metformin do they report? (choice=Retching)"
label(data$cv_sideeffectsreported___5) <- "What GI side effects from metformin do they report? (choice=Bloating)"
label(data$cv_sideeffectsreported___6) <- "What GI side effects from metformin do they report? (choice=Constipation)"
label(data$cv_sideeffectsreported___7) <- "What GI side effects from metformin do they report? (choice=Not specified)"
label(data$cv_sideeffectsreported___60) <- "What GI side effects from metformin do they report? (choice=Other {cv_other_symptoms})"
label(data$cv_sideeffectsreported___unk) <- "What GI side effects from metformin do they report? (choice=Unknown/Not recorded)"
label(data$cv_sideeffectsreported___na) <- "What GI side effects from metformin do they report? (choice=Not applicable)"
label(data$cv_sideeffectsreported___oth) <- "What GI side effects from metformin do they report? (choice=Other)"
label(data$cv_sideeffectsreported___pm) <- "What GI side effects from metformin do they report? (choice=Premenarchal)"
label(data$cv_other_symptoms) <- "Other GI side effects from metformin"
label(data$cv_meds_adhd_name) <- "Name of ADHD medication "
label(data$cv_weight) <- "Weight (kg)"
label(data$cv_height) <- "Height (cm)"
label(data$cv_bmi) <- "BMI"
label(data$cv_bmi_percentile) <- "BMI percentile"
label(data$cv_bmi_z) <- "BMI z-score"
label(data$cv_waist) <- "Waist circumference (cm) "
label(data$cv_hip) <- "Hip circumference (cm) "
label(data$cv_wh_ratio) <- "Waist-to-hip ratio"
label(data$cv_sbp) <- "Systolic BP (mmHg) "
label(data$cv_dbp) <- "Diastolic BP (mmHg) "
label(data$cv_hr) <- "Heart rate (bpm)"
label(data$cv_menfreq) <- "Menses in last 6 months"
label(data$cv_sexuallyactive) <- "Sexually active at any time prior to this visit"
label(data$cv_genderid) <- "Patients gender identity"
label(data$cv_genderid_other) <- "Other"
label(data$cv_sexualpref) <- "Sexual preference"
label(data$cv_sexualpref_other) <- "Other"
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
label(data$labs_this_visit___1) <- "Labs/imaging collected this visit: (choice=Pregnancy test)"
label(data$labs_this_visit___2) <- "Labs/imaging collected this visit: (choice=Free testosterone)"
label(data$labs_this_visit___3) <- "Labs/imaging collected this visit: (choice=Total testosterone)"
label(data$labs_this_visit___4) <- "Labs/imaging collected this visit: (choice=SHBG)"
label(data$labs_this_visit___5) <- "Labs/imaging collected this visit: (choice=DHEAS)"
label(data$labs_this_visit___6) <- "Labs/imaging collected this visit: (choice=Androstendione)"
label(data$labs_this_visit___7) <- "Labs/imaging collected this visit: (choice=AMH)"
label(data$labs_this_visit___8) <- "Labs/imaging collected this visit: (choice=17OHP)"
label(data$labs_this_visit___9) <- "Labs/imaging collected this visit: (choice=Estradiol)"
label(data$labs_this_visit___10) <- "Labs/imaging collected this visit: (choice=LH)"
label(data$labs_this_visit___11) <- "Labs/imaging collected this visit: (choice=FSH)"
label(data$labs_this_visit___12) <- "Labs/imaging collected this visit: (choice=Prolactin)"
label(data$labs_this_visit___13) <- "Labs/imaging collected this visit: (choice=TSH)"
label(data$labs_this_visit___14) <- "Labs/imaging collected this visit: (choice=Total T4)"
label(data$labs_this_visit___15) <- "Labs/imaging collected this visit: (choice=Free T4)"
label(data$labs_this_visit___16) <- "Labs/imaging collected this visit: (choice=HbA1C)"
label(data$labs_this_visit___17) <- "Labs/imaging collected this visit: (choice=Fasting glucose)"
label(data$labs_this_visit___18) <- "Labs/imaging collected this visit: (choice=2 hour glucose from OGTT)"
label(data$labs_this_visit___19) <- "Labs/imaging collected this visit: (choice=Fasting insulin)"
label(data$labs_this_visit___20) <- "Labs/imaging collected this visit: (choice=Hemoglobin)"
label(data$labs_this_visit___21) <- "Labs/imaging collected this visit: (choice=Platelets)"
label(data$labs_this_visit___22) <- "Labs/imaging collected this visit: (choice=Ferritin)"
label(data$labs_this_visit___23) <- "Labs/imaging collected this visit: (choice=25-Hydroxy Vitamin D)"
label(data$labs_this_visit___24) <- "Labs/imaging collected this visit: (choice=Triglycerides)"
label(data$labs_this_visit___25) <- "Labs/imaging collected this visit: (choice=HDL)"
label(data$labs_this_visit___26) <- "Labs/imaging collected this visit: (choice=LDL)"
label(data$labs_this_visit___27) <- "Labs/imaging collected this visit: (choice=Total cholesterol)"
label(data$labs_this_visit___28) <- "Labs/imaging collected this visit: (choice=ALT)"
label(data$labs_this_visit___29) <- "Labs/imaging collected this visit: (choice=AST)"
label(data$labs_this_visit___31) <- "Labs/imaging collected this visit: (choice=BUN)"
label(data$labs_this_visit___32) <- "Labs/imaging collected this visit: (choice=Creatinine)"
label(data$labs_this_visit___30) <- "Labs/imaging collected this visit: (choice=Liver ultrasound/MRI results)"
label(data$labs_this_visit___33) <- "Labs/imaging collected this visit: (choice=Pelvic ultrasound)"
label(data$labs_this_visit___na) <- "Labs/imaging collected this visit: (choice=No lab ordered this visit)"
label(data$labs_this_visit___unk) <- "Labs/imaging collected this visit: (choice=Unknown)"
label(data$labs_this_visit___oth) <- "Labs/imaging collected this visit: (choice=Other)"
label(data$labs_this_visit___pm) <- "Labs/imaging collected this visit: (choice=Premenarchal)"
label(data$cv_upt_result) <- "Pregnancy test result"
label(data$cv_ft) <- "Free testosterone (pg/mL) "
label(data$cv_ft_uln) <- "Free Testosterone (pg/mL) value for upper limit of normal for sex, age/tanner stage  for assay"
label(data$cv_ft_perc) <- "Free testosterone (pg/mL) percentage of upper limit of normal"
label(data$cv_tt) <- "Total testosterone (ng/dL) "
label(data$cv_tt_assay) <- "Assay for total testosterone"
label(data$cv_tt_assay_oth) <- "Assay other"
label(data$cv_tt_uln) <- "Total testosterone (ng/dL) value for upper limit of normal for assay"
label(data$cv_tt_perc) <- "Total testosterone percentage of upper limit of normal"
label(data$cv_shbg) <- "SHBG (nmol/L)"
label(data$cv_dheas) <- "DHEAS (µg/dL)"
label(data$cv_dheas_uln) <- "DHEAS value for upper limit of normal for assay (µg/dL)"
label(data$cv_dheas_perc) <- "DHEAS percentage of upper limit of normal"
label(data$cv_androstendione) <- "Androstendione (ng/dL)"
label(data$cv_androstendione_uln) <- "Androstendione upper limit of normal for assay (ng/dL)"
label(data$cv_androstendione_perc) <- "Androstendione percentage of upper limit of normal"
label(data$cv_amh) <- "AMH (ng/mL)"
label(data$cv_amh_uln) <- "AMH (ng/mL) value for upper limit of normal for assay"
label(data$cv_amh_perc) <- "AMH (ng/mL) percentage of upper limit of normal"
label(data$cv_17ohp) <- "17OHP (ng/dL)"
label(data$cv_17ohp_uln) <- "17OHP for upper limit of normal for assay (ng/dL)"
label(data$cv_17ohp_perc) <- "17OHP percentage of upper limit of normal"
label(data$cv_estradiol) <- "Estradiol"
label(data$cv_estradiol_uln) <- "Estradiol upper limit of normal for assay "
label(data$cv_estradiol_perc) <- "Estradiol percentage of upper limit of normal"
label(data$cv_lh) <- "LH (mIU/mL)"
label(data$cv_fsh) <- "FSH (mIU/mL)"
label(data$cv_prolactin) <- "Prolactin (ng/mL)"
label(data$cv_tsh) <- "TSH (mIU/L)"
label(data$cv_total_t4) <- "Total T4 (µg/dL)"
label(data$cv_free_t4) <- "Free T4 (ng/dL)"
label(data$cv_a1c) <- "HbA1C (%)"
label(data$cv_fbg) <- "Fasting glucose (mg/dL)"
label(data$cv_2hrglucoseogtt) <- "2 hour glucose from OGTT (mg/dL)"
label(data$cv_fastinsulin) <- "Fasting insulin (mIU/mL)"
label(data$cv_hgb) <- "Hemoglobin (g/dL)"
label(data$cv_platelets) <- "Platelets (10^9 cells/L)"
label(data$cv_ferritin) <- "Ferritin (µg/L)"
label(data$cv_25ohd) <- "25-Hydroxy Vitamin D (ng/mL)"
label(data$cv_tg_fasting) <- "Conditions of triglyceride draw"
label(data$cv_tg) <- "Triglycerides (mg/dL)"
label(data$cv_hdl) <- "HDL (mg/dL)"
label(data$cv_ldl) <- "LDL (mg/dL)"
label(data$cv_tc) <- "Total cholesterol (mg/dL)"
label(data$cv_alt) <- "ALT (U/L)"
label(data$cv_ast) <- "AST (U/L)"
label(data$cv_bun) <- "BUN (mg/dL)"
label(data$cv_creatinine) <- "Creatinine (mg/dL)"
label(data$cv_liverimaging_results) <- "Liver ultrasound/MRI results"
label(data$cv_pelvicus_age) <- "Age during pelvic ultrasound (years)"
label(data$cv_pelvicus_rtcount) <- "Right ovary follicle count"
label(data$cv_pelvicus_rtvol) <- "Right ovary volume (mL)"
label(data$cv_pelvicus_ltcount) <- "Left ovary follicle count"
label(data$cv_pelvicus_ltvol) <- "Left ovary volume (mL)"
label(data$cv_pelvicus_uterussize) <- "Uterus volume (mL)"
label(data$cv_pelvicus_uterussize3d) <- "Uterus size (cm x cm x cm)"
label(data$cv_pelvicus_stripe) <- "Thickness of endometrial stripe (mm)"
label(data$cv_pelvicus_comments) <- "Ultrasound Comments and Clinical Indication"
label(data$lab_comments) <- "Lab work/imaging comments"
label(data$cv_osa_sx) <- "Symptoms of sleep apnea (snoring, frequent awakening, AM headaches, daytime fatigue, napping during day)"
label(data$cv_osa_comments) <- "Sleep apnea symptoms comments"
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
label(data$cv_newmeds_pill) <- "Brand of pill"
label(data$cv_newmeds_diabetes) <- "Secondary diabetes meds"
label(data$cv_reflux_newmeds) <- "Reflux meds"
label(data$cv_sleep_newmeds) <- "Sleep aids"
label(data$cv_newcont_meds_othercomment) <- "New or continued medication/supplement comments (list the other medications prescribed and provide dose)"
label(data$cv_newcont_meds_hormonecomments) <- "Hormone therapy comments"
label(data$cv_newcont_meds_discontinued) <- "If any medications were discontinued, please specify which medication and the reason for discontinuation"
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
label(data$cv_fitnesstesting_yn) <- "Was fitness testing performed?"
label(data$cv_fitness_sitreach) <- "Sit and reach score"
label(data$cv_fitness_step) <- "3-minute Step Test score"
label(data$cv_fitness_dynamometer) <- "Dynamometer score"
label(data$cv_exerciseplan) <- "Received an exercise plan at visit"
label(data$clinical_visit_complete) <- "Complete?"
# Setting Factors(will create new variable for factors)
data$redcap_repeat_instrument.factor <- factor(data$redcap_repeat_instrument, levels = c("clinical_visit"))
data$redcap_data_access_group.factor <- factor(data$redcap_data_access_group, levels = c("bch", "chi", "chla", "chop", "csmc", "dc", "dch", "den", "kc", "lc", "nch", "nycol", "penn", "pitt", "tamu", "uab", "ucm", "uf", "uva", "van", "vcu"))
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
data$pcosdx_famhx_parent___27.factor <- factor(data$pcosdx_famhx_parent___27, levels = c("0", "1"))
data$pcosdx_famhx_parent___28.factor <- factor(data$pcosdx_famhx_parent___28, levels = c("0", "1"))
data$pcosdx_famhx_parent___30.factor <- factor(data$pcosdx_famhx_parent___30, levels = c("0", "1"))
data$pcosdx_famhx_parent___5.factor <- factor(data$pcosdx_famhx_parent___5, levels = c("0", "1"))
data$pcosdx_famhx_parent___29.factor <- factor(data$pcosdx_famhx_parent___29, levels = c("0", "1"))
data$pcosdx_famhx_parent___7.factor <- factor(data$pcosdx_famhx_parent___7, levels = c("0", "1"))
data$pcosdx_famhx_parent___24.factor <- factor(data$pcosdx_famhx_parent___24, levels = c("0", "1"))
data$pcosdx_famhx_parent___8.factor <- factor(data$pcosdx_famhx_parent___8, levels = c("0", "1"))
data$pcosdx_famhx_parent___25.factor <- factor(data$pcosdx_famhx_parent___25, levels = c("0", "1"))
data$pcosdx_famhx_parent___12.factor <- factor(data$pcosdx_famhx_parent___12, levels = c("0", "1"))
data$pcosdx_famhx_parent___22.factor <- factor(data$pcosdx_famhx_parent___22, levels = c("0", "1"))
data$pcosdx_famhx_parent___13.factor <- factor(data$pcosdx_famhx_parent___13, levels = c("0", "1"))
data$pcosdx_famhx_parent___14.factor <- factor(data$pcosdx_famhx_parent___14, levels = c("0", "1"))
data$pcosdx_famhx_parent___26.factor <- factor(data$pcosdx_famhx_parent___26, levels = c("0", "1"))
data$pcosdx_famhx_parent___15.factor <- factor(data$pcosdx_famhx_parent___15, levels = c("0", "1"))
data$pcosdx_famhx_parent___16.factor <- factor(data$pcosdx_famhx_parent___16, levels = c("0", "1"))
data$pcosdx_famhx_parent___17.factor <- factor(data$pcosdx_famhx_parent___17, levels = c("0", "1"))
data$pcosdx_famhx_parent___18.factor <- factor(data$pcosdx_famhx_parent___18, levels = c("0", "1"))
data$pcosdx_famhx_parent___31.factor <- factor(data$pcosdx_famhx_parent___31, levels = c("0", "1"))
data$pcosdx_famhx_parent___32.factor <- factor(data$pcosdx_famhx_parent___32, levels = c("0", "1"))
data$pcosdx_famhx_parent___33.factor <- factor(data$pcosdx_famhx_parent___33, levels = c("0", "1"))
data$pcosdx_famhx_parent___60.factor <- factor(data$pcosdx_famhx_parent___60, levels = c("0", "1"))
data$pcosdx_famhx_parent___0.factor <- factor(data$pcosdx_famhx_parent___0, levels = c("0", "1"))
data$pcosdx_famhx_parent___unk.factor <- factor(data$pcosdx_famhx_parent___unk, levels = c("0", "1"))
data$pcosdx_famhx_parent___na.factor <- factor(data$pcosdx_famhx_parent___na, levels = c("0", "1"))
data$pcosdx_famhx_parent___oth.factor <- factor(data$pcosdx_famhx_parent___oth, levels = c("0", "1"))
data$pcosdx_famhx_parent___pm.factor <- factor(data$pcosdx_famhx_parent___pm, levels = c("0", "1"))
data$pcosdx_famhx___27.factor <- factor(data$pcosdx_famhx___27, levels = c("0", "1"))
data$pcosdx_famhx___28.factor <- factor(data$pcosdx_famhx___28, levels = c("0", "1"))
data$pcosdx_famhx___30.factor <- factor(data$pcosdx_famhx___30, levels = c("0", "1"))
data$pcosdx_famhx___5.factor <- factor(data$pcosdx_famhx___5, levels = c("0", "1"))
data$pcosdx_famhx___29.factor <- factor(data$pcosdx_famhx___29, levels = c("0", "1"))
data$pcosdx_famhx___7.factor <- factor(data$pcosdx_famhx___7, levels = c("0", "1"))
data$pcosdx_famhx___24.factor <- factor(data$pcosdx_famhx___24, levels = c("0", "1"))
data$pcosdx_famhx___8.factor <- factor(data$pcosdx_famhx___8, levels = c("0", "1"))
data$pcosdx_famhx___25.factor <- factor(data$pcosdx_famhx___25, levels = c("0", "1"))
data$pcosdx_famhx___12.factor <- factor(data$pcosdx_famhx___12, levels = c("0", "1"))
data$pcosdx_famhx___22.factor <- factor(data$pcosdx_famhx___22, levels = c("0", "1"))
data$pcosdx_famhx___13.factor <- factor(data$pcosdx_famhx___13, levels = c("0", "1"))
data$pcosdx_famhx___14.factor <- factor(data$pcosdx_famhx___14, levels = c("0", "1"))
data$pcosdx_famhx___26.factor <- factor(data$pcosdx_famhx___26, levels = c("0", "1"))
data$pcosdx_famhx___15.factor <- factor(data$pcosdx_famhx___15, levels = c("0", "1"))
data$pcosdx_famhx___16.factor <- factor(data$pcosdx_famhx___16, levels = c("0", "1"))
data$pcosdx_famhx___17.factor <- factor(data$pcosdx_famhx___17, levels = c("0", "1"))
data$pcosdx_famhx___18.factor <- factor(data$pcosdx_famhx___18, levels = c("0", "1"))
data$pcosdx_famhx___31.factor <- factor(data$pcosdx_famhx___31, levels = c("0", "1"))
data$pcosdx_famhx___32.factor <- factor(data$pcosdx_famhx___32, levels = c("0", "1"))
data$pcosdx_famhx___33.factor <- factor(data$pcosdx_famhx___33, levels = c("0", "1"))
data$pcosdx_famhx___60.factor <- factor(data$pcosdx_famhx___60, levels = c("0", "1"))
data$pcosdx_famhx___0.factor <- factor(data$pcosdx_famhx___0, levels = c("0", "1"))
data$pcosdx_famhx___unk.factor <- factor(data$pcosdx_famhx___unk, levels = c("0", "1"))
data$pcosdx_famhx___na.factor <- factor(data$pcosdx_famhx___na, levels = c("0", "1"))
data$pcosdx_famhx___oth.factor <- factor(data$pcosdx_famhx___oth, levels = c("0", "1"))
data$pcosdx_famhx___pm.factor <- factor(data$pcosdx_famhx___pm, levels = c("0", "1"))
data$pcosdx_famhx_pcos_specify___1.factor <- factor(data$pcosdx_famhx_pcos_specify___1, levels = c("0", "1"))
data$pcosdx_famhx_pcos_specify___2.factor <- factor(data$pcosdx_famhx_pcos_specify___2, levels = c("0", "1"))
data$pcosdx_famhx_pcos_specify___3.factor <- factor(data$pcosdx_famhx_pcos_specify___3, levels = c("0", "1"))
data$pcosdx_famhx_pcos_specify___unk.factor <- factor(data$pcosdx_famhx_pcos_specify___unk, levels = c("0", "1"))
data$pcosdx_famhx_pcos_specify___na.factor <- factor(data$pcosdx_famhx_pcos_specify___na, levels = c("0", "1"))
data$pcosdx_famhx_pcos_specify___oth.factor <- factor(data$pcosdx_famhx_pcos_specify___oth, levels = c("0", "1"))
data$pcosdx_famhx_pcos_specify___pm.factor <- factor(data$pcosdx_famhx_pcos_specify___pm, levels = c("0", "1"))
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
data$pcosdx_padx_body_odor.factor <- factor(data$pcosdx_padx_body_odor, levels = c("1", "0", "UNK"))
data$pcosdx_padx_axillary_hair.factor <- factor(data$pcosdx_padx_axillary_hair, levels = c("1", "0", "UNK"))
data$pcosdx_padx_tanner_pubic.factor <- factor(data$pcosdx_padx_tanner_pubic, levels = c("1", "2", "3", "4", "5", "UNK"))
data$pcosdx_padx_tanner_breast.factor <- factor(data$pcosdx_padx_tanner_breast, levels = c("1", "2", "3", "4", "5", "UNK"))
data$pcosdx_padx_acne.factor <- factor(data$pcosdx_padx_acne, levels = c("1", "0", "UNK"))
data$pcosdx_padx_acne_severity.factor <- factor(data$pcosdx_padx_acne_severity, levels = c("1", "2", "3", "0", "UNK"))
data$pcosdx_surghx___1.factor <- factor(data$pcosdx_surghx___1, levels = c("0", "1"))
data$pcosdx_surghx___2.factor <- factor(data$pcosdx_surghx___2, levels = c("0", "1"))
data$pcosdx_surghx___3.factor <- factor(data$pcosdx_surghx___3, levels = c("0", "1"))
data$pcosdx_surghx___4.factor <- factor(data$pcosdx_surghx___4, levels = c("0", "1"))
data$pcosdx_surghx___5.factor <- factor(data$pcosdx_surghx___5, levels = c("0", "1"))
data$pcosdx_surghx___6.factor <- factor(data$pcosdx_surghx___6, levels = c("0", "1"))
data$pcosdx_surghx___7.factor <- factor(data$pcosdx_surghx___7, levels = c("0", "1"))
data$pcosdx_surghx___8.factor <- factor(data$pcosdx_surghx___8, levels = c("0", "1"))
data$pcosdx_surghx___9.factor <- factor(data$pcosdx_surghx___9, levels = c("0", "1"))
data$pcosdx_surghx___10.factor <- factor(data$pcosdx_surghx___10, levels = c("0", "1"))
data$pcosdx_surghx___60.factor <- factor(data$pcosdx_surghx___60, levels = c("0", "1"))
data$pcosdx_surghx___0.factor <- factor(data$pcosdx_surghx___0, levels = c("0", "1"))
data$pcosdx_surghx___unk.factor <- factor(data$pcosdx_surghx___unk, levels = c("0", "1"))
data$pcosdx_surghx___na.factor <- factor(data$pcosdx_surghx___na, levels = c("0", "1"))
data$pcosdx_surghx___oth.factor <- factor(data$pcosdx_surghx___oth, levels = c("0", "1"))
data$pcosdx_surghx___pm.factor <- factor(data$pcosdx_surghx___pm, levels = c("0", "1"))
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
data$pcosdx_schooltype.factor <- factor(data$pcosdx_schooltype, levels = c("1", "2", "60", "UNK"))
data$pcosdx_schoolability.factor <- factor(data$pcosdx_schoolability, levels = c("1", "2", "UNK"))
data$pcosdx_learning.factor <- factor(data$pcosdx_learning, levels = c("1", "0", "UNK"))
data$pcosdx_specialty.factor <- factor(data$pcosdx_specialty, levels = c("1", "2", "3", "4", "5", "60", "UNK"))
data$pcosdx_irregular_menses.factor <- factor(data$pcosdx_irregular_menses, levels = c("1", "2", "4", "5", "6", "UNK"))
data$pcosdx_hyperandrogenism.factor <- factor(data$pcosdx_hyperandrogenism, levels = c("1", "2", "3", "UNK"))
data$pcosdx_lmp_monthsago.factor <- factor(data$pcosdx_lmp_monthsago, levels = c("1", "2", "3", "4", "5", "UNK"))
data$pcosdx_menfreq.factor <- factor(data$pcosdx_menfreq, levels = c("0", "1", "2", "3", "4", "5", "UNK"))
data$pcosdx_us.factor <- factor(data$pcosdx_us, levels = c("1", "0", "UNK"))
data$history_complete.factor <- factor(data$history_complete, levels = c("0", "1", "2"))
data$cv_visittype.factor <- factor(data$cv_visittype, levels = c("1", "2"))
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
data$cv_newdx___3.factor <- factor(data$cv_newdx___3, levels = c("0", "1"))
data$cv_newdx___4.factor <- factor(data$cv_newdx___4, levels = c("0", "1"))
data$cv_newdx___5.factor <- factor(data$cv_newdx___5, levels = c("0", "1"))
data$cv_newdx___6.factor <- factor(data$cv_newdx___6, levels = c("0", "1"))
data$cv_newdx___29.factor <- factor(data$cv_newdx___29, levels = c("0", "1"))
data$cv_newdx___7.factor <- factor(data$cv_newdx___7, levels = c("0", "1"))
data$cv_newdx___8.factor <- factor(data$cv_newdx___8, levels = c("0", "1"))
data$cv_newdx___9.factor <- factor(data$cv_newdx___9, levels = c("0", "1"))
data$cv_newdx___10.factor <- factor(data$cv_newdx___10, levels = c("0", "1"))
data$cv_newdx___11.factor <- factor(data$cv_newdx___11, levels = c("0", "1"))
data$cv_newdx___12.factor <- factor(data$cv_newdx___12, levels = c("0", "1"))
data$cv_newdx___22.factor <- factor(data$cv_newdx___22, levels = c("0", "1"))
data$cv_newdx___13.factor <- factor(data$cv_newdx___13, levels = c("0", "1"))
data$cv_newdx___14.factor <- factor(data$cv_newdx___14, levels = c("0", "1"))
data$cv_newdx___15.factor <- factor(data$cv_newdx___15, levels = c("0", "1"))
data$cv_newdx___16.factor <- factor(data$cv_newdx___16, levels = c("0", "1"))
data$cv_newdx___27.factor <- factor(data$cv_newdx___27, levels = c("0", "1"))
data$cv_newdx___28.factor <- factor(data$cv_newdx___28, levels = c("0", "1"))
data$cv_newdx___17.factor <- factor(data$cv_newdx___17, levels = c("0", "1"))
data$cv_newdx___18.factor <- factor(data$cv_newdx___18, levels = c("0", "1"))
data$cv_newdx___19.factor <- factor(data$cv_newdx___19, levels = c("0", "1"))
data$cv_newdx___20.factor <- factor(data$cv_newdx___20, levels = c("0", "1"))
data$cv_newdx___21.factor <- factor(data$cv_newdx___21, levels = c("0", "1"))
data$cv_newdx___23.factor <- factor(data$cv_newdx___23, levels = c("0", "1"))
data$cv_newdx___24.factor <- factor(data$cv_newdx___24, levels = c("0", "1"))
data$cv_newdx___25.factor <- factor(data$cv_newdx___25, levels = c("0", "1"))
data$cv_newdx___26.factor <- factor(data$cv_newdx___26, levels = c("0", "1"))
data$cv_newdx___60.factor <- factor(data$cv_newdx___60, levels = c("0", "1"))
data$cv_newdx___0.factor <- factor(data$cv_newdx___0, levels = c("0", "1"))
data$cv_newdx___unk.factor <- factor(data$cv_newdx___unk, levels = c("0", "1"))
data$cv_newdx___na.factor <- factor(data$cv_newdx___na, levels = c("0", "1"))
data$cv_newdx___oth.factor <- factor(data$cv_newdx___oth, levels = c("0", "1"))
data$cv_newdx___pm.factor <- factor(data$cv_newdx___pm, levels = c("0", "1"))
data$cv_surghx___1.factor <- factor(data$cv_surghx___1, levels = c("0", "1"))
data$cv_surghx___2.factor <- factor(data$cv_surghx___2, levels = c("0", "1"))
data$cv_surghx___3.factor <- factor(data$cv_surghx___3, levels = c("0", "1"))
data$cv_surghx___4.factor <- factor(data$cv_surghx___4, levels = c("0", "1"))
data$cv_surghx___5.factor <- factor(data$cv_surghx___5, levels = c("0", "1"))
data$cv_surghx___6.factor <- factor(data$cv_surghx___6, levels = c("0", "1"))
data$cv_surghx___7.factor <- factor(data$cv_surghx___7, levels = c("0", "1"))
data$cv_surghx___8.factor <- factor(data$cv_surghx___8, levels = c("0", "1"))
data$cv_surghx___9.factor <- factor(data$cv_surghx___9, levels = c("0", "1"))
data$cv_surghx___10.factor <- factor(data$cv_surghx___10, levels = c("0", "1"))
data$cv_surghx___60.factor <- factor(data$cv_surghx___60, levels = c("0", "1"))
data$cv_surghx___0.factor <- factor(data$cv_surghx___0, levels = c("0", "1"))
data$cv_surghx___unk.factor <- factor(data$cv_surghx___unk, levels = c("0", "1"))
data$cv_surghx___na.factor <- factor(data$cv_surghx___na, levels = c("0", "1"))
data$cv_surghx___oth.factor <- factor(data$cv_surghx___oth, levels = c("0", "1"))
data$cv_surghx___pm.factor <- factor(data$cv_surghx___pm, levels = c("0", "1"))
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
data$cv_reflux_meds.factor <- factor(data$cv_reflux_meds, levels = c("1", "2", "3", "4"))
data$cv_sleep_meds.factor <- factor(data$cv_sleep_meds, levels = c("1", "2", "3", "4", "5"))
data$cv_metforminformulation.factor <- factor(data$cv_metforminformulation, levels = c("1", "2", "3", "60"))
data$cv_metformintdd.factor <- factor(data$cv_metformintdd, levels = c("500", "750", "850", "1000", "1500", "1700", "2000", "60"))
data$cv_metforminfreq.factor <- factor(data$cv_metforminfreq, levels = c("1", "2", "3"))
data$cv_misseddoses.factor <- factor(data$cv_misseddoses, levels = c("1", "2", "3"))
data$cv_gisymptoms.factor <- factor(data$cv_gisymptoms, levels = c("1", "0", "UNK"))
data$cv_sideeffectsreported___1.factor <- factor(data$cv_sideeffectsreported___1, levels = c("0", "1"))
data$cv_sideeffectsreported___2.factor <- factor(data$cv_sideeffectsreported___2, levels = c("0", "1"))
data$cv_sideeffectsreported___3.factor <- factor(data$cv_sideeffectsreported___3, levels = c("0", "1"))
data$cv_sideeffectsreported___4.factor <- factor(data$cv_sideeffectsreported___4, levels = c("0", "1"))
data$cv_sideeffectsreported___5.factor <- factor(data$cv_sideeffectsreported___5, levels = c("0", "1"))
data$cv_sideeffectsreported___6.factor <- factor(data$cv_sideeffectsreported___6, levels = c("0", "1"))
data$cv_sideeffectsreported___7.factor <- factor(data$cv_sideeffectsreported___7, levels = c("0", "1"))
data$cv_sideeffectsreported___60.factor <- factor(data$cv_sideeffectsreported___60, levels = c("0", "1"))
data$cv_sideeffectsreported___unk.factor <- factor(data$cv_sideeffectsreported___unk, levels = c("0", "1"))
data$cv_sideeffectsreported___na.factor <- factor(data$cv_sideeffectsreported___na, levels = c("0", "1"))
data$cv_sideeffectsreported___oth.factor <- factor(data$cv_sideeffectsreported___oth, levels = c("0", "1"))
data$cv_sideeffectsreported___pm.factor <- factor(data$cv_sideeffectsreported___pm, levels = c("0", "1"))
data$cv_menfreq.factor <- factor(data$cv_menfreq, levels = c("0", "1", "2", "3", "4"))
data$cv_sexuallyactive.factor <- factor(data$cv_sexuallyactive, levels = c("1", "0", "UNK"))
data$cv_genderid.factor <- factor(data$cv_genderid, levels = c("1", "2", "3", "60", "UNK"))
data$cv_sexualpref.factor <- factor(data$cv_sexualpref, levels = c("1", "2", "3", "4", "5", "60", "UNK"))
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
data$labs_this_visit___1.factor <- factor(data$labs_this_visit___1, levels = c("0", "1"))
data$labs_this_visit___2.factor <- factor(data$labs_this_visit___2, levels = c("0", "1"))
data$labs_this_visit___3.factor <- factor(data$labs_this_visit___3, levels = c("0", "1"))
data$labs_this_visit___4.factor <- factor(data$labs_this_visit___4, levels = c("0", "1"))
data$labs_this_visit___5.factor <- factor(data$labs_this_visit___5, levels = c("0", "1"))
data$labs_this_visit___6.factor <- factor(data$labs_this_visit___6, levels = c("0", "1"))
data$labs_this_visit___7.factor <- factor(data$labs_this_visit___7, levels = c("0", "1"))
data$labs_this_visit___8.factor <- factor(data$labs_this_visit___8, levels = c("0", "1"))
data$labs_this_visit___9.factor <- factor(data$labs_this_visit___9, levels = c("0", "1"))
data$labs_this_visit___10.factor <- factor(data$labs_this_visit___10, levels = c("0", "1"))
data$labs_this_visit___11.factor <- factor(data$labs_this_visit___11, levels = c("0", "1"))
data$labs_this_visit___12.factor <- factor(data$labs_this_visit___12, levels = c("0", "1"))
data$labs_this_visit___13.factor <- factor(data$labs_this_visit___13, levels = c("0", "1"))
data$labs_this_visit___14.factor <- factor(data$labs_this_visit___14, levels = c("0", "1"))
data$labs_this_visit___15.factor <- factor(data$labs_this_visit___15, levels = c("0", "1"))
data$labs_this_visit___16.factor <- factor(data$labs_this_visit___16, levels = c("0", "1"))
data$labs_this_visit___17.factor <- factor(data$labs_this_visit___17, levels = c("0", "1"))
data$labs_this_visit___18.factor <- factor(data$labs_this_visit___18, levels = c("0", "1"))
data$labs_this_visit___19.factor <- factor(data$labs_this_visit___19, levels = c("0", "1"))
data$labs_this_visit___20.factor <- factor(data$labs_this_visit___20, levels = c("0", "1"))
data$labs_this_visit___21.factor <- factor(data$labs_this_visit___21, levels = c("0", "1"))
data$labs_this_visit___22.factor <- factor(data$labs_this_visit___22, levels = c("0", "1"))
data$labs_this_visit___23.factor <- factor(data$labs_this_visit___23, levels = c("0", "1"))
data$labs_this_visit___24.factor <- factor(data$labs_this_visit___24, levels = c("0", "1"))
data$labs_this_visit___25.factor <- factor(data$labs_this_visit___25, levels = c("0", "1"))
data$labs_this_visit___26.factor <- factor(data$labs_this_visit___26, levels = c("0", "1"))
data$labs_this_visit___27.factor <- factor(data$labs_this_visit___27, levels = c("0", "1"))
data$labs_this_visit___28.factor <- factor(data$labs_this_visit___28, levels = c("0", "1"))
data$labs_this_visit___29.factor <- factor(data$labs_this_visit___29, levels = c("0", "1"))
data$labs_this_visit___31.factor <- factor(data$labs_this_visit___31, levels = c("0", "1"))
data$labs_this_visit___32.factor <- factor(data$labs_this_visit___32, levels = c("0", "1"))
data$labs_this_visit___30.factor <- factor(data$labs_this_visit___30, levels = c("0", "1"))
data$labs_this_visit___33.factor <- factor(data$labs_this_visit___33, levels = c("0", "1"))
data$labs_this_visit___na.factor <- factor(data$labs_this_visit___na, levels = c("0", "1"))
data$labs_this_visit___unk.factor <- factor(data$labs_this_visit___unk, levels = c("0", "1"))
data$labs_this_visit___oth.factor <- factor(data$labs_this_visit___oth, levels = c("0", "1"))
data$labs_this_visit___pm.factor <- factor(data$labs_this_visit___pm, levels = c("0", "1"))
data$cv_upt_result.factor <- factor(data$cv_upt_result, levels = c("1", "0"))
data$cv_tt_assay.factor <- factor(data$cv_tt_assay, levels = c("1", "OTH", "UNK"))
data$cv_tg_fasting.factor <- factor(data$cv_tg_fasting, levels = c("1", "2"))
data$cv_liverimaging_results.factor <- factor(data$cv_liverimaging_results, levels = c("1", "0"))
data$cv_osa_sx.factor <- factor(data$cv_osa_sx, levels = c("1", "2", "4", "3", "0", "UNK"))
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
data$cv_reflux_newmeds.factor <- factor(data$cv_reflux_newmeds, levels = c("1", "2", "3", "4"))
data$cv_sleep_newmeds.factor <- factor(data$cv_sleep_newmeds, levels = c("1", "2", "3", "4", "5"))
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
data$cv_fitnesstesting_yn.factor <- factor(data$cv_fitnesstesting_yn, levels = c("1", "0"))
data$cv_exerciseplan.factor <- factor(data$cv_exerciseplan, levels = c("1", "0", "2", "UNK"))
data$clinical_visit_complete.factor <- factor(data$clinical_visit_complete, levels = c("0", "1", "2"))
levels(data$redcap_repeat_instrument.factor) <- c("Clinical Visit")
levels(data$redcap_data_access_group.factor) <- c("BCH", "CHI", "CHLA", "CHOP", "CSMC", "DC", "DCH", "DEN", "KC", "LC", "NCH", "NY/COL", "PENN", "PITT", "TAMU", "UAB", "UCM", "UF", "UVA", "VAN", "VCU")
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
levels(data$pcosdx_famhx_parent___27.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___28.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___30.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___5.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___29.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___7.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___24.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___8.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___25.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___12.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___22.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___13.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___14.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___26.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___15.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___16.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___17.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___18.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___31.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___32.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___33.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___60.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___0.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___unk.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___na.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___oth.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_parent___pm.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___27.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___28.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___30.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___5.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___29.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___7.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___24.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___8.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___25.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___12.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___22.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___13.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___14.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___26.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___15.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___16.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___17.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___18.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___31.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___32.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___33.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___60.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___0.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___unk.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___na.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___oth.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx___pm.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_pcos_specify___1.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_pcos_specify___2.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_pcos_specify___3.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_pcos_specify___unk.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_pcos_specify___na.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_pcos_specify___oth.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_famhx_pcos_specify___pm.factor) <- c("Unchecked", "Checked")
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
levels(data$pcosdx_padx_body_odor.factor) <- c("Yes", "No", "Unknown/Not recorded")
levels(data$pcosdx_padx_axillary_hair.factor) <- c("Yes", "No", "Unknown/Not recorded")
levels(data$pcosdx_padx_tanner_pubic.factor) <- c("Tanner 1", "Tanner 2", "Tanner 3", "Tanner 4", "Tanner 5", "Unknown/Not recorded")
levels(data$pcosdx_padx_tanner_breast.factor) <- c("Tanner 1", "Tanner 2", "Tanner 3", "Tanner 4", "Tanner 5", "Unknown/Not recorded")
levels(data$pcosdx_padx_acne.factor) <- c("Yes", "No", "Unknown/Not recorded")
levels(data$pcosdx_padx_acne_severity.factor) <- c("Mild", "Moderate (pustular)", "Severe (nodular, cystic, scarring)", "None", "Unknown/Not recorded")
levels(data$pcosdx_surghx___1.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___2.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___3.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___4.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___5.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___6.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___7.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___8.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___9.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___10.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___60.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___0.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___unk.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___na.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___oth.factor) <- c("Unchecked", "Checked")
levels(data$pcosdx_surghx___pm.factor) <- c("Unchecked", "Checked")
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
levels(data$pcosdx_schooltype.factor) <- c("Traditional school", "Home/online school", "Other {schooltype_other}", "Unknown/Not recorded")
levels(data$pcosdx_schoolability.factor) <- c("Acceptable", "Poor or family concerned", "Unknown/Not recorded")
levels(data$pcosdx_learning.factor) <- c("Yes", "No", "Unknown/Not recorded")
levels(data$pcosdx_specialty.factor) <- c("Endocrinology", "Primary Care (Pediatrics or Family Medicine)", "Gynecology", "Dermatology", "Adolescent Medicine", "Other {pcosdx_specialty_other}", "Unknown/Not recorded")
levels(data$pcosdx_irregular_menses.factor) <- c(">1 year post-menarche: > 90 days for any one cycle", ">1 year post-menarche: > or = 6 week cycles (< or = 8 cycles a year)", "Primary amenorrhea by age 15", "Primary amenorrhea > 3 years post-thelarche", "Excess non-ovulatory menses", "Unknown/Not recorded")
levels(data$pcosdx_hyperandrogenism.factor) <- c("Only clinical", "Only biochemical", "Both clinical and biochemical", "Unknown/Not recorded")
levels(data$pcosdx_lmp_monthsago.factor) <- c("1 month ago", "1-3 months ago", "3-6 months ago", "6-9 months ago", ">9 months ago", "Unknown/Not recorded")
levels(data$pcosdx_menfreq.factor) <- c("0 menses", "1 or 2 menses", "3-4 menses", "5-8 menses", "8-12 menses", ">12 menses", "Unknown/Not recorded")
levels(data$pcosdx_us.factor) <- c("Yes", "No", "Unknown/Not recorded")
levels(data$history_complete.factor) <- c("Incomplete", "Unverified", "Complete")
levels(data$cv_visittype.factor) <- c("Diagnostic Visit", "Follow-up")
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
levels(data$cv_newdx___3.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___4.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___5.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___6.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___29.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___7.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___8.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___9.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___10.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___11.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___12.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___22.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___13.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___14.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___15.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___16.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___27.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___28.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___17.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___18.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___19.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___20.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___21.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___23.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___24.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___25.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___26.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___60.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___0.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___unk.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___na.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___oth.factor) <- c("Unchecked", "Checked")
levels(data$cv_newdx___pm.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___1.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___2.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___3.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___4.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___5.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___6.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___7.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___8.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___9.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___10.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___60.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___0.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___unk.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___na.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___oth.factor) <- c("Unchecked", "Checked")
levels(data$cv_surghx___pm.factor) <- c("Unchecked", "Checked")
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
levels(data$cv_reflux_meds.factor) <- c("H2 Blocker", "PPI", "Carafate", "Tums")
levels(data$cv_sleep_meds.factor) <- c("Ambien", "Melatonin", "Ciproheptadine", "Clonidine", "Trazadone")
levels(data$cv_metforminformulation.factor) <- c("Standard release", "Extended release", "Riomet (liquid)", "Other {metformin_form_other}")
levels(data$cv_metformintdd.factor) <- c("500", "750", "850", "1000", "1500", "1700", "2000", "Other dose {cv_meds_metformin_otherdose}")
levels(data$cv_metforminfreq.factor) <- c("1", "2", "3")
levels(data$cv_misseddoses.factor) <- c("0-1 doses", "2-4 doses", "5+ doses")
levels(data$cv_gisymptoms.factor) <- c("Yes", "No", "Unknown/not recorded")
levels(data$cv_sideeffectsreported___1.factor) <- c("Unchecked", "Checked")
levels(data$cv_sideeffectsreported___2.factor) <- c("Unchecked", "Checked")
levels(data$cv_sideeffectsreported___3.factor) <- c("Unchecked", "Checked")
levels(data$cv_sideeffectsreported___4.factor) <- c("Unchecked", "Checked")
levels(data$cv_sideeffectsreported___5.factor) <- c("Unchecked", "Checked")
levels(data$cv_sideeffectsreported___6.factor) <- c("Unchecked", "Checked")
levels(data$cv_sideeffectsreported___7.factor) <- c("Unchecked", "Checked")
levels(data$cv_sideeffectsreported___60.factor) <- c("Unchecked", "Checked")
levels(data$cv_sideeffectsreported___unk.factor) <- c("Unchecked", "Checked")
levels(data$cv_sideeffectsreported___na.factor) <- c("Unchecked", "Checked")
levels(data$cv_sideeffectsreported___oth.factor) <- c("Unchecked", "Checked")
levels(data$cv_sideeffectsreported___pm.factor) <- c("Unchecked", "Checked")
levels(data$cv_menfreq.factor) <- c("0 menses", "1-2 menses", "3-4 menses", "5-6 menses", ">6 menses")
levels(data$cv_sexuallyactive.factor) <- c("Yes", "No", "Unknown/Not recorded")
levels(data$cv_genderid.factor) <- c("Female (cisgender)", "Male (transgender)", "Nonbinary", "Other {cv_genderid_other}", "Unknown/Not recorded")
levels(data$cv_sexualpref.factor) <- c("Interested in male partners (heterosexual)", "Interested in female partners (gay or lesbian)", "Bisexual", "Pansexual", "Asexual", "Other {cv_sexualpref_other}", "Unknown/Not recorded")
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
levels(data$labs_this_visit___1.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___2.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___3.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___4.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___5.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___6.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___7.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___8.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___9.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___10.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___11.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___12.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___13.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___14.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___15.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___16.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___17.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___18.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___19.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___20.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___21.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___22.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___23.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___24.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___25.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___26.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___27.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___28.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___29.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___31.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___32.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___30.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___33.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___na.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___unk.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___oth.factor) <- c("Unchecked", "Checked")
levels(data$labs_this_visit___pm.factor) <- c("Unchecked", "Checked")
levels(data$cv_upt_result.factor) <- c("Positive", "Negative")
levels(data$cv_tt_assay.factor) <- c("LCMSMS", "Other {cv_tt_assay_oth}", "Not recorded/unknown")
levels(data$cv_tg_fasting.factor) <- c("Fasting", "Non-fasting")
levels(data$cv_liverimaging_results.factor) <- c("Nonalcoholic fatty liver disease", "No nonalcoholic fatty liver disease")
levels(data$cv_osa_sx.factor) <- c("Yes, not addressed", "Yes, prescribed CPAP but not using it", "Yes, referral placed for sleep study", "No, treated with CPAP", "No, not currently treated", "Unknown/not recorded")
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
levels(data$cv_reflux_newmeds.factor) <- c("H2 Blocker", "PPI", "Carafate", "Tums")
levels(data$cv_sleep_newmeds.factor) <- c("Ambien", "Melatonin", "Ciproheptadine", "Clonidine", "Trazadone")
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
levels(data$cv_fitnesstesting_yn.factor) <- c("Yes", "No")
levels(data$cv_exerciseplan.factor) <- c("Yes", "No", "No, but referred to a exercise specialist", "Unknown/not recorded")
levels(data$clinical_visit_complete.factor) <- c("Incomplete", "Unverified", "Complete")
# Label factor variables as well
factors <- which(duplicated(sub("\\.factor", "", colnames(data))))
# Create lists of variables for analysis aims
aim1_vars <- c(
  "record_number", "redcap_repeat_instrument", "redcap_repeat_instance",
  "cv_bmi", "cv_bmi_percentile", "cv_bmi_z", "cv_hirsutism_num",
  "cv_hirsutism_cat", "cv_acneface", "cv_acneother___1", "cv_acneother___2",
  "cv_acneother___0", "cv_acneother___unk", "cv_acneother___na",
  "cv_acneother___oth", "cv_acneother___pm", "cv_ft", "cv_ft_perc", "cv_tt",
  "cv_tt_perc", "cv_dheas", "cv_dheas_perc", "cv_androstendione",
  "cv_androstendione_perc", "cv_lh", "cv_fsh", "cv_amh", "cv_a1c", "cv_tg",
  "cv_hdl", "cv_shbg", "cv_alt", "cv_fbg", "cv_fastinsulin",
  "cv_2hrglucoseogtt", "cv_acanthosisneck", "cv_waist", "cv_osa_sx",
  "cv_medications___1", "cv_medications___2", "cv_medications___3",
  "cv_medications___4", "cv_medications___5", "cv_medications___6",
  "cv_medications___7", "cv_medications___8", "cv_medications___9",
  "cv_medications___10", "cv_medications___11", "cv_medications___12",
  "cv_medications___13", "cv_medications___14", "cv_medications___15",
  "cv_medications___16", "cv_medications___17", "cv_medications___18",
  "cv_medications___19", "cv_medications___20", "cv_medications___21",
  "cv_medications___22", "cv_medications___23", "cv_medications___32",
  "cv_medications___24", "cv_medications___25", "cv_medications___26",
  "cv_medications___27", "cv_medications___28", "cv_medications___29",
  "cv_medications___30", "cv_medications___31", "cv_medications___60",
  "cv_medications___0", "cv_medications___unk", "cv_medications___na",
  "cv_medications___oth", "cv_medications___pm",
  "redcap_repeat_instrument.factor", "cv_hirsutism_cat.factor",
  "cv_acneface.factor", "cv_acneother___1.factor", "cv_acneother___2.factor",
  "cv_acneother___0.factor", "cv_acneother___unk.factor",
  "cv_acneother___na.factor", "cv_acneother___oth.factor",
  "cv_acneother___pm.factor", "cv_acanthosisneck.factor", "cv_osa_sx.factor",
  "cv_medications___1.factor", "cv_medications___2.factor",
  "cv_medications___3.factor", "cv_medications___4.factor",
  "cv_medications___5.factor", "cv_medications___6.factor",
  "cv_medications___7.factor", "cv_medications___8.factor",
  "cv_medications___9.factor", "cv_medications___10.factor",
  "cv_medications___11.factor", "cv_medications___12.factor",
  "cv_medications___13.factor", "cv_medications___14.factor",
  "cv_medications___15.factor", "cv_medications___16.factor",
  "cv_medications___17.factor", "cv_medications___18.factor",
  "cv_medications___19.factor", "cv_medications___20.factor",
  "cv_medications___21.factor", "cv_medications___22.factor",
  "cv_medications___23.factor", "cv_medications___32.factor",
  "cv_medications___24.factor", "cv_medications___25.factor",
  "cv_medications___26.factor", "cv_medications___27.factor",
  "cv_medications___28.factor", "cv_medications___29.factor",
  "cv_medications___30.factor", "cv_medications___31.factor",
  "cv_medications___60.factor", "cv_medications___0.factor",
  "cv_medications___unk.factor", "cv_medications___na.factor",
  "cv_medications___oth.factor", "cv_medications___pm.factor"
)
