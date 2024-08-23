# Read Data
data <- read.csv("/home/tim/OneDrive/Vigers/BDC/Janet Snell-Bergeon/CALICO/Data_Raw/CALICOMulticenterPCO-Aim1_DATA_2024-08-23_1534.csv",na.strings = "")
# Setting Labels
label(data$record_number) <- "Record number"
label(data$redcap_repeat_instrument) <- "Repeat Instrument"
label(data$redcap_repeat_instance) <- "Repeat Instance"
label(data$cv_bmi) <- "BMI"
label(data$cv_hirsutism_num) <- "Numerical FGS score "
label(data$cv_hirsutism_cat) <- "FGS score category"
label(data$cv_acneface) <- "Acne severity on the face"
label(data$cv_acneother) <- "Acne elsewhere on the body"
label(data$cv_ft) <- "Free testosterone (pg/mL) "
label(data$cv_ft_perc) <- "Free testosterone (pg/mL) percentage of upper limit of normal"
label(data$cv_tt) <- "Total testosterone (ng/dL) "
label(data$cv_tt_perc) <- "Total testosterone percentage of upper limit of normal"
label(data$cv_dheas) <- "DHEAS (µg/dL)"
label(data$cv_dheas_perc) <- "DHEAS percentage of upper limit of normal"
label(data$cv_androstendione) <- "Androstendione (ng/dL)"
label(data$cv_androstendione_perc) <- "Androstendione percentage of upper limit of normal"
label(data$cv_lh) <- "LH (mIU/mL)"
label(data$cv_fsh) <- "FSH (mIU/mL)"
label(data$cv_amh) <- "AMH (ng/mL)"
label(data$cv_a1c) <- "HbA1C (%)"
label(data$cv_tg) <- "Triglycerides (mg/dL)"
label(data$cv_hdl) <- "HDL (mg/dL)"
label(data$cv_shbg) <- "SHBG (nmol/L)"
label(data$cv_alt) <- "ALT (U/L)"
label(data$cv_fbg) <- "Fasting glucose (mg/dL)"
label(data$cv_fastinsulin) <- "Fasting insulin (mIU/mL)"
label(data$cv_2hrglucoseogtt) <- "2 hour glucose from OGTT (mg/dL)"
label(data$cv_acanthosisneck) <- "Acanthosis at the neck"
label(data$cv_waist) <- "Waist circumference (cm) "
label(data$cv_osa_sx) <- "Symptoms of sleep apnea (snoring, frequent awakening, AM headaches, daytime fatigue, napping during day)"
label(data$cv_medications) <- "Current treatment (not prescribed this visit)"
# Setting Factors(will create new variable for factors)
data$redcap_repeat_instrument.factor <- factor(data$redcap_repeat_instrument, levels = c("clinical_visit"))
data$cv_hirsutism_cat.factor <- factor(data$cv_hirsutism_cat, levels = c("0", "1", "2", "3", "4", "NA", "UNK"))
data$cv_acneface.factor <- factor(data$cv_acneface, levels = c("1", "2", "3", "0", "4", "UNK"))
data$cv_acanthosisneck.factor <- factor(data$cv_acanthosisneck, levels = c("1", "2", "3", "4", "0", "5", "UNK"))
data$cv_osa_sx.factor <- factor(data$cv_osa_sx, levels = c("1", "2", "4", "3", "0", "UNK"))
levels(data$redcap_repeat_instrument.factor) <- c("Clinical Visit")
levels(data$cv_hirsutism_cat.factor) <- c("None (0-6)", "Mild (7-9)", "Moderate (10-15)", "Severe (>15)", "Noted by clinician, severity unknown.", "Not applicable, patient had permanent hair removal (i.e. laser hair removal)", "Unknown/Not recorded")
levels(data$cv_acneface.factor) <- c("Mild", "Moderate (pustular)", "Severe (nodular, cystic, scarring)", "None", "Acne noted by clinician, severity unknown.", "Unknown/Not recorded")
levels(data$cv_acanthosisneck.factor) <- c("Back of neck, barely visible (minimal)", "Back of neck, obvious (mild)", "Lateral sides of neck (moderate)", "Circumferential (severe)", "None", "Noted by clinician, severity unknown.", "Unknown/Not recorded")
levels(data$cv_osa_sx.factor) <- c("Yes, not addressed", "Yes, prescribed CPAP but not using it", "Yes, referral placed for sleep study", "No, treated with CPAP", "No, not currently treated", "Unknown/not recorded")
