# Read Data
data <- read.csv("/home/tim/OneDrive/Vigers/BDC/Janet Snell-Bergeon/CALICO/Data_Raw/CALICOMulticenterPCO-Demographics_DATA_2024-08-23_1532.csv",na.strings = "")
# Setting Labels
label(data$record_number) <- "Record number"
label(data$redcap_repeat_instrument) <- "Repeat Instrument"
label(data$redcap_repeat_instance) <- "Repeat Instance"
label(data$site) <- "Study Site"
label(data$pcosdx_age) <- "Age at time of PCOS diagnosis"
label(data$race) <- "Race"
label(data$ethnicity) <- "Ethnicity "
label(data$insur_type) <- "Insurance Type"
label(data$pcosdx_irregular_menses) <- "Irregular menses defined"
label(data$pcosdx_anxietydx_age) <- "Age when anxiety diagnosed (years)"
# Setting Factors(will create new variable for factors)
data$redcap_repeat_instrument.factor <- factor(data$redcap_repeat_instrument, levels = c("clinical_visit"))
data$site.factor <- factor(data$site, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"))
data$ethnicity.factor <- factor(data$ethnicity, levels = c("0", "1", "UNK"))
data$insur_type.factor <- factor(data$insur_type, levels = c("1", "2", "3", "0", "60", "UNK"))
data$pcosdx_irregular_menses.factor <- factor(data$pcosdx_irregular_menses, levels = c("1", "2", "4", "5", "6", "UNK"))
levels(data$redcap_repeat_instrument.factor) <- c("Clinical Visit")
levels(data$site.factor) <- c("Childrens Hospital Colorado, Denver (DEN)", "NIH/Childrens National (DC)", "Childrens Hospital of Philadelphia (CHOP)", "PENN", "Childrens Hospital of Pittsburgh (PIT)", "University of Florida, Gainesville (UFPCOS)", "Mercy Childrens Hospital-Kansas City (KC)", "Childrens Hospital Los Angeles (L)", "Cook County Hospital, Chicago (CHI)", "New York/Columbia (COL)", "Northwestern/Lurie Childrens (LC)", "Boston Childrens (BCH)", "University of Alabama (UAB)", "Los Angeles /Cedars Sinai (CSMC)", "Virginia Commonwealth University (VCU)", "University of Virginia (UVA)", "Nationwide Childrens Hospital, Ohio (NCH)", "Vanderbilt (VAN)", "UChicago Medicine (UCM)", "Texas A&M (TAMU)", "Driscoll Childrens Hospital (DCH)")
levels(data$ethnicity.factor) <- c("Non-Hispanic", "Hispanic", "Unknown")
levels(data$insur_type.factor) <- c("Public", "Private", "Military", "None", "Other {insur_other}", "Unknown/Not recorded")
levels(data$pcosdx_irregular_menses.factor) <- c(">1 year post-menarche: > 90 days for any one cycle", ">1 year post-menarche: > or = 6 week cycles (< or = 8 cycles a year)", "Primary amenorrhea by age 15", "Primary amenorrhea > 3 years post-thelarche", "Excess non-ovulatory menses", "Unknown/Not recorded")
