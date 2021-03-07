# Read in
data=read.csv('./NelsonVitalStatusDEC_DATA_2021-02-04_1616.csv',na.strings = "")
#Setting Labels
# label(data$record_id)="Patient ID"
# label(data$redcap_event_name)="Event Name"
# label(data$redcap_repeat_instrument)="Repeat Instrument"
# label(data$redcap_repeat_instance)="Repeat Instance"
# label(data$esrd_start_date)="Start date of chronic renal replacement therapy:"
# label(data$esrd_start_date_est)="Date of chronic renal replacement therapy is estimated?"
# label(data$esrd_type)="Type of chronic renal replacement therapy at initiation of therapy:"
# label(data$esrd_is_transplanted)="Ever transplanted:"
# label(data$esrd_diagnosis)="Primary diagnosis of kidney failure:"
# label(data$esrd_complete_date)="Date form completed:"
# label(data$esrd_initials)="Initials of person completing form:"
# label(data$end_stage_renal_disease_form_complete)="Complete?"
# label(data$transplant_date)="Date of transplant:"
# label(data$transplant_complete_date)="Date form completed:"
# label(data$transplant_initials)="Initials of person completing form:"
# label(data$transplant_form_complete)="Complete?"
# label(data$dod)="Date of death"
# label(data$codunerlying)="Cause of death from certificate (underlying) - ICD code"
# label(data$coddrunerlying)="Cause of death from medical review (DR review underlying) - ICD code"
# label(data$physiciancr)="Physician performing chart review"
# label(data$initials)="Initials of person who completed form"
# label(data$death_comments)="Comments"
# label(data$death_notice_form_complete)="Complete?"
# label(data$comment)=""
# label(data$comment_form_complete)="Complete?"
#Setting Factors(will create new variable for factors)
data$redcap_event_name.factor = factor(data$redcap_event_name,levels=c("esrd_arm_1","transplant_arm_1","death_arm_1"))
data$redcap_repeat_instrument.factor = factor(data$redcap_repeat_instrument,levels=c("comment_form"))
data$esrd_start_date_est.factor = factor(data$esrd_start_date_est,levels=c("1","0"))
data$esrd_type.factor = factor(data$esrd_type,levels=c("1","2","3"))
data$esrd_is_transplanted.factor = factor(data$esrd_is_transplanted,levels=c("0","1","9"))
data$esrd_diagnosis.factor = factor(data$esrd_diagnosis,levels=c("1","2","3","4","5","6","9"))
data$end_stage_renal_disease_form_complete.factor = factor(data$end_stage_renal_disease_form_complete,levels=c("0","1","2"))
data$transplant_form_complete.factor = factor(data$transplant_form_complete,levels=c("0","1","2"))
data$physiciancr.factor = factor(data$physiciancr,levels=c("1","2"))
data$death_notice_form_complete.factor = factor(data$death_notice_form_complete,levels=c("0","1","2"))
data$comment_form_complete.factor = factor(data$comment_form_complete,levels=c("0","1","2"))

levels(data$redcap_event_name.factor)=c("ESRD","Transplant","Death")
levels(data$redcap_repeat_instrument.factor)=c("Comment Form")
levels(data$esrd_start_date_est.factor)=c("Yes","No")
levels(data$esrd_type.factor)=c("Hemodialysis","Peritoneal Dialysis","Viable Transplant")
levels(data$esrd_is_transplanted.factor)=c("No","Yes","Unknown")
levels(data$esrd_diagnosis.factor)=c("Diabetes","Glomerulonephritis","Hypertension","Cystic Kidney Disease","Obstructive Uropathy","Other","Unknown")
levels(data$end_stage_renal_disease_form_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$transplant_form_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$physiciancr.factor)=c("Nelson","Sievers")
levels(data$death_notice_form_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$comment_form_complete.factor)=c("Incomplete","Unverified","Complete")

vital_status = data
rm(data)