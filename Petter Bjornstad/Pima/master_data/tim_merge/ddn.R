#Read in
data=read.csv('./Nelson13DKN151Determ_DATA_2021-02-04_1610.csv',na.strings = "")
#Setting Labels
# label(data$record_id)="Patient ID"
# label(data$redcap_event_name)="Event Name"
# label(data$redcap_repeat_instrument)="Repeat Instrument"
# label(data$redcap_repeat_instance)="Repeat Instance"
# label(data$lastname_mv)="Missing value verified for: Lastname"
# label(data$firstname_mv)="Missing value verified for: Firstname"
# label(data$stratum)="Patient stratum"
# label(data$stratum_mv)="Missing value verified for: Patient stratum"
# label(data$randomnum)="Patient radomization number"
# label(data$randomnum_mv)="Missing value verified for: Patient radomization number"
# label(data$phx_number)="PIMC number"
# label(data$phx_number_mv)="Missing value verified for: PIMC number"
# label(data$dob_mv)="Missing value verified for: DOB"
# label(data$gender)="Gender"
# label(data$gender_mv)="Missing value verified for: Gender"
# label(data$height)="Height (cm)"
# label(data$height_mv)="Missing value verified for: Height (cm)"
# label(data$allergies)="Allergies?"
# label(data$allergies_mv)="Missing value verified for: Allergies?"
# label(data$allerglist_mv)="Missing value verified for: Allergies (specify)"
# label(data$district_mv)="Missing value verified for: District"
# label(data$phonenum_mv)="Missing value verified for: Phone number"
# label(data$bldprsr)="Blood pressure arm"
# label(data$bldprsr_mv)="Missing value verified for: Blood pressure arm"
# label(data$dtdiabonset)="Date of Diabetes Onset?"
# label(data$dtdiabonset_mv)="Missing value verified for: Date of Diabetes Onset?"
# label(data$typerecruit)="Type of Recruitment?"
# label(data$typerecruit_mv)="Missing value verified for: Type of Recruitment?"
# label(data$typerecruit_other_mv)="Missing value verified for: Type of Recruitment (Other)"
# label(data$familyrelation)="Family Relation"
# label(data$familyrelation_mv)="Missing value verified for: Family Relation"
# label(data$familyrelation_other_mv)="Missing value verified for: Family Relation (Other)"
# label(data$nih_number_mv)="Missing value verified for: NIH Number"
# label(data$demoform_complete)="Complete?"
# label(data$sex_code)="Sex Code"
# label(data$sex_code_mv)="Missing value verified for: Sex Code"
# label(data$start_date)="Admission date"
# label(data$start_date_mv)="Missing value verified for: Admission date"
# label(data$end_date)="Discharge date"
# label(data$end_date_mv)="Missing value verified for: Discharge date"
# label(data$study_status)="Study status"
# label(data$study_status_mv)="Missing value verified for: Study status"
# label(data$clrncvisittype)="Type of Visit"
# label(data$clrncvisittype_mv)="Missing value verified for: Type of Visit"
# label(data$clrnctestintvrl)="Test Interval"
# label(data$clrnctestintvrl_mv)="Missing value verified for: Test Interval"
# label(data$comments_mv)="Missing value verified for: Comments"
# label(data$study_code)="Study code"
# label(data$study_code_mv)="Missing value verified for: Study code"
# label(data$encounter_id_mv)="Missing value verified for: Encounter ID"
# label(data$clearancevisit_complete)="Complete?"
# label(data$brainmri_done)="Brain MRI done?"
# label(data$brainmri_done_mv)="Missing value verified for: Brain MRI done?"
# label(data$brainmri_date)="Brain MRI date"
# label(data$brainmri_date_mv)="Missing value verified for: Brain MRI date"
# label(data$cognitivetest_done)="Cognitive test done?"
# label(data$cognitivetest_done_mv)="Missing value verified for: Cognitive test done?"
# label(data$cognitivetest_date)="Cognitive test date"
# label(data$cognitivetest_date_mv)="Missing value verified for: Cognitive test date"
# label(data$kidneybiopsy_done)="Kidney biopsy done?"
# label(data$kidneybiopsy_done_mv)="Missing value verified for: Kidney biopsy done?"
# label(data$kidneybiopsy_date)="Kidney biopsy date"
# label(data$kidneybiopsy_date_mv)="Missing value verified for: Kidney biopsy date"
# label(data$kidneybiopsy_site)="Kidney biopsy site"
# label(data$kidneybiopsy_site_mv)="Missing value verified for: Kidney biopsy site"
# label(data$kidneylivermre_done)="Kidney/liver MRE done?"
# label(data$kidneylivermre_done_mv)="Missing value verified for: Kidney/liver MRE done?"
# label(data$kidneylivermre_date)="Kidney/liver MRE Date"
# label(data$kidneylivermre_date_mv)="Missing value verified for: Kidney/liver MRE Date"
# label(data$kidneymri_done)="Kidney/liver MRI done?"
# label(data$kidneymri_done_mv)="Missing value verified for: Kidney/liver MRI done?"
# label(data$kidneymri_date)="Kidney/liver MRI date"
# label(data$kidneymri_date_mv)="Missing value verified for: Kidney/liver MRI date"
# label(data$skinbiopsy_done)="Skin biopsy done?"
# label(data$skinbiopsy_done_mv)="Missing value verified for: Skin biopsy done?"
# label(data$skinbiopsy_date)="Skin biopsy date"
# label(data$skinbiopsy_date_mv)="Missing value verified for: Skin biopsy date"
# label(data$skinbiopsy_site)="Skin biopsy site?"
# label(data$skinbiopsy_site_mv)="Missing value verified for: Skin biopsy site?"
# label(data$stemcells_done)="Stem cells collected?"
# label(data$stemcells_done_mv)="Missing value verified for: Stem cells collected?"
# label(data$stemcells_date)="Stem cells collection date"
# label(data$stemcells_date_mv)="Missing value verified for: Stem cells collection date"
# label(data$stool_done)="Stool sample collection done?"
# label(data$stool_done_mv)="Missing value verified for: Stool sample collection done?"
# label(data$stool_date)="Stool sample collection date"
# label(data$stool_date_mv)="Missing value verified for: Stool sample collection date"
# label(data$neuropathy_done)="Was a neuropathy test performed?"
# label(data$neuropathy_done_mv)="Missing value verified for: Was a neuropathy test performed?"
# label(data$neuropathy_date)="Neuropathy test date"
# label(data$neuropathy_date_mv)="Missing value verified for: Neuropathy test date"
# label(data$retinopathy_done)="Was a retinopathy test performed?"
# label(data$retinopathy_done_mv)="Missing value verified for: Was a retinopathy test performed?"
# label(data$retinopathy_date)="Retinopathy test date"
# label(data$retinopathy_date_mv)="Missing value verified for: Retinopathy test date"
# label(data$sp_comments_mv)="Missing value verified for: Miscellaneous comments for procedures above"
# label(data$sp_encounter_id)="Encounter ID"
# label(data$study_procedures_status_complete)="Complete?"
# label(data$f4bldprsr)="Blood pressure arm"
# label(data$f4bldprsr_mv)="Missing value verified for: Blood pressure arm"
# label(data$f4frstmsrsys)="1st measure - Systolic"
# label(data$f4frstmsrsys_mv)="Missing value verified for: 1st measure - Systolic"
# label(data$f4frstmsrdia)="1st measure - Diastolic"
# label(data$f4frstmsrdia_mv)="Missing value verified for: 1st measure - Diastolic"
# label(data$f4scndmsrsys)="2nd measure - Systolic"
# label(data$f4scndmsrsys_mv)="Missing value verified for: 2nd measure - Systolic"
# label(data$f4scndmsrdia)="2nd measure - Diastolic"
# label(data$f4scndmsrdia_mv)="Missing value verified for: 2nd measure - Diastolic"
# label(data$f4heartrate)="Heart rate (beats/min)"
# label(data$f4heartrate_mv)="Missing value verified for: Heart rate (beats/min)"
# label(data$f4pregnant)="Pregnancy Test?"
# label(data$f4pregnant_mv)="Missing value verified for: Pregnancy Test?"
# label(data$f4medsnow)="Is the subect currently taking any other medicine?"
# label(data$f4medsnow_mv)="Missing value verified for: Is the subect currently taking any other medicine?"
# label(data$f4srmpotasm)="Results of serum potassium (mmol/L)"
# label(data$f4srmpotasm_mv)="Missing value verified for: Results of serum potassium (mmol/L)"
# label(data$f4srmglu)="Results of serum glucose (mg/dL)"
# label(data$f4srmglu_mv)="Missing value verified for: Results of serum glucose (mg/dL)"
# label(data$f4forminit_mv)="Missing value verified for: Initials of person who completed form"
# label(data$f4bpam_encounter_id)="Encounter ID"
# label(data$f4bldprsrandmeds_complete)="Complete?"
# label(data$f4drugcode1)="Drug Code"
# label(data$f4drugcode1_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage1)="Dosage"
# label(data$f4dosage1_mv)="Missing value verified for: Dosage"
# label(data$f4timesday1)="Number of times taken"
# label(data$f4timesday1_mv)="Missing value verified for: Number of times taken"
# label(data$f4dwm1)="Times/day, week or month"
# label(data$f4dwm1_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode2)="Drug Code"
# label(data$f4drugcode2_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage2)="Dosage"
# label(data$f4dosage2_mv)="Missing value verified for: Dosage"
# label(data$f4timesday2)="Number of time taken"
# label(data$f4timesday2_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm2)="Times/day, week or month"
# label(data$f4dwm2_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode3)="Drug Code"
# label(data$f4drugcode3_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage3)="Dosage"
# label(data$f4dosage3_mv)="Missing value verified for: Dosage"
# label(data$f4timesday3)="Number of time taken"
# label(data$f4timesday3_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm3)="Times/day, week or month"
# label(data$f4dwm3_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode4)="Drug Code"
# label(data$f4drugcode4_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage4)="Dosage"
# label(data$f4dosage4_mv)="Missing value verified for: Dosage"
# label(data$f4timesday4)="Number of time taken"
# label(data$f4timesday4_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm4)="Times/day, week or month"
# label(data$f4dwm4_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode5)="Drug Code"
# label(data$f4drugcode5_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage5)="Dosage"
# label(data$f4dosage5_mv)="Missing value verified for: Dosage"
# label(data$f4timesday5)="Number of time taken"
# label(data$f4timesday5_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm5)="Times/day, week or month"
# label(data$f4dwm5_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode6)="Drug Code"
# label(data$f4drugcode6_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage6)="Dosage"
# label(data$f4dosage6_mv)="Missing value verified for: Dosage"
# label(data$f4timesday6)="Number of time taken"
# label(data$f4timesday6_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm6)="Times/day, week or month"
# label(data$f4dwm6_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode7)="Drug Code"
# label(data$f4drugcode7_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage7)="Dosage"
# label(data$f4dosage7_mv)="Missing value verified for: Dosage"
# label(data$f4timesday7)="Number of time taken"
# label(data$f4timesday7_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm7)="Times/day, week or month"
# label(data$f4dwm7_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode8)="Drug Code"
# label(data$f4drugcode8_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage8)="Dosage"
# label(data$f4dosage8_mv)="Missing value verified for: Dosage"
# label(data$f4timesday8)="Number of time taken"
# label(data$f4timesday8_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm8)="Times/day, week or month"
# label(data$f4dwm8_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode9)="Drug Code"
# label(data$f4drugcode9_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage9)="Dosage"
# label(data$f4dosage9_mv)="Missing value verified for: Dosage"
# label(data$f4timesday9)="Number of time taken"
# label(data$f4timesday9_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm9)="Times/day, week or month"
# label(data$f4dwm9_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode10)="Drug Code"
# label(data$f4drugcode10_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage10)="Dosage"
# label(data$f4dosage10_mv)="Missing value verified for: Dosage"
# label(data$f4timesday10)="Number of time taken"
# label(data$f4timesday10_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm10)="Times/day, week or month"
# label(data$f4dwm10_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode11)="Drug Code"
# label(data$f4drugcode11_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage11)="Dosage"
# label(data$f4dosage11_mv)="Missing value verified for: Dosage"
# label(data$f4timesday11)="Number of time taken"
# label(data$f4timesday11_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm11)="Times/day, week or month"
# label(data$f4dwm11_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode12)="Drug Code"
# label(data$f4drugcode12_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage12)="Dosage"
# label(data$f4dosage12_mv)="Missing value verified for: Dosage"
# label(data$f4timesday12)="Number of time taken"
# label(data$f4timesday12_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm12)="Times/day, week or month"
# label(data$f4dwm12_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode13)="Drug Code"
# label(data$f4drugcode13_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage13)="Dosage"
# label(data$f4dosage13_mv)="Missing value verified for: Dosage"
# label(data$f4timesday13)="Number of time taken"
# label(data$f4timesday13_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm13)="Times/day, week or month"
# label(data$f4dwm13_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode14)="Drug Code"
# label(data$f4drugcode14_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage14)="Dosage"
# label(data$f4dosage14_mv)="Missing value verified for: Dosage"
# label(data$f4timesday14)="Number of time taken"
# label(data$f4timesday14_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm14)="Times/day, week or month"
# label(data$f4dwm14_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode15)="Drug Code"
# label(data$f4drugcode15_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage15)="Dosage"
# label(data$f4dosage15_mv)="Missing value verified for: Dosage"
# label(data$f4timesday15)="Number of time taken"
# label(data$f4timesday15_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm15)="Times/day, week or month"
# label(data$f4dwm15_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode16)="Drug Code"
# label(data$f4drugcode16_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage16)="Dosage"
# label(data$f4dosage16_mv)="Missing value verified for: Dosage"
# label(data$f4timesday16)="Number of time taken"
# label(data$f4timesday16_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm16)="Times/day, week or month"
# label(data$f4dwm16_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode17)="Drug Code"
# label(data$f4drugcode17_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage17)="Dosage"
# label(data$f4dosage17_mv)="Missing value verified for: Dosage"
# label(data$f4timesday17)="Number of time taken"
# label(data$f4timesday17_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm17)="Times/day, week or month"
# label(data$f4dwm17_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode18)="Drug Code"
# label(data$f4drugcode18_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage18)="Dosage"
# label(data$f4dosage18_mv)="Missing value verified for: Dosage"
# label(data$f4timesday18)="Number of time taken"
# label(data$f4timesday18_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm18)="Times/day, week or month"
# label(data$f4dwm18_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode19)="Drug Code"
# label(data$f4drugcode19_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage19)="Dosage"
# label(data$f4dosage19_mv)="Missing value verified for: Dosage"
# label(data$f4timesday19)="Number of time taken"
# label(data$f4timesday19_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm19)="Times/day, week or month"
# label(data$f4dwm19_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode20)="Drug Code"
# label(data$f4drugcode20_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage20)="Dosage"
# label(data$f4dosage20_mv)="Missing value verified for: Dosage"
# label(data$f4timesday20)="Number of time taken"
# label(data$f4timesday20_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm20)="Times/day, week or month"
# label(data$f4dwm20_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode21)="Drug Code"
# label(data$f4drugcode21_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage21)="Dosage"
# label(data$f4dosage21_mv)="Missing value verified for: Dosage"
# label(data$f4timesday21)="Number of time taken"
# label(data$f4timesday21_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm21)="Times/day, week or month"
# label(data$f4dwm21_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode22)="Drug Code"
# label(data$f4drugcode22_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage22)="Dosage"
# label(data$f4dosage22_mv)="Missing value verified for: Dosage"
# label(data$f4timesday22)="Number of time taken"
# label(data$f4timesday22_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm22)="Times/day, week or month"
# label(data$f4dwm22_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode23)="Drug Code"
# label(data$f4drugcode23_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage23)="Dosage"
# label(data$f4dosage23_mv)="Missing value verified for: Dosage"
# label(data$f4timesday23)="Number of time taken"
# label(data$f4timesday23_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm23)="Times/day, week or month"
# label(data$f4dwm23_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode24)="Drug Code"
# label(data$f4drugcode24_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage24)="Dosage"
# label(data$f4dosage24_mv)="Missing value verified for: Dosage"
# label(data$f4timesday24)="Number of time taken"
# label(data$f4timesday24_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm24)="Times/day, week or month"
# label(data$f4dwm24_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode25)="Drug Code"
# label(data$f4drugcode25_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage25)="Dosage"
# label(data$f4dosage25_mv)="Missing value verified for: Dosage"
# label(data$f4timesday25)="Number of time taken"
# label(data$f4timesday25_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm25)="Times/day, week or month"
# label(data$f4dwm25_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode26)="Drug Code"
# label(data$f4drugcode26_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage26)="Dosage"
# label(data$f4dosage26_mv)="Missing value verified for: Dosage"
# label(data$f4timesday26)="Number of time taken"
# label(data$f4timesday26_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm26)="Times/day, week or month"
# label(data$f4dwm26_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode27)="Drug Code"
# label(data$f4drugcode27_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage27)="Dosage"
# label(data$f4dosage27_mv)="Missing value verified for: Dosage"
# label(data$f4timesday27)="Number of time taken"
# label(data$f4timesday27_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm27)="Times/day, week or month"
# label(data$f4dwm27_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode28)="Drug Code"
# label(data$f4drugcode28_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage28)="Dosage"
# label(data$f4dosage28_mv)="Missing value verified for: Dosage"
# label(data$f4timesday28)="Number of time taken"
# label(data$f4timesday28_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm28)="Times/day, week or month"
# label(data$f4dwm28_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode29)="Drug Code"
# label(data$f4drugcode29_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage29)="Dosage"
# label(data$f4dosage29_mv)="Missing value verified for: Dosage"
# label(data$f4timesday29)="Number of time taken"
# label(data$f4timesday29_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm29)="Times/day, week or month"
# label(data$f4dwm29_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4drugcode30)="Drug Code"
# label(data$f4drugcode30_mv)="Missing value verified for: Drug Code"
# label(data$f4dosage30)="Dosage"
# label(data$f4dosage30_mv)="Missing value verified for: Dosage"
# label(data$f4timesday30)="Number of time taken"
# label(data$f4timesday30_mv)="Missing value verified for: Number of time taken"
# label(data$f4dwm30)="Times/day, week or month"
# label(data$f4dwm30_mv)="Missing value verified for: Times/day, week or month"
# label(data$f4meds_encounter_id)="Encounter ID"
# label(data$f4meds_complete)="Complete?"
# label(data$f5rcntillns)="Recent Illnesses Rqring Med Att."
# label(data$f5rcntillns_mv)="Missing value verified for: Recent Illnesses Rqring Med Att."
# label(data$f5outpat)="Out Patient?"
# label(data$f5outpat_mv)="Missing value verified for: Out Patient?"
# label(data$f5outpatspec_mv)="Missing value verified for: Out Patient - Specify"
# label(data$f5hospwosrg)="Hosp w/o Surgery"
# label(data$f5hospwosrg_mv)="Missing value verified for: Hosp w/o Surgery"
# label(data$f5hospwoexp_mv)="Missing value verified for: Hosp w/o Surgery - Exp"
# label(data$f5hospwsrg)="Hosp w/Surgery"
# label(data$f5hospwsrg_mv)="Missing value verified for: Hosp w/Surgery"
# label(data$f5hospwexp_mv)="Missing value verified for: Hosp w/Surgery - Exp"
# label(data$f5medprbsdoc)="Any Docmntd Meds Prob since Last Phys?"
# label(data$f5medprbsdoc_mv)="Missing value verified for: Any Docmntd Meds Prob since Last Phys?"
# label(data$f5crnryartds)="Coronary artery disease?"
# label(data$f5crnryartds_mv)="Missing value verified for: Coronary artery disease?"
# label(data$f5cancer)="Cancer?"
# label(data$f5cancer_mv)="Missing value verified for: Cancer?"
# label(data$f5cancerexp_mv)="Missing value verified for: Cancer - Explain"
# label(data$f5crbrlvasds)="Cerebral Vascular disease?"
# label(data$f5crbrlvasds_mv)="Missing value verified for: Cerebral Vascular disease?"
# label(data$f5peripvasds)="Peripheral vascular disease?"
# label(data$f5peripvasds_mv)="Missing value verified for: Peripheral vascular disease?"
# label(data$f5hypertensn)="Hypertension?"
# label(data$f5hypertensn_mv)="Missing value verified for: Hypertension?"
# label(data$f5seizures)="Seizures?"
# label(data$f5seizures_mv)="Missing value verified for: Seizures?"
# label(data$f5gnitrnryds)="Genitourinary disease?"
# label(data$f5gnitrnryds_mv)="Missing value verified for: Genitourinary disease?"
# label(data$f5gnitrnyexp_mv)="Missing value verified for: GD Explain"
# label(data$f5lungdsz)="Lung Disease?"
# label(data$f5lungdsz_mv)="Missing value verified for: Lung Disease?"
# label(data$f5majsurg)="Major surgery?"
# label(data$f5majsurg_mv)="Missing value verified for: Major surgery?"
# label(data$f5majsurgexp_mv)="Missing value verified for: MS Explain"
# label(data$f5othrmeddia)="Other medical diagnosis?"
# label(data$f5othrmeddia_mv)="Missing value verified for: Other medical diagnosis?"
# label(data$f5othrmedexp_mv)="Missing value verified for: OMD Explain"
# label(data$f5lungs)="Lungs"
# label(data$f5lungs_mv)="Missing value verified for: Lungs"
# label(data$f5lungsexp_mv)="Missing value verified for: Explain"
# label(data$f5heart)="Heart"
# label(data$f5heart_mv)="Missing value verified for: Heart"
# label(data$f5heartexp_mv)="Missing value verified for: Explain"
# label(data$f5skin)="Skin"
# label(data$f5skin_mv)="Missing value verified for: Skin"
# label(data$f5skinexp_mv)="Missing value verified for: Explain"
# label(data$f5edema)="Edema"
# label(data$f5edema_mv)="Missing value verified for: Edema"
# label(data$f5edemaexp_mv)="Missing value verified for: Explain"
# label(data$f5hape_encounter_id)="Encounter ID"
# label(data$f5histandphysexm_complete)="Complete?"
# label(data$f6ivsites)="IV Sites"
# label(data$f6ivsites_mv)="Missing value verified for: IV Sites"
# label(data$f6weight)="Weight (kg)"
# label(data$f6weight_mv)="Missing value verified for: Weight (kg)"
# label(data$f6u0strttime)="U0 start time (24 hr. clock)"
# label(data$f6u0strttime_mv)="Missing value verified for: U0 start time (24 hr. clock)"
# label(data$f6u0volume)="U0 volume (ml)"
# label(data$f6u0volume_mv)="Missing value verified for: U0 volume (ml)"
# label(data$f6bldprsrarm)="Bld Prsr (arm used)"
# label(data$f6bldprsrarm_mv)="Missing value verified for: Bld Prsr (arm used)"
# label(data$f6frstmsrsty)="1st measure (systolic)"
# label(data$f6frstmsrsty_mv)="Missing value verified for: 1st measure (systolic)"
# label(data$f6frstmsrdia)="1st measure (diastolic)"
# label(data$f6frstmsrdia_mv)="Missing value verified for: 1st measure (diastolic)"
# label(data$f6scndmsrsys)="2nd measure (systolic)"
# label(data$f6scndmsrsys_mv)="Missing value verified for: 2nd measure (systolic)"
# label(data$f6scndmsrdia)="2nd measure (diastolic)"
# label(data$f6scndmsrdia_mv)="Missing value verified for: 2nd measure (diastolic)"
# label(data$f6heartrate)="Heart rate (beats/min)"
# label(data$f6heartrate_mv)="Missing value verified for: Heart rate (beats/min)"
# label(data$f6u0h20)="U0 H20 + 10cc/kg load (ml)"
# label(data$f6u0h20_mv)="Missing value verified for: U0 H20 + 10cc/kg load (ml)"
# label(data$f6fnleqlurtm)="Final equilibration urine time (min)"
# label(data$f6fnleqlurtm_mv)="Missing value verified for: Final equilibration urine time (min)"
# label(data$f6p1time)="P1 time (min)"
# label(data$f6p1time_mv)="Missing value verified for: P1 time (min)"
# label(data$f6urnvol)="Total equilibration period urine volume (ml)"
# label(data$f6urnvol_mv)="Missing value verified for: Total equilibration period urine volume (ml)"
# label(data$f6h20vol)="Total equilibration H20 volume (ml)"
# label(data$f6h20vol_mv)="Missing value verified for: Total equilibration H20 volume (ml)"
# label(data$f6u1time)="U1 time (min)"
# label(data$f6u1time_mv)="Missing value verified for: U1 time (min)"
# label(data$f6p2time)="P2 time (min)"
# label(data$f6p2time_mv)="Missing value verified for: P2 time (min)"
# label(data$f6u1vol)="U1 volume (ml)"
# label(data$f6u1vol_mv)="Missing value verified for: U1 volume (ml)"
# label(data$f6u1h20vol)="U1 H20 volume (ml)"
# label(data$f6u1h20vol_mv)="Missing value verified for: U1 H20 volume (ml)"
# label(data$f6u2time)="U2 time (min)"
# label(data$f6u2time_mv)="Missing value verified for: U2 time (min)"
# label(data$f6p3time)="P3 time (min)"
# label(data$f6p3time_mv)="Missing value verified for: P3 time (min)"
# label(data$f6u2vol)="U2 volume (ml)"
# label(data$f6u2vol_mv)="Missing value verified for: U2 volume (ml)"
# label(data$f6u2h20vol)="U2 H20 volume (ml)"
# label(data$f6u2h20vol_mv)="Missing value verified for: U2 H20 volume (ml)"
# label(data$f6u3time)="U3 time (min)"
# label(data$f6u3time_mv)="Missing value verified for: U3 time (min)"
# label(data$f6p4time)="P4 time (min)"
# label(data$f6p4time_mv)="Missing value verified for: P4 time (min)"
# label(data$f6u3vol)="U3 volume (ml)"
# label(data$f6u3vol_mv)="Missing value verified for: U3 volume (ml)"
# label(data$f6u3h20vol)="U3 H20 volume (ml)"
# label(data$f6u3h20vol_mv)="Missing value verified for: U3 H20 volume (ml)"
# label(data$f6u4time)="U4 time (min)"
# label(data$f6u4time_mv)="Missing value verified for: U4 time (min)"
# label(data$f6p5time)="P5 time (min)"
# label(data$f6p5time_mv)="Missing value verified for: P5 time (min)"
# label(data$f6u4vol)="U4 volume (ml)"
# label(data$f6u4vol_mv)="Missing value verified for: U4 volume (ml)"
# label(data$f6u4h20vol)="U4 H20 volume (ml)"
# label(data$f6u4h20vol_mv)="Missing value verified for: U4 H20 volume (ml)"
# label(data$f6u5time)="U5 time (min)"
# label(data$f6u5time_mv)="Missing value verified for: U5 time (min)"
# label(data$f6p6time)="P6 time (min)"
# label(data$f6p6time_mv)="Missing value verified for: P6 time (min)"
# label(data$f6u5vol)="U5 volume (ml)"
# label(data$f6u5vol_mv)="Missing value verified for: U5 volume (ml)"
# label(data$f6u5h20vol)="U5 H2O volume (ml)"
# label(data$f6u5h20vol_mv)="Missing value verified for: U5 H2O volume (ml)"
# label(data$f6urndate)="Urinalysis date"
# label(data$f6urndate_mv)="Missing value verified for: Urinalysis date"
# label(data$f6glucose)="Glucose"
# label(data$f6glucose_mv)="Missing value verified for: Glucose"
# label(data$f6bilirubin)="Bilirubin"
# label(data$f6bilirubin_mv)="Missing value verified for: Bilirubin"
# label(data$f6ketones)="Ketones"
# label(data$f6ketones_mv)="Missing value verified for: Ketones"
# label(data$f6specgrvty)="Specific gravity"
# label(data$f6specgrvty_mv)="Missing value verified for: Specific gravity"
# label(data$f6bloodocult)="Blood, occult"
# label(data$f6bloodocult_mv)="Missing value verified for: Blood, occult"
# label(data$f6ph)="pH"
# label(data$f6ph_mv)="Missing value verified for: pH"
# label(data$f6protein)="Protein"
# label(data$f6protein_mv)="Missing value verified for: Protein"
# label(data$f6casts)="Casts"
# label(data$f6casts_mv)="Missing value verified for: Casts"
# label(data$f6rbccast)="RBC"
# label(data$f6rbccast_mv)="Missing value verified for: RBC"
# label(data$f6wbccast)="WBC"
# label(data$f6wbccast_mv)="Missing value verified for: WBC"
# label(data$f6hyalinecast)="Hyaline"
# label(data$f6hyalinecast_mv)="Missing value verified for: Hyaline"
# label(data$f6granularcast)="Granular"
# label(data$f6granularcast_mv)="Missing value verified for: Granular"
# label(data$f6othercast)="Other "
# label(data$f6othercast_mv)="Missing value verified for: Other "
# label(data$f6othercast_exp_mv)="Missing value verified for: Other (specify)"
# label(data$f6wbc)="WBC"
# label(data$f6wbc_mv)="Missing value verified for: WBC"
# label(data$f6rbc)="RBC"
# label(data$f6rbc_mv)="Missing value verified for: RBC"
# label(data$f6epithcells)="Epithelial cells"
# label(data$f6epithcells_mv)="Missing value verified for: Epithelial cells"
# label(data$f6bacteria)="Bacteria"
# label(data$f6bacteria_mv)="Missing value verified for: Bacteria"
# label(data$f6rc_encounter_id)="Encounter ID"
# label(data$f6renalclrnc_complete)="Complete?"
# label(data$f8visitdate)="Date sample collected"
# label(data$f8visitdate_mv)="Missing value verified for: Date sample collected"
# label(data$f8analdate)="Date sample analyzed"
# label(data$f8analdate_mv)="Missing value verified for: Date sample analyzed"
# label(data$f8alti)="ALTI (U/L)"
# label(data$f8alti_mv)="Missing value verified for: ALTI (U/L)"
# label(data$f8glucose)="Glucose (mg/dL)"
# label(data$f8glucose_mv)="Missing value verified for: Glucose (mg/dL)"
# label(data$f8bun)="BUN (mg/dL)"
# label(data$f8bun_mv)="Missing value verified for: BUN (mg/dL)"
# label(data$f8creatinine)="Creatinine (mg/dL)"
# label(data$f8creatinine_mv)="Missing value verified for: Creatinine (mg/dL)"
# label(data$f8uricacid)="Uric acid (mg/dL)"
# label(data$f8uricacid_mv)="Missing value verified for: Uric acid (mg/dL)"
# label(data$f8sodium)="Sodium (mmol/L)"
# label(data$f8sodium_mv)="Missing value verified for: Sodium (mmol/L)"
# label(data$f8potassium)="Potassium (mmol/L)"
# label(data$f8potassium_mv)="Missing value verified for: Potassium (mmol/L)"
# label(data$f8chloride)="Chloride (mmol/L)"
# label(data$f8chloride_mv)="Missing value verified for: Chloride (mmol/L)"
# label(data$f8calcium)="Calcium (mg/dL)"
# label(data$f8calcium_mv)="Missing value verified for: Calcium (mg/dL)"
# label(data$f8totlprotn)="Total protein (gm/dL)"
# label(data$f8totlprotn_mv)="Missing value verified for: Total protein (gm/dL)"
# label(data$f8albumin)="Albumin (gm/dL)"
# label(data$f8albumin_mv)="Missing value verified for: Albumin (gm/dL)"
# label(data$f8cholestero)="Cholesterol (mg/dL)"
# label(data$f8cholestero_mv)="Missing value verified for: Cholesterol (mg/dL)"
# label(data$f8triglyceri)="Triglycerides (mg/dL)"
# label(data$f8triglyceri_mv)="Missing value verified for: Triglycerides (mg/dL)"
# label(data$f8hdl)="HDL Cholesterol (mg/dl)"
# label(data$f8hdl_mv)="Missing value verified for: HDL Cholesterol (mg/dl)"
# label(data$f8ldl)="LDL Cholesterol (mg/dl)"
# label(data$f8ldl_mv)="Missing value verified for: LDL Cholesterol (mg/dl)"
# label(data$f8totlbilrbn)="Total bilirubin (mg/dL)"
# label(data$f8totlbilrbn_mv)="Missing value verified for: Total bilirubin (mg/dL)"
# label(data$f8dirctbilrbn)="Direct bilirubin (mg/dl)"
# label(data$f8dirctbilrbn_mv)="Missing value verified for: Direct bilirubin (mg/dl)"
# label(data$f8alkp04)="ALK P04 (alkaline phosphatase) (U/L)"
# label(data$f8alkp04_mv)="Missing value verified for: ALK P04 (alkaline phosphatase) (U/L)"
# label(data$f8ast)="AST (SGOT) (U/L)"
# label(data$f8ast_mv)="Missing value verified for: AST (SGOT) (U/L)"
# label(data$f8sc_encounter_id)="Encounter ID"
# label(data$f8smacchempanl_complete)="Complete?"
# label(data$f9dateanal)="Date of sample collection"
# label(data$f9dateanal_mv)="Missing value verified for: Date of sample collection"
# label(data$f9analyzed_date)="Date specimen analyzed"
# label(data$f9analyzed_date_mv)="Missing value verified for: Date specimen analyzed"
# label(data$f9wbc)="WBC (x 10^3)"
# label(data$f9wbc_mv)="Missing value verified for: WBC (x 10^3)"
# label(data$f9rbc)="RBC (x 10^6)"
# label(data$f9rbc_mv)="Missing value verified for: RBC (x 10^6)"
# label(data$f9hemoglob)="Hemoglobin (g/dL)"
# label(data$f9hemoglob_mv)="Missing value verified for: Hemoglobin (g/dL)"
# label(data$f9hematocrit)="Hematocrit (%)"
# label(data$f9hematocrit_mv)="Missing value verified for: Hematocrit (%)"
# label(data$f9mcv)="MCV (fL)"
# label(data$f9mcv_mv)="Missing value verified for: MCV (fL)"
# label(data$f9mch)="MCH (pg)"
# label(data$f9mch_mv)="Missing value verified for: MCH (pg)"
# label(data$f9mchc)="MCHC (%)"
# label(data$f9mchc_mv)="Missing value verified for: MCHC (%)"
# label(data$f9platelets)="Platelets (x 10^3)"
# label(data$f9platelets_mv)="Missing value verified for: Platelets (x 10^3)"
# label(data$f9lymphinst)="Lymph, Inst %"
# label(data$f9lymphinst_mv)="Missing value verified for: Lymph, Inst %"
# label(data$f9monoinst)="Mono, Inst %"
# label(data$f9monoinst_mv)="Missing value verified for: Mono, Inst %"
# label(data$f9neutinst)="Neut, Inst %"
# label(data$f9neutinst_mv)="Missing value verified for: Neut, Inst %"
# label(data$f9eosinst)="Eos, Inst %"
# label(data$f9eosinst_mv)="Missing value verified for: Eos, Inst %"
# label(data$f9basoinst)="Baso, Inst %"
# label(data$f9basoinst_mv)="Missing value verified for: Baso, Inst %"
# label(data$f9gran)="Gran, %"
# label(data$f9gran_mv)="Missing value verified for: Gran, %"
# label(data$f9pt)="PROTIME (PT)"
# label(data$f9pt_mv)="Missing value verified for: PROTIME (PT)"
# label(data$f9inr)="INR"
# label(data$f9inr_mv)="Missing value verified for: INR"
# label(data$f9ptt)="PTT"
# label(data$f9ptt_mv)="Missing value verified for: PTT"
# label(data$f9cbc_encounter_id)="Encounter ID"
# label(data$f9cbc_complete)="Complete?"
# label(data$f16qcnumber)="Quality control number"
# label(data$f16qcnumber_mv)="Missing value verified for: Quality control number"
# label(data$f16analdate)="Date sample analyzed"
# label(data$f16analdate_mv)="Missing value verified for: Date sample analyzed"
# label(data$f16glucose)="Glucose (mg/dL)"
# label(data$f16glucose_mv)="Missing value verified for: Glucose (mg/dL)"
# label(data$f16bun)="BUN (mg/dL)"
# label(data$f16bun_mv)="Missing value verified for: BUN (mg/dL)"
# label(data$f16creatinine)="Creatinine (mg/dL)"
# label(data$f16creatinine_mv)="Missing value verified for: Creatinine (mg/dL)"
# label(data$f16sodium)="Sodium (mmol/L)"
# label(data$f16sodium_mv)="Missing value verified for: Sodium (mmol/L)"
# label(data$f16potassium)="Potassium (mmol/L)"
# label(data$f16potassium_mv)="Missing value verified for: Potassium (mmol/L)"
# label(data$f16chloride)="Chloride (mmol/L)"
# label(data$f16chloride_mv)="Missing value verified for: Chloride (mmol/L)"
# label(data$f16calcium)="Calcium (mg/dL)"
# label(data$f16calcium_mv)="Missing value verified for: Calcium (mg/dL)"
# label(data$f16phosphorus)="Phosphorus (mg/dL)"
# label(data$f16phosphorus_mv)="Missing value verified for: Phosphorus (mg/dL)"
# label(data$f16totlprotn)="Total protein (gm/dL)"
# label(data$f16totlprotn_mv)="Missing value verified for: Total protein (gm/dL)"
# label(data$f16albumin)="Albumin (gm/dL)"
# label(data$f16albumin_mv)="Missing value verified for: Albumin (gm/dL)"
# label(data$f16cholestero)="Cholesterol (mg/dL)"
# label(data$f16cholestero_mv)="Missing value verified for: Cholesterol (mg/dL)"
# label(data$f16triglyceri)="Triglycerides (mg/dL)"
# label(data$f16triglyceri_mv)="Missing value verified for: Triglycerides (mg/dL)"
# label(data$f16totlbilrbn)="Total bilirubin (mg/dL)"
# label(data$f16totlbilrbn_mv)="Missing value verified for: Total bilirubin (mg/dL)"
# label(data$f16alkp04)="ALK P04 (alkaline phosphatase) (U/L)"
# label(data$f16alkp04_mv)="Missing value verified for: ALK P04 (alkaline phosphatase) (U/L)"
# label(data$f16ast)="AST (SGOT) (U/L)"
# label(data$f16ast_mv)="Missing value verified for: AST (SGOT) (U/L)"
# label(data$f16sqc_encounter_id)="Encounter ID"
# label(data$f16smacqc_complete)="Complete?"
# label(data$f17qcnumber)="Quality control number"
# label(data$f17qcnumber_mv)="Missing value verified for: Quality control number"
# label(data$f17analdate)="Date specimen analyzed"
# label(data$f17analdate_mv)="Missing value verified for: Date specimen analyzed"
# label(data$f17wbc)="WBC (x 10^3)"
# label(data$f17wbc_mv)="Missing value verified for: WBC (x 10^3)"
# label(data$f17rbc)="RBC (x 10^6)"
# label(data$f17rbc_mv)="Missing value verified for: RBC (x 10^6)"
# label(data$f17hemoglob)="Hemoglobin (g/dL)"
# label(data$f17hemoglob_mv)="Missing value verified for: Hemoglobin (g/dL)"
# label(data$f17hematocrit)="Hematocrit (%)"
# label(data$f17hematocrit_mv)="Missing value verified for: Hematocrit (%)"
# label(data$f17mcv)="MCV (femtoliter)"
# label(data$f17mcv_mv)="Missing value verified for: MCV (femtoliter)"
# label(data$f17mch)="MCH (pg)"
# label(data$f17mch_mv)="Missing value verified for: MCH (pg)"
# label(data$f17mchc)="MCHC (%)"
# label(data$f17mchc_mv)="Missing value verified for: MCHC (%)"
# label(data$f17platelets)="Platelets (x 10^3)"
# label(data$f17platelets_mv)="Missing value verified for: Platelets (x 10^3)"
# label(data$f17cbcqc_encounter_id)="Encounter ID"
# label(data$f17cbcqc_complete)="Complete?"
# label(data$f13dtlastvst)="Date of Last Visit"
# label(data$f13dtlastvst_mv)="Missing value verified for: Date of Last Visit"
# label(data$f13typevist)="Type visit"
# label(data$f13typevist_mv)="Missing value verified for: Type visit"
# label(data$f13tstintrvl)="Test interval"
# label(data$f13tstintrvl_mv)="Missing value verified for: Test interval"
# label(data$f13endstgrnl)="End stage renal disease"
# label(data$f13endstgrnl_mv)="Missing value verified for: End stage renal disease"
# label(data$f13esrddate)="Date of diagnosis"
# label(data$f13esrddate_mv)="Missing value verified for: Date of diagnosis"
# label(data$f13nondiabkid)="Non-diabetic kidney disease"
# label(data$f13nondiabkid_mv)="Missing value verified for: Non-diabetic kidney disease"
# label(data$f13ndkddate)="Date of diagnosis"
# label(data$f13ndkddate_mv)="Missing value verified for: Date of diagnosis"
# label(data$f13imprdbldr)="Impaired bladder functioning"
# label(data$f13imprdbldr_mv)="Missing value verified for: Impaired bladder functioning"
# label(data$f13ibfdate)="Date of diagnosis"
# label(data$f13ibfdate_mv)="Missing value verified for: Date of diagnosis"
# label(data$f13cngstvhrt)="Congestive heart failure"
# label(data$f13cngstvhrt_mv)="Missing value verified for: Congestive heart failure"
# label(data$f13chfdate)="Date of diagnosis"
# label(data$f13chfdate_mv)="Missing value verified for: Date of diagnosis"
# label(data$f13ascites)="Ascites"
# label(data$f13ascites_mv)="Missing value verified for: Ascites"
# label(data$f13ascitesdt)="Date of diagnosis"
# label(data$f13ascitesdt_mv)="Missing value verified for: Date of diagnosis"
# label(data$f13aceinhib)="Potential reactions to medication"
# label(data$f13aceinhib_mv)="Missing value verified for: Potential reactions to medication"
# label(data$f13typeofreaction)="Type of reaction"
# label(data$f13typeofreaction_mv)="Missing value verified for: Type of reaction"
# label(data$f13otherdiag)="Other diagnosis"
# label(data$f13otherdiag_mv)="Missing value verified for: Other diagnosis"
# label(data$f13reactionspecify_mv)="Missing value verified for: Other (specify)"
# label(data$f13datereactn)="Date of reaction"
# label(data$f13datereactn_mv)="Missing value verified for: Date of reaction"
# label(data$f13prfrmwdclrnc)="Should the withdrawal clearance study be performed"
# label(data$f13prfrmwdclrnc_mv)="Missing value verified for: Should the withdrawal clearance study be performed"
# label(data$f13prfrmwdclrncdate)="Date scheduled"
# label(data$f13prfrmwdclrncdate_mv)="Missing value verified for: Date scheduled"
# label(data$f13dtstoppntdeclrd)="Date stop point declared"
# label(data$f13dtstoppntdeclrd_mv)="Missing value verified for: Date stop point declared"
# label(data$f13dtfrmclmpltd)="Date form completed"
# label(data$f13dtfrmclmpltd_mv)="Missing value verified for: Date form completed"
# label(data$f13frmcmpltdby_mv)="Missing value verified for: Initials of person who completed form"
# label(data$f13stoppnt_complete)="Complete?"
# label(data$f7qcnumber)="Quality control number"
# label(data$f7qcnumber_mv)="Missing value verified for: Quality control number"
# label(data$f7visitdate)="Date of sample collection"
# label(data$f7visitdate_mv)="Missing value verified for: Date of sample collection"
# label(data$f7tstintrvl)="Test interval"
# label(data$f7tstintrvl_mv)="Missing value verified for: Test interval"
# label(data$f7typesample)="Type of sample"
# label(data$f7typesample_mv)="Missing value verified for: Type of sample"
# label(data$f7labqc_encounter_id)="Encounter ID"
# label(data$f7labqc_complete)="Complete?"
# label(data$f12expdatevis)="Expected date of visit"
# label(data$f12expdatevis_mv)="Missing value verified for: Expected date of visit"
# label(data$f12typevisit)="Type of visit"
# label(data$f12typevisit_mv)="Missing value verified for: Type of visit"
# label(data$f12testintrvl)="Test interval"
# label(data$f12testintrvl_mv)="Missing value verified for: Test interval"
# label(data$f12reasnmissd)="Reason visit was missed"
# label(data$f12reasnmissd_mv)="Missing value verified for: Reason visit was missed"
# label(data$f12otherspec_mv)="Missing value verified for: Reason visit was missed, other (specify)"
# label(data$f12missedvisit_complete)="Complete?"
# label(data$f18visitdate)="Visit date"
# label(data$f18visitdate_mv)="Missing value verified for: Visit date"
# label(data$f18qcnumber)="Quality control number"
# label(data$f18qcnumber_mv)="Missing value verified for: Quality control number"
# label(data$f18hba1c)="HbA1c (%)"
# label(data$f18hba1c_mv)="Missing value verified for: HbA1c (%)"
# label(data$f18p0creatn)="P0 serum creatinine (mg/dL)"
# label(data$f18p0creatn_mv)="Missing value verified for: P0 serum creatinine (mg/dL)"
# label(data$f18p0albumin)="P0 albumin (mg/L)"
# label(data$f18p0albumin_mv)="Missing value verified for: P0 albumin (mg/L)"
# label(data$f18p0igg)="P0 IgG (mg/L)"
# label(data$f18p0igg_mv)="Missing value verified for: P0 IgG (mg/L)"
# label(data$f18u0creatn)="U0 creatinine (g/L)"
# label(data$f18u0creatn_mv)="Missing value verified for: U0 creatinine (g/L)"
# label(data$f18u0albumin)="U0 albumin (mg/L)"
# label(data$f18u0albumin_mv)="Missing value verified for: U0 albumin (mg/L)"
# label(data$f18u0igg)="U0 IgG (mg/L)"
# label(data$f18u0igg_mv)="Missing value verified for: U0 IgG (mg/L)"
# label(data$daescqc_encounter_id)="Encounter ID"
# label(data$f18daesclrncqc_complete)="Complete?"
# label(data$f19qcnumber)="Quality control number"
# label(data$f19qcnumber_mv)="Missing value verified for: Quality control number"
# label(data$f19urnflowrate)="Urine flow rate >= 10 ml/min "
# label(data$f19urnflowrate_mv)="Missing value verified for: Urine flow rate >= 10 ml/min "
# label(data$f19assaydate)="Date of assay"
# label(data$f19assaydate_mv)="Missing value verified for: Date of assay"
# label(data$f19urnioth3qc)="Urine Ioth. 3 (mg/dl) QC Result"
# label(data$f19urnioth3qc_mv)="Missing value verified for: Urine Ioth. 3 (mg/dl) QC Result"
# label(data$f19urnioth3df)="Urine Ioth. 3 (mg/dl) Dilution Factor"
# label(data$f19urnioth3df_mv)="Missing value verified for: Urine Ioth. 3 (mg/dl) Dilution Factor"
# label(data$f19serumioth3qc)="Serum Ioth. 3 (mg/dl) QC Result"
# label(data$f19serumioth3qc_mv)="Missing value verified for: Serum Ioth. 3 (mg/dl) QC Result"
# label(data$f19serumioth3df)="Serum Ioth. 3 (mg/dl) Dilution Factor"
# label(data$f19serumioth3df_mv)="Missing value verified for: Serum Ioth. 3 (mg/dl) Dilution Factor"
# label(data$f19serumioth4qc)="Serum Ioth. 4 (mg/dl) QC Result"
# label(data$f19serumioth4qc_mv)="Missing value verified for: Serum Ioth. 4 (mg/dl) QC Result"
# label(data$f19serumioth4df)="Serum Ioth. 4 (mg/dl) Dilution Factor"
# label(data$f19serumioth4df_mv)="Missing value verified for: Serum Ioth. 4 (mg/dl) Dilution Factor"
# label(data$f19urnpah3qc)="Urine PAH 3 (mg/dl) QC Result"
# label(data$f19urnpah3qc_mv)="Missing value verified for: Urine PAH 3 (mg/dl) QC Result"
# label(data$f19urnpah3df)="Urine PAH 3 (mg/dl) Dilution Factor"
# label(data$f19urnpah3df_mv)="Missing value verified for: Urine PAH 3 (mg/dl) Dilution Factor"
# label(data$f19serumpah3qc)="Serum PAH 3 (mg/dl) QC Result"
# label(data$f19serumpah3qc_mv)="Missing value verified for: Serum PAH 3 (mg/dl) QC Result"
# label(data$f19serumpah3df)="Serum PAH 3 (mg/dl) Dilution Factor"
# label(data$f19serumpah3df_mv)="Missing value verified for: Serum PAH 3 (mg/dl) Dilution Factor"
# label(data$f19serumpah4qc)="Serum PAH 4 (mg/dl) QC Result"
# label(data$f19serumpah4qc_mv)="Missing value verified for: Serum PAH 4 (mg/dl) QC Result"
# label(data$f19serumpah4df)="Serum PAH 4 (mg/dl) Dilution Factor"
# label(data$f19serumpah4df_mv)="Missing value verified for: Serum PAH 4 (mg/dl) Dilution Factor"
# label(data$f19rfqcr_encounter_id)="Encounter ID"
# label(data$f19renlfuncqcrslt_complete)="Complete?"
# label(data$phxlabs_visitdate)="Visit Date"
# label(data$phxlabs_visitdate_mv)="Missing value verified for: Visit Date"
# label(data$scr)="Serum Creatinine"
# label(data$scr_mv)="Missing value verified for: Serum Creatinine"
# label(data$gfr)="Glomerular Filtration Rate (GFR)"
# label(data$gfr_mv)="Missing value verified for: Glomerular Filtration Rate (GFR)"
# label(data$uricacid)="Uric acid (mg/dL)"
# label(data$uricacid_mv)="Missing value verified for: Uric acid (mg/dL)"
# label(data$hba1a)="Hemoglobin A1a"
# label(data$hba1a_mv)="Missing value verified for: Hemoglobin A1a"
# label(data$hba1b)="Hemoglobin A1b"
# label(data$hba1b_mv)="Missing value verified for: Hemoglobin A1b"
# label(data$hba1c)="Hemoglobin A1c"
# label(data$hba1c_mv)="Missing value verified for: Hemoglobin A1c"
# label(data$hba1o)="Hemoglobin A1o"
# label(data$hba1o_mv)="Missing value verified for: Hemoglobin A1o"
# label(data$hbf)="Hemoglobin F"
# label(data$hbf_mv)="Missing value verified for: Hemoglobin F"
# label(data$p0albumin)="P0 Albumin"
# label(data$p0albumin_mv)="Missing value verified for: P0 Albumin"
# label(data$p0gluc)="P0 Glucose"
# label(data$p0gluc_mv)="Missing value verified for: P0 Glucose"
# label(data$p0igg)="P0 IgG"
# label(data$p0igg_mv)="Missing value verified for: P0 IgG"
# label(data$p0oncoticprsr)="P0 Oncotic Pressure"
# label(data$p0oncoticprsr_mv)="Missing value verified for: P0 Oncotic Pressure"
# label(data$p3albumin)="P3 Albumin"
# label(data$p3albumin_mv)="Missing value verified for: P3 Albumin"
# label(data$p3igg)="P3 IgG"
# label(data$p3igg_mv)="Missing value verified for: P3 IgG"
# label(data$p3oncoticprsr)="P3 Oncotic Pressure"
# label(data$p3oncoticprsr_mv)="Missing value verified for: P3 Oncotic Pressure"
# label(data$pahclrnc)="PAH Clearance"
# label(data$pahclrnc_mv)="Missing value verified for: PAH Clearance"
# label(data$u0igg)="U0 IgG"
# label(data$u0igg_is_below_limit)="U0 IgG is below detection limit of assay?"
# label(data$u0igg_is_below_limit_mv)="Missing value verified for: U0 IgG is below detection limit of assay?"
# label(data$u0igg_mv)="Missing value verified for: U0 IgG"
# label(data$u3flow_gfr)="U3 Flow"
# label(data$u3flow_gfr_mv)="Missing value verified for: U3 Flow"
# label(data$u3gfr)="U3 GFR"
# label(data$u3gfr_mv)="Missing value verified for: U3 GFR"
# label(data$u3use_gfr)="U3 Use"
# label(data$u3use_gfr_mv)="Missing value verified for: U3 Use"
# label(data$u3albumin)="U3 Albumin"
# label(data$u3albumin_is_below_limit)="U3 Albumin is below detection limit of assay?"
# label(data$u3albumin_is_below_limit_mv)="Missing value verified for: U3 Albumin is below detection limit of assay?"
# label(data$u3albumin_mv)="Missing value verified for: U3 Albumin"
# label(data$u3igg)="U3 IgG"
# label(data$u3igg_is_below_limit)="U3 IgG is below detection limit of assay?"
# label(data$u3igg_is_below_limit_mv)="Missing value verified for: U3 IgG is below detection limit of assay?"
# label(data$u3igg_mv)="Missing value verified for: U3 IgG"
# label(data$ualb)="Urine Albumin"
# label(data$ualb_is_below_limit)="Urine Albumin is below detection limit of assay?"
# label(data$ualb_is_below_limit_mv)="Missing value verified for: Urine Albumin is below detection limit of assay?"
# label(data$ualb_mv)="Missing value verified for: Urine Albumin"
# label(data$ucr)="Urine Creatinine"
# label(data$ucr_mv)="Missing value verified for: Urine Creatinine"
# label(data$phxlabs_encounter_id)="Encounter ID"
# label(data$phoenix_labs_complete)="Complete?"
# label(data$cantpv_2)="Date of test"
# label(data$cantpv_2_mv)="Missing value verified for: Date of test"
# label(data$cantpv_3)="Has the patient had anything to eat or drink (except water) in the past 8 hours If yes, subject is not eligible for testing and should be rescheduled."
# label(data$cantpv_3_mv)="Missing value verified for: Has the patient had anything to eat or drink (except water) in the past 8 hours If yes, subject is not eligible for testing and should be rescheduled."
# label(data$cantpv_4)="Any tobacco products in the past 8 hours?  If yes, subject should be rescheduled."
# label(data$cantpv_4_mv)="Missing value verified for: Any tobacco products in the past 8 hours?  If yes, subject should be rescheduled."
# label(data$cantpv_5)="Any vigorous exercise in the last 24 hours? (Any exercise not part of patients daily routine, i.e., routine jogging O.K., but marathon running is not.  NO exercise morning of test) If yes, subject should be rescheduled."
# label(data$cantpv_5_mv)="Missing value verified for: Any vigorous exercise in the last 24 hours? (Any exercise not part of patients daily routine, i.e., routine jogging O.K., but marathon running is not.  NO exercise morning of test) If yes, subject should be rescheduled."
# label(data$cantpv_6)="Any significant emotional upset in the last 24 hours? (Depression, crying episodes, anxiety from personal trauma) If yes, subject should be rescheduled."
# label(data$cantpv_6_mv)="Missing value verified for: Any significant emotional upset in the last 24 hours? (Depression, crying episodes, anxiety from personal trauma) If yes, subject should be rescheduled."
# label(data$cantpv_7)="Acute illness in last 48 hours? (cold, flu, fever, measles, etc). If yes, subject should be rescheduled."
# label(data$cantpv_7_mv)="Missing value verified for: Acute illness in last 48 hours? (cold, flu, fever, measles, etc). If yes, subject should be rescheduled."
# label(data$cantpv_8)="Any hypoglycemic episodes (blood sugar less than 70 or symptoms of low blood sugar) in the past 8 hours?   If yes, subject is not eligible for testing.  Reschedule for another day."
# label(data$cantpv_8_mv)="Missing value verified for: Any hypoglycemic episodes (blood sugar less than 70 or symptoms of low blood sugar) in the past 8 hours?   If yes, subject is not eligible for testing.  Reschedule for another day."
# label(data$cantpv_9)="Fasting blood sugar value (finger stick method O.K.)"
# label(data$cantpv_9_mv)="Missing value verified for: Fasting blood sugar value (finger stick method O.K.)"
# label(data$cantpv_10)="Fasting blood sugar value below 50 and/or signs or symptoms of hypoglycemia? If yes, subject is not eligible for testing.  Reschedule for another day."
# label(data$cantpv_10_mv)="Missing value verified for: Fasting blood sugar value below 50 and/or signs or symptoms of hypoglycemia? If yes, subject is not eligible for testing.  Reschedule for another day."
# label(data$cantpv_11)="Are any items above (A1 to A8) answered Yes?"
# label(data$cantpv_12_mv)="Missing value verified for: If yes, and you plan to proceed with testing, please explain why:"
# label(data$cantpv_13)="Were any medications (other than insulin or other glucose lowering medications) taken in the last 8 hours?"
# label(data$cantpv_13_mv)="Missing value verified for: Were any medications (other than insulin or other glucose lowering medications) taken in the last 8 hours?"
# label(data$cantpv_14_mv)="Missing value verified for: If yes, list any medications taken in the last 8 hours.  Include prescribed and OTC medications:"
# label(data$cantpv_15)="Were any glucose lowering medications taken in the last 8 hours?"
# label(data$cantpv_15_mv)="Missing value verified for: Were any glucose lowering medications taken in the last 8 hours?"
# label(data$cantpv_16_mv)="Missing value verified for: If yes, list any glucose lowering medications taken in the last 8 hours:"
# label(data$cantpv_17)="Was any insulin (aside from basal insulin via insulin pump) taken in the 1 hour prior to planned start of testing?   If yes, wait at least 1 hour from time of last insulin dose before beginning CAN test."
# label(data$cantpv_17_mv)="Missing value verified for: Was any insulin (aside from basal insulin via insulin pump) taken in the 1 hour prior to planned start of testing?   If yes, wait at least 1 hour from time of last insulin dose before beginning CAN test."
# label(data$cantpv_18)="Does the patient take a beta-blocker?"
# label(data$cantpv_18_mv)="Missing value verified for: Does the patient take a beta-blocker?"
# label(data$cantpv_19_mv)="Missing value verified for: Completed by"
# label(data$cantpv_20)="Resting Heart Rate (5 minute) blood pressure - systolic"
# label(data$cantpv_20_mv)="Missing value verified for: Resting Heart Rate (5 minute) blood pressure - systolic"
# label(data$cantpv_21)="Resting Heart Rate (5 minute) blood pressure - diastolic"
# label(data$cantpv_21_mv)="Missing value verified for: Resting Heart Rate (5 minute) blood pressure - diastolic"
# label(data$cantpv_22_mv)="Missing value verified for: Comments"
# label(data$cantpv_23)="HRV (6 minute paced breathing) - systolic start"
# label(data$cantpv_23_mv)="Missing value verified for: HRV (6 minute paced breathing) - systolic start"
# label(data$cantpv_24)="HRV (6 minute paced breathing) - diatolic start"
# label(data$cantpv_24_mv)="Missing value verified for: HRV (6 minute paced breathing) - diatolic start"
# label(data$cantpv_25)="HRV (6 minute paced breathing) - systolic end"
# label(data$cantpv_25_mv)="Missing value verified for: HRV (6 minute paced breathing) - systolic end"
# label(data$cantpv_26)="HRV (6 minute paced breathing) - diastolic end"
# label(data$cantpv_26_mv)="Missing value verified for: HRV (6 minute paced breathing) - diastolic end"
# label(data$cantpv_27_mv)="Missing value verified for: Comments"
# label(data$cantpv_28_mv)="Missing value verified for: Test completed by:"
# label(data$cantpv_29_mv)="Missing value verified for: Certification number:"
# label(data$cantpv_encounter_id)="Encounter ID"
# label(data$cardiac_autonomic_neuropathy_test_preparation_veri_complete)="Complete?"
# label(data$mnsi_cl_2)="Date of assessment"
# label(data$mnsi_cl_2_mv)="Missing value verified for: Date of assessment"
# label(data$mnsi_cl_3)="Normal - right"
# label(data$mnsi_cl_3_mv)="Missing value verified for: Normal - right"
# label(data$mnsi_cl_4___1)="If no, check all that apply (choice=Deformities)"
# label(data$mnsi_cl_4___2)="If no, check all that apply (choice=Dry skin callus)"
# label(data$mnsi_cl_4___3)="If no, check all that apply (choice=Infection)"
# label(data$mnsi_cl_4___4)="If no, check all that apply (choice=Fissure)"
# label(data$mnsi_cl_4___0)="If no, check all that apply (choice=Other)"
# label(data$mnsi_cl_5_mv)="Missing value verified for: Specify:"
# label(data$mnsi_cl_6)="Ulceration - right"
# label(data$mnsi_cl_6_mv)="Missing value verified for: Ulceration - right"
# label(data$mnsi_cl_7)="Ankle reflex - right"
# label(data$mnsi_cl_7_mv)="Missing value verified for: Ankle reflex - right"
# label(data$mnsi_cl_8)="Vibration perception at big toe - right"
# label(data$mnsi_cl_8_mv)="Missing value verified for: Vibration perception at big toe - right"
# label(data$mnsi_cl_9)="Monofilament - right"
# label(data$mnsi_cl_9_mv)="Missing value verified for: Monofilament - right"
# label(data$mnsi_cl_10)="Normal - left"
# label(data$mnsi_cl_10_mv)="Missing value verified for: Normal - left"
# label(data$mnsi_cl_11___1)="If no, check all that apply (choice=Deformities)"
# label(data$mnsi_cl_11___2)="If no, check all that apply (choice=Dry skin callus)"
# label(data$mnsi_cl_11___3)="If no, check all that apply (choice=Infection)"
# label(data$mnsi_cl_11___4)="If no, check all that apply (choice=Fissure)"
# label(data$mnsi_cl_11___0)="If no, check all that apply (choice=Other)"
# label(data$mnsi_cl_12_mv)="Missing value verified for: Specify:"
# label(data$mnsi_cl_13)="Ulceration - left"
# label(data$mnsi_cl_13_mv)="Missing value verified for: Ulceration - left"
# label(data$mnsi_cl_14)="Ankle reflex - left"
# label(data$mnsi_cl_14_mv)="Missing value verified for: Ankle reflex - left"
# label(data$mnsi_cl_15)="Vibration perception at big toe - left"
# label(data$mnsi_cl_15_mv)="Missing value verified for: Vibration perception at big toe - left"
# label(data$mnsi_cl_16)="Monofilament - left"
# label(data$mnsi_cl_16_mv)="Missing value verified for: Monofilament - left"
# label(data$mnsi_cl_17_mv)="Missing value verified for: Notes"
# label(data$mnsi_cl_18)="Score"
# label(data$mnsifc_encounter_id)="Encounter ID"
# label(data$michigan_neuropathy_screening_instrument_for_clini_complete)="Complete?"
# label(data$mnsi_subj_2)="Date"
# label(data$mnsi_subj_2_mv)="Missing value verified for: Date"
# label(data$mnsi_subj_3)="Are your legs and/or feet numb?"
# label(data$mnsi_subj_3_mv)="Missing value verified for: Are your legs and/or feet numb?"
# label(data$mnsi_subj_4)="Do you ever have any burning pain in your legs and/or feet?"
# label(data$mnsi_subj_4_mv)="Missing value verified for: Do you ever have any burning pain in your legs and/or feet?"
# label(data$mnsi_subj_5)="Are your feet too sensitive to touch?"
# label(data$mnsi_subj_5_mv)="Missing value verified for: Are your feet too sensitive to touch?"
# label(data$mnsi_subj_6)="Do you get muscle cramps in your legs and/or feet?"
# label(data$mnsi_subj_6_mv)="Missing value verified for: Do you get muscle cramps in your legs and/or feet?"
# label(data$mnsi_subj_7)="Do you ever have any prickling feelings in your legs or feet?"
# label(data$mnsi_subj_7_mv)="Missing value verified for: Do you ever have any prickling feelings in your legs or feet?"
# label(data$mnsi_subj_8)="Does it hurt when the bed covers touch your skin?"
# label(data$mnsi_subj_8_mv)="Missing value verified for: Does it hurt when the bed covers touch your skin?"
# label(data$mnsi_subj_9)="When you get into the tub or shower, are you able to tell the hot water from the cold water?"
# label(data$mnsi_subj_9_mv)="Missing value verified for: When you get into the tub or shower, are you able to tell the hot water from the cold water?"
# label(data$mnsi_subj_10)="Have you ever had an open sore on your foot?"
# label(data$mnsi_subj_10_mv)="Missing value verified for: Have you ever had an open sore on your foot?"
# label(data$mnsi_subj_11)="Has your doctor ever told you that you have diabetic neuropathy?"
# label(data$mnsi_subj_11_mv)="Missing value verified for: Has your doctor ever told you that you have diabetic neuropathy?"
# label(data$mnsi_subj_12)="Do you feel weak all over most of the time?"
# label(data$mnsi_subj_12_mv)="Missing value verified for: Do you feel weak all over most of the time?"
# label(data$mnsi_subj_13)="Are your symptoms worse at night?"
# label(data$mnsi_subj_13_mv)="Missing value verified for: Are your symptoms worse at night?"
# label(data$mnsi_subj_14)="Do your legs hurt when you walk?"
# label(data$mnsi_subj_14_mv)="Missing value verified for: Do your legs hurt when you walk?"
# label(data$mnsi_subj_15)="Are you able to sense your feet when you walk?"
# label(data$mnsi_subj_15_mv)="Missing value verified for: Are you able to sense your feet when you walk?"
# label(data$mnsi_subj_16)="Is the skin on your feet so dry that it cracks open?"
# label(data$mnsi_subj_16_mv)="Missing value verified for: Is the skin on your feet so dry that it cracks open?"
# label(data$mnsi_subj_17)="Have you ever had an amputation?"
# label(data$mnsi_subj_17_mv)="Missing value verified for: Have you ever had an amputation?"
# label(data$mnsi_subj_18_mv)="Missing value verified for: Notes"
# label(data$mnsi_subj_19)="Total"
# label(data$mnsifs_encounter_id)="Encounter ID"
# label(data$michigan_neuropathy_screening_instrument_for_subje_complete)="Complete?"
# label(data$aspq_1)="Date"
# label(data$aspq_1_mv)="Missing value verified for: Date"
# label(data$aspq_2_mv)="Missing value verified for: Patient initials"
# label(data$aspq_3)="Visit"
# label(data$aspq_3_mv)="Missing value verified for: Visit"
# label(data$aspq_4)="In the past year, have you ever felt faint, dizzy, or goofy or had difficulty thinking soon after standing up from a sitting or lying down position?"
# label(data$aspq_4_mv)="Missing value verified for: In the past year, have you ever felt faint, dizzy, or goofy or had difficulty thinking soon after standing up from a sitting or lying down position?"
# label(data$aspq_5)="When standing up, how frequently do you get these feelings or symptoms?"
# label(data$aspq_5_mv)="Missing value verified for: When standing up, how frequently do you get these feelings or symptoms?"
# label(data$aspq_6)="How would you rate the severity of these feelings or symptoms?"
# label(data$aspq_6_mv)="Missing value verified for: How would you rate the severity of these feelings or symptoms?"
# label(data$aspq_7)="For how long have you been experiencing these feelings or symptoms?"
# label(data$aspq_7_mv)="Missing value verified for: For how long have you been experiencing these feelings or symptoms?"
# label(data$aspq_8)="In the past year, how often have you ended up fainting soon after standing up from a sitting or lying down position?"
# label(data$aspq_8_mv)="Missing value verified for: In the past year, how often have you ended up fainting soon after standing up from a sitting or lying down position?"
# label(data$aspq_9)="How cautious are you about standing up from a sitting or lying down position?"
# label(data$aspq_9_mv)="Missing value verified for: How cautious are you about standing up from a sitting or lying down position?"
# label(data$aspq_10)="What part of the day are these feelings worst?"
# label(data$aspq_10_mv)="Missing value verified for: What part of the day are these feelings worst?"
# label(data$aspq_11_mv)="Missing value verified for: Please specify"
# label(data$aspq_12)="In the past year, have these feelings or symptoms that you have experienced:"
# label(data$aspq_12_mv)="Missing value verified for: In the past year, have these feelings or symptoms that you have experienced:"
# label(data$aspq_14)="Rapid or increased heart rate? (palpitations)?"
# label(data$aspq_14_mv)="Missing value verified for: Rapid or increased heart rate? (palpitations)?"
# label(data$aspq_15)="Sickness to your stomach (nausea) or vomiting?"
# label(data$aspq_15_mv)="Missing value verified for: Sickness to your stomach (nausea) or vomiting?"
# label(data$aspq_16)="A spinning or swimming sensation?"
# label(data$aspq_16_mv)="Missing value verified for: A spinning or swimming sensation?"
# label(data$aspq_17)="Dizziness?"
# label(data$aspq_17_mv)="Missing value verified for: Dizziness?"
# label(data$aspq_18)="Blurred vision?"
# label(data$aspq_18_mv)="Missing value verified for: Blurred vision?"
# label(data$aspq_19)="Feeling of weakness?"
# label(data$aspq_19_mv)="Missing value verified for: Feeling of weakness?"
# label(data$aspq_20)="Feeling shaky or shaking sensation?"
# label(data$aspq_20_mv)="Missing value verified for: Feeling shaky or shaking sensation?"
# label(data$aspq_21)="Feeling anxious or nervous?"
# label(data$aspq_21_mv)="Missing value verified for: Feeling anxious or nervous?"
# label(data$aspq_22)="Turning pale?"
# label(data$aspq_22_mv)="Missing value verified for: Turning pale?"
# label(data$aspq_23)="Clammy feeling to your skin?"
# label(data$aspq_23_mv)="Missing value verified for: Clammy feeling to your skin?"
# label(data$aspq_24)="Do you have any biologic (blood, natural) relatives among your parents, grandparents, brothers, sisters, or children who have frequent dizziness after standing from a sitting or lying down position?"
# label(data$aspq_24_mv)="Missing value verified for: Do you have any biologic (blood, natural) relatives among your parents, grandparents, brothers, sisters, or children who have frequent dizziness after standing from a sitting or lying down position?"
# label(data$aspq_25_mv)="Missing value verified for: If Yes, please list their names and relationship to you."
# label(data$aspq_26)="In the past year, have you ever felt faint, dizzy, or goofy or had difficulty thinking soon after a mean?"
# label(data$aspq_26_mv)="Missing value verified for: In the past year, have you ever felt faint, dizzy, or goofy or had difficulty thinking soon after a mean?"
# label(data$aspq_27)="In the past year, have you ever felt faint, dizzy, or goofy or had difficulty thinking after standing for a long time?"
# label(data$aspq_27_mv)="Missing value verified for: In the past year, have you ever felt faint, dizzy, or goofy or had difficulty thinking after standing for a long time?"
# label(data$aspq_28)="In the past year, have you ever felt faint, dizzy, or goofy or had difficulty thinking during or soon after physical activity or exercise?"
# label(data$aspq_28_mv)="Missing value verified for: In the past year, have you ever felt faint, dizzy, or goofy or had difficulty thinking during or soon after physical activity or exercise?"
# label(data$aspq_29)="In the past year, have you ever felt faint, dizzy, or goofy or had difficulty thinking during or soon after being in a hot bath, shower, tub, or sauna?"
# label(data$aspq_29_mv)="Missing value verified for: In the past year, have you ever felt faint, dizzy, or goofy or had difficulty thinking during or soon after being in a hot bath, shower, tub, or sauna?"
# label(data$aspq_30)="Have you ever felt dizzy or faint or actually fainted when you saw blood or had a blood sample taken?"
# label(data$aspq_30_mv)="Missing value verified for: Have you ever felt dizzy or faint or actually fainted when you saw blood or had a blood sample taken?"
# label(data$aspq_31)="In the past year, have you fainted while passing urine?"
# label(data$aspq_31_mv)="Missing value verified for: In the past year, have you fainted while passing urine?"
# label(data$aspq_32)="In the past year, have you fainted while coughing?"
# label(data$aspq_32_mv)="Missing value verified for: In the past year, have you fainted while coughing?"
# label(data$aspq_34)="In the past year, have you fainted while side of neck?"
# label(data$aspq_34_mv)="Missing value verified for: In the past year, have you fainted while side of neck?"
# label(data$aspq_35)="In the past year, have you fainted before a public speech?"
# label(data$aspq_35_mv)="Missing value verified for: In the past year, have you fainted before a public speech?"
# label(data$aspq_36)="In the past year, have you fainted any other time?"
# label(data$aspq_36_mv)="Missing value verified for: In the past year, have you fainted any other time?"
# label(data$aspq_37_mv)="Missing value verified for: If you checked Yes to any of these questions on fainting, please describe circumstances."
# label(data$aspq_38)="In the past year, have you ever completely lost consciousness after a spell of dizziness?"
# label(data$aspq_38_mv)="Missing value verified for: In the past year, have you ever completely lost consciousness after a spell of dizziness?"
# label(data$aspq_39)="In the past year, have you had any seizures or convulsions?"
# label(data$aspq_39_mv)="Missing value verified for: In the past year, have you had any seizures or convulsions?"
# label(data$aspq_40_mv)="Missing value verified for: If yes, please describe the circumstances."
# label(data$aspq_encounter_id)="Encounter ID"
# label(data$autonomic_symptoms_profile_questionnaire_complete)="Complete?"
# label(data$nsqlq_2_mv)="Missing value verified for: Date"
# label(data$nsqlq_3_mv)="Missing value verified for: Initials"
# label(data$nsqlq_4)="Visit"
# label(data$nsqlq_5)="In the past 4 weeks how often have you experienced burning in your legs or feet?"
# label(data$nsqlq_5_mv)="Missing value verified for: In the past 4 weeks how often have you experienced burning in your legs or feet?"
# label(data$nsqlq_6)="How much bother did this cause you?"
# label(data$nsqlq_6_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_7)="In the past 4 weeks how often have you experienced excessive heat or cold in your legs or feet?"
# label(data$nsqlq_7_mv)="Missing value verified for: In the past 4 weeks how often have you experienced excessive heat or cold in your legs or feet?"
# label(data$nsqlq_8)="How much bother did this cause you?"
# label(data$nsqlq_8_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_9)="In the past 4 weeks how often have you experienced pins and needles in your legs or feet?"
# label(data$nsqlq_9_mv)="Missing value verified for: In the past 4 weeks how often have you experienced pins and needles in your legs or feet?"
# label(data$nsqlq_10)="How much bother did this cause you?"
# label(data$nsqlq_10_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_11)="In the past 4 weeks how often have you experienced shooting or stabbing pain in your legs or feet?"
# label(data$nsqlq_11_mv)="Missing value verified for: In the past 4 weeks how often have you experienced shooting or stabbing pain in your legs or feet?"
# label(data$nsqlq_12)="How much bother did this cause you?"
# label(data$nsqlq_12_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_13)="In the past 4 weeks how often have you experienced throbbing in your legs or feet?"
# label(data$nsqlq_13_mv)="Missing value verified for: In the past 4 weeks how often have you experienced throbbing in your legs or feet?"
# label(data$nsqlq_14)="How much bother did this cause you?"
# label(data$nsqlq_14_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_15)="In the past 4 weeks how often have you experienced sensations in your legs or feet that make them jump?"
# label(data$nsqlq_15_mv)="Missing value verified for: In the past 4 weeks how often have you experienced sensations in your legs or feet that make them jump?"
# label(data$nsqlq_16)="How much bother did this cause you?"
# label(data$nsqlq_16_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_17)="In the past 4 weeks how often have you experienced irritation of the skin caused by something touching your feet, such as bedsheets or socks?"
# label(data$nsqlq_17_mv)="Missing value verified for: In the past 4 weeks how often have you experienced irritation of the skin caused by something touching your feet, such as bedsheets or socks?"
# label(data$nsqlq_18)="How much bother did this cause you?"
# label(data$nsqlq_18_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_19)="Have these painful symptoms reduced your quality of life?"
# label(data$nsqlq_19_mv)="Missing value verified for: Have these painful symptoms reduced your quality of life?"
# label(data$nsqlq_20)="In the past 4 weeks how often have you experienced numbness in your feet?"
# label(data$nsqlq_20_mv)="Missing value verified for: In the past 4 weeks how often have you experienced numbness in your feet?"
# label(data$nsqlq_21)="How much bother did this cause you?"
# label(data$nsqlq_21_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_22)="In the past 4 weeks how often have you experienced inability to feel the difference between hot and cold with your feet?"
# label(data$nsqlq_22_mv)="Missing value verified for: In the past 4 weeks how often have you experienced inability to feel the difference between hot and cold with your feet?"
# label(data$nsqlq_23)="How much bother did this cause you?"
# label(data$nsqlq_23_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_24)="In the past 4 weeks how often have you experienced inability to feel objects with your feet?"
# label(data$nsqlq_24_mv)="Missing value verified for: In the past 4 weeks how often have you experienced inability to feel objects with your feet?"
# label(data$nsqlq_25)="How much bother did this cause you?"
# label(data$nsqlq_25_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_26)="Have these last three symptoms reduced your quality of life?"
# label(data$nsqlq_26_mv)="Missing value verified for: Have these last three symptoms reduced your quality of life?"
# label(data$nsqlq_27)="In the past 4 weeks how often have you experienced weakness in your hands?"
# label(data$nsqlq_27_mv)="Missing value verified for: In the past 4 weeks how often have you experienced weakness in your hands?"
# label(data$nsqlq_28)="How much bother did this cause you?"
# label(data$nsqlq_28_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_29)="In the past 4 weeks how often have you experienced problems with balance or unsteadiness while walking?"
# label(data$nsqlq_29_mv)="Missing value verified for: In the past 4 weeks how often have you experienced problems with balance or unsteadiness while walking?"
# label(data$nsqlq_30)="How much bother did this cause you?"
# label(data$nsqlq_30_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_31)="In the past 4 weeks how often have you experienced problems with balance or unsteadiness while standing?"
# label(data$nsqlq_31_mv)="Missing value verified for: In the past 4 weeks how often have you experienced problems with balance or unsteadiness while standing?"
# label(data$nsqlq_32)="How much bother did this cause you?"
# label(data$nsqlq_32_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_33)="Have these last three symptoms reduced your quality of life?"
# label(data$nsqlq_33_mv)="Missing value verified for: Have these last three symptoms reduced your quality of life?"
# label(data$nsqlq_34)="Are you in PAID WORK?"
# label(data$nsqlq_34_mv)="Missing value verified for: Are you in PAID WORK?"
# label(data$nsqlq_35)="In the past 4 weeks HOW MUCH have your foot problems interferred with your ability to perform your paid work?"
# label(data$nsqlq_35_mv)="Missing value verified for: In the past 4 weeks HOW MUCH have your foot problems interferred with your ability to perform your paid work?"
# label(data$nsqlq_36)="How important is this aspect of your life to you?"
# label(data$nsqlq_36_mv)="Missing value verified for: How important is this aspect of your life to you?"
# label(data$nsqlq_37)="In the past 4 weeks HOW MUCH have your foot problems interferred with your ability to perform tasks around the house or garden?"
# label(data$nsqlq_37_mv)="Missing value verified for: In the past 4 weeks HOW MUCH have your foot problems interferred with your ability to perform tasks around the house or garden?"
# label(data$nsqlq_38)="How important is this aspect of your life to you?"
# label(data$nsqlq_38_mv)="Missing value verified for: How important is this aspect of your life to you?"
# label(data$nsqlq_39)="In the past 4 weeks HOW MUCH have your foot problems interferred with your ability to take part in leisure activities?"
# label(data$nsqlq_39_mv)="Missing value verified for: In the past 4 weeks HOW MUCH have your foot problems interferred with your ability to take part in leisure activities?"
# label(data$nsqlq_40)="How important is this aspect of your life to you?"
# label(data$nsqlq_40_mv)="Missing value verified for: How important is this aspect of your life to you?"
# label(data$nsqlq_41)="Have these changes in daily activities as a result of your foot problems reduced your quality of life?"
# label(data$nsqlq_41_mv)="Missing value verified for: Have these changes in daily activities as a result of your foot problems reduced your quality of life?"
# label(data$nsqlq_42)="In the past 4 weeks how much have your foot problems interferred with your relationships with people close to you?"
# label(data$nsqlq_42_mv)="Missing value verified for: In the past 4 weeks how much have your foot problems interferred with your relationships with people close to you?"
# label(data$nsqlq_43)="How important is this aspect of your life to you?"
# label(data$nsqlq_43_mv)="Missing value verified for: How important is this aspect of your life to you?"
# label(data$nsqlq_44)="In the past 4 weeks have you felt more physically dependent than you would like to be on people close to you as a result of your foot problems?"
# label(data$nsqlq_44_mv)="Missing value verified for: In the past 4 weeks have you felt more physically dependent than you would like to be on people close to you as a result of your foot problems?"
# label(data$nsqlq_45)="How important is this aspect of your life to you?"
# label(data$nsqlq_45_mv)="Missing value verified for: How important is this aspect of your life to you?"
# label(data$nsqlq_46)="In the past 4 weeks have you felt more emotionally dependent than you would like to be on people close to you as a result of your foot problems?"
# label(data$nsqlq_46_mv)="Missing value verified for: In the past 4 weeks have you felt more emotionally dependent than you would like to be on people close to you as a result of your foot problems?"
# label(data$nsqlq_47)="How important is this aspect of your life to you?"
# label(data$nsqlq_47_mv)="Missing value verified for: How important is this aspect of your life to you?"
# label(data$nsqlq_48)="In the past 4 weeks has your role in the family changed as a result of your foot problems?"
# label(data$nsqlq_48_mv)="Missing value verified for: In the past 4 weeks has your role in the family changed as a result of your foot problems?"
# label(data$nsqlq_49)="How important is this aspect of your life to you?"
# label(data$nsqlq_49_mv)="Missing value verified for: How important is this aspect of your life to you?"
# label(data$nsqlq_50)="Have these changes in relationships with other people as a result of your foot problems reduced your quality of life?"
# label(data$nsqlq_50_mv)="Missing value verified for: Have these changes in relationships with other people as a result of your foot problems reduced your quality of life?"
# label(data$nsqlq_51)="People treat me differently from other people as a result of my foot problems."
# label(data$nsqlq_51_mv)="Missing value verified for: People treat me differently from other people as a result of my foot problems."
# label(data$nsqlq_52)="How much bother did this cause you?"
# label(data$nsqlq_52_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_53)="I feel older than my years as a result of my foot problems."
# label(data$nsqlq_53_mv)="Missing value verified for: I feel older than my years as a result of my foot problems."
# label(data$nsqlq_54)="How much bother did this cause you?"
# label(data$nsqlq_54_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_55)="My self-confidence is affected as a result of my foot problems."
# label(data$nsqlq_55_mv)="Missing value verified for: My self-confidence is affected as a result of my foot problems."
# label(data$nsqlq_56)="How much bother did this cause you?"
# label(data$nsqlq_56_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_57)="My foot problems make my life a struggle."
# label(data$nsqlq_57_mv)="Missing value verified for: My foot problems make my life a struggle."
# label(data$nsqlq_58)="How much bother did this cause you?"
# label(data$nsqlq_58_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_59)="I generally feel frustrated because of my foot problems."
# label(data$nsqlq_59_mv)="Missing value verified for: I generally feel frustrated because of my foot problems."
# label(data$nsqlq_60)="How much bother did this cause you?"
# label(data$nsqlq_60_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_61)="My foot problems cause me embarrassment."
# label(data$nsqlq_61_mv)="Missing value verified for: My foot problems cause me embarrassment."
# label(data$nsqlq_62)="How much bother did this cause you?"
# label(data$nsqlq_62_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_63)="How depressed have you felt because of your foot problems?"
# label(data$nsqlq_63_mv)="Missing value verified for: How depressed have you felt because of your foot problems?"
# label(data$nsqlq_64)="How much bother did this cause you?"
# label(data$nsqlq_64_mv)="Missing value verified for: How much bother did this cause you?"
# label(data$nsqlq_65)="Have these feelings about yourself as a  result of your foot problems reduced your quality of life?"
# label(data$nsqlq_65_mv)="Missing value verified for: Have these feelings about yourself as a  result of your foot problems reduced your quality of life?"
# label(data$nsqlq_66)="Overall, I would say problems with my feet reduced my quality of life:"
# label(data$nsqlq_66_mv)="Missing value verified for: Overall, I would say problems with my feet reduced my quality of life:"
# label(data$nsqlq_67)="Overall, I would rate my quality of life as:"
# label(data$nsqlq_67_mv)="Missing value verified for: Overall, I would rate my quality of life as:"
# label(data$nsqolq_encounter_id)="Encounter ID"
# label(data$neuropathy_specific_quality_of_life_questionnaire_complete)="Complete?"
# label(data$jvn_collecdate)="Date of NIH Study"
# label(data$jvn_collecdate_mv)="Missing value verified for: Date of NIH Study"
# label(data$jvn_date)="Date of JVN Image"
# label(data$jvn_date_mv)="Missing value verified for: Date of JVN Image"
# label(data$jvn_npdr_od)="Right Eye - (OD)"
# label(data$jvn_npdr_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_npdr_os)="Left Eye - (OS)"
# label(data$jvn_npdr_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_npdr_notes_mv)="Missing value verified for: Notes"
# label(data$jvn_pdr_od)="Right Eye - (OD)"
# label(data$jvn_pdr_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_pdr_os)="Left Eye - (OS)"
# label(data$jvn_pdr_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_pdr_notes_mv)="Missing value verified for: Notes"
# label(data$jvn_me_od)="Right Eye - (OD)"
# label(data$jvn_me_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_me_os)="Left Eye - (OS)"
# label(data$jvn_me_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_me_notes_mv)="Missing value verified for: Notes"
# label(data$jvn_amd_od)="Right Eye - (OD)"
# label(data$jvn_amd_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_amd_os)="Left Eye - (OS)"
# label(data$jvn_amd_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_asterhy_od)="Right Eye - (OD)"
# label(data$jvn_asterhy_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_asterhy_os)="Left Eye - (OS)"
# label(data$jvn_asterhy_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_cat_od)="Right Eye - (OD)"
# label(data$jvn_cat_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_cat_os)="Left Eye - (OS)"
# label(data$jvn_cat_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_chorles_od)="Right Eye - (OD)"
# label(data$jvn_chorles_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_chorles_os)="Left Eye - (OS)"
# label(data$jvn_chorles_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_arcus_od)="Right Eye - (OD)"
# label(data$jvn_arcus_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_arcus_os)="Left Eye - (OS)"
# label(data$jvn_arcus_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_cratscar_od)="Right Eye - (OD)"
# label(data$jvn_cratscar_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_cratscar_os)="Left Eye - (OS)"
# label(data$jvn_cratscar_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_epiret_od)="Right Eye - (OD)"
# label(data$jvn_epiret_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_epiret_os)="Left Eye - (OS)"
# label(data$jvn_epiret_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_iol_od)="Right Eye - (OD)"
# label(data$jvn_iol_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_iol_os)="Left Eye - (OS)"
# label(data$jvn_iol_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_nevus_od)="Right Eye - (OD)"
# label(data$jvn_nevus_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_nevus_os)="Left Eye - (OS)"
# label(data$jvn_nevus_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_hardex_od)="Right Eye - (OD)"
# label(data$jvn_hardex_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_hardex_os)="Left Eye - (OS)"
# label(data$jvn_hardex_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_lidles_od)="Right Eye - (OD)"
# label(data$jvn_lidles_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_lidles_os)="Left Eye - (OS)"
# label(data$jvn_lidles_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_drusen_od)="Right Eye - (OD)"
# label(data$jvn_drusen_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_drusen_os)="Left Eye - (OS)"
# label(data$jvn_drusen_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_ondrus_od)="Right Eye - (OD)"
# label(data$jvn_ondrus_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_ondrus_os)="Left Eye - (OS)"
# label(data$jvn_ondrus_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_miscas_od)="Right Eye - (OD)"
# label(data$jvn_miscas_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_miscas_os)="Left Eye - (OS)"
# label(data$jvn_miscas_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_miscvr_od)="Right Eye - (OD)"
# label(data$jvn_miscvr_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_miscvr_os)="Left Eye - (OS)"
# label(data$jvn_miscvr_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_ppa_od)="Right Eye - (OD)"
# label(data$jvn_ppa_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_ppa_os)="Left Eye - (OS)"
# label(data$jvn_ppa_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_prfib_od)="Right Eye - (OD)"
# label(data$jvn_prfib_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_prfib_os)="Left Eye - (OS)"
# label(data$jvn_prfib_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_htnret_od)="Right Eye - (OD)"
# label(data$jvn_htnret_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_htnret_os)="Left Eye - (OS)"
# label(data$jvn_htnret_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvp_rpe_od)="Right Eye - (OD)"
# label(data$jvp_rpe_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_rpe_os)="Left Eye - (OS)"
# label(data$jvn_rpe_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_rrd_od)="Right Eye - (OD)"
# label(data$jvn_rrd_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_rrd_os)="Left Eye - (OS)"
# label(data$jvn_rrd_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_trd_od)="Right Eye - (OD)"
# label(data$jvn_trd_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_trd_os)="Left Eye - (OS)"
# label(data$jvn_trd_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_vithem_od)="Right Eye - (OD)"
# label(data$jvn_vithem_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_vithem_os)="Left Eye - (OS)"
# label(data$jvn_vithem_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_other_od)="Right Eye - (OD)"
# label(data$jvn_other_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_other_os)="Left Eye - (OS)"
# label(data$jvn_other_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_laser_od)="Right Eye - (OD)"
# label(data$jvn_laser_od_mv)="Missing value verified for: Right Eye - (OD)"
# label(data$jvn_laser_os)="Left Eye - (OS)"
# label(data$jvn_laser_os_mv)="Missing value verified for: Left Eye - (OS)"
# label(data$jvn_notes_mv)="Missing value verified for: Notes"
# label(data$jvn_deck_8_mv)="Missing value verified for: Deck"
# label(data$jvnrif_encounter_id)="Encounter ID"
# label(data$jvn_retinal_imaging_form_complete)="Complete?"
#Setting Factors(will create new variable for factors)
data$redcap_event_name.factor = factor(data$redcap_event_name,levels=c("interval_0_arm_1","interval_12_arm_1","interval_24_arm_1","interval_36_arm_1","interval_48_arm_1","interval_60_arm_1","interval_72_arm_1","interval_84_arm_1","interval_96_arm_1","interval_108_arm_1","interval_120_arm_1","interval_132_arm_1","interval_144_arm_1","interval_156_arm_1","interval_168_arm_1","interval_180_arm_1","interval_192_arm_1","interval_204_arm_1","interval_216_arm_1","interval_228_arm_1"))
data$redcap_repeat_instrument.factor = factor(data$redcap_repeat_instrument,levels=c(""))
data$lastname_mv.factor = factor(data$lastname_mv,levels=c("1","0"))
data$firstname_mv.factor = factor(data$firstname_mv,levels=c("1","0"))
data$stratum.factor = factor(data$stratum,levels=c("1","2","3"))
data$stratum_mv.factor = factor(data$stratum_mv,levels=c("1","0"))
data$randomnum_mv.factor = factor(data$randomnum_mv,levels=c("1","0"))
data$phx_number_mv.factor = factor(data$phx_number_mv,levels=c("1","0"))
data$dob_mv.factor = factor(data$dob_mv,levels=c("1","0"))
data$gender.factor = factor(data$gender,levels=c("1","2","3"))
data$gender_mv.factor = factor(data$gender_mv,levels=c("1","0"))
data$height_mv.factor = factor(data$height_mv,levels=c("1","0"))
data$allergies.factor = factor(data$allergies,levels=c("1","2"))
data$allergies_mv.factor = factor(data$allergies_mv,levels=c("1","0"))
data$allerglist_mv.factor = factor(data$allerglist_mv,levels=c("1","0"))
data$district_mv.factor = factor(data$district_mv,levels=c("1","0"))
data$phonenum_mv.factor = factor(data$phonenum_mv,levels=c("1","0"))
data$bldprsr.factor = factor(data$bldprsr,levels=c("1","2"))
data$bldprsr_mv.factor = factor(data$bldprsr_mv,levels=c("1","0"))
data$dtdiabonset_mv.factor = factor(data$dtdiabonset_mv,levels=c("1","0"))
data$typerecruit.factor = factor(data$typerecruit,levels=c("1","2","3","4","5"))
data$typerecruit_mv.factor = factor(data$typerecruit_mv,levels=c("1","0"))
data$typerecruit_other_mv.factor = factor(data$typerecruit_other_mv,levels=c("1","0"))
data$familyrelation.factor = factor(data$familyrelation,levels=c("1","2","3","4"))
data$familyrelation_mv.factor = factor(data$familyrelation_mv,levels=c("1","0"))
data$familyrelation_other_mv.factor = factor(data$familyrelation_other_mv,levels=c("1","0"))
data$nih_number_mv.factor = factor(data$nih_number_mv,levels=c("1","0"))
data$demoform_complete.factor = factor(data$demoform_complete,levels=c("0","1","2"))
data$sex_code.factor = factor(data$sex_code,levels=c("1","2"))
data$sex_code_mv.factor = factor(data$sex_code_mv,levels=c("1","0"))
data$start_date_mv.factor = factor(data$start_date_mv,levels=c("1","0"))
data$end_date_mv.factor = factor(data$end_date_mv,levels=c("1","0"))
data$study_status.factor = factor(data$study_status,levels=c("C","I"))
data$study_status_mv.factor = factor(data$study_status_mv,levels=c("1","0"))
data$clrncvisittype.factor = factor(data$clrncvisittype,levels=c("3","5"))
data$clrncvisittype_mv.factor = factor(data$clrncvisittype_mv,levels=c("1","0"))
data$clrnctestintvrl_mv.factor = factor(data$clrnctestintvrl_mv,levels=c("1","0"))
data$comments_mv.factor = factor(data$comments_mv,levels=c("1","0"))
data$study_code.factor = factor(data$study_code,levels=c("FTP","EIR","61","RET","EEP","81","82","5YR","DBD","PAB","EE1","EE2","EE3","EE4","EE5","EEA","WP5","EE6","EE7","DPP","LP3","T4SS","T4SW1","36"))
data$study_code_mv.factor = factor(data$study_code_mv,levels=c("1","0"))
data$encounter_id_mv.factor = factor(data$encounter_id_mv,levels=c("1","0"))
data$clearancevisit_complete.factor = factor(data$clearancevisit_complete,levels=c("0","1","2"))
data$brainmri_done.factor = factor(data$brainmri_done,levels=c("1","0"))
data$brainmri_done_mv.factor = factor(data$brainmri_done_mv,levels=c("1","0"))
data$brainmri_date_mv.factor = factor(data$brainmri_date_mv,levels=c("1","0"))
data$cognitivetest_done.factor = factor(data$cognitivetest_done,levels=c("1","0"))
data$cognitivetest_done_mv.factor = factor(data$cognitivetest_done_mv,levels=c("1","0"))
data$cognitivetest_date_mv.factor = factor(data$cognitivetest_date_mv,levels=c("1","0"))
data$kidneybiopsy_done.factor = factor(data$kidneybiopsy_done,levels=c("1","0"))
data$kidneybiopsy_done_mv.factor = factor(data$kidneybiopsy_done_mv,levels=c("1","0"))
data$kidneybiopsy_date_mv.factor = factor(data$kidneybiopsy_date_mv,levels=c("1","0"))
data$kidneybiopsy_site.factor = factor(data$kidneybiopsy_site,levels=c("1","2"))
data$kidneybiopsy_site_mv.factor = factor(data$kidneybiopsy_site_mv,levels=c("1","0"))
data$kidneylivermre_done.factor = factor(data$kidneylivermre_done,levels=c("1","0"))
data$kidneylivermre_done_mv.factor = factor(data$kidneylivermre_done_mv,levels=c("1","0"))
data$kidneylivermre_date_mv.factor = factor(data$kidneylivermre_date_mv,levels=c("1","0"))
data$kidneymri_done.factor = factor(data$kidneymri_done,levels=c("1","0"))
data$kidneymri_done_mv.factor = factor(data$kidneymri_done_mv,levels=c("1","0"))
data$kidneymri_date_mv.factor = factor(data$kidneymri_date_mv,levels=c("1","0"))
data$skinbiopsy_done.factor = factor(data$skinbiopsy_done,levels=c("1","0"))
data$skinbiopsy_done_mv.factor = factor(data$skinbiopsy_done_mv,levels=c("1","0"))
data$skinbiopsy_date_mv.factor = factor(data$skinbiopsy_date_mv,levels=c("1","0"))
data$skinbiopsy_site.factor = factor(data$skinbiopsy_site,levels=c("1","2"))
data$skinbiopsy_site_mv.factor = factor(data$skinbiopsy_site_mv,levels=c("1","0"))
data$stemcells_done.factor = factor(data$stemcells_done,levels=c("1","0"))
data$stemcells_done_mv.factor = factor(data$stemcells_done_mv,levels=c("1","0"))
data$stemcells_date_mv.factor = factor(data$stemcells_date_mv,levels=c("1","0"))
data$stool_done.factor = factor(data$stool_done,levels=c("1","0"))
data$stool_done_mv.factor = factor(data$stool_done_mv,levels=c("1","0"))
data$stool_date_mv.factor = factor(data$stool_date_mv,levels=c("1","0"))
data$neuropathy_done.factor = factor(data$neuropathy_done,levels=c("1","0"))
data$neuropathy_done_mv.factor = factor(data$neuropathy_done_mv,levels=c("1","0"))
data$neuropathy_date_mv.factor = factor(data$neuropathy_date_mv,levels=c("1","0"))
data$retinopathy_done.factor = factor(data$retinopathy_done,levels=c("1","0"))
data$retinopathy_done_mv.factor = factor(data$retinopathy_done_mv,levels=c("1","0"))
data$retinopathy_date_mv.factor = factor(data$retinopathy_date_mv,levels=c("1","0"))
data$sp_comments_mv.factor = factor(data$sp_comments_mv,levels=c("1","0"))
data$study_procedures_status_complete.factor = factor(data$study_procedures_status_complete,levels=c("0","1","2"))
data$f4bldprsr.factor = factor(data$f4bldprsr,levels=c("1","2"))
data$f4bldprsr_mv.factor = factor(data$f4bldprsr_mv,levels=c("1","0"))
data$f4frstmsrsys_mv.factor = factor(data$f4frstmsrsys_mv,levels=c("1","0"))
data$f4frstmsrdia_mv.factor = factor(data$f4frstmsrdia_mv,levels=c("1","0"))
data$f4scndmsrsys_mv.factor = factor(data$f4scndmsrsys_mv,levels=c("1","0"))
data$f4scndmsrdia_mv.factor = factor(data$f4scndmsrdia_mv,levels=c("1","0"))
data$f4heartrate_mv.factor = factor(data$f4heartrate_mv,levels=c("1","0"))
data$f4pregnant.factor = factor(data$f4pregnant,levels=c("1","2"))
data$f4pregnant_mv.factor = factor(data$f4pregnant_mv,levels=c("1","0"))
data$f4medsnow.factor = factor(data$f4medsnow,levels=c("1","2"))
data$f4medsnow_mv.factor = factor(data$f4medsnow_mv,levels=c("1","0"))
data$f4srmpotasm_mv.factor = factor(data$f4srmpotasm_mv,levels=c("1","0"))
data$f4srmglu_mv.factor = factor(data$f4srmglu_mv,levels=c("1","0"))
data$f4forminit_mv.factor = factor(data$f4forminit_mv,levels=c("1","0"))
data$f4bldprsrandmeds_complete.factor = factor(data$f4bldprsrandmeds_complete,levels=c("0","1","2"))
data$f4drugcode1_mv.factor = factor(data$f4drugcode1_mv,levels=c("1","0"))
data$f4dosage1_mv.factor = factor(data$f4dosage1_mv,levels=c("1","0"))
data$f4timesday1_mv.factor = factor(data$f4timesday1_mv,levels=c("1","0"))
data$f4dwm1.factor = factor(data$f4dwm1,levels=c("D","W","M"))
data$f4dwm1_mv.factor = factor(data$f4dwm1_mv,levels=c("1","0"))
data$f4drugcode2_mv.factor = factor(data$f4drugcode2_mv,levels=c("1","0"))
data$f4dosage2_mv.factor = factor(data$f4dosage2_mv,levels=c("1","0"))
data$f4timesday2_mv.factor = factor(data$f4timesday2_mv,levels=c("1","0"))
data$f4dwm2.factor = factor(data$f4dwm2,levels=c("D","W","M"))
data$f4dwm2_mv.factor = factor(data$f4dwm2_mv,levels=c("1","0"))
data$f4drugcode3_mv.factor = factor(data$f4drugcode3_mv,levels=c("1","0"))
data$f4dosage3_mv.factor = factor(data$f4dosage3_mv,levels=c("1","0"))
data$f4timesday3_mv.factor = factor(data$f4timesday3_mv,levels=c("1","0"))
data$f4dwm3.factor = factor(data$f4dwm3,levels=c("D","W","M"))
data$f4dwm3_mv.factor = factor(data$f4dwm3_mv,levels=c("1","0"))
data$f4drugcode4_mv.factor = factor(data$f4drugcode4_mv,levels=c("1","0"))
data$f4dosage4_mv.factor = factor(data$f4dosage4_mv,levels=c("1","0"))
data$f4timesday4_mv.factor = factor(data$f4timesday4_mv,levels=c("1","0"))
data$f4dwm4.factor = factor(data$f4dwm4,levels=c("D","W","M"))
data$f4dwm4_mv.factor = factor(data$f4dwm4_mv,levels=c("1","0"))
data$f4drugcode5_mv.factor = factor(data$f4drugcode5_mv,levels=c("1","0"))
data$f4dosage5_mv.factor = factor(data$f4dosage5_mv,levels=c("1","0"))
data$f4timesday5_mv.factor = factor(data$f4timesday5_mv,levels=c("1","0"))
data$f4dwm5.factor = factor(data$f4dwm5,levels=c("D","W","M"))
data$f4dwm5_mv.factor = factor(data$f4dwm5_mv,levels=c("1","0"))
data$f4drugcode6_mv.factor = factor(data$f4drugcode6_mv,levels=c("1","0"))
data$f4dosage6_mv.factor = factor(data$f4dosage6_mv,levels=c("1","0"))
data$f4timesday6_mv.factor = factor(data$f4timesday6_mv,levels=c("1","0"))
data$f4dwm6.factor = factor(data$f4dwm6,levels=c("D","W","M"))
data$f4dwm6_mv.factor = factor(data$f4dwm6_mv,levels=c("1","0"))
data$f4drugcode7_mv.factor = factor(data$f4drugcode7_mv,levels=c("1","0"))
data$f4dosage7_mv.factor = factor(data$f4dosage7_mv,levels=c("1","0"))
data$f4timesday7_mv.factor = factor(data$f4timesday7_mv,levels=c("1","0"))
data$f4dwm7.factor = factor(data$f4dwm7,levels=c("D","W","M"))
data$f4dwm7_mv.factor = factor(data$f4dwm7_mv,levels=c("1","0"))
data$f4drugcode8_mv.factor = factor(data$f4drugcode8_mv,levels=c("1","0"))
data$f4dosage8_mv.factor = factor(data$f4dosage8_mv,levels=c("1","0"))
data$f4timesday8_mv.factor = factor(data$f4timesday8_mv,levels=c("1","0"))
data$f4dwm8.factor = factor(data$f4dwm8,levels=c("D","W","M"))
data$f4dwm8_mv.factor = factor(data$f4dwm8_mv,levels=c("1","0"))
data$f4drugcode9_mv.factor = factor(data$f4drugcode9_mv,levels=c("1","0"))
data$f4dosage9_mv.factor = factor(data$f4dosage9_mv,levels=c("1","0"))
data$f4timesday9_mv.factor = factor(data$f4timesday9_mv,levels=c("1","0"))
data$f4dwm9.factor = factor(data$f4dwm9,levels=c("D","W","M"))
data$f4dwm9_mv.factor = factor(data$f4dwm9_mv,levels=c("1","0"))
data$f4drugcode10_mv.factor = factor(data$f4drugcode10_mv,levels=c("1","0"))
data$f4dosage10_mv.factor = factor(data$f4dosage10_mv,levels=c("1","0"))
data$f4timesday10_mv.factor = factor(data$f4timesday10_mv,levels=c("1","0"))
data$f4dwm10.factor = factor(data$f4dwm10,levels=c("D","W","M"))
data$f4dwm10_mv.factor = factor(data$f4dwm10_mv,levels=c("1","0"))
data$f4drugcode11_mv.factor = factor(data$f4drugcode11_mv,levels=c("1","0"))
data$f4dosage11_mv.factor = factor(data$f4dosage11_mv,levels=c("1","0"))
data$f4timesday11_mv.factor = factor(data$f4timesday11_mv,levels=c("1","0"))
data$f4dwm11.factor = factor(data$f4dwm11,levels=c("D","W","M"))
data$f4dwm11_mv.factor = factor(data$f4dwm11_mv,levels=c("1","0"))
data$f4drugcode12_mv.factor = factor(data$f4drugcode12_mv,levels=c("1","0"))
data$f4dosage12_mv.factor = factor(data$f4dosage12_mv,levels=c("1","0"))
data$f4timesday12_mv.factor = factor(data$f4timesday12_mv,levels=c("1","0"))
data$f4dwm12.factor = factor(data$f4dwm12,levels=c("D","W","M"))
data$f4dwm12_mv.factor = factor(data$f4dwm12_mv,levels=c("1","0"))
data$f4drugcode13_mv.factor = factor(data$f4drugcode13_mv,levels=c("1","0"))
data$f4dosage13_mv.factor = factor(data$f4dosage13_mv,levels=c("1","0"))
data$f4timesday13_mv.factor = factor(data$f4timesday13_mv,levels=c("1","0"))
data$f4dwm13.factor = factor(data$f4dwm13,levels=c("D","W","M"))
data$f4dwm13_mv.factor = factor(data$f4dwm13_mv,levels=c("1","0"))
data$f4drugcode14_mv.factor = factor(data$f4drugcode14_mv,levels=c("1","0"))
data$f4dosage14_mv.factor = factor(data$f4dosage14_mv,levels=c("1","0"))
data$f4timesday14_mv.factor = factor(data$f4timesday14_mv,levels=c("1","0"))
data$f4dwm14.factor = factor(data$f4dwm14,levels=c("D","W","M"))
data$f4dwm14_mv.factor = factor(data$f4dwm14_mv,levels=c("1","0"))
data$f4drugcode15_mv.factor = factor(data$f4drugcode15_mv,levels=c("1","0"))
data$f4dosage15_mv.factor = factor(data$f4dosage15_mv,levels=c("1","0"))
data$f4timesday15_mv.factor = factor(data$f4timesday15_mv,levels=c("1","0"))
data$f4dwm15.factor = factor(data$f4dwm15,levels=c("D","W","M"))
data$f4dwm15_mv.factor = factor(data$f4dwm15_mv,levels=c("1","0"))
data$f4drugcode16_mv.factor = factor(data$f4drugcode16_mv,levels=c("1","0"))
data$f4dosage16_mv.factor = factor(data$f4dosage16_mv,levels=c("1","0"))
data$f4timesday16_mv.factor = factor(data$f4timesday16_mv,levels=c("1","0"))
data$f4dwm16.factor = factor(data$f4dwm16,levels=c("D","W","M"))
data$f4dwm16_mv.factor = factor(data$f4dwm16_mv,levels=c("1","0"))
data$f4drugcode17_mv.factor = factor(data$f4drugcode17_mv,levels=c("1","0"))
data$f4dosage17_mv.factor = factor(data$f4dosage17_mv,levels=c("1","0"))
data$f4timesday17_mv.factor = factor(data$f4timesday17_mv,levels=c("1","0"))
data$f4dwm17.factor = factor(data$f4dwm17,levels=c("D","W","M"))
data$f4dwm17_mv.factor = factor(data$f4dwm17_mv,levels=c("1","0"))
data$f4drugcode18_mv.factor = factor(data$f4drugcode18_mv,levels=c("1","0"))
data$f4dosage18_mv.factor = factor(data$f4dosage18_mv,levels=c("1","0"))
data$f4timesday18_mv.factor = factor(data$f4timesday18_mv,levels=c("1","0"))
data$f4dwm18.factor = factor(data$f4dwm18,levels=c("D","W","M"))
data$f4dwm18_mv.factor = factor(data$f4dwm18_mv,levels=c("1","0"))
data$f4drugcode19_mv.factor = factor(data$f4drugcode19_mv,levels=c("1","0"))
data$f4dosage19_mv.factor = factor(data$f4dosage19_mv,levels=c("1","0"))
data$f4timesday19_mv.factor = factor(data$f4timesday19_mv,levels=c("1","0"))
data$f4dwm19.factor = factor(data$f4dwm19,levels=c("D","W","M"))
data$f4dwm19_mv.factor = factor(data$f4dwm19_mv,levels=c("1","0"))
data$f4drugcode20_mv.factor = factor(data$f4drugcode20_mv,levels=c("1","0"))
data$f4dosage20_mv.factor = factor(data$f4dosage20_mv,levels=c("1","0"))
data$f4timesday20_mv.factor = factor(data$f4timesday20_mv,levels=c("1","0"))
data$f4dwm20.factor = factor(data$f4dwm20,levels=c("D","W","M"))
data$f4dwm20_mv.factor = factor(data$f4dwm20_mv,levels=c("1","0"))
data$f4drugcode21_mv.factor = factor(data$f4drugcode21_mv,levels=c("1","0"))
data$f4dosage21_mv.factor = factor(data$f4dosage21_mv,levels=c("1","0"))
data$f4timesday21_mv.factor = factor(data$f4timesday21_mv,levels=c("1","0"))
data$f4dwm21.factor = factor(data$f4dwm21,levels=c("D","W","M"))
data$f4dwm21_mv.factor = factor(data$f4dwm21_mv,levels=c("1","0"))
data$f4drugcode22_mv.factor = factor(data$f4drugcode22_mv,levels=c("1","0"))
data$f4dosage22_mv.factor = factor(data$f4dosage22_mv,levels=c("1","0"))
data$f4timesday22_mv.factor = factor(data$f4timesday22_mv,levels=c("1","0"))
data$f4dwm22.factor = factor(data$f4dwm22,levels=c("D","W","M"))
data$f4dwm22_mv.factor = factor(data$f4dwm22_mv,levels=c("1","0"))
data$f4drugcode23_mv.factor = factor(data$f4drugcode23_mv,levels=c("1","0"))
data$f4dosage23_mv.factor = factor(data$f4dosage23_mv,levels=c("1","0"))
data$f4timesday23_mv.factor = factor(data$f4timesday23_mv,levels=c("1","0"))
data$f4dwm23.factor = factor(data$f4dwm23,levels=c("D","W","M"))
data$f4dwm23_mv.factor = factor(data$f4dwm23_mv,levels=c("1","0"))
data$f4drugcode24_mv.factor = factor(data$f4drugcode24_mv,levels=c("1","0"))
data$f4dosage24_mv.factor = factor(data$f4dosage24_mv,levels=c("1","0"))
data$f4timesday24_mv.factor = factor(data$f4timesday24_mv,levels=c("1","0"))
data$f4dwm24.factor = factor(data$f4dwm24,levels=c("D","W","M"))
data$f4dwm24_mv.factor = factor(data$f4dwm24_mv,levels=c("1","0"))
data$f4drugcode25_mv.factor = factor(data$f4drugcode25_mv,levels=c("1","0"))
data$f4dosage25_mv.factor = factor(data$f4dosage25_mv,levels=c("1","0"))
data$f4timesday25_mv.factor = factor(data$f4timesday25_mv,levels=c("1","0"))
data$f4dwm25.factor = factor(data$f4dwm25,levels=c("D","W","M"))
data$f4dwm25_mv.factor = factor(data$f4dwm25_mv,levels=c("1","0"))
data$f4drugcode26_mv.factor = factor(data$f4drugcode26_mv,levels=c("1","0"))
data$f4dosage26_mv.factor = factor(data$f4dosage26_mv,levels=c("1","0"))
data$f4timesday26_mv.factor = factor(data$f4timesday26_mv,levels=c("1","0"))
data$f4dwm26.factor = factor(data$f4dwm26,levels=c("D","W","M"))
data$f4dwm26_mv.factor = factor(data$f4dwm26_mv,levels=c("1","0"))
data$f4drugcode27_mv.factor = factor(data$f4drugcode27_mv,levels=c("1","0"))
data$f4dosage27_mv.factor = factor(data$f4dosage27_mv,levels=c("1","0"))
data$f4timesday27_mv.factor = factor(data$f4timesday27_mv,levels=c("1","0"))
data$f4dwm27.factor = factor(data$f4dwm27,levels=c("D","W","M"))
data$f4dwm27_mv.factor = factor(data$f4dwm27_mv,levels=c("1","0"))
data$f4drugcode28_mv.factor = factor(data$f4drugcode28_mv,levels=c("1","0"))
data$f4dosage28_mv.factor = factor(data$f4dosage28_mv,levels=c("1","0"))
data$f4timesday28_mv.factor = factor(data$f4timesday28_mv,levels=c("1","0"))
data$f4dwm28.factor = factor(data$f4dwm28,levels=c("D","W","M"))
data$f4dwm28_mv.factor = factor(data$f4dwm28_mv,levels=c("1","0"))
data$f4drugcode29_mv.factor = factor(data$f4drugcode29_mv,levels=c("1","0"))
data$f4dosage29_mv.factor = factor(data$f4dosage29_mv,levels=c("1","0"))
data$f4timesday29_mv.factor = factor(data$f4timesday29_mv,levels=c("1","0"))
data$f4dwm29.factor = factor(data$f4dwm29,levels=c("D","W","M"))
data$f4dwm29_mv.factor = factor(data$f4dwm29_mv,levels=c("1","0"))
data$f4drugcode30_mv.factor = factor(data$f4drugcode30_mv,levels=c("1","0"))
data$f4dosage30_mv.factor = factor(data$f4dosage30_mv,levels=c("1","0"))
data$f4timesday30_mv.factor = factor(data$f4timesday30_mv,levels=c("1","0"))
data$f4dwm30.factor = factor(data$f4dwm30,levels=c("D","W","M"))
data$f4dwm30_mv.factor = factor(data$f4dwm30_mv,levels=c("1","0"))
data$f4meds_complete.factor = factor(data$f4meds_complete,levels=c("0","1","2"))
data$f5rcntillns.factor = factor(data$f5rcntillns,levels=c("1","2"))
data$f5rcntillns_mv.factor = factor(data$f5rcntillns_mv,levels=c("1","0"))
data$f5outpat.factor = factor(data$f5outpat,levels=c("1","2"))
data$f5outpat_mv.factor = factor(data$f5outpat_mv,levels=c("1","0"))
data$f5outpatspec_mv.factor = factor(data$f5outpatspec_mv,levels=c("1","0"))
data$f5hospwosrg.factor = factor(data$f5hospwosrg,levels=c("1","2"))
data$f5hospwosrg_mv.factor = factor(data$f5hospwosrg_mv,levels=c("1","0"))
data$f5hospwoexp_mv.factor = factor(data$f5hospwoexp_mv,levels=c("1","0"))
data$f5hospwsrg.factor = factor(data$f5hospwsrg,levels=c("1","2"))
data$f5hospwsrg_mv.factor = factor(data$f5hospwsrg_mv,levels=c("1","0"))
data$f5hospwexp_mv.factor = factor(data$f5hospwexp_mv,levels=c("1","0"))
data$f5medprbsdoc.factor = factor(data$f5medprbsdoc,levels=c("1","2"))
data$f5medprbsdoc_mv.factor = factor(data$f5medprbsdoc_mv,levels=c("1","0"))
data$f5crnryartds.factor = factor(data$f5crnryartds,levels=c("1","2"))
data$f5crnryartds_mv.factor = factor(data$f5crnryartds_mv,levels=c("1","0"))
data$f5cancer.factor = factor(data$f5cancer,levels=c("1","2"))
data$f5cancer_mv.factor = factor(data$f5cancer_mv,levels=c("1","0"))
data$f5cancerexp_mv.factor = factor(data$f5cancerexp_mv,levels=c("1","0"))
data$f5crbrlvasds.factor = factor(data$f5crbrlvasds,levels=c("1","2"))
data$f5crbrlvasds_mv.factor = factor(data$f5crbrlvasds_mv,levels=c("1","0"))
data$f5peripvasds.factor = factor(data$f5peripvasds,levels=c("1","2"))
data$f5peripvasds_mv.factor = factor(data$f5peripvasds_mv,levels=c("1","0"))
data$f5hypertensn.factor = factor(data$f5hypertensn,levels=c("1","2"))
data$f5hypertensn_mv.factor = factor(data$f5hypertensn_mv,levels=c("1","0"))
data$f5seizures.factor = factor(data$f5seizures,levels=c("1","2"))
data$f5seizures_mv.factor = factor(data$f5seizures_mv,levels=c("1","0"))
data$f5gnitrnryds.factor = factor(data$f5gnitrnryds,levels=c("1","2"))
data$f5gnitrnryds_mv.factor = factor(data$f5gnitrnryds_mv,levels=c("1","0"))
data$f5gnitrnyexp_mv.factor = factor(data$f5gnitrnyexp_mv,levels=c("1","0"))
data$f5lungdsz.factor = factor(data$f5lungdsz,levels=c("1","2"))
data$f5lungdsz_mv.factor = factor(data$f5lungdsz_mv,levels=c("1","0"))
data$f5majsurg.factor = factor(data$f5majsurg,levels=c("1","2"))
data$f5majsurg_mv.factor = factor(data$f5majsurg_mv,levels=c("1","0"))
data$f5majsurgexp_mv.factor = factor(data$f5majsurgexp_mv,levels=c("1","0"))
data$f5othrmeddia.factor = factor(data$f5othrmeddia,levels=c("1","2"))
data$f5othrmeddia_mv.factor = factor(data$f5othrmeddia_mv,levels=c("1","0"))
data$f5othrmedexp_mv.factor = factor(data$f5othrmedexp_mv,levels=c("1","0"))
data$f5lungs.factor = factor(data$f5lungs,levels=c("1","2"))
data$f5lungs_mv.factor = factor(data$f5lungs_mv,levels=c("1","0"))
data$f5lungsexp_mv.factor = factor(data$f5lungsexp_mv,levels=c("1","0"))
data$f5heart.factor = factor(data$f5heart,levels=c("1","2"))
data$f5heart_mv.factor = factor(data$f5heart_mv,levels=c("1","0"))
data$f5heartexp_mv.factor = factor(data$f5heartexp_mv,levels=c("1","0"))
data$f5skin.factor = factor(data$f5skin,levels=c("1","2"))
data$f5skin_mv.factor = factor(data$f5skin_mv,levels=c("1","0"))
data$f5skinexp_mv.factor = factor(data$f5skinexp_mv,levels=c("1","0"))
data$f5edema.factor = factor(data$f5edema,levels=c("1","2"))
data$f5edema_mv.factor = factor(data$f5edema_mv,levels=c("1","0"))
data$f5edemaexp_mv.factor = factor(data$f5edemaexp_mv,levels=c("1","0"))
data$f5histandphysexm_complete.factor = factor(data$f5histandphysexm_complete,levels=c("0","1","2"))
data$f6ivsites_mv.factor = factor(data$f6ivsites_mv,levels=c("1","0"))
data$f6weight_mv.factor = factor(data$f6weight_mv,levels=c("1","0"))
data$f6u0strttime_mv.factor = factor(data$f6u0strttime_mv,levels=c("1","0"))
data$f6u0volume_mv.factor = factor(data$f6u0volume_mv,levels=c("1","0"))
data$f6bldprsrarm.factor = factor(data$f6bldprsrarm,levels=c("1","2"))
data$f6bldprsrarm_mv.factor = factor(data$f6bldprsrarm_mv,levels=c("1","0"))
data$f6frstmsrsty_mv.factor = factor(data$f6frstmsrsty_mv,levels=c("1","0"))
data$f6frstmsrdia_mv.factor = factor(data$f6frstmsrdia_mv,levels=c("1","0"))
data$f6scndmsrsys_mv.factor = factor(data$f6scndmsrsys_mv,levels=c("1","0"))
data$f6scndmsrdia_mv.factor = factor(data$f6scndmsrdia_mv,levels=c("1","0"))
data$f6heartrate_mv.factor = factor(data$f6heartrate_mv,levels=c("1","0"))
data$f6u0h20_mv.factor = factor(data$f6u0h20_mv,levels=c("1","0"))
data$f6fnleqlurtm_mv.factor = factor(data$f6fnleqlurtm_mv,levels=c("1","0"))
data$f6p1time_mv.factor = factor(data$f6p1time_mv,levels=c("1","0"))
data$f6urnvol_mv.factor = factor(data$f6urnvol_mv,levels=c("1","0"))
data$f6h20vol_mv.factor = factor(data$f6h20vol_mv,levels=c("1","0"))
data$f6u1time_mv.factor = factor(data$f6u1time_mv,levels=c("1","0"))
data$f6p2time_mv.factor = factor(data$f6p2time_mv,levels=c("1","0"))
data$f6u1vol_mv.factor = factor(data$f6u1vol_mv,levels=c("1","0"))
data$f6u1h20vol_mv.factor = factor(data$f6u1h20vol_mv,levels=c("1","0"))
data$f6u2time_mv.factor = factor(data$f6u2time_mv,levels=c("1","0"))
data$f6p3time_mv.factor = factor(data$f6p3time_mv,levels=c("1","0"))
data$f6u2vol_mv.factor = factor(data$f6u2vol_mv,levels=c("1","0"))
data$f6u2h20vol_mv.factor = factor(data$f6u2h20vol_mv,levels=c("1","0"))
data$f6u3time_mv.factor = factor(data$f6u3time_mv,levels=c("1","0"))
data$f6p4time_mv.factor = factor(data$f6p4time_mv,levels=c("1","0"))
data$f6u3vol_mv.factor = factor(data$f6u3vol_mv,levels=c("1","0"))
data$f6u3h20vol_mv.factor = factor(data$f6u3h20vol_mv,levels=c("1","0"))
data$f6u4time_mv.factor = factor(data$f6u4time_mv,levels=c("1","0"))
data$f6p5time_mv.factor = factor(data$f6p5time_mv,levels=c("1","0"))
data$f6u4vol_mv.factor = factor(data$f6u4vol_mv,levels=c("1","0"))
data$f6u4h20vol_mv.factor = factor(data$f6u4h20vol_mv,levels=c("1","0"))
data$f6u5time_mv.factor = factor(data$f6u5time_mv,levels=c("1","0"))
data$f6p6time_mv.factor = factor(data$f6p6time_mv,levels=c("1","0"))
data$f6u5vol_mv.factor = factor(data$f6u5vol_mv,levels=c("1","0"))
data$f6u5h20vol_mv.factor = factor(data$f6u5h20vol_mv,levels=c("1","0"))
data$f6urndate_mv.factor = factor(data$f6urndate_mv,levels=c("1","0"))
data$f6glucose.factor = factor(data$f6glucose,levels=c("0","T","1","2","3","4"))
data$f6glucose_mv.factor = factor(data$f6glucose_mv,levels=c("1","0"))
data$f6bilirubin.factor = factor(data$f6bilirubin,levels=c("0","1","2","3"))
data$f6bilirubin_mv.factor = factor(data$f6bilirubin_mv,levels=c("1","0"))
data$f6ketones.factor = factor(data$f6ketones,levels=c("0","T","1","2","3","4"))
data$f6ketones_mv.factor = factor(data$f6ketones_mv,levels=c("1","0"))
data$f6specgrvty_mv.factor = factor(data$f6specgrvty_mv,levels=c("1","0"))
data$f6bloodocult.factor = factor(data$f6bloodocult,levels=c("1","2"))
data$f6bloodocult_mv.factor = factor(data$f6bloodocult_mv,levels=c("1","0"))
data$f6ph_mv.factor = factor(data$f6ph_mv,levels=c("1","0"))
data$f6protein.factor = factor(data$f6protein,levels=c("0","T","1","2","3","4"))
data$f6protein_mv.factor = factor(data$f6protein_mv,levels=c("1","0"))
data$f6casts_mv.factor = factor(data$f6casts_mv,levels=c("1","0"))
data$f6rbccast.factor = factor(data$f6rbccast,levels=c("1","2"))
data$f6rbccast_mv.factor = factor(data$f6rbccast_mv,levels=c("1","0"))
data$f6wbccast.factor = factor(data$f6wbccast,levels=c("1","2"))
data$f6wbccast_mv.factor = factor(data$f6wbccast_mv,levels=c("1","0"))
data$f6hyalinecast.factor = factor(data$f6hyalinecast,levels=c("1","2"))
data$f6hyalinecast_mv.factor = factor(data$f6hyalinecast_mv,levels=c("1","0"))
data$f6granularcast.factor = factor(data$f6granularcast,levels=c("1","2"))
data$f6granularcast_mv.factor = factor(data$f6granularcast_mv,levels=c("1","0"))
data$f6othercast.factor = factor(data$f6othercast,levels=c("1","2"))
data$f6othercast_mv.factor = factor(data$f6othercast_mv,levels=c("1","0"))
data$f6othercast_exp_mv.factor = factor(data$f6othercast_exp_mv,levels=c("1","0"))
data$f6wbc.factor = factor(data$f6wbc,levels=c("0","1","2"))
data$f6wbc_mv.factor = factor(data$f6wbc_mv,levels=c("1","0"))
data$f6rbc.factor = factor(data$f6rbc,levels=c("0","1","2"))
data$f6rbc_mv.factor = factor(data$f6rbc_mv,levels=c("1","0"))
data$f6epithcells.factor = factor(data$f6epithcells,levels=c("0","1","2"))
data$f6epithcells_mv.factor = factor(data$f6epithcells_mv,levels=c("1","0"))
data$f6bacteria.factor = factor(data$f6bacteria,levels=c("0","1","2"))
data$f6bacteria_mv.factor = factor(data$f6bacteria_mv,levels=c("1","0"))
data$f6renalclrnc_complete.factor = factor(data$f6renalclrnc_complete,levels=c("0","1","2"))
data$f8visitdate_mv.factor = factor(data$f8visitdate_mv,levels=c("1","0"))
data$f8analdate_mv.factor = factor(data$f8analdate_mv,levels=c("1","0"))
data$f8alti_mv.factor = factor(data$f8alti_mv,levels=c("1","0"))
data$f8glucose_mv.factor = factor(data$f8glucose_mv,levels=c("1","0"))
data$f8bun_mv.factor = factor(data$f8bun_mv,levels=c("1","0"))
data$f8creatinine_mv.factor = factor(data$f8creatinine_mv,levels=c("1","0"))
data$f8uricacid_mv.factor = factor(data$f8uricacid_mv,levels=c("1","0"))
data$f8sodium_mv.factor = factor(data$f8sodium_mv,levels=c("1","0"))
data$f8potassium_mv.factor = factor(data$f8potassium_mv,levels=c("1","0"))
data$f8chloride_mv.factor = factor(data$f8chloride_mv,levels=c("1","0"))
data$f8calcium_mv.factor = factor(data$f8calcium_mv,levels=c("1","0"))
data$f8totlprotn_mv.factor = factor(data$f8totlprotn_mv,levels=c("1","0"))
data$f8albumin_mv.factor = factor(data$f8albumin_mv,levels=c("1","0"))
data$f8cholestero_mv.factor = factor(data$f8cholestero_mv,levels=c("1","0"))
data$f8triglyceri_mv.factor = factor(data$f8triglyceri_mv,levels=c("1","0"))
data$f8hdl_mv.factor = factor(data$f8hdl_mv,levels=c("1","0"))
data$f8ldl_mv.factor = factor(data$f8ldl_mv,levels=c("1","0"))
data$f8totlbilrbn_mv.factor = factor(data$f8totlbilrbn_mv,levels=c("1","0"))
data$f8dirctbilrbn_mv.factor = factor(data$f8dirctbilrbn_mv,levels=c("1","0"))
data$f8alkp04_mv.factor = factor(data$f8alkp04_mv,levels=c("1","0"))
data$f8ast_mv.factor = factor(data$f8ast_mv,levels=c("1","0"))
data$f8smacchempanl_complete.factor = factor(data$f8smacchempanl_complete,levels=c("0","1","2"))
data$f9dateanal_mv.factor = factor(data$f9dateanal_mv,levels=c("1","0"))
data$f9analyzed_date_mv.factor = factor(data$f9analyzed_date_mv,levels=c("1","0"))
data$f9wbc_mv.factor = factor(data$f9wbc_mv,levels=c("1","0"))
data$f9rbc_mv.factor = factor(data$f9rbc_mv,levels=c("1","0"))
data$f9hemoglob_mv.factor = factor(data$f9hemoglob_mv,levels=c("1","0"))
data$f9hematocrit_mv.factor = factor(data$f9hematocrit_mv,levels=c("1","0"))
data$f9mcv_mv.factor = factor(data$f9mcv_mv,levels=c("1","0"))
data$f9mch_mv.factor = factor(data$f9mch_mv,levels=c("1","0"))
data$f9mchc_mv.factor = factor(data$f9mchc_mv,levels=c("1","0"))
data$f9platelets_mv.factor = factor(data$f9platelets_mv,levels=c("1","0"))
data$f9lymphinst_mv.factor = factor(data$f9lymphinst_mv,levels=c("1","0"))
data$f9monoinst_mv.factor = factor(data$f9monoinst_mv,levels=c("1","0"))
data$f9neutinst_mv.factor = factor(data$f9neutinst_mv,levels=c("1","0"))
data$f9eosinst_mv.factor = factor(data$f9eosinst_mv,levels=c("1","0"))
data$f9basoinst_mv.factor = factor(data$f9basoinst_mv,levels=c("1","0"))
data$f9gran_mv.factor = factor(data$f9gran_mv,levels=c("1","0"))
data$f9pt_mv.factor = factor(data$f9pt_mv,levels=c("1","0"))
data$f9inr_mv.factor = factor(data$f9inr_mv,levels=c("1","0"))
data$f9ptt_mv.factor = factor(data$f9ptt_mv,levels=c("1","0"))
data$f9cbc_complete.factor = factor(data$f9cbc_complete,levels=c("0","1","2"))
data$f16qcnumber_mv.factor = factor(data$f16qcnumber_mv,levels=c("1","0"))
data$f16analdate_mv.factor = factor(data$f16analdate_mv,levels=c("1","0"))
data$f16glucose_mv.factor = factor(data$f16glucose_mv,levels=c("1","0"))
data$f16bun_mv.factor = factor(data$f16bun_mv,levels=c("1","0"))
data$f16creatinine_mv.factor = factor(data$f16creatinine_mv,levels=c("1","0"))
data$f16sodium_mv.factor = factor(data$f16sodium_mv,levels=c("1","0"))
data$f16potassium_mv.factor = factor(data$f16potassium_mv,levels=c("1","0"))
data$f16chloride_mv.factor = factor(data$f16chloride_mv,levels=c("1","0"))
data$f16calcium_mv.factor = factor(data$f16calcium_mv,levels=c("1","0"))
data$f16phosphorus_mv.factor = factor(data$f16phosphorus_mv,levels=c("1","0"))
data$f16totlprotn_mv.factor = factor(data$f16totlprotn_mv,levels=c("1","0"))
data$f16albumin_mv.factor = factor(data$f16albumin_mv,levels=c("1","0"))
data$f16cholestero_mv.factor = factor(data$f16cholestero_mv,levels=c("1","0"))
data$f16triglyceri_mv.factor = factor(data$f16triglyceri_mv,levels=c("1","0"))
data$f16totlbilrbn_mv.factor = factor(data$f16totlbilrbn_mv,levels=c("1","0"))
data$f16alkp04_mv.factor = factor(data$f16alkp04_mv,levels=c("1","0"))
data$f16ast_mv.factor = factor(data$f16ast_mv,levels=c("1","0"))
data$f16smacqc_complete.factor = factor(data$f16smacqc_complete,levels=c("0","1","2"))
data$f17qcnumber_mv.factor = factor(data$f17qcnumber_mv,levels=c("1","0"))
data$f17analdate_mv.factor = factor(data$f17analdate_mv,levels=c("1","0"))
data$f17wbc_mv.factor = factor(data$f17wbc_mv,levels=c("1","0"))
data$f17rbc_mv.factor = factor(data$f17rbc_mv,levels=c("1","0"))
data$f17hemoglob_mv.factor = factor(data$f17hemoglob_mv,levels=c("1","0"))
data$f17hematocrit_mv.factor = factor(data$f17hematocrit_mv,levels=c("1","0"))
data$f17mcv_mv.factor = factor(data$f17mcv_mv,levels=c("1","0"))
data$f17mch_mv.factor = factor(data$f17mch_mv,levels=c("1","0"))
data$f17mchc_mv.factor = factor(data$f17mchc_mv,levels=c("1","0"))
data$f17platelets_mv.factor = factor(data$f17platelets_mv,levels=c("1","0"))
data$f17cbcqc_complete.factor = factor(data$f17cbcqc_complete,levels=c("0","1","2"))
data$f13dtlastvst_mv.factor = factor(data$f13dtlastvst_mv,levels=c("1","0"))
data$f13typevist.factor = factor(data$f13typevist,levels=c("3","5"))
data$f13typevist_mv.factor = factor(data$f13typevist_mv,levels=c("1","0"))
data$f13tstintrvl_mv.factor = factor(data$f13tstintrvl_mv,levels=c("1","0"))
data$f13endstgrnl.factor = factor(data$f13endstgrnl,levels=c("1","2"))
data$f13endstgrnl_mv.factor = factor(data$f13endstgrnl_mv,levels=c("1","0"))
data$f13esrddate_mv.factor = factor(data$f13esrddate_mv,levels=c("1","0"))
data$f13nondiabkid.factor = factor(data$f13nondiabkid,levels=c("1","2"))
data$f13nondiabkid_mv.factor = factor(data$f13nondiabkid_mv,levels=c("1","0"))
data$f13ndkddate_mv.factor = factor(data$f13ndkddate_mv,levels=c("1","0"))
data$f13imprdbldr.factor = factor(data$f13imprdbldr,levels=c("1","2"))
data$f13imprdbldr_mv.factor = factor(data$f13imprdbldr_mv,levels=c("1","0"))
data$f13ibfdate_mv.factor = factor(data$f13ibfdate_mv,levels=c("1","0"))
data$f13cngstvhrt.factor = factor(data$f13cngstvhrt,levels=c("1","2"))
data$f13cngstvhrt_mv.factor = factor(data$f13cngstvhrt_mv,levels=c("1","0"))
data$f13chfdate_mv.factor = factor(data$f13chfdate_mv,levels=c("1","0"))
data$f13ascites.factor = factor(data$f13ascites,levels=c("1","2"))
data$f13ascites_mv.factor = factor(data$f13ascites_mv,levels=c("1","0"))
data$f13ascitesdt_mv.factor = factor(data$f13ascitesdt_mv,levels=c("1","0"))
data$f13aceinhib.factor = factor(data$f13aceinhib,levels=c("1","2"))
data$f13aceinhib_mv.factor = factor(data$f13aceinhib_mv,levels=c("1","0"))
data$f13typeofreaction.factor = factor(data$f13typeofreaction,levels=c("1","2","3","4","5","6","7"))
data$f13typeofreaction_mv.factor = factor(data$f13typeofreaction_mv,levels=c("1","0"))
data$f13otherdiag.factor = factor(data$f13otherdiag,levels=c("1","2"))
data$f13otherdiag_mv.factor = factor(data$f13otherdiag_mv,levels=c("1","0"))
data$f13reactionspecify_mv.factor = factor(data$f13reactionspecify_mv,levels=c("1","0"))
data$f13datereactn_mv.factor = factor(data$f13datereactn_mv,levels=c("1","0"))
data$f13prfrmwdclrnc.factor = factor(data$f13prfrmwdclrnc,levels=c("1","2"))
data$f13prfrmwdclrnc_mv.factor = factor(data$f13prfrmwdclrnc_mv,levels=c("1","0"))
data$f13prfrmwdclrncdate_mv.factor = factor(data$f13prfrmwdclrncdate_mv,levels=c("1","0"))
data$f13dtstoppntdeclrd_mv.factor = factor(data$f13dtstoppntdeclrd_mv,levels=c("1","0"))
data$f13dtfrmclmpltd_mv.factor = factor(data$f13dtfrmclmpltd_mv,levels=c("1","0"))
data$f13frmcmpltdby_mv.factor = factor(data$f13frmcmpltdby_mv,levels=c("1","0"))
data$f13stoppnt_complete.factor = factor(data$f13stoppnt_complete,levels=c("0","1","2"))
data$f7qcnumber_mv.factor = factor(data$f7qcnumber_mv,levels=c("1","0"))
data$f7visitdate_mv.factor = factor(data$f7visitdate_mv,levels=c("1","0"))
data$f7tstintrvl_mv.factor = factor(data$f7tstintrvl_mv,levels=c("1","0"))
data$f7typesample.factor = factor(data$f7typesample,levels=c("1","2","3","4","5","6"))
data$f7typesample_mv.factor = factor(data$f7typesample_mv,levels=c("1","0"))
data$f7labqc_complete.factor = factor(data$f7labqc_complete,levels=c("0","1","2"))
data$f12expdatevis_mv.factor = factor(data$f12expdatevis_mv,levels=c("1","0"))
data$f12typevisit.factor = factor(data$f12typevisit,levels=c("3","5"))
data$f12typevisit_mv.factor = factor(data$f12typevisit_mv,levels=c("1","0"))
data$f12testintrvl_mv.factor = factor(data$f12testintrvl_mv,levels=c("1","0"))
data$f12reasnmissd.factor = factor(data$f12reasnmissd,levels=c("1","2","3","4","5","6","7","8","9"))
data$f12reasnmissd_mv.factor = factor(data$f12reasnmissd_mv,levels=c("1","0"))
data$f12otherspec_mv.factor = factor(data$f12otherspec_mv,levels=c("1","0"))
data$f12missedvisit_complete.factor = factor(data$f12missedvisit_complete,levels=c("0","1","2"))
data$f18visitdate_mv.factor = factor(data$f18visitdate_mv,levels=c("1","0"))
data$f18qcnumber_mv.factor = factor(data$f18qcnumber_mv,levels=c("1","0"))
data$f18hba1c_mv.factor = factor(data$f18hba1c_mv,levels=c("1","0"))
data$f18p0creatn_mv.factor = factor(data$f18p0creatn_mv,levels=c("1","0"))
data$f18p0albumin_mv.factor = factor(data$f18p0albumin_mv,levels=c("1","0"))
data$f18p0igg_mv.factor = factor(data$f18p0igg_mv,levels=c("1","0"))
data$f18u0creatn_mv.factor = factor(data$f18u0creatn_mv,levels=c("1","0"))
data$f18u0albumin_mv.factor = factor(data$f18u0albumin_mv,levels=c("1","0"))
data$f18u0igg_mv.factor = factor(data$f18u0igg_mv,levels=c("1","0"))
data$f18daesclrncqc_complete.factor = factor(data$f18daesclrncqc_complete,levels=c("0","1","2"))
data$f19qcnumber_mv.factor = factor(data$f19qcnumber_mv,levels=c("1","0"))
data$f19urnflowrate.factor = factor(data$f19urnflowrate,levels=c("1","2"))
data$f19urnflowrate_mv.factor = factor(data$f19urnflowrate_mv,levels=c("1","0"))
data$f19assaydate_mv.factor = factor(data$f19assaydate_mv,levels=c("1","0"))
data$f19urnioth3qc_mv.factor = factor(data$f19urnioth3qc_mv,levels=c("1","0"))
data$f19urnioth3df_mv.factor = factor(data$f19urnioth3df_mv,levels=c("1","0"))
data$f19serumioth3qc_mv.factor = factor(data$f19serumioth3qc_mv,levels=c("1","0"))
data$f19serumioth3df_mv.factor = factor(data$f19serumioth3df_mv,levels=c("1","0"))
data$f19serumioth4qc_mv.factor = factor(data$f19serumioth4qc_mv,levels=c("1","0"))
data$f19serumioth4df_mv.factor = factor(data$f19serumioth4df_mv,levels=c("1","0"))
data$f19urnpah3qc_mv.factor = factor(data$f19urnpah3qc_mv,levels=c("1","0"))
data$f19urnpah3df_mv.factor = factor(data$f19urnpah3df_mv,levels=c("1","0"))
data$f19serumpah3qc_mv.factor = factor(data$f19serumpah3qc_mv,levels=c("1","0"))
data$f19serumpah3df_mv.factor = factor(data$f19serumpah3df_mv,levels=c("1","0"))
data$f19serumpah4qc_mv.factor = factor(data$f19serumpah4qc_mv,levels=c("1","0"))
data$f19serumpah4df_mv.factor = factor(data$f19serumpah4df_mv,levels=c("1","0"))
data$f19renlfuncqcrslt_complete.factor = factor(data$f19renlfuncqcrslt_complete,levels=c("0","1","2"))
data$phxlabs_visitdate_mv.factor = factor(data$phxlabs_visitdate_mv,levels=c("1","0"))
data$scr_mv.factor = factor(data$scr_mv,levels=c("1","0"))
data$gfr_mv.factor = factor(data$gfr_mv,levels=c("1","0"))
data$uricacid_mv.factor = factor(data$uricacid_mv,levels=c("1","0"))
data$hba1a_mv.factor = factor(data$hba1a_mv,levels=c("1","0"))
data$hba1b_mv.factor = factor(data$hba1b_mv,levels=c("1","0"))
data$hba1c_mv.factor = factor(data$hba1c_mv,levels=c("1","0"))
data$hba1o_mv.factor = factor(data$hba1o_mv,levels=c("1","0"))
data$hbf_mv.factor = factor(data$hbf_mv,levels=c("1","0"))
data$p0albumin_mv.factor = factor(data$p0albumin_mv,levels=c("1","0"))
data$p0gluc_mv.factor = factor(data$p0gluc_mv,levels=c("1","0"))
data$p0igg_mv.factor = factor(data$p0igg_mv,levels=c("1","0"))
data$p0oncoticprsr_mv.factor = factor(data$p0oncoticprsr_mv,levels=c("1","0"))
data$p3albumin_mv.factor = factor(data$p3albumin_mv,levels=c("1","0"))
data$p3igg_mv.factor = factor(data$p3igg_mv,levels=c("1","0"))
data$p3oncoticprsr_mv.factor = factor(data$p3oncoticprsr_mv,levels=c("1","0"))
data$pahclrnc_mv.factor = factor(data$pahclrnc_mv,levels=c("1","0"))
data$u0igg_is_below_limit.factor = factor(data$u0igg_is_below_limit,levels=c("1","0"))
data$u0igg_is_below_limit_mv.factor = factor(data$u0igg_is_below_limit_mv,levels=c("1","0"))
data$u0igg_mv.factor = factor(data$u0igg_mv,levels=c("1","0"))
data$u3flow_gfr_mv.factor = factor(data$u3flow_gfr_mv,levels=c("1","0"))
data$u3gfr_mv.factor = factor(data$u3gfr_mv,levels=c("1","0"))
data$u3use_gfr_mv.factor = factor(data$u3use_gfr_mv,levels=c("1","0"))
data$u3albumin_is_below_limit.factor = factor(data$u3albumin_is_below_limit,levels=c("1","0"))
data$u3albumin_is_below_limit_mv.factor = factor(data$u3albumin_is_below_limit_mv,levels=c("1","0"))
data$u3albumin_mv.factor = factor(data$u3albumin_mv,levels=c("1","0"))
data$u3igg_is_below_limit.factor = factor(data$u3igg_is_below_limit,levels=c("1","0"))
data$u3igg_is_below_limit_mv.factor = factor(data$u3igg_is_below_limit_mv,levels=c("1","0"))
data$u3igg_mv.factor = factor(data$u3igg_mv,levels=c("1","0"))
data$ualb_is_below_limit.factor = factor(data$ualb_is_below_limit,levels=c("1","0"))
data$ualb_is_below_limit_mv.factor = factor(data$ualb_is_below_limit_mv,levels=c("1","0"))
data$ualb_mv.factor = factor(data$ualb_mv,levels=c("1","0"))
data$ucr_mv.factor = factor(data$ucr_mv,levels=c("1","0"))
data$phoenix_labs_complete.factor = factor(data$phoenix_labs_complete,levels=c("0","1","2"))
data$cantpv_2_mv.factor = factor(data$cantpv_2_mv,levels=c("1","0"))
data$cantpv_3.factor = factor(data$cantpv_3,levels=c("1","0"))
data$cantpv_3_mv.factor = factor(data$cantpv_3_mv,levels=c("1","0"))
data$cantpv_4.factor = factor(data$cantpv_4,levels=c("1","0"))
data$cantpv_4_mv.factor = factor(data$cantpv_4_mv,levels=c("1","0"))
data$cantpv_5.factor = factor(data$cantpv_5,levels=c("1","0"))
data$cantpv_5_mv.factor = factor(data$cantpv_5_mv,levels=c("1","0"))
data$cantpv_6.factor = factor(data$cantpv_6,levels=c("1","0"))
data$cantpv_6_mv.factor = factor(data$cantpv_6_mv,levels=c("1","0"))
data$cantpv_7.factor = factor(data$cantpv_7,levels=c("1","0"))
data$cantpv_7_mv.factor = factor(data$cantpv_7_mv,levels=c("1","0"))
data$cantpv_8.factor = factor(data$cantpv_8,levels=c("1","0"))
data$cantpv_8_mv.factor = factor(data$cantpv_8_mv,levels=c("1","0"))
data$cantpv_9_mv.factor = factor(data$cantpv_9_mv,levels=c("1","0"))
data$cantpv_10.factor = factor(data$cantpv_10,levels=c("1","0"))
data$cantpv_10_mv.factor = factor(data$cantpv_10_mv,levels=c("1","0"))
data$cantpv_12_mv.factor = factor(data$cantpv_12_mv,levels=c("1","0"))
data$cantpv_13.factor = factor(data$cantpv_13,levels=c("1","0"))
data$cantpv_13_mv.factor = factor(data$cantpv_13_mv,levels=c("1","0"))
data$cantpv_14_mv.factor = factor(data$cantpv_14_mv,levels=c("1","0"))
data$cantpv_15.factor = factor(data$cantpv_15,levels=c("1","0"))
data$cantpv_15_mv.factor = factor(data$cantpv_15_mv,levels=c("1","0"))
data$cantpv_16_mv.factor = factor(data$cantpv_16_mv,levels=c("1","0"))
data$cantpv_17.factor = factor(data$cantpv_17,levels=c("1","0"))
data$cantpv_17_mv.factor = factor(data$cantpv_17_mv,levels=c("1","0"))
data$cantpv_18.factor = factor(data$cantpv_18,levels=c("1","0"))
data$cantpv_18_mv.factor = factor(data$cantpv_18_mv,levels=c("1","0"))
data$cantpv_19_mv.factor = factor(data$cantpv_19_mv,levels=c("1","0"))
data$cantpv_20_mv.factor = factor(data$cantpv_20_mv,levels=c("1","0"))
data$cantpv_21_mv.factor = factor(data$cantpv_21_mv,levels=c("1","0"))
data$cantpv_22_mv.factor = factor(data$cantpv_22_mv,levels=c("1","0"))
data$cantpv_23_mv.factor = factor(data$cantpv_23_mv,levels=c("1","0"))
data$cantpv_24_mv.factor = factor(data$cantpv_24_mv,levels=c("1","0"))
data$cantpv_25_mv.factor = factor(data$cantpv_25_mv,levels=c("1","0"))
data$cantpv_26_mv.factor = factor(data$cantpv_26_mv,levels=c("1","0"))
data$cantpv_27_mv.factor = factor(data$cantpv_27_mv,levels=c("1","0"))
data$cantpv_28_mv.factor = factor(data$cantpv_28_mv,levels=c("1","0"))
data$cantpv_29_mv.factor = factor(data$cantpv_29_mv,levels=c("1","0"))
data$cardiac_autonomic_neuropathy_test_preparation_veri_complete.factor = factor(data$cardiac_autonomic_neuropathy_test_preparation_veri_complete,levels=c("0","1","2"))
data$mnsi_cl_2_mv.factor = factor(data$mnsi_cl_2_mv,levels=c("1","0"))
data$mnsi_cl_3.factor = factor(data$mnsi_cl_3,levels=c("0","1"))
data$mnsi_cl_3_mv.factor = factor(data$mnsi_cl_3_mv,levels=c("1","0"))
data$mnsi_cl_4___1.factor = factor(data$mnsi_cl_4___1,levels=c("0","1"))
data$mnsi_cl_4___2.factor = factor(data$mnsi_cl_4___2,levels=c("0","1"))
data$mnsi_cl_4___3.factor = factor(data$mnsi_cl_4___3,levels=c("0","1"))
data$mnsi_cl_4___4.factor = factor(data$mnsi_cl_4___4,levels=c("0","1"))
data$mnsi_cl_4___0.factor = factor(data$mnsi_cl_4___0,levels=c("0","1"))
data$mnsi_cl_5_mv.factor = factor(data$mnsi_cl_5_mv,levels=c("1","0"))
data$mnsi_cl_6.factor = factor(data$mnsi_cl_6,levels=c("0","1"))
data$mnsi_cl_6_mv.factor = factor(data$mnsi_cl_6_mv,levels=c("1","0"))
data$mnsi_cl_7.factor = factor(data$mnsi_cl_7,levels=c("0","0.5","1"))
data$mnsi_cl_7_mv.factor = factor(data$mnsi_cl_7_mv,levels=c("1","0"))
data$mnsi_cl_8.factor = factor(data$mnsi_cl_8,levels=c("0","0.5","1"))
data$mnsi_cl_8_mv.factor = factor(data$mnsi_cl_8_mv,levels=c("1","0"))
data$mnsi_cl_9.factor = factor(data$mnsi_cl_9,levels=c("0","0.5","1"))
data$mnsi_cl_9_mv.factor = factor(data$mnsi_cl_9_mv,levels=c("1","0"))
data$mnsi_cl_10.factor = factor(data$mnsi_cl_10,levels=c("0","1"))
data$mnsi_cl_10_mv.factor = factor(data$mnsi_cl_10_mv,levels=c("1","0"))
data$mnsi_cl_11___1.factor = factor(data$mnsi_cl_11___1,levels=c("0","1"))
data$mnsi_cl_11___2.factor = factor(data$mnsi_cl_11___2,levels=c("0","1"))
data$mnsi_cl_11___3.factor = factor(data$mnsi_cl_11___3,levels=c("0","1"))
data$mnsi_cl_11___4.factor = factor(data$mnsi_cl_11___4,levels=c("0","1"))
data$mnsi_cl_11___0.factor = factor(data$mnsi_cl_11___0,levels=c("0","1"))
data$mnsi_cl_12_mv.factor = factor(data$mnsi_cl_12_mv,levels=c("1","0"))
data$mnsi_cl_13.factor = factor(data$mnsi_cl_13,levels=c("0","1"))
data$mnsi_cl_13_mv.factor = factor(data$mnsi_cl_13_mv,levels=c("1","0"))
data$mnsi_cl_14.factor = factor(data$mnsi_cl_14,levels=c("0","0.5","1"))
data$mnsi_cl_14_mv.factor = factor(data$mnsi_cl_14_mv,levels=c("1","0"))
data$mnsi_cl_15.factor = factor(data$mnsi_cl_15,levels=c("0","0.5","1"))
data$mnsi_cl_15_mv.factor = factor(data$mnsi_cl_15_mv,levels=c("1","0"))
data$mnsi_cl_16.factor = factor(data$mnsi_cl_16,levels=c("0","0.5","1"))
data$mnsi_cl_16_mv.factor = factor(data$mnsi_cl_16_mv,levels=c("1","0"))
data$mnsi_cl_17_mv.factor = factor(data$mnsi_cl_17_mv,levels=c("1","0"))
data$michigan_neuropathy_screening_instrument_for_clini_complete.factor = factor(data$michigan_neuropathy_screening_instrument_for_clini_complete,levels=c("0","1","2"))
data$mnsi_subj_2_mv.factor = factor(data$mnsi_subj_2_mv,levels=c("1","0"))
data$mnsi_subj_3.factor = factor(data$mnsi_subj_3,levels=c("1","0"))
data$mnsi_subj_3_mv.factor = factor(data$mnsi_subj_3_mv,levels=c("1","0"))
data$mnsi_subj_4.factor = factor(data$mnsi_subj_4,levels=c("1","0"))
data$mnsi_subj_4_mv.factor = factor(data$mnsi_subj_4_mv,levels=c("1","0"))
data$mnsi_subj_5.factor = factor(data$mnsi_subj_5,levels=c("1","0"))
data$mnsi_subj_5_mv.factor = factor(data$mnsi_subj_5_mv,levels=c("1","0"))
data$mnsi_subj_6.factor = factor(data$mnsi_subj_6,levels=c("1","0"))
data$mnsi_subj_6_mv.factor = factor(data$mnsi_subj_6_mv,levels=c("1","0"))
data$mnsi_subj_7.factor = factor(data$mnsi_subj_7,levels=c("1","0"))
data$mnsi_subj_7_mv.factor = factor(data$mnsi_subj_7_mv,levels=c("1","0"))
data$mnsi_subj_8.factor = factor(data$mnsi_subj_8,levels=c("1","0"))
data$mnsi_subj_8_mv.factor = factor(data$mnsi_subj_8_mv,levels=c("1","0"))
data$mnsi_subj_9.factor = factor(data$mnsi_subj_9,levels=c("0","1"))
data$mnsi_subj_9_mv.factor = factor(data$mnsi_subj_9_mv,levels=c("1","0"))
data$mnsi_subj_10.factor = factor(data$mnsi_subj_10,levels=c("1","0"))
data$mnsi_subj_10_mv.factor = factor(data$mnsi_subj_10_mv,levels=c("1","0"))
data$mnsi_subj_11.factor = factor(data$mnsi_subj_11,levels=c("1","0"))
data$mnsi_subj_11_mv.factor = factor(data$mnsi_subj_11_mv,levels=c("1","0"))
data$mnsi_subj_12.factor = factor(data$mnsi_subj_12,levels=c("1","0"))
data$mnsi_subj_12_mv.factor = factor(data$mnsi_subj_12_mv,levels=c("1","0"))
data$mnsi_subj_13.factor = factor(data$mnsi_subj_13,levels=c("1","0"))
data$mnsi_subj_13_mv.factor = factor(data$mnsi_subj_13_mv,levels=c("1","0"))
data$mnsi_subj_14.factor = factor(data$mnsi_subj_14,levels=c("1","0"))
data$mnsi_subj_14_mv.factor = factor(data$mnsi_subj_14_mv,levels=c("1","0"))
data$mnsi_subj_15.factor = factor(data$mnsi_subj_15,levels=c("0","1"))
data$mnsi_subj_15_mv.factor = factor(data$mnsi_subj_15_mv,levels=c("1","0"))
data$mnsi_subj_16.factor = factor(data$mnsi_subj_16,levels=c("1","0"))
data$mnsi_subj_16_mv.factor = factor(data$mnsi_subj_16_mv,levels=c("1","0"))
data$mnsi_subj_17.factor = factor(data$mnsi_subj_17,levels=c("1","0"))
data$mnsi_subj_17_mv.factor = factor(data$mnsi_subj_17_mv,levels=c("1","0"))
data$mnsi_subj_18_mv.factor = factor(data$mnsi_subj_18_mv,levels=c("1","0"))
data$michigan_neuropathy_screening_instrument_for_subje_complete.factor = factor(data$michigan_neuropathy_screening_instrument_for_subje_complete,levels=c("0","1","2"))
data$aspq_1_mv.factor = factor(data$aspq_1_mv,levels=c("1","0"))
data$aspq_2_mv.factor = factor(data$aspq_2_mv,levels=c("1","0"))
data$aspq_3_mv.factor = factor(data$aspq_3_mv,levels=c("1","0"))
data$aspq_4.factor = factor(data$aspq_4,levels=c("1","0"))
data$aspq_4_mv.factor = factor(data$aspq_4_mv,levels=c("1","0"))
data$aspq_5.factor = factor(data$aspq_5,levels=c("1","2","3","4"))
data$aspq_5_mv.factor = factor(data$aspq_5_mv,levels=c("1","0"))
data$aspq_6.factor = factor(data$aspq_6,levels=c("1","2","3"))
data$aspq_6_mv.factor = factor(data$aspq_6_mv,levels=c("1","0"))
data$aspq_7.factor = factor(data$aspq_7,levels=c("1","2","3","4","5","6"))
data$aspq_7_mv.factor = factor(data$aspq_7_mv,levels=c("1","0"))
data$aspq_8.factor = factor(data$aspq_8,levels=c("0","1","2","3","4","5"))
data$aspq_8_mv.factor = factor(data$aspq_8_mv,levels=c("1","0"))
data$aspq_9.factor = factor(data$aspq_9,levels=c("1","2","3"))
data$aspq_9_mv.factor = factor(data$aspq_9_mv,levels=c("1","0"))
data$aspq_10.factor = factor(data$aspq_10,levels=c("1","2","3","4","5","6","7"))
data$aspq_10_mv.factor = factor(data$aspq_10_mv,levels=c("1","0"))
data$aspq_11_mv.factor = factor(data$aspq_11_mv,levels=c("1","0"))
data$aspq_12.factor = factor(data$aspq_12,levels=c("1","2","3","4","5","6"))
data$aspq_12_mv.factor = factor(data$aspq_12_mv,levels=c("1","0"))
data$aspq_14.factor = factor(data$aspq_14,levels=c("1","2","3","4"))
data$aspq_14_mv.factor = factor(data$aspq_14_mv,levels=c("1","0"))
data$aspq_15.factor = factor(data$aspq_15,levels=c("1","2","3","4"))
data$aspq_15_mv.factor = factor(data$aspq_15_mv,levels=c("1","0"))
data$aspq_16.factor = factor(data$aspq_16,levels=c("1","2","3","4"))
data$aspq_16_mv.factor = factor(data$aspq_16_mv,levels=c("1","0"))
data$aspq_17.factor = factor(data$aspq_17,levels=c("1","2","3","4"))
data$aspq_17_mv.factor = factor(data$aspq_17_mv,levels=c("1","0"))
data$aspq_18.factor = factor(data$aspq_18,levels=c("1","2","3","4"))
data$aspq_18_mv.factor = factor(data$aspq_18_mv,levels=c("1","0"))
data$aspq_19.factor = factor(data$aspq_19,levels=c("1","2","3","4"))
data$aspq_19_mv.factor = factor(data$aspq_19_mv,levels=c("1","0"))
data$aspq_20.factor = factor(data$aspq_20,levels=c("1","2","3","4"))
data$aspq_20_mv.factor = factor(data$aspq_20_mv,levels=c("1","0"))
data$aspq_21.factor = factor(data$aspq_21,levels=c("1","2","3","4"))
data$aspq_21_mv.factor = factor(data$aspq_21_mv,levels=c("1","0"))
data$aspq_22.factor = factor(data$aspq_22,levels=c("1","2","3","4"))
data$aspq_22_mv.factor = factor(data$aspq_22_mv,levels=c("1","0"))
data$aspq_23.factor = factor(data$aspq_23,levels=c("1","2","3","4"))
data$aspq_23_mv.factor = factor(data$aspq_23_mv,levels=c("1","0"))
data$aspq_24.factor = factor(data$aspq_24,levels=c("1","0"))
data$aspq_24_mv.factor = factor(data$aspq_24_mv,levels=c("1","0"))
data$aspq_25_mv.factor = factor(data$aspq_25_mv,levels=c("1","0"))
data$aspq_26.factor = factor(data$aspq_26,levels=c("1","0"))
data$aspq_26_mv.factor = factor(data$aspq_26_mv,levels=c("1","0"))
data$aspq_27.factor = factor(data$aspq_27,levels=c("1","0"))
data$aspq_27_mv.factor = factor(data$aspq_27_mv,levels=c("1","0"))
data$aspq_28.factor = factor(data$aspq_28,levels=c("1","0"))
data$aspq_28_mv.factor = factor(data$aspq_28_mv,levels=c("1","0"))
data$aspq_29.factor = factor(data$aspq_29,levels=c("1","0"))
data$aspq_29_mv.factor = factor(data$aspq_29_mv,levels=c("1","0"))
data$aspq_30.factor = factor(data$aspq_30,levels=c("1","0"))
data$aspq_30_mv.factor = factor(data$aspq_30_mv,levels=c("1","0"))
data$aspq_31.factor = factor(data$aspq_31,levels=c("1","0"))
data$aspq_31_mv.factor = factor(data$aspq_31_mv,levels=c("1","0"))
data$aspq_32.factor = factor(data$aspq_32,levels=c("1","0"))
data$aspq_32_mv.factor = factor(data$aspq_32_mv,levels=c("1","0"))
data$aspq_34.factor = factor(data$aspq_34,levels=c("1","0"))
data$aspq_34_mv.factor = factor(data$aspq_34_mv,levels=c("1","0"))
data$aspq_35.factor = factor(data$aspq_35,levels=c("1","0"))
data$aspq_35_mv.factor = factor(data$aspq_35_mv,levels=c("1","0"))
data$aspq_36.factor = factor(data$aspq_36,levels=c("1","0"))
data$aspq_36_mv.factor = factor(data$aspq_36_mv,levels=c("1","0"))
data$aspq_37_mv.factor = factor(data$aspq_37_mv,levels=c("1","0"))
data$aspq_38.factor = factor(data$aspq_38,levels=c("1","0"))
data$aspq_38_mv.factor = factor(data$aspq_38_mv,levels=c("1","0"))
data$aspq_39.factor = factor(data$aspq_39,levels=c("1","0"))
data$aspq_39_mv.factor = factor(data$aspq_39_mv,levels=c("1","0"))
data$aspq_40_mv.factor = factor(data$aspq_40_mv,levels=c("1","0"))
data$autonomic_symptoms_profile_questionnaire_complete.factor = factor(data$autonomic_symptoms_profile_questionnaire_complete,levels=c("0","1","2"))
data$nsqlq_2_mv.factor = factor(data$nsqlq_2_mv,levels=c("1","0"))
data$nsqlq_3_mv.factor = factor(data$nsqlq_3_mv,levels=c("1","0"))
data$nsqlq_5.factor = factor(data$nsqlq_5,levels=c("1","2","3","4","0"))
data$nsqlq_5_mv.factor = factor(data$nsqlq_5_mv,levels=c("1","0"))
data$nsqlq_6.factor = factor(data$nsqlq_6,levels=c("1","2","0"))
data$nsqlq_6_mv.factor = factor(data$nsqlq_6_mv,levels=c("1","0"))
data$nsqlq_7.factor = factor(data$nsqlq_7,levels=c("1","2","3","4","0"))
data$nsqlq_7_mv.factor = factor(data$nsqlq_7_mv,levels=c("1","0"))
data$nsqlq_8.factor = factor(data$nsqlq_8,levels=c("1","2","0"))
data$nsqlq_8_mv.factor = factor(data$nsqlq_8_mv,levels=c("1","0"))
data$nsqlq_9.factor = factor(data$nsqlq_9,levels=c("1","2","3","4","0"))
data$nsqlq_9_mv.factor = factor(data$nsqlq_9_mv,levels=c("1","0"))
data$nsqlq_10.factor = factor(data$nsqlq_10,levels=c("1","2","0"))
data$nsqlq_10_mv.factor = factor(data$nsqlq_10_mv,levels=c("1","0"))
data$nsqlq_11.factor = factor(data$nsqlq_11,levels=c("1","2","3","4","0"))
data$nsqlq_11_mv.factor = factor(data$nsqlq_11_mv,levels=c("1","0"))
data$nsqlq_12.factor = factor(data$nsqlq_12,levels=c("1","2","0"))
data$nsqlq_12_mv.factor = factor(data$nsqlq_12_mv,levels=c("1","0"))
data$nsqlq_13.factor = factor(data$nsqlq_13,levels=c("1","2","3","4","0"))
data$nsqlq_13_mv.factor = factor(data$nsqlq_13_mv,levels=c("1","0"))
data$nsqlq_14.factor = factor(data$nsqlq_14,levels=c("1","2","0"))
data$nsqlq_14_mv.factor = factor(data$nsqlq_14_mv,levels=c("1","0"))
data$nsqlq_15.factor = factor(data$nsqlq_15,levels=c("1","2","3","4","0"))
data$nsqlq_15_mv.factor = factor(data$nsqlq_15_mv,levels=c("1","0"))
data$nsqlq_16.factor = factor(data$nsqlq_16,levels=c("1","2","0"))
data$nsqlq_16_mv.factor = factor(data$nsqlq_16_mv,levels=c("1","0"))
data$nsqlq_17.factor = factor(data$nsqlq_17,levels=c("1","2","3","4","0"))
data$nsqlq_17_mv.factor = factor(data$nsqlq_17_mv,levels=c("1","0"))
data$nsqlq_18.factor = factor(data$nsqlq_18,levels=c("1","2","0"))
data$nsqlq_18_mv.factor = factor(data$nsqlq_18_mv,levels=c("1","0"))
data$nsqlq_19.factor = factor(data$nsqlq_19,levels=c("1","2","3","4","0"))
data$nsqlq_19_mv.factor = factor(data$nsqlq_19_mv,levels=c("1","0"))
data$nsqlq_20.factor = factor(data$nsqlq_20,levels=c("1","2","3","4","0"))
data$nsqlq_20_mv.factor = factor(data$nsqlq_20_mv,levels=c("1","0"))
data$nsqlq_21.factor = factor(data$nsqlq_21,levels=c("1","2","0"))
data$nsqlq_21_mv.factor = factor(data$nsqlq_21_mv,levels=c("1","0"))
data$nsqlq_22.factor = factor(data$nsqlq_22,levels=c("1","2","3","4","0"))
data$nsqlq_22_mv.factor = factor(data$nsqlq_22_mv,levels=c("1","0"))
data$nsqlq_23.factor = factor(data$nsqlq_23,levels=c("1","2","0"))
data$nsqlq_23_mv.factor = factor(data$nsqlq_23_mv,levels=c("1","0"))
data$nsqlq_24.factor = factor(data$nsqlq_24,levels=c("1","2","3","4","0"))
data$nsqlq_24_mv.factor = factor(data$nsqlq_24_mv,levels=c("1","0"))
data$nsqlq_25.factor = factor(data$nsqlq_25,levels=c("1","2","0"))
data$nsqlq_25_mv.factor = factor(data$nsqlq_25_mv,levels=c("1","0"))
data$nsqlq_26.factor = factor(data$nsqlq_26,levels=c("1","2","3","4","0"))
data$nsqlq_26_mv.factor = factor(data$nsqlq_26_mv,levels=c("1","0"))
data$nsqlq_27.factor = factor(data$nsqlq_27,levels=c("1","2","3","4","0"))
data$nsqlq_27_mv.factor = factor(data$nsqlq_27_mv,levels=c("1","0"))
data$nsqlq_28.factor = factor(data$nsqlq_28,levels=c("1","2","0"))
data$nsqlq_28_mv.factor = factor(data$nsqlq_28_mv,levels=c("1","0"))
data$nsqlq_29.factor = factor(data$nsqlq_29,levels=c("1","2","3","4","0"))
data$nsqlq_29_mv.factor = factor(data$nsqlq_29_mv,levels=c("1","0"))
data$nsqlq_30.factor = factor(data$nsqlq_30,levels=c("1","2","0"))
data$nsqlq_30_mv.factor = factor(data$nsqlq_30_mv,levels=c("1","0"))
data$nsqlq_31.factor = factor(data$nsqlq_31,levels=c("1","2","3","4","0"))
data$nsqlq_31_mv.factor = factor(data$nsqlq_31_mv,levels=c("1","0"))
data$nsqlq_32.factor = factor(data$nsqlq_32,levels=c("1","2","0"))
data$nsqlq_32_mv.factor = factor(data$nsqlq_32_mv,levels=c("1","0"))
data$nsqlq_33.factor = factor(data$nsqlq_33,levels=c("1","2","3","4","0"))
data$nsqlq_33_mv.factor = factor(data$nsqlq_33_mv,levels=c("1","0"))
data$nsqlq_34.factor = factor(data$nsqlq_34,levels=c("1","0"))
data$nsqlq_34_mv.factor = factor(data$nsqlq_34_mv,levels=c("1","0"))
data$nsqlq_35.factor = factor(data$nsqlq_35,levels=c("1","2","3","4","0"))
data$nsqlq_35_mv.factor = factor(data$nsqlq_35_mv,levels=c("1","0"))
data$nsqlq_36.factor = factor(data$nsqlq_36,levels=c("1","2","0"))
data$nsqlq_36_mv.factor = factor(data$nsqlq_36_mv,levels=c("1","0"))
data$nsqlq_37.factor = factor(data$nsqlq_37,levels=c("1","2","3","4","0"))
data$nsqlq_37_mv.factor = factor(data$nsqlq_37_mv,levels=c("1","0"))
data$nsqlq_38.factor = factor(data$nsqlq_38,levels=c("1","2","0"))
data$nsqlq_38_mv.factor = factor(data$nsqlq_38_mv,levels=c("1","0"))
data$nsqlq_39.factor = factor(data$nsqlq_39,levels=c("1","2","3","4","0"))
data$nsqlq_39_mv.factor = factor(data$nsqlq_39_mv,levels=c("1","0"))
data$nsqlq_40.factor = factor(data$nsqlq_40,levels=c("1","2","0"))
data$nsqlq_40_mv.factor = factor(data$nsqlq_40_mv,levels=c("1","0"))
data$nsqlq_41.factor = factor(data$nsqlq_41,levels=c("1","2","3","4","0"))
data$nsqlq_41_mv.factor = factor(data$nsqlq_41_mv,levels=c("1","0"))
data$nsqlq_42.factor = factor(data$nsqlq_42,levels=c("1","2","3","4","0"))
data$nsqlq_42_mv.factor = factor(data$nsqlq_42_mv,levels=c("1","0"))
data$nsqlq_43.factor = factor(data$nsqlq_43,levels=c("1","2","0"))
data$nsqlq_43_mv.factor = factor(data$nsqlq_43_mv,levels=c("1","0"))
data$nsqlq_44.factor = factor(data$nsqlq_44,levels=c("1","2","3","4","0"))
data$nsqlq_44_mv.factor = factor(data$nsqlq_44_mv,levels=c("1","0"))
data$nsqlq_45.factor = factor(data$nsqlq_45,levels=c("1","2","0"))
data$nsqlq_45_mv.factor = factor(data$nsqlq_45_mv,levels=c("1","0"))
data$nsqlq_46.factor = factor(data$nsqlq_46,levels=c("1","2","3","4","0"))
data$nsqlq_46_mv.factor = factor(data$nsqlq_46_mv,levels=c("1","0"))
data$nsqlq_47.factor = factor(data$nsqlq_47,levels=c("1","2","0"))
data$nsqlq_47_mv.factor = factor(data$nsqlq_47_mv,levels=c("1","0"))
data$nsqlq_48.factor = factor(data$nsqlq_48,levels=c("1","2","3","4","0"))
data$nsqlq_48_mv.factor = factor(data$nsqlq_48_mv,levels=c("1","0"))
data$nsqlq_49.factor = factor(data$nsqlq_49,levels=c("1","2","0"))
data$nsqlq_49_mv.factor = factor(data$nsqlq_49_mv,levels=c("1","0"))
data$nsqlq_50.factor = factor(data$nsqlq_50,levels=c("1","2","3","4","0"))
data$nsqlq_50_mv.factor = factor(data$nsqlq_50_mv,levels=c("1","0"))
data$nsqlq_51.factor = factor(data$nsqlq_51,levels=c("1","2","3","4","5"))
data$nsqlq_51_mv.factor = factor(data$nsqlq_51_mv,levels=c("1","0"))
data$nsqlq_52.factor = factor(data$nsqlq_52,levels=c("1","2","0"))
data$nsqlq_52_mv.factor = factor(data$nsqlq_52_mv,levels=c("1","0"))
data$nsqlq_53.factor = factor(data$nsqlq_53,levels=c("1","2","3","4","5"))
data$nsqlq_53_mv.factor = factor(data$nsqlq_53_mv,levels=c("1","0"))
data$nsqlq_54.factor = factor(data$nsqlq_54,levels=c("1","2","0"))
data$nsqlq_54_mv.factor = factor(data$nsqlq_54_mv,levels=c("1","0"))
data$nsqlq_55.factor = factor(data$nsqlq_55,levels=c("1","2","3","4","5"))
data$nsqlq_55_mv.factor = factor(data$nsqlq_55_mv,levels=c("1","0"))
data$nsqlq_56.factor = factor(data$nsqlq_56,levels=c("1","2","0"))
data$nsqlq_56_mv.factor = factor(data$nsqlq_56_mv,levels=c("1","0"))
data$nsqlq_57.factor = factor(data$nsqlq_57,levels=c("1","2","3","4","5"))
data$nsqlq_57_mv.factor = factor(data$nsqlq_57_mv,levels=c("1","0"))
data$nsqlq_58.factor = factor(data$nsqlq_58,levels=c("1","2","0"))
data$nsqlq_58_mv.factor = factor(data$nsqlq_58_mv,levels=c("1","0"))
data$nsqlq_59.factor = factor(data$nsqlq_59,levels=c("1","2","3","4","5"))
data$nsqlq_59_mv.factor = factor(data$nsqlq_59_mv,levels=c("1","0"))
data$nsqlq_60.factor = factor(data$nsqlq_60,levels=c("1","2","0"))
data$nsqlq_60_mv.factor = factor(data$nsqlq_60_mv,levels=c("1","0"))
data$nsqlq_61.factor = factor(data$nsqlq_61,levels=c("1","2","3","4","5"))
data$nsqlq_61_mv.factor = factor(data$nsqlq_61_mv,levels=c("1","0"))
data$nsqlq_62.factor = factor(data$nsqlq_62,levels=c("1","2","0"))
data$nsqlq_62_mv.factor = factor(data$nsqlq_62_mv,levels=c("1","0"))
data$nsqlq_63.factor = factor(data$nsqlq_63,levels=c("1","2","3","4","5"))
data$nsqlq_63_mv.factor = factor(data$nsqlq_63_mv,levels=c("1","0"))
data$nsqlq_64.factor = factor(data$nsqlq_64,levels=c("1","2","0"))
data$nsqlq_64_mv.factor = factor(data$nsqlq_64_mv,levels=c("1","0"))
data$nsqlq_65.factor = factor(data$nsqlq_65,levels=c("1","2","3","4","0"))
data$nsqlq_65_mv.factor = factor(data$nsqlq_65_mv,levels=c("1","0"))
data$nsqlq_66.factor = factor(data$nsqlq_66,levels=c("1","2","3","4","0"))
data$nsqlq_66_mv.factor = factor(data$nsqlq_66_mv,levels=c("1","0"))
data$nsqlq_67.factor = factor(data$nsqlq_67,levels=c("1","2","3","4","5"))
data$nsqlq_67_mv.factor = factor(data$nsqlq_67_mv,levels=c("1","0"))
data$neuropathy_specific_quality_of_life_questionnaire_complete.factor = factor(data$neuropathy_specific_quality_of_life_questionnaire_complete,levels=c("0","1","2"))
data$jvn_collecdate_mv.factor = factor(data$jvn_collecdate_mv,levels=c("1","0"))
data$jvn_date_mv.factor = factor(data$jvn_date_mv,levels=c("1","0"))
data$jvn_npdr_od.factor = factor(data$jvn_npdr_od,levels=c("0","1","2","3","4","5","9"))
data$jvn_npdr_od_mv.factor = factor(data$jvn_npdr_od_mv,levels=c("1","0"))
data$jvn_npdr_os.factor = factor(data$jvn_npdr_os,levels=c("0","1","2","3","4","5","9"))
data$jvn_npdr_os_mv.factor = factor(data$jvn_npdr_os_mv,levels=c("1","0"))
data$jvn_npdr_notes_mv.factor = factor(data$jvn_npdr_notes_mv,levels=c("1","0"))
data$jvn_pdr_od.factor = factor(data$jvn_pdr_od,levels=c("0","1","2","3","9"))
data$jvn_pdr_od_mv.factor = factor(data$jvn_pdr_od_mv,levels=c("1","0"))
data$jvn_pdr_os.factor = factor(data$jvn_pdr_os,levels=c("0","1","2","3","9"))
data$jvn_pdr_os_mv.factor = factor(data$jvn_pdr_os_mv,levels=c("1","0"))
data$jvn_pdr_notes_mv.factor = factor(data$jvn_pdr_notes_mv,levels=c("1","0"))
data$jvn_me_od.factor = factor(data$jvn_me_od,levels=c("0","1","2","9"))
data$jvn_me_od_mv.factor = factor(data$jvn_me_od_mv,levels=c("1","0"))
data$jvn_me_os.factor = factor(data$jvn_me_os,levels=c("0","1","2","9"))
data$jvn_me_os_mv.factor = factor(data$jvn_me_os_mv,levels=c("1","0"))
data$jvn_me_notes_mv.factor = factor(data$jvn_me_notes_mv,levels=c("1","0"))
data$jvn_amd_od.factor = factor(data$jvn_amd_od,levels=c("1","0"))
data$jvn_amd_od_mv.factor = factor(data$jvn_amd_od_mv,levels=c("1","0"))
data$jvn_amd_os.factor = factor(data$jvn_amd_os,levels=c("1","0"))
data$jvn_amd_os_mv.factor = factor(data$jvn_amd_os_mv,levels=c("1","0"))
data$jvn_asterhy_od.factor = factor(data$jvn_asterhy_od,levels=c("1","0"))
data$jvn_asterhy_od_mv.factor = factor(data$jvn_asterhy_od_mv,levels=c("1","0"))
data$jvn_asterhy_os.factor = factor(data$jvn_asterhy_os,levels=c("1","0"))
data$jvn_asterhy_os_mv.factor = factor(data$jvn_asterhy_os_mv,levels=c("1","0"))
data$jvn_cat_od.factor = factor(data$jvn_cat_od,levels=c("1","0"))
data$jvn_cat_od_mv.factor = factor(data$jvn_cat_od_mv,levels=c("1","0"))
data$jvn_cat_os.factor = factor(data$jvn_cat_os,levels=c("1","0"))
data$jvn_cat_os_mv.factor = factor(data$jvn_cat_os_mv,levels=c("1","0"))
data$jvn_chorles_od.factor = factor(data$jvn_chorles_od,levels=c("1","0"))
data$jvn_chorles_od_mv.factor = factor(data$jvn_chorles_od_mv,levels=c("1","0"))
data$jvn_chorles_os.factor = factor(data$jvn_chorles_os,levels=c("1","0"))
data$jvn_chorles_os_mv.factor = factor(data$jvn_chorles_os_mv,levels=c("1","0"))
data$jvn_arcus_od.factor = factor(data$jvn_arcus_od,levels=c("1","0"))
data$jvn_arcus_od_mv.factor = factor(data$jvn_arcus_od_mv,levels=c("1","0"))
data$jvn_arcus_os.factor = factor(data$jvn_arcus_os,levels=c("1","0"))
data$jvn_arcus_os_mv.factor = factor(data$jvn_arcus_os_mv,levels=c("1","0"))
data$jvn_cratscar_od.factor = factor(data$jvn_cratscar_od,levels=c("1","0"))
data$jvn_cratscar_od_mv.factor = factor(data$jvn_cratscar_od_mv,levels=c("1","0"))
data$jvn_cratscar_os.factor = factor(data$jvn_cratscar_os,levels=c("1","0"))
data$jvn_cratscar_os_mv.factor = factor(data$jvn_cratscar_os_mv,levels=c("1","0"))
data$jvn_epiret_od.factor = factor(data$jvn_epiret_od,levels=c("1","0"))
data$jvn_epiret_od_mv.factor = factor(data$jvn_epiret_od_mv,levels=c("1","0"))
data$jvn_epiret_os.factor = factor(data$jvn_epiret_os,levels=c("1","0"))
data$jvn_epiret_os_mv.factor = factor(data$jvn_epiret_os_mv,levels=c("1","0"))
data$jvn_iol_od.factor = factor(data$jvn_iol_od,levels=c("1","0"))
data$jvn_iol_od_mv.factor = factor(data$jvn_iol_od_mv,levels=c("1","0"))
data$jvn_iol_os.factor = factor(data$jvn_iol_os,levels=c("1","0"))
data$jvn_iol_os_mv.factor = factor(data$jvn_iol_os_mv,levels=c("1","0"))
data$jvn_nevus_od.factor = factor(data$jvn_nevus_od,levels=c("1","0"))
data$jvn_nevus_od_mv.factor = factor(data$jvn_nevus_od_mv,levels=c("1","0"))
data$jvn_nevus_os.factor = factor(data$jvn_nevus_os,levels=c("1","0"))
data$jvn_nevus_os_mv.factor = factor(data$jvn_nevus_os_mv,levels=c("1","0"))
data$jvn_hardex_od.factor = factor(data$jvn_hardex_od,levels=c("1","0"))
data$jvn_hardex_od_mv.factor = factor(data$jvn_hardex_od_mv,levels=c("1","0"))
data$jvn_hardex_os.factor = factor(data$jvn_hardex_os,levels=c("1","0"))
data$jvn_hardex_os_mv.factor = factor(data$jvn_hardex_os_mv,levels=c("1","0"))
data$jvn_lidles_od.factor = factor(data$jvn_lidles_od,levels=c("1","0"))
data$jvn_lidles_od_mv.factor = factor(data$jvn_lidles_od_mv,levels=c("1","0"))
data$jvn_lidles_os.factor = factor(data$jvn_lidles_os,levels=c("1","0"))
data$jvn_lidles_os_mv.factor = factor(data$jvn_lidles_os_mv,levels=c("1","0"))
data$jvn_drusen_od.factor = factor(data$jvn_drusen_od,levels=c("1","0"))
data$jvn_drusen_od_mv.factor = factor(data$jvn_drusen_od_mv,levels=c("1","0"))
data$jvn_drusen_os.factor = factor(data$jvn_drusen_os,levels=c("1","0"))
data$jvn_drusen_os_mv.factor = factor(data$jvn_drusen_os_mv,levels=c("1","0"))
data$jvn_ondrus_od.factor = factor(data$jvn_ondrus_od,levels=c("1","0"))
data$jvn_ondrus_od_mv.factor = factor(data$jvn_ondrus_od_mv,levels=c("1","0"))
data$jvn_ondrus_os.factor = factor(data$jvn_ondrus_os,levels=c("1","0"))
data$jvn_ondrus_os_mv.factor = factor(data$jvn_ondrus_os_mv,levels=c("1","0"))
data$jvn_miscas_od.factor = factor(data$jvn_miscas_od,levels=c("1","0"))
data$jvn_miscas_od_mv.factor = factor(data$jvn_miscas_od_mv,levels=c("1","0"))
data$jvn_miscas_os.factor = factor(data$jvn_miscas_os,levels=c("1","0"))
data$jvn_miscas_os_mv.factor = factor(data$jvn_miscas_os_mv,levels=c("1","0"))
data$jvn_miscvr_od.factor = factor(data$jvn_miscvr_od,levels=c("1","0"))
data$jvn_miscvr_od_mv.factor = factor(data$jvn_miscvr_od_mv,levels=c("1","0"))
data$jvn_miscvr_os.factor = factor(data$jvn_miscvr_os,levels=c("1","0"))
data$jvn_miscvr_os_mv.factor = factor(data$jvn_miscvr_os_mv,levels=c("1","0"))
data$jvn_ppa_od.factor = factor(data$jvn_ppa_od,levels=c("1","0"))
data$jvn_ppa_od_mv.factor = factor(data$jvn_ppa_od_mv,levels=c("1","0"))
data$jvn_ppa_os.factor = factor(data$jvn_ppa_os,levels=c("1","0"))
data$jvn_ppa_os_mv.factor = factor(data$jvn_ppa_os_mv,levels=c("1","0"))
data$jvn_prfib_od.factor = factor(data$jvn_prfib_od,levels=c("1","0"))
data$jvn_prfib_od_mv.factor = factor(data$jvn_prfib_od_mv,levels=c("1","0"))
data$jvn_prfib_os.factor = factor(data$jvn_prfib_os,levels=c("1","0"))
data$jvn_prfib_os_mv.factor = factor(data$jvn_prfib_os_mv,levels=c("1","0"))
data$jvn_htnret_od.factor = factor(data$jvn_htnret_od,levels=c("1","0"))
data$jvn_htnret_od_mv.factor = factor(data$jvn_htnret_od_mv,levels=c("1","0"))
data$jvn_htnret_os.factor = factor(data$jvn_htnret_os,levels=c("1","0"))
data$jvn_htnret_os_mv.factor = factor(data$jvn_htnret_os_mv,levels=c("1","0"))
data$jvp_rpe_od.factor = factor(data$jvp_rpe_od,levels=c("1","0"))
data$jvp_rpe_od_mv.factor = factor(data$jvp_rpe_od_mv,levels=c("1","0"))
data$jvn_rpe_os.factor = factor(data$jvn_rpe_os,levels=c("1","0"))
data$jvn_rpe_os_mv.factor = factor(data$jvn_rpe_os_mv,levels=c("1","0"))
data$jvn_rrd_od.factor = factor(data$jvn_rrd_od,levels=c("1","0"))
data$jvn_rrd_od_mv.factor = factor(data$jvn_rrd_od_mv,levels=c("1","0"))
data$jvn_rrd_os.factor = factor(data$jvn_rrd_os,levels=c("1","0"))
data$jvn_rrd_os_mv.factor = factor(data$jvn_rrd_os_mv,levels=c("1","0"))
data$jvn_trd_od.factor = factor(data$jvn_trd_od,levels=c("1","0"))
data$jvn_trd_od_mv.factor = factor(data$jvn_trd_od_mv,levels=c("1","0"))
data$jvn_trd_os.factor = factor(data$jvn_trd_os,levels=c("1","0"))
data$jvn_trd_os_mv.factor = factor(data$jvn_trd_os_mv,levels=c("1","0"))
data$jvn_vithem_od.factor = factor(data$jvn_vithem_od,levels=c("1","0"))
data$jvn_vithem_od_mv.factor = factor(data$jvn_vithem_od_mv,levels=c("1","0"))
data$jvn_vithem_os.factor = factor(data$jvn_vithem_os,levels=c("1","0"))
data$jvn_vithem_os_mv.factor = factor(data$jvn_vithem_os_mv,levels=c("1","0"))
data$jvn_other_od.factor = factor(data$jvn_other_od,levels=c("1","0"))
data$jvn_other_od_mv.factor = factor(data$jvn_other_od_mv,levels=c("1","0"))
data$jvn_other_os.factor = factor(data$jvn_other_os,levels=c("1","0"))
data$jvn_other_os_mv.factor = factor(data$jvn_other_os_mv,levels=c("1","0"))
data$jvn_laser_od.factor = factor(data$jvn_laser_od,levels=c("1","0"))
data$jvn_laser_od_mv.factor = factor(data$jvn_laser_od_mv,levels=c("1","0"))
data$jvn_laser_os.factor = factor(data$jvn_laser_os,levels=c("1","0"))
data$jvn_laser_os_mv.factor = factor(data$jvn_laser_os_mv,levels=c("1","0"))
data$jvn_notes_mv.factor = factor(data$jvn_notes_mv,levels=c("1","0"))
data$jvn_deck_8_mv.factor = factor(data$jvn_deck_8_mv,levels=c("1","0"))
data$jvn_retinal_imaging_form_complete.factor = factor(data$jvn_retinal_imaging_form_complete,levels=c("0","1","2"))

levels(data$redcap_event_name.factor)=c("Interval 0","Interval 12","Interval 24","Interval 36","Interval 48","Interval 60","Interval 72","Interval 84","Interval 96","Interval 108","Interval 120","Interval 132","Interval 144","Interval 156","Interval 168","Interval 180","Interval 192","Interval 204","Interval 216","Interval 228")
levels(data$redcap_repeat_instrument.factor)=c("")
levels(data$lastname_mv.factor)=c("Yes","No")
levels(data$firstname_mv.factor)=c("Yes","No")
levels(data$stratum.factor)=c("Normal UAE","Micro UAE","Other")
levels(data$stratum_mv.factor)=c("Yes","No")
levels(data$randomnum_mv.factor)=c("Yes","No")
levels(data$phx_number_mv.factor)=c("Yes","No")
levels(data$dob_mv.factor)=c("Yes","No")
levels(data$gender.factor)=c("male","non-fertile female","fertile female")
levels(data$gender_mv.factor)=c("Yes","No")
levels(data$height_mv.factor)=c("Yes","No")
levels(data$allergies.factor)=c("No","Yes")
levels(data$allergies_mv.factor)=c("Yes","No")
levels(data$allerglist_mv.factor)=c("Yes","No")
levels(data$district_mv.factor)=c("Yes","No")
levels(data$phonenum_mv.factor)=c("Yes","No")
levels(data$bldprsr.factor)=c("Right","Left")
levels(data$bldprsr_mv.factor)=c("Yes","No")
levels(data$dtdiabonset_mv.factor)=c("Yes","No")
levels(data$typerecruit.factor)=c("Walk-in/Called us","Word of mouth","Referral from current study participant","Community event","Other")
levels(data$typerecruit_mv.factor)=c("Yes","No")
levels(data$typerecruit_other_mv.factor)=c("Yes","No")
levels(data$familyrelation.factor)=c("Parent","Sibling","Child","Other (Specify)")
levels(data$familyrelation_mv.factor)=c("Yes","No")
levels(data$familyrelation_other_mv.factor)=c("Yes","No")
levels(data$nih_number_mv.factor)=c("Yes","No")
levels(data$demoform_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$sex_code.factor)=c("Male","Female")
levels(data$sex_code_mv.factor)=c("Yes","No")
levels(data$start_date_mv.factor)=c("Yes","No")
levels(data$end_date_mv.factor)=c("Yes","No")
levels(data$study_status.factor)=c("Complete","Incomplete")
levels(data$study_status_mv.factor)=c("Yes","No")
levels(data$clrncvisittype.factor)=c("Clearance","Follow-up")
levels(data$clrncvisittype_mv.factor)=c("Yes","No")
levels(data$clrnctestintvrl_mv.factor)=c("Yes","No")
levels(data$comments_mv.factor)=c("Yes","No")
levels(data$study_code.factor)=c("FTP","EIR","61","RET","EEP","81","82","5YR","DBD","PAB","EE1","EE2","EE3","EE4","EE5","EEA","WP5","EE6","EE7","DPP","LP3","T4SS","T4SW1","36")
levels(data$study_code_mv.factor)=c("Yes","No")
levels(data$encounter_id_mv.factor)=c("Yes","No")
levels(data$clearancevisit_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$brainmri_done.factor)=c("Yes","No")
levels(data$brainmri_done_mv.factor)=c("Yes","No")
levels(data$brainmri_date_mv.factor)=c("Yes","No")
levels(data$cognitivetest_done.factor)=c("Yes","No")
levels(data$cognitivetest_done_mv.factor)=c("Yes","No")
levels(data$cognitivetest_date_mv.factor)=c("Yes","No")
levels(data$kidneybiopsy_done.factor)=c("Yes","No")
levels(data$kidneybiopsy_done_mv.factor)=c("Yes","No")
levels(data$kidneybiopsy_date_mv.factor)=c("Yes","No")
levels(data$kidneybiopsy_site.factor)=c("Right","Left")
levels(data$kidneybiopsy_site_mv.factor)=c("Yes","No")
levels(data$kidneylivermre_done.factor)=c("Yes","No")
levels(data$kidneylivermre_done_mv.factor)=c("Yes","No")
levels(data$kidneylivermre_date_mv.factor)=c("Yes","No")
levels(data$kidneymri_done.factor)=c("Yes","No")
levels(data$kidneymri_done_mv.factor)=c("Yes","No")
levels(data$kidneymri_date_mv.factor)=c("Yes","No")
levels(data$skinbiopsy_done.factor)=c("Yes","No")
levels(data$skinbiopsy_done_mv.factor)=c("Yes","No")
levels(data$skinbiopsy_date_mv.factor)=c("Yes","No")
levels(data$skinbiopsy_site.factor)=c("Right proximal calf","Left proximal calf")
levels(data$skinbiopsy_site_mv.factor)=c("Yes","No")
levels(data$stemcells_done.factor)=c("Yes","No")
levels(data$stemcells_done_mv.factor)=c("Yes","No")
levels(data$stemcells_date_mv.factor)=c("Yes","No")
levels(data$stool_done.factor)=c("Yes","No")
levels(data$stool_done_mv.factor)=c("Yes","No")
levels(data$stool_date_mv.factor)=c("Yes","No")
levels(data$neuropathy_done.factor)=c("Yes","No")
levels(data$neuropathy_done_mv.factor)=c("Yes","No")
levels(data$neuropathy_date_mv.factor)=c("Yes","No")
levels(data$retinopathy_done.factor)=c("Yes","No")
levels(data$retinopathy_done_mv.factor)=c("Yes","No")
levels(data$retinopathy_date_mv.factor)=c("Yes","No")
levels(data$sp_comments_mv.factor)=c("Yes","No")
levels(data$study_procedures_status_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f4bldprsr.factor)=c("Right","Left")
levels(data$f4bldprsr_mv.factor)=c("Yes","No")
levels(data$f4frstmsrsys_mv.factor)=c("Yes","No")
levels(data$f4frstmsrdia_mv.factor)=c("Yes","No")
levels(data$f4scndmsrsys_mv.factor)=c("Yes","No")
levels(data$f4scndmsrdia_mv.factor)=c("Yes","No")
levels(data$f4heartrate_mv.factor)=c("Yes","No")
levels(data$f4pregnant.factor)=c("No","Yes")
levels(data$f4pregnant_mv.factor)=c("Yes","No")
levels(data$f4medsnow.factor)=c("No","Yes")
levels(data$f4medsnow_mv.factor)=c("Yes","No")
levels(data$f4srmpotasm_mv.factor)=c("Yes","No")
levels(data$f4srmglu_mv.factor)=c("Yes","No")
levels(data$f4forminit_mv.factor)=c("Yes","No")
levels(data$f4bldprsrandmeds_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f4drugcode1_mv.factor)=c("Yes","No")
levels(data$f4dosage1_mv.factor)=c("Yes","No")
levels(data$f4timesday1_mv.factor)=c("Yes","No")
levels(data$f4dwm1.factor)=c("Day","Week","Month")
levels(data$f4dwm1_mv.factor)=c("Yes","No")
levels(data$f4drugcode2_mv.factor)=c("Yes","No")
levels(data$f4dosage2_mv.factor)=c("Yes","No")
levels(data$f4timesday2_mv.factor)=c("Yes","No")
levels(data$f4dwm2.factor)=c("Day","Week","Month")
levels(data$f4dwm2_mv.factor)=c("Yes","No")
levels(data$f4drugcode3_mv.factor)=c("Yes","No")
levels(data$f4dosage3_mv.factor)=c("Yes","No")
levels(data$f4timesday3_mv.factor)=c("Yes","No")
levels(data$f4dwm3.factor)=c("Day","Week","Month")
levels(data$f4dwm3_mv.factor)=c("Yes","No")
levels(data$f4drugcode4_mv.factor)=c("Yes","No")
levels(data$f4dosage4_mv.factor)=c("Yes","No")
levels(data$f4timesday4_mv.factor)=c("Yes","No")
levels(data$f4dwm4.factor)=c("Day","Week","Month")
levels(data$f4dwm4_mv.factor)=c("Yes","No")
levels(data$f4drugcode5_mv.factor)=c("Yes","No")
levels(data$f4dosage5_mv.factor)=c("Yes","No")
levels(data$f4timesday5_mv.factor)=c("Yes","No")
levels(data$f4dwm5.factor)=c("Day","Week","Month")
levels(data$f4dwm5_mv.factor)=c("Yes","No")
levels(data$f4drugcode6_mv.factor)=c("Yes","No")
levels(data$f4dosage6_mv.factor)=c("Yes","No")
levels(data$f4timesday6_mv.factor)=c("Yes","No")
levels(data$f4dwm6.factor)=c("Day","Week","Month")
levels(data$f4dwm6_mv.factor)=c("Yes","No")
levels(data$f4drugcode7_mv.factor)=c("Yes","No")
levels(data$f4dosage7_mv.factor)=c("Yes","No")
levels(data$f4timesday7_mv.factor)=c("Yes","No")
levels(data$f4dwm7.factor)=c("Day","Week","Month")
levels(data$f4dwm7_mv.factor)=c("Yes","No")
levels(data$f4drugcode8_mv.factor)=c("Yes","No")
levels(data$f4dosage8_mv.factor)=c("Yes","No")
levels(data$f4timesday8_mv.factor)=c("Yes","No")
levels(data$f4dwm8.factor)=c("Day","Week","Month")
levels(data$f4dwm8_mv.factor)=c("Yes","No")
levels(data$f4drugcode9_mv.factor)=c("Yes","No")
levels(data$f4dosage9_mv.factor)=c("Yes","No")
levels(data$f4timesday9_mv.factor)=c("Yes","No")
levels(data$f4dwm9.factor)=c("Day","Week","Month")
levels(data$f4dwm9_mv.factor)=c("Yes","No")
levels(data$f4drugcode10_mv.factor)=c("Yes","No")
levels(data$f4dosage10_mv.factor)=c("Yes","No")
levels(data$f4timesday10_mv.factor)=c("Yes","No")
levels(data$f4dwm10.factor)=c("Day","Week","Month")
levels(data$f4dwm10_mv.factor)=c("Yes","No")
levels(data$f4drugcode11_mv.factor)=c("Yes","No")
levels(data$f4dosage11_mv.factor)=c("Yes","No")
levels(data$f4timesday11_mv.factor)=c("Yes","No")
levels(data$f4dwm11.factor)=c("Day","Week","Month")
levels(data$f4dwm11_mv.factor)=c("Yes","No")
levels(data$f4drugcode12_mv.factor)=c("Yes","No")
levels(data$f4dosage12_mv.factor)=c("Yes","No")
levels(data$f4timesday12_mv.factor)=c("Yes","No")
levels(data$f4dwm12.factor)=c("Day","Week","Month")
levels(data$f4dwm12_mv.factor)=c("Yes","No")
levels(data$f4drugcode13_mv.factor)=c("Yes","No")
levels(data$f4dosage13_mv.factor)=c("Yes","No")
levels(data$f4timesday13_mv.factor)=c("Yes","No")
levels(data$f4dwm13.factor)=c("Day","Week","Month")
levels(data$f4dwm13_mv.factor)=c("Yes","No")
levels(data$f4drugcode14_mv.factor)=c("Yes","No")
levels(data$f4dosage14_mv.factor)=c("Yes","No")
levels(data$f4timesday14_mv.factor)=c("Yes","No")
levels(data$f4dwm14.factor)=c("Day","Week","Month")
levels(data$f4dwm14_mv.factor)=c("Yes","No")
levels(data$f4drugcode15_mv.factor)=c("Yes","No")
levels(data$f4dosage15_mv.factor)=c("Yes","No")
levels(data$f4timesday15_mv.factor)=c("Yes","No")
levels(data$f4dwm15.factor)=c("Day","Week","Month")
levels(data$f4dwm15_mv.factor)=c("Yes","No")
levels(data$f4drugcode16_mv.factor)=c("Yes","No")
levels(data$f4dosage16_mv.factor)=c("Yes","No")
levels(data$f4timesday16_mv.factor)=c("Yes","No")
levels(data$f4dwm16.factor)=c("Day","Week","Month")
levels(data$f4dwm16_mv.factor)=c("Yes","No")
levels(data$f4drugcode17_mv.factor)=c("Yes","No")
levels(data$f4dosage17_mv.factor)=c("Yes","No")
levels(data$f4timesday17_mv.factor)=c("Yes","No")
levels(data$f4dwm17.factor)=c("Day","Week","Month")
levels(data$f4dwm17_mv.factor)=c("Yes","No")
levels(data$f4drugcode18_mv.factor)=c("Yes","No")
levels(data$f4dosage18_mv.factor)=c("Yes","No")
levels(data$f4timesday18_mv.factor)=c("Yes","No")
levels(data$f4dwm18.factor)=c("Day","Week","Month")
levels(data$f4dwm18_mv.factor)=c("Yes","No")
levels(data$f4drugcode19_mv.factor)=c("Yes","No")
levels(data$f4dosage19_mv.factor)=c("Yes","No")
levels(data$f4timesday19_mv.factor)=c("Yes","No")
levels(data$f4dwm19.factor)=c("Day","Week","Month")
levels(data$f4dwm19_mv.factor)=c("Yes","No")
levels(data$f4drugcode20_mv.factor)=c("Yes","No")
levels(data$f4dosage20_mv.factor)=c("Yes","No")
levels(data$f4timesday20_mv.factor)=c("Yes","No")
levels(data$f4dwm20.factor)=c("Day","Week","Month")
levels(data$f4dwm20_mv.factor)=c("Yes","No")
levels(data$f4drugcode21_mv.factor)=c("Yes","No")
levels(data$f4dosage21_mv.factor)=c("Yes","No")
levels(data$f4timesday21_mv.factor)=c("Yes","No")
levels(data$f4dwm21.factor)=c("Day","Week","Month")
levels(data$f4dwm21_mv.factor)=c("Yes","No")
levels(data$f4drugcode22_mv.factor)=c("Yes","No")
levels(data$f4dosage22_mv.factor)=c("Yes","No")
levels(data$f4timesday22_mv.factor)=c("Yes","No")
levels(data$f4dwm22.factor)=c("Day","Week","Month")
levels(data$f4dwm22_mv.factor)=c("Yes","No")
levels(data$f4drugcode23_mv.factor)=c("Yes","No")
levels(data$f4dosage23_mv.factor)=c("Yes","No")
levels(data$f4timesday23_mv.factor)=c("Yes","No")
levels(data$f4dwm23.factor)=c("Day","Week","Month")
levels(data$f4dwm23_mv.factor)=c("Yes","No")
levels(data$f4drugcode24_mv.factor)=c("Yes","No")
levels(data$f4dosage24_mv.factor)=c("Yes","No")
levels(data$f4timesday24_mv.factor)=c("Yes","No")
levels(data$f4dwm24.factor)=c("Day","Week","Month")
levels(data$f4dwm24_mv.factor)=c("Yes","No")
levels(data$f4drugcode25_mv.factor)=c("Yes","No")
levels(data$f4dosage25_mv.factor)=c("Yes","No")
levels(data$f4timesday25_mv.factor)=c("Yes","No")
levels(data$f4dwm25.factor)=c("Day","Week","Month")
levels(data$f4dwm25_mv.factor)=c("Yes","No")
levels(data$f4drugcode26_mv.factor)=c("Yes","No")
levels(data$f4dosage26_mv.factor)=c("Yes","No")
levels(data$f4timesday26_mv.factor)=c("Yes","No")
levels(data$f4dwm26.factor)=c("Day","Week","Month")
levels(data$f4dwm26_mv.factor)=c("Yes","No")
levels(data$f4drugcode27_mv.factor)=c("Yes","No")
levels(data$f4dosage27_mv.factor)=c("Yes","No")
levels(data$f4timesday27_mv.factor)=c("Yes","No")
levels(data$f4dwm27.factor)=c("Day","Week","Month")
levels(data$f4dwm27_mv.factor)=c("Yes","No")
levels(data$f4drugcode28_mv.factor)=c("Yes","No")
levels(data$f4dosage28_mv.factor)=c("Yes","No")
levels(data$f4timesday28_mv.factor)=c("Yes","No")
levels(data$f4dwm28.factor)=c("Day","Week","Month")
levels(data$f4dwm28_mv.factor)=c("Yes","No")
levels(data$f4drugcode29_mv.factor)=c("Yes","No")
levels(data$f4dosage29_mv.factor)=c("Yes","No")
levels(data$f4timesday29_mv.factor)=c("Yes","No")
levels(data$f4dwm29.factor)=c("Day","Week","Month")
levels(data$f4dwm29_mv.factor)=c("Yes","No")
levels(data$f4drugcode30_mv.factor)=c("Yes","No")
levels(data$f4dosage30_mv.factor)=c("Yes","No")
levels(data$f4timesday30_mv.factor)=c("Yes","No")
levels(data$f4dwm30.factor)=c("Day","Week","Month")
levels(data$f4dwm30_mv.factor)=c("Yes","No")
levels(data$f4meds_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f5rcntillns.factor)=c("No","Yes")
levels(data$f5rcntillns_mv.factor)=c("Yes","No")
levels(data$f5outpat.factor)=c("No","Yes")
levels(data$f5outpat_mv.factor)=c("Yes","No")
levels(data$f5outpatspec_mv.factor)=c("Yes","No")
levels(data$f5hospwosrg.factor)=c("No","Yes")
levels(data$f5hospwosrg_mv.factor)=c("Yes","No")
levels(data$f5hospwoexp_mv.factor)=c("Yes","No")
levels(data$f5hospwsrg.factor)=c("No","Yes")
levels(data$f5hospwsrg_mv.factor)=c("Yes","No")
levels(data$f5hospwexp_mv.factor)=c("Yes","No")
levels(data$f5medprbsdoc.factor)=c("No","Yes")
levels(data$f5medprbsdoc_mv.factor)=c("Yes","No")
levels(data$f5crnryartds.factor)=c("No","Yes")
levels(data$f5crnryartds_mv.factor)=c("Yes","No")
levels(data$f5cancer.factor)=c("No","Yes")
levels(data$f5cancer_mv.factor)=c("Yes","No")
levels(data$f5cancerexp_mv.factor)=c("Yes","No")
levels(data$f5crbrlvasds.factor)=c("No","Yes")
levels(data$f5crbrlvasds_mv.factor)=c("Yes","No")
levels(data$f5peripvasds.factor)=c("No","Yes")
levels(data$f5peripvasds_mv.factor)=c("Yes","No")
levels(data$f5hypertensn.factor)=c("No","Yes")
levels(data$f5hypertensn_mv.factor)=c("Yes","No")
levels(data$f5seizures.factor)=c("No","Yes")
levels(data$f5seizures_mv.factor)=c("Yes","No")
levels(data$f5gnitrnryds.factor)=c("No","Yes")
levels(data$f5gnitrnryds_mv.factor)=c("Yes","No")
levels(data$f5gnitrnyexp_mv.factor)=c("Yes","No")
levels(data$f5lungdsz.factor)=c("No","Yes")
levels(data$f5lungdsz_mv.factor)=c("Yes","No")
levels(data$f5majsurg.factor)=c("No","Yes")
levels(data$f5majsurg_mv.factor)=c("Yes","No")
levels(data$f5majsurgexp_mv.factor)=c("Yes","No")
levels(data$f5othrmeddia.factor)=c("No","Yes")
levels(data$f5othrmeddia_mv.factor)=c("Yes","No")
levels(data$f5othrmedexp_mv.factor)=c("Yes","No")
levels(data$f5lungs.factor)=c("Normal","Abnormal")
levels(data$f5lungs_mv.factor)=c("Yes","No")
levels(data$f5lungsexp_mv.factor)=c("Yes","No")
levels(data$f5heart.factor)=c("Normal","Abnormal")
levels(data$f5heart_mv.factor)=c("Yes","No")
levels(data$f5heartexp_mv.factor)=c("Yes","No")
levels(data$f5skin.factor)=c("Normal","Abnormal")
levels(data$f5skin_mv.factor)=c("Yes","No")
levels(data$f5skinexp_mv.factor)=c("Yes","No")
levels(data$f5edema.factor)=c("Normal","Abnormal")
levels(data$f5edema_mv.factor)=c("Yes","No")
levels(data$f5edemaexp_mv.factor)=c("Yes","No")
levels(data$f5histandphysexm_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f6ivsites_mv.factor)=c("Yes","No")
levels(data$f6weight_mv.factor)=c("Yes","No")
levels(data$f6u0strttime_mv.factor)=c("Yes","No")
levels(data$f6u0volume_mv.factor)=c("Yes","No")
levels(data$f6bldprsrarm.factor)=c("Right","Left")
levels(data$f6bldprsrarm_mv.factor)=c("Yes","No")
levels(data$f6frstmsrsty_mv.factor)=c("Yes","No")
levels(data$f6frstmsrdia_mv.factor)=c("Yes","No")
levels(data$f6scndmsrsys_mv.factor)=c("Yes","No")
levels(data$f6scndmsrdia_mv.factor)=c("Yes","No")
levels(data$f6heartrate_mv.factor)=c("Yes","No")
levels(data$f6u0h20_mv.factor)=c("Yes","No")
levels(data$f6fnleqlurtm_mv.factor)=c("Yes","No")
levels(data$f6p1time_mv.factor)=c("Yes","No")
levels(data$f6urnvol_mv.factor)=c("Yes","No")
levels(data$f6h20vol_mv.factor)=c("Yes","No")
levels(data$f6u1time_mv.factor)=c("Yes","No")
levels(data$f6p2time_mv.factor)=c("Yes","No")
levels(data$f6u1vol_mv.factor)=c("Yes","No")
levels(data$f6u1h20vol_mv.factor)=c("Yes","No")
levels(data$f6u2time_mv.factor)=c("Yes","No")
levels(data$f6p3time_mv.factor)=c("Yes","No")
levels(data$f6u2vol_mv.factor)=c("Yes","No")
levels(data$f6u2h20vol_mv.factor)=c("Yes","No")
levels(data$f6u3time_mv.factor)=c("Yes","No")
levels(data$f6p4time_mv.factor)=c("Yes","No")
levels(data$f6u3vol_mv.factor)=c("Yes","No")
levels(data$f6u3h20vol_mv.factor)=c("Yes","No")
levels(data$f6u4time_mv.factor)=c("Yes","No")
levels(data$f6p5time_mv.factor)=c("Yes","No")
levels(data$f6u4vol_mv.factor)=c("Yes","No")
levels(data$f6u4h20vol_mv.factor)=c("Yes","No")
levels(data$f6u5time_mv.factor)=c("Yes","No")
levels(data$f6p6time_mv.factor)=c("Yes","No")
levels(data$f6u5vol_mv.factor)=c("Yes","No")
levels(data$f6u5h20vol_mv.factor)=c("Yes","No")
levels(data$f6urndate_mv.factor)=c("Yes","No")
levels(data$f6glucose.factor)=c("Neg","Trace (100 mg/dl)","250 mg/dl","500 mg/dl","1000 mg/dl",">2000 mg/dl")
levels(data$f6glucose_mv.factor)=c("Yes","No")
levels(data$f6bilirubin.factor)=c("Neg","1+ (small)","2+ (moderate)","3+ (large)")
levels(data$f6bilirubin_mv.factor)=c("Yes","No")
levels(data$f6ketones.factor)=c("Neg","Trace (5)","Small (15)","Moderate (40)","80",">160")
levels(data$f6ketones_mv.factor)=c("Yes","No")
levels(data$f6specgrvty_mv.factor)=c("Yes","No")
levels(data$f6bloodocult.factor)=c("No","Yes")
levels(data$f6bloodocult_mv.factor)=c("Yes","No")
levels(data$f6ph_mv.factor)=c("Yes","No")
levels(data$f6protein.factor)=c("Neg","Trace","1+ 30 mg/dl","2+ 100 mg/dl","3+ 300 mg/dl","4+ > 2000 mg/dl")
levels(data$f6protein_mv.factor)=c("Yes","No")
levels(data$f6casts_mv.factor)=c("Yes","No")
levels(data$f6rbccast.factor)=c("No","Yes")
levels(data$f6rbccast_mv.factor)=c("Yes","No")
levels(data$f6wbccast.factor)=c("No","Yes")
levels(data$f6wbccast_mv.factor)=c("Yes","No")
levels(data$f6hyalinecast.factor)=c("No","Yes")
levels(data$f6hyalinecast_mv.factor)=c("Yes","No")
levels(data$f6granularcast.factor)=c("No","Yes")
levels(data$f6granularcast_mv.factor)=c("Yes","No")
levels(data$f6othercast.factor)=c("No","Yes")
levels(data$f6othercast_mv.factor)=c("Yes","No")
levels(data$f6othercast_exp_mv.factor)=c("Yes","No")
levels(data$f6wbc.factor)=c("0-5/HPF,Negative,Trace,Rare,Occassional, Few","6-20/HPF, 1+, 2+, Some, Moderate",">20/HPF,3+,4+,Great,Many,TNTC")
levels(data$f6wbc_mv.factor)=c("Yes","No")
levels(data$f6rbc.factor)=c("0-5/HPF,Negative,Trace,Rare,Occassional, Few","6-20/HPF, 1+, 2+, Some, Moderate",">20/HPF,3+,4+,Great,Many,TNTC")
levels(data$f6rbc_mv.factor)=c("Yes","No")
levels(data$f6epithcells.factor)=c("0-5/HPF,Negative,Trace,Rare,Occassional, Few","6-20/HPF, 1+, 2+, Some, Moderate",">20/HPF,3+,4+,Great,Many,TNTC")
levels(data$f6epithcells_mv.factor)=c("Yes","No")
levels(data$f6bacteria.factor)=c("0-5/HPF,Negative,Trace,Rare,Occassional, Few","6-20/HPF, 1+, 2+, Some, Moderate",">20/HPF,3+,4+,Great,Many,TNTC")
levels(data$f6bacteria_mv.factor)=c("Yes","No")
levels(data$f6renalclrnc_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f8visitdate_mv.factor)=c("Yes","No")
levels(data$f8analdate_mv.factor)=c("Yes","No")
levels(data$f8alti_mv.factor)=c("Yes","No")
levels(data$f8glucose_mv.factor)=c("Yes","No")
levels(data$f8bun_mv.factor)=c("Yes","No")
levels(data$f8creatinine_mv.factor)=c("Yes","No")
levels(data$f8uricacid_mv.factor)=c("Yes","No")
levels(data$f8sodium_mv.factor)=c("Yes","No")
levels(data$f8potassium_mv.factor)=c("Yes","No")
levels(data$f8chloride_mv.factor)=c("Yes","No")
levels(data$f8calcium_mv.factor)=c("Yes","No")
levels(data$f8totlprotn_mv.factor)=c("Yes","No")
levels(data$f8albumin_mv.factor)=c("Yes","No")
levels(data$f8cholestero_mv.factor)=c("Yes","No")
levels(data$f8triglyceri_mv.factor)=c("Yes","No")
levels(data$f8hdl_mv.factor)=c("Yes","No")
levels(data$f8ldl_mv.factor)=c("Yes","No")
levels(data$f8totlbilrbn_mv.factor)=c("Yes","No")
levels(data$f8dirctbilrbn_mv.factor)=c("Yes","No")
levels(data$f8alkp04_mv.factor)=c("Yes","No")
levels(data$f8ast_mv.factor)=c("Yes","No")
levels(data$f8smacchempanl_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f9dateanal_mv.factor)=c("Yes","No")
levels(data$f9analyzed_date_mv.factor)=c("Yes","No")
levels(data$f9wbc_mv.factor)=c("Yes","No")
levels(data$f9rbc_mv.factor)=c("Yes","No")
levels(data$f9hemoglob_mv.factor)=c("Yes","No")
levels(data$f9hematocrit_mv.factor)=c("Yes","No")
levels(data$f9mcv_mv.factor)=c("Yes","No")
levels(data$f9mch_mv.factor)=c("Yes","No")
levels(data$f9mchc_mv.factor)=c("Yes","No")
levels(data$f9platelets_mv.factor)=c("Yes","No")
levels(data$f9lymphinst_mv.factor)=c("Yes","No")
levels(data$f9monoinst_mv.factor)=c("Yes","No")
levels(data$f9neutinst_mv.factor)=c("Yes","No")
levels(data$f9eosinst_mv.factor)=c("Yes","No")
levels(data$f9basoinst_mv.factor)=c("Yes","No")
levels(data$f9gran_mv.factor)=c("Yes","No")
levels(data$f9pt_mv.factor)=c("Yes","No")
levels(data$f9inr_mv.factor)=c("Yes","No")
levels(data$f9ptt_mv.factor)=c("Yes","No")
levels(data$f9cbc_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f16qcnumber_mv.factor)=c("Yes","No")
levels(data$f16analdate_mv.factor)=c("Yes","No")
levels(data$f16glucose_mv.factor)=c("Yes","No")
levels(data$f16bun_mv.factor)=c("Yes","No")
levels(data$f16creatinine_mv.factor)=c("Yes","No")
levels(data$f16sodium_mv.factor)=c("Yes","No")
levels(data$f16potassium_mv.factor)=c("Yes","No")
levels(data$f16chloride_mv.factor)=c("Yes","No")
levels(data$f16calcium_mv.factor)=c("Yes","No")
levels(data$f16phosphorus_mv.factor)=c("Yes","No")
levels(data$f16totlprotn_mv.factor)=c("Yes","No")
levels(data$f16albumin_mv.factor)=c("Yes","No")
levels(data$f16cholestero_mv.factor)=c("Yes","No")
levels(data$f16triglyceri_mv.factor)=c("Yes","No")
levels(data$f16totlbilrbn_mv.factor)=c("Yes","No")
levels(data$f16alkp04_mv.factor)=c("Yes","No")
levels(data$f16ast_mv.factor)=c("Yes","No")
levels(data$f16smacqc_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f17qcnumber_mv.factor)=c("Yes","No")
levels(data$f17analdate_mv.factor)=c("Yes","No")
levels(data$f17wbc_mv.factor)=c("Yes","No")
levels(data$f17rbc_mv.factor)=c("Yes","No")
levels(data$f17hemoglob_mv.factor)=c("Yes","No")
levels(data$f17hematocrit_mv.factor)=c("Yes","No")
levels(data$f17mcv_mv.factor)=c("Yes","No")
levels(data$f17mch_mv.factor)=c("Yes","No")
levels(data$f17mchc_mv.factor)=c("Yes","No")
levels(data$f17platelets_mv.factor)=c("Yes","No")
levels(data$f17cbcqc_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f13dtlastvst_mv.factor)=c("Yes","No")
levels(data$f13typevist.factor)=c("Clearance","Follow-up")
levels(data$f13typevist_mv.factor)=c("Yes","No")
levels(data$f13tstintrvl_mv.factor)=c("Yes","No")
levels(data$f13endstgrnl.factor)=c("No","Yes")
levels(data$f13endstgrnl_mv.factor)=c("Yes","No")
levels(data$f13esrddate_mv.factor)=c("Yes","No")
levels(data$f13nondiabkid.factor)=c("No","Yes")
levels(data$f13nondiabkid_mv.factor)=c("Yes","No")
levels(data$f13ndkddate_mv.factor)=c("Yes","No")
levels(data$f13imprdbldr.factor)=c("No","Yes")
levels(data$f13imprdbldr_mv.factor)=c("Yes","No")
levels(data$f13ibfdate_mv.factor)=c("Yes","No")
levels(data$f13cngstvhrt.factor)=c("No","Yes")
levels(data$f13cngstvhrt_mv.factor)=c("Yes","No")
levels(data$f13chfdate_mv.factor)=c("Yes","No")
levels(data$f13ascites.factor)=c("No","Yes")
levels(data$f13ascites_mv.factor)=c("Yes","No")
levels(data$f13ascitesdt_mv.factor)=c("Yes","No")
levels(data$f13aceinhib.factor)=c("No","Yes")
levels(data$f13aceinhib_mv.factor)=c("Yes","No")
levels(data$f13typeofreaction.factor)=c("Hyperkalemia","Angioedema","Neuthropenia","Thrombocytopenia","Symptomatic hypotension","Intractable cough","Other")
levels(data$f13typeofreaction_mv.factor)=c("Yes","No")
levels(data$f13otherdiag.factor)=c("No","Yes")
levels(data$f13otherdiag_mv.factor)=c("Yes","No")
levels(data$f13reactionspecify_mv.factor)=c("Yes","No")
levels(data$f13datereactn_mv.factor)=c("Yes","No")
levels(data$f13prfrmwdclrnc.factor)=c("No","Yes")
levels(data$f13prfrmwdclrnc_mv.factor)=c("Yes","No")
levels(data$f13prfrmwdclrncdate_mv.factor)=c("Yes","No")
levels(data$f13dtstoppntdeclrd_mv.factor)=c("Yes","No")
levels(data$f13dtfrmclmpltd_mv.factor)=c("Yes","No")
levels(data$f13frmcmpltdby_mv.factor)=c("Yes","No")
levels(data$f13stoppnt_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f7qcnumber_mv.factor)=c("Yes","No")
levels(data$f7visitdate_mv.factor)=c("Yes","No")
levels(data$f7tstintrvl_mv.factor)=c("Yes","No")
levels(data$f7typesample.factor)=c("Stanford Lab GFR sample (iothalamate and PAH)","NIH Lab P0 creatinine, albumin, IgG - DAES clearance results","SMAC","CBC","U0 albumin, IgG, and creatinine (spot urine) - DAES clearance result","NIH Lab glycosylated hemoglobin (HbA1c)")
levels(data$f7typesample_mv.factor)=c("Yes","No")
levels(data$f7labqc_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f12expdatevis_mv.factor)=c("Yes","No")
levels(data$f12typevisit.factor)=c("Clearance","Follow-up")
levels(data$f12typevisit_mv.factor)=c("Yes","No")
levels(data$f12testintrvl_mv.factor)=c("Yes","No")
levels(data$f12reasnmissd.factor)=c("illness","hospitalization","personal/family business","employment","weather","unknown","pregnancy","convalescence","Other")
levels(data$f12reasnmissd_mv.factor)=c("Yes","No")
levels(data$f12otherspec_mv.factor)=c("Yes","No")
levels(data$f12missedvisit_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f18visitdate_mv.factor)=c("Yes","No")
levels(data$f18qcnumber_mv.factor)=c("Yes","No")
levels(data$f18hba1c_mv.factor)=c("Yes","No")
levels(data$f18p0creatn_mv.factor)=c("Yes","No")
levels(data$f18p0albumin_mv.factor)=c("Yes","No")
levels(data$f18p0igg_mv.factor)=c("Yes","No")
levels(data$f18u0creatn_mv.factor)=c("Yes","No")
levels(data$f18u0albumin_mv.factor)=c("Yes","No")
levels(data$f18u0igg_mv.factor)=c("Yes","No")
levels(data$f18daesclrncqc_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f19qcnumber_mv.factor)=c("Yes","No")
levels(data$f19urnflowrate.factor)=c("No","Yes")
levels(data$f19urnflowrate_mv.factor)=c("Yes","No")
levels(data$f19assaydate_mv.factor)=c("Yes","No")
levels(data$f19urnioth3qc_mv.factor)=c("Yes","No")
levels(data$f19urnioth3df_mv.factor)=c("Yes","No")
levels(data$f19serumioth3qc_mv.factor)=c("Yes","No")
levels(data$f19serumioth3df_mv.factor)=c("Yes","No")
levels(data$f19serumioth4qc_mv.factor)=c("Yes","No")
levels(data$f19serumioth4df_mv.factor)=c("Yes","No")
levels(data$f19urnpah3qc_mv.factor)=c("Yes","No")
levels(data$f19urnpah3df_mv.factor)=c("Yes","No")
levels(data$f19serumpah3qc_mv.factor)=c("Yes","No")
levels(data$f19serumpah3df_mv.factor)=c("Yes","No")
levels(data$f19serumpah4qc_mv.factor)=c("Yes","No")
levels(data$f19serumpah4df_mv.factor)=c("Yes","No")
levels(data$f19renlfuncqcrslt_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$phxlabs_visitdate_mv.factor)=c("Yes","No")
levels(data$scr_mv.factor)=c("Yes","No")
levels(data$gfr_mv.factor)=c("Yes","No")
levels(data$uricacid_mv.factor)=c("Yes","No")
levels(data$hba1a_mv.factor)=c("Yes","No")
levels(data$hba1b_mv.factor)=c("Yes","No")
levels(data$hba1c_mv.factor)=c("Yes","No")
levels(data$hba1o_mv.factor)=c("Yes","No")
levels(data$hbf_mv.factor)=c("Yes","No")
levels(data$p0albumin_mv.factor)=c("Yes","No")
levels(data$p0gluc_mv.factor)=c("Yes","No")
levels(data$p0igg_mv.factor)=c("Yes","No")
levels(data$p0oncoticprsr_mv.factor)=c("Yes","No")
levels(data$p3albumin_mv.factor)=c("Yes","No")
levels(data$p3igg_mv.factor)=c("Yes","No")
levels(data$p3oncoticprsr_mv.factor)=c("Yes","No")
levels(data$pahclrnc_mv.factor)=c("Yes","No")
levels(data$u0igg_is_below_limit.factor)=c("Yes","No")
levels(data$u0igg_is_below_limit_mv.factor)=c("Yes","No")
levels(data$u0igg_mv.factor)=c("Yes","No")
levels(data$u3flow_gfr_mv.factor)=c("Yes","No")
levels(data$u3gfr_mv.factor)=c("Yes","No")
levels(data$u3use_gfr_mv.factor)=c("Yes","No")
levels(data$u3albumin_is_below_limit.factor)=c("Yes","No")
levels(data$u3albumin_is_below_limit_mv.factor)=c("Yes","No")
levels(data$u3albumin_mv.factor)=c("Yes","No")
levels(data$u3igg_is_below_limit.factor)=c("Yes","No")
levels(data$u3igg_is_below_limit_mv.factor)=c("Yes","No")
levels(data$u3igg_mv.factor)=c("Yes","No")
levels(data$ualb_is_below_limit.factor)=c("Yes","No")
levels(data$ualb_is_below_limit_mv.factor)=c("Yes","No")
levels(data$ualb_mv.factor)=c("Yes","No")
levels(data$ucr_mv.factor)=c("Yes","No")
levels(data$phoenix_labs_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$cantpv_2_mv.factor)=c("Yes","No")
levels(data$cantpv_3.factor)=c("Yes","No")
levels(data$cantpv_3_mv.factor)=c("Yes","No")
levels(data$cantpv_4.factor)=c("Yes","No")
levels(data$cantpv_4_mv.factor)=c("Yes","No")
levels(data$cantpv_5.factor)=c("Yes","No")
levels(data$cantpv_5_mv.factor)=c("Yes","No")
levels(data$cantpv_6.factor)=c("Yes","No")
levels(data$cantpv_6_mv.factor)=c("Yes","No")
levels(data$cantpv_7.factor)=c("Yes","No")
levels(data$cantpv_7_mv.factor)=c("Yes","No")
levels(data$cantpv_8.factor)=c("Yes","No")
levels(data$cantpv_8_mv.factor)=c("Yes","No")
levels(data$cantpv_9_mv.factor)=c("Yes","No")
levels(data$cantpv_10.factor)=c("Yes","No")
levels(data$cantpv_10_mv.factor)=c("Yes","No")
levels(data$cantpv_12_mv.factor)=c("Yes","No")
levels(data$cantpv_13.factor)=c("Yes","No")
levels(data$cantpv_13_mv.factor)=c("Yes","No")
levels(data$cantpv_14_mv.factor)=c("Yes","No")
levels(data$cantpv_15.factor)=c("Yes","No")
levels(data$cantpv_15_mv.factor)=c("Yes","No")
levels(data$cantpv_16_mv.factor)=c("Yes","No")
levels(data$cantpv_17.factor)=c("Yes","No")
levels(data$cantpv_17_mv.factor)=c("Yes","No")
levels(data$cantpv_18.factor)=c("Yes","No")
levels(data$cantpv_18_mv.factor)=c("Yes","No")
levels(data$cantpv_19_mv.factor)=c("Yes","No")
levels(data$cantpv_20_mv.factor)=c("Yes","No")
levels(data$cantpv_21_mv.factor)=c("Yes","No")
levels(data$cantpv_22_mv.factor)=c("Yes","No")
levels(data$cantpv_23_mv.factor)=c("Yes","No")
levels(data$cantpv_24_mv.factor)=c("Yes","No")
levels(data$cantpv_25_mv.factor)=c("Yes","No")
levels(data$cantpv_26_mv.factor)=c("Yes","No")
levels(data$cantpv_27_mv.factor)=c("Yes","No")
levels(data$cantpv_28_mv.factor)=c("Yes","No")
levels(data$cantpv_29_mv.factor)=c("Yes","No")
levels(data$cardiac_autonomic_neuropathy_test_preparation_veri_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$mnsi_cl_2_mv.factor)=c("Yes","No")
levels(data$mnsi_cl_3.factor)=c("Yes","No")
levels(data$mnsi_cl_3_mv.factor)=c("Yes","No")
levels(data$mnsi_cl_4___1.factor)=c("Unchecked","Checked")
levels(data$mnsi_cl_4___2.factor)=c("Unchecked","Checked")
levels(data$mnsi_cl_4___3.factor)=c("Unchecked","Checked")
levels(data$mnsi_cl_4___4.factor)=c("Unchecked","Checked")
levels(data$mnsi_cl_4___0.factor)=c("Unchecked","Checked")
levels(data$mnsi_cl_5_mv.factor)=c("Yes","No")
levels(data$mnsi_cl_6.factor)=c("Absent","Present")
levels(data$mnsi_cl_6_mv.factor)=c("Yes","No")
levels(data$mnsi_cl_7.factor)=c("Present","Present/Reinforcement","Absent")
levels(data$mnsi_cl_7_mv.factor)=c("Yes","No")
levels(data$mnsi_cl_8.factor)=c("Present","Decreased","Absent")
levels(data$mnsi_cl_8_mv.factor)=c("Yes","No")
levels(data$mnsi_cl_9.factor)=c("Normal","Reduced","Absent")
levels(data$mnsi_cl_9_mv.factor)=c("Yes","No")
levels(data$mnsi_cl_10.factor)=c("Yes","No")
levels(data$mnsi_cl_10_mv.factor)=c("Yes","No")
levels(data$mnsi_cl_11___1.factor)=c("Unchecked","Checked")
levels(data$mnsi_cl_11___2.factor)=c("Unchecked","Checked")
levels(data$mnsi_cl_11___3.factor)=c("Unchecked","Checked")
levels(data$mnsi_cl_11___4.factor)=c("Unchecked","Checked")
levels(data$mnsi_cl_11___0.factor)=c("Unchecked","Checked")
levels(data$mnsi_cl_12_mv.factor)=c("Yes","No")
levels(data$mnsi_cl_13.factor)=c("Absent","Present")
levels(data$mnsi_cl_13_mv.factor)=c("Yes","No")
levels(data$mnsi_cl_14.factor)=c("Present","Present/Reinforcement","Absent")
levels(data$mnsi_cl_14_mv.factor)=c("Yes","No")
levels(data$mnsi_cl_15.factor)=c("Present","Decreased","Absent")
levels(data$mnsi_cl_15_mv.factor)=c("Yes","No")
levels(data$mnsi_cl_16.factor)=c("Normal","Reduced","Absent")
levels(data$mnsi_cl_16_mv.factor)=c("Yes","No")
levels(data$mnsi_cl_17_mv.factor)=c("Yes","No")
levels(data$michigan_neuropathy_screening_instrument_for_clini_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$mnsi_subj_2_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_3.factor)=c("Yes","No")
levels(data$mnsi_subj_3_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_4.factor)=c("Yes","No")
levels(data$mnsi_subj_4_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_5.factor)=c("Yes","No")
levels(data$mnsi_subj_5_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_6.factor)=c("Yes","No")
levels(data$mnsi_subj_6_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_7.factor)=c("Yes","No")
levels(data$mnsi_subj_7_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_8.factor)=c("Yes","No")
levels(data$mnsi_subj_8_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_9.factor)=c("Yes","No")
levels(data$mnsi_subj_9_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_10.factor)=c("Yes","No")
levels(data$mnsi_subj_10_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_11.factor)=c("Yes","No")
levels(data$mnsi_subj_11_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_12.factor)=c("Yes","No")
levels(data$mnsi_subj_12_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_13.factor)=c("Yes","No")
levels(data$mnsi_subj_13_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_14.factor)=c("Yes","No")
levels(data$mnsi_subj_14_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_15.factor)=c("Yes","No")
levels(data$mnsi_subj_15_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_16.factor)=c("Yes","No")
levels(data$mnsi_subj_16_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_17.factor)=c("Yes","No")
levels(data$mnsi_subj_17_mv.factor)=c("Yes","No")
levels(data$mnsi_subj_18_mv.factor)=c("Yes","No")
levels(data$michigan_neuropathy_screening_instrument_for_subje_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$aspq_1_mv.factor)=c("Yes","No")
levels(data$aspq_2_mv.factor)=c("Yes","No")
levels(data$aspq_3_mv.factor)=c("Yes","No")
levels(data$aspq_4.factor)=c("Yes","No")
levels(data$aspq_4_mv.factor)=c("Yes","No")
levels(data$aspq_5.factor)=c("Rarely","Occasionally","Frequently","Almost always")
levels(data$aspq_5_mv.factor)=c("Yes","No")
levels(data$aspq_6.factor)=c("Mild","Moderate","Severe")
levels(data$aspq_6_mv.factor)=c("Yes","No")
levels(data$aspq_7.factor)=c("Less than 3 months","3 to 6 months","7 to 12 months","13 months to 5 years","More than 5 years","As long as I can remember")
levels(data$aspq_7_mv.factor)=c("Yes","No")
levels(data$aspq_8.factor)=c("Never","Once","Twice","Three times","Four times","Five or more times")
levels(data$aspq_8_mv.factor)=c("Yes","No")
levels(data$aspq_9.factor)=c("Not cautious at all","Somewhat cautious","Extremely cautious")
levels(data$aspq_9_mv.factor)=c("Yes","No")
levels(data$aspq_10.factor)=c("Early morning","Rest of morning","Afternoon","Evening","At night, when I get up after Ive been asleep","No particular time is worst","Other time")
levels(data$aspq_10_mv.factor)=c("Yes","No")
levels(data$aspq_11_mv.factor)=c("Yes","No")
levels(data$aspq_12.factor)=c("Gotten much worse","Gotten somewhat worse","Stayed about the same","Gotten somewhat better","Gotten much better","Completely gone")
levels(data$aspq_12_mv.factor)=c("Yes","No")
levels(data$aspq_14.factor)=c("Have not had","Mild","Moderate","Severe")
levels(data$aspq_14_mv.factor)=c("Yes","No")
levels(data$aspq_15.factor)=c("Have not had","Mild","Moderate","Severe")
levels(data$aspq_15_mv.factor)=c("Yes","No")
levels(data$aspq_16.factor)=c("Have not had","Mild","Moderate","Severe")
levels(data$aspq_16_mv.factor)=c("Yes","No")
levels(data$aspq_17.factor)=c("Have not had","Mild","Moderate","Severe")
levels(data$aspq_17_mv.factor)=c("Yes","No")
levels(data$aspq_18.factor)=c("Have not had","Mild","Moderate","Severe")
levels(data$aspq_18_mv.factor)=c("Yes","No")
levels(data$aspq_19.factor)=c("Have not had","Mild","Moderate","Severe")
levels(data$aspq_19_mv.factor)=c("Yes","No")
levels(data$aspq_20.factor)=c("Have not had","Mild","Moderate","Severe")
levels(data$aspq_20_mv.factor)=c("Yes","No")
levels(data$aspq_21.factor)=c("Have not had","Mild","Moderate","Severe")
levels(data$aspq_21_mv.factor)=c("Yes","No")
levels(data$aspq_22.factor)=c("Have not had","Mild","Moderate","Severe")
levels(data$aspq_22_mv.factor)=c("Yes","No")
levels(data$aspq_23.factor)=c("Have not had","Mild","Moderate","Severe")
levels(data$aspq_23_mv.factor)=c("Yes","No")
levels(data$aspq_24.factor)=c("Yes","No")
levels(data$aspq_24_mv.factor)=c("Yes","No")
levels(data$aspq_25_mv.factor)=c("Yes","No")
levels(data$aspq_26.factor)=c("Yes","No")
levels(data$aspq_26_mv.factor)=c("Yes","No")
levels(data$aspq_27.factor)=c("Yes","No")
levels(data$aspq_27_mv.factor)=c("Yes","No")
levels(data$aspq_28.factor)=c("Yes","No")
levels(data$aspq_28_mv.factor)=c("Yes","No")
levels(data$aspq_29.factor)=c("Yes","No")
levels(data$aspq_29_mv.factor)=c("Yes","No")
levels(data$aspq_30.factor)=c("Yes","No")
levels(data$aspq_30_mv.factor)=c("Yes","No")
levels(data$aspq_31.factor)=c("Yes","No")
levels(data$aspq_31_mv.factor)=c("Yes","No")
levels(data$aspq_32.factor)=c("Yes","No")
levels(data$aspq_32_mv.factor)=c("Yes","No")
levels(data$aspq_34.factor)=c("Yes","No")
levels(data$aspq_34_mv.factor)=c("Yes","No")
levels(data$aspq_35.factor)=c("Yes","No")
levels(data$aspq_35_mv.factor)=c("Yes","No")
levels(data$aspq_36.factor)=c("Yes","No")
levels(data$aspq_36_mv.factor)=c("Yes","No")
levels(data$aspq_37_mv.factor)=c("Yes","No")
levels(data$aspq_38.factor)=c("Yes","No")
levels(data$aspq_38_mv.factor)=c("Yes","No")
levels(data$aspq_39.factor)=c("Yes","No")
levels(data$aspq_39_mv.factor)=c("Yes","No")
levels(data$aspq_40_mv.factor)=c("Yes","No")
levels(data$autonomic_symptoms_profile_questionnaire_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$nsqlq_2_mv.factor)=c("Yes","No")
levels(data$nsqlq_3_mv.factor)=c("Yes","No")
levels(data$nsqlq_5.factor)=c("All the time","Most of the time","Some of the time","Occasionally","Never")
levels(data$nsqlq_5_mv.factor)=c("Yes","No")
levels(data$nsqlq_6.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_6_mv.factor)=c("Yes","No")
levels(data$nsqlq_7.factor)=c("All the time","Most of the time","Some of the time","Occasionally","Never")
levels(data$nsqlq_7_mv.factor)=c("Yes","No")
levels(data$nsqlq_8.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_8_mv.factor)=c("Yes","No")
levels(data$nsqlq_9.factor)=c("All the time","Most of the time","Some of the time","Occasionally","Never")
levels(data$nsqlq_9_mv.factor)=c("Yes","No")
levels(data$nsqlq_10.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_10_mv.factor)=c("Yes","No")
levels(data$nsqlq_11.factor)=c("All the time","Most of the time","Some of the time","Occasionally","Never")
levels(data$nsqlq_11_mv.factor)=c("Yes","No")
levels(data$nsqlq_12.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_12_mv.factor)=c("Yes","No")
levels(data$nsqlq_13.factor)=c("All the time","Most of the time","Some of the time","Occasionally","Never")
levels(data$nsqlq_13_mv.factor)=c("Yes","No")
levels(data$nsqlq_14.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_14_mv.factor)=c("Yes","No")
levels(data$nsqlq_15.factor)=c("All the time","Most of the time","Some of the time","Occasionally","Never")
levels(data$nsqlq_15_mv.factor)=c("Yes","No")
levels(data$nsqlq_16.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_16_mv.factor)=c("Yes","No")
levels(data$nsqlq_17.factor)=c("All the time","Most of the time","Some of the time","Occasionally","Never")
levels(data$nsqlq_17_mv.factor)=c("Yes","No")
levels(data$nsqlq_18.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_18_mv.factor)=c("Yes","No")
levels(data$nsqlq_19.factor)=c("Very much","Quite a lot","Somewhat","A little","Not at all")
levels(data$nsqlq_19_mv.factor)=c("Yes","No")
levels(data$nsqlq_20.factor)=c("All the time","Most of the time","Some of the time","Occasionally","Never")
levels(data$nsqlq_20_mv.factor)=c("Yes","No")
levels(data$nsqlq_21.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_21_mv.factor)=c("Yes","No")
levels(data$nsqlq_22.factor)=c("All the time","Most of the time","Some of the time","Occasionally","Never")
levels(data$nsqlq_22_mv.factor)=c("Yes","No")
levels(data$nsqlq_23.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_23_mv.factor)=c("Yes","No")
levels(data$nsqlq_24.factor)=c("All the time","Most of the time","Some of the time","Occasionally","Never")
levels(data$nsqlq_24_mv.factor)=c("Yes","No")
levels(data$nsqlq_25.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_25_mv.factor)=c("Yes","No")
levels(data$nsqlq_26.factor)=c("Very much","Quite a lot","Somewhat","A little","Not at all")
levels(data$nsqlq_26_mv.factor)=c("Yes","No")
levels(data$nsqlq_27.factor)=c("All the time","Most of the time","Some of the time","Occasionally","Never")
levels(data$nsqlq_27_mv.factor)=c("Yes","No")
levels(data$nsqlq_28.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_28_mv.factor)=c("Yes","No")
levels(data$nsqlq_29.factor)=c("All the time","Most of the time","Some of the time","Occasionally","Never")
levels(data$nsqlq_29_mv.factor)=c("Yes","No")
levels(data$nsqlq_30.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_30_mv.factor)=c("Yes","No")
levels(data$nsqlq_31.factor)=c("All the time","Most of the time","Some of the time","Occasionally","Never")
levels(data$nsqlq_31_mv.factor)=c("Yes","No")
levels(data$nsqlq_32.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_32_mv.factor)=c("Yes","No")
levels(data$nsqlq_33.factor)=c("Very much","Quite a lot","Somewhat","A little","Not at all")
levels(data$nsqlq_33_mv.factor)=c("Yes","No")
levels(data$nsqlq_34.factor)=c("Yes","No")
levels(data$nsqlq_34_mv.factor)=c("Yes","No")
levels(data$nsqlq_35.factor)=c("Very much","Quite a lot","Somewhat","A little","Not at all")
levels(data$nsqlq_35_mv.factor)=c("Yes","No")
levels(data$nsqlq_36.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_36_mv.factor)=c("Yes","No")
levels(data$nsqlq_37.factor)=c("Very much","Quite a lot","Somewhat","A little","Not at all")
levels(data$nsqlq_37_mv.factor)=c("Yes","No")
levels(data$nsqlq_38.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_38_mv.factor)=c("Yes","No")
levels(data$nsqlq_39.factor)=c("Very much","Quite a lot","Somewhat","A little","Not at all")
levels(data$nsqlq_39_mv.factor)=c("Yes","No")
levels(data$nsqlq_40.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_40_mv.factor)=c("Yes","No")
levels(data$nsqlq_41.factor)=c("Very much","Quite a lot","Somewhat","A little","Not at all")
levels(data$nsqlq_41_mv.factor)=c("Yes","No")
levels(data$nsqlq_42.factor)=c("Very much","Quite a lot","Somewhat","A little","Not at all")
levels(data$nsqlq_42_mv.factor)=c("Yes","No")
levels(data$nsqlq_43.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_43_mv.factor)=c("Yes","No")
levels(data$nsqlq_44.factor)=c("Very much","Quite a lot","Somewhat","A little","Not at all")
levels(data$nsqlq_44_mv.factor)=c("Yes","No")
levels(data$nsqlq_45.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_45_mv.factor)=c("Yes","No")
levels(data$nsqlq_46.factor)=c("Very much","Quite a lot","Somewhat","A little","Not at all")
levels(data$nsqlq_46_mv.factor)=c("Yes","No")
levels(data$nsqlq_47.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_47_mv.factor)=c("Yes","No")
levels(data$nsqlq_48.factor)=c("Very much","Quite a lot","Somewhat","A little","Not at all")
levels(data$nsqlq_48_mv.factor)=c("Yes","No")
levels(data$nsqlq_49.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_49_mv.factor)=c("Yes","No")
levels(data$nsqlq_50.factor)=c("Very much","Quite a lot","Somewhat","A little","Not at all")
levels(data$nsqlq_50_mv.factor)=c("Yes","No")
levels(data$nsqlq_51.factor)=c("Completely agree","Partly agree","Neither agree nor disagree","Partly disagree","Completely disagree")
levels(data$nsqlq_51_mv.factor)=c("Yes","No")
levels(data$nsqlq_52.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_52_mv.factor)=c("Yes","No")
levels(data$nsqlq_53.factor)=c("Completely agree","Partly agree","Neither agree nor disagree","Partly disagree","Completely disagree")
levels(data$nsqlq_53_mv.factor)=c("Yes","No")
levels(data$nsqlq_54.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_54_mv.factor)=c("Yes","No")
levels(data$nsqlq_55.factor)=c("Completely agree","Partly agree","Neither agree nor disagree","Partly disagree","Completely disagree")
levels(data$nsqlq_55_mv.factor)=c("Yes","No")
levels(data$nsqlq_56.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_56_mv.factor)=c("Yes","No")
levels(data$nsqlq_57.factor)=c("Completely agree","Partly agree","Neither agree nor disagree","Partly disagree","Completely disagree")
levels(data$nsqlq_57_mv.factor)=c("Yes","No")
levels(data$nsqlq_58.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_58_mv.factor)=c("Yes","No")
levels(data$nsqlq_59.factor)=c("Completely agree","Partly agree","Neither agree nor disagree","Partly disagree","Completely disagree")
levels(data$nsqlq_59_mv.factor)=c("Yes","No")
levels(data$nsqlq_60.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_60_mv.factor)=c("Yes","No")
levels(data$nsqlq_61.factor)=c("Completely agree","Partly agree","Neither agree nor disagree","Partly disagree","Completely disagree")
levels(data$nsqlq_61_mv.factor)=c("Yes","No")
levels(data$nsqlq_62.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_62_mv.factor)=c("Yes","No")
levels(data$nsqlq_63.factor)=c("Very depressed","Quite depressed","Somewhat depressed","A little depressed","Not at all depressed")
levels(data$nsqlq_63_mv.factor)=c("Yes","No")
levels(data$nsqlq_64.factor)=c("Very much","Some bother","None")
levels(data$nsqlq_64_mv.factor)=c("Yes","No")
levels(data$nsqlq_65.factor)=c("Very much","Quite a lot","Somewhat","A little","Not at all")
levels(data$nsqlq_65_mv.factor)=c("Yes","No")
levels(data$nsqlq_66.factor)=c("Very much","Quite a lot","Somewhat","A little","Not at all")
levels(data$nsqlq_66_mv.factor)=c("Yes","No")
levels(data$nsqlq_67.factor)=c("Excellent","Very good","Good","Fair","Poor")
levels(data$nsqlq_67_mv.factor)=c("Yes","No")
levels(data$neuropathy_specific_quality_of_life_questionnaire_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$jvn_collecdate_mv.factor)=c("Yes","No")
levels(data$jvn_date_mv.factor)=c("Yes","No")
levels(data$jvn_npdr_od.factor)=c("No Evidence","Very Mild","Mild","Moderate","Severe","Very Severe","Ungradeable")
levels(data$jvn_npdr_od_mv.factor)=c("Yes","No")
levels(data$jvn_npdr_os.factor)=c("No Evidence","Very Mild","Mild","Moderate","Severe","Very Severe","Ungradeable")
levels(data$jvn_npdr_os_mv.factor)=c("Yes","No")
levels(data$jvn_npdr_notes_mv.factor)=c("Yes","No")
levels(data$jvn_pdr_od.factor)=c("No Evidence","Less than High Risk","High Risk","Quiescent","Ungradeable")
levels(data$jvn_pdr_od_mv.factor)=c("Yes","No")
levels(data$jvn_pdr_os.factor)=c("No Evidence","Less than High Risk","High Risk","Quiescent","Ungradeable")
levels(data$jvn_pdr_os_mv.factor)=c("Yes","No")
levels(data$jvn_pdr_notes_mv.factor)=c("Yes","No")
levels(data$jvn_me_od.factor)=c("No Evidence","Macular Edema, Not Clinically significant","Macular Edema, Clinically significant","Ungradeable")
levels(data$jvn_me_od_mv.factor)=c("Yes","No")
levels(data$jvn_me_os.factor)=c("No Evidence","Macular Edema, Not Clinically significant","Macular Edema, Clinically significant","Ungradeable")
levels(data$jvn_me_os_mv.factor)=c("Yes","No")
levels(data$jvn_me_notes_mv.factor)=c("Yes","No")
levels(data$jvn_amd_od.factor)=c("Yes","No")
levels(data$jvn_amd_od_mv.factor)=c("Yes","No")
levels(data$jvn_amd_os.factor)=c("Yes","No")
levels(data$jvn_amd_os_mv.factor)=c("Yes","No")
levels(data$jvn_asterhy_od.factor)=c("Yes","No")
levels(data$jvn_asterhy_od_mv.factor)=c("Yes","No")
levels(data$jvn_asterhy_os.factor)=c("Yes","No")
levels(data$jvn_asterhy_os_mv.factor)=c("Yes","No")
levels(data$jvn_cat_od.factor)=c("Yes","No")
levels(data$jvn_cat_od_mv.factor)=c("Yes","No")
levels(data$jvn_cat_os.factor)=c("Yes","No")
levels(data$jvn_cat_os_mv.factor)=c("Yes","No")
levels(data$jvn_chorles_od.factor)=c("Yes","No")
levels(data$jvn_chorles_od_mv.factor)=c("Yes","No")
levels(data$jvn_chorles_os.factor)=c("Yes","No")
levels(data$jvn_chorles_os_mv.factor)=c("Yes","No")
levels(data$jvn_arcus_od.factor)=c("Yes","No")
levels(data$jvn_arcus_od_mv.factor)=c("Yes","No")
levels(data$jvn_arcus_os.factor)=c("Yes","No")
levels(data$jvn_arcus_os_mv.factor)=c("Yes","No")
levels(data$jvn_cratscar_od.factor)=c("Yes","No")
levels(data$jvn_cratscar_od_mv.factor)=c("Yes","No")
levels(data$jvn_cratscar_os.factor)=c("Yes","No")
levels(data$jvn_cratscar_os_mv.factor)=c("Yes","No")
levels(data$jvn_epiret_od.factor)=c("Yes","No")
levels(data$jvn_epiret_od_mv.factor)=c("Yes","No")
levels(data$jvn_epiret_os.factor)=c("Yes","No")
levels(data$jvn_epiret_os_mv.factor)=c("Yes","No")
levels(data$jvn_iol_od.factor)=c("Yes","No")
levels(data$jvn_iol_od_mv.factor)=c("Yes","No")
levels(data$jvn_iol_os.factor)=c("Yes","No")
levels(data$jvn_iol_os_mv.factor)=c("Yes","No")
levels(data$jvn_nevus_od.factor)=c("Yes","No")
levels(data$jvn_nevus_od_mv.factor)=c("Yes","No")
levels(data$jvn_nevus_os.factor)=c("Yes","No")
levels(data$jvn_nevus_os_mv.factor)=c("Yes","No")
levels(data$jvn_hardex_od.factor)=c("Yes","No")
levels(data$jvn_hardex_od_mv.factor)=c("Yes","No")
levels(data$jvn_hardex_os.factor)=c("Yes","No")
levels(data$jvn_hardex_os_mv.factor)=c("Yes","No")
levels(data$jvn_lidles_od.factor)=c("Yes","No")
levels(data$jvn_lidles_od_mv.factor)=c("Yes","No")
levels(data$jvn_lidles_os.factor)=c("Yes","No")
levels(data$jvn_lidles_os_mv.factor)=c("Yes","No")
levels(data$jvn_drusen_od.factor)=c("Yes","No")
levels(data$jvn_drusen_od_mv.factor)=c("Yes","No")
levels(data$jvn_drusen_os.factor)=c("Yes","No")
levels(data$jvn_drusen_os_mv.factor)=c("Yes","No")
levels(data$jvn_ondrus_od.factor)=c("Yes","No")
levels(data$jvn_ondrus_od_mv.factor)=c("Yes","No")
levels(data$jvn_ondrus_os.factor)=c("Yes","No")
levels(data$jvn_ondrus_os_mv.factor)=c("Yes","No")
levels(data$jvn_miscas_od.factor)=c("Yes","No")
levels(data$jvn_miscas_od_mv.factor)=c("Yes","No")
levels(data$jvn_miscas_os.factor)=c("Yes","No")
levels(data$jvn_miscas_os_mv.factor)=c("Yes","No")
levels(data$jvn_miscvr_od.factor)=c("Yes","No")
levels(data$jvn_miscvr_od_mv.factor)=c("Yes","No")
levels(data$jvn_miscvr_os.factor)=c("Yes","No")
levels(data$jvn_miscvr_os_mv.factor)=c("Yes","No")
levels(data$jvn_ppa_od.factor)=c("Yes","No")
levels(data$jvn_ppa_od_mv.factor)=c("Yes","No")
levels(data$jvn_ppa_os.factor)=c("Yes","No")
levels(data$jvn_ppa_os_mv.factor)=c("Yes","No")
levels(data$jvn_prfib_od.factor)=c("Yes","No")
levels(data$jvn_prfib_od_mv.factor)=c("Yes","No")
levels(data$jvn_prfib_os.factor)=c("Yes","No")
levels(data$jvn_prfib_os_mv.factor)=c("Yes","No")
levels(data$jvn_htnret_od.factor)=c("Yes","No")
levels(data$jvn_htnret_od_mv.factor)=c("Yes","No")
levels(data$jvn_htnret_os.factor)=c("Yes","No")
levels(data$jvn_htnret_os_mv.factor)=c("Yes","No")
levels(data$jvp_rpe_od.factor)=c("Yes","No")
levels(data$jvp_rpe_od_mv.factor)=c("Yes","No")
levels(data$jvn_rpe_os.factor)=c("Yes","No")
levels(data$jvn_rpe_os_mv.factor)=c("Yes","No")
levels(data$jvn_rrd_od.factor)=c("Yes","No")
levels(data$jvn_rrd_od_mv.factor)=c("Yes","No")
levels(data$jvn_rrd_os.factor)=c("Yes","No")
levels(data$jvn_rrd_os_mv.factor)=c("Yes","No")
levels(data$jvn_trd_od.factor)=c("Yes","No")
levels(data$jvn_trd_od_mv.factor)=c("Yes","No")
levels(data$jvn_trd_os.factor)=c("Yes","No")
levels(data$jvn_trd_os_mv.factor)=c("Yes","No")
levels(data$jvn_vithem_od.factor)=c("Yes","No")
levels(data$jvn_vithem_od_mv.factor)=c("Yes","No")
levels(data$jvn_vithem_os.factor)=c("Yes","No")
levels(data$jvn_vithem_os_mv.factor)=c("Yes","No")
levels(data$jvn_other_od.factor)=c("Yes","No")
levels(data$jvn_other_od_mv.factor)=c("Yes","No")
levels(data$jvn_other_os.factor)=c("Yes","No")
levels(data$jvn_other_os_mv.factor)=c("Yes","No")
levels(data$jvn_laser_od.factor)=c("Yes","No")
levels(data$jvn_laser_od_mv.factor)=c("Yes","No")
levels(data$jvn_laser_os.factor)=c("Yes","No")
levels(data$jvn_laser_os_mv.factor)=c("Yes","No")
levels(data$jvn_notes_mv.factor)=c("Yes","No")
levels(data$jvn_deck_8_mv.factor)=c("Yes","No")
levels(data$jvn_retinal_imaging_form_complete.factor)=c("Incomplete","Unverified","Complete")

ddn = data
rm(data)