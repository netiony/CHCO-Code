#Load Hmisc library
library(Hmisc)
#Read Data
data=read.csv('E:/Petter Bjornstad/Pima/Master data/Raw data/NelsonPECRBRenoprote_DATA_2021-02-04_1619.csv',na.strings = "")
#Setting Labels

label(data$record_id)="Patient ID"
label(data$redcap_event_name)="Event Name"
label(data$redcap_repeat_instrument)="Repeat Instrument"
label(data$redcap_repeat_instance)="Repeat Instance"
label(data$lastname)="Lastname"
label(data$lastname_mv)="Missing value verified for: Lastname"
label(data$firstname)="Firstname"
label(data$firstname_mv)="Missing value verified for: Firstname"
label(data$stratum)="Patient stratum"
label(data$stratum_mv)="Missing value verified for: Patient stratum"
label(data$randomnum)="Patient radomization number"
label(data$randomnum_mv)="Missing value verified for: Patient radomization number"
label(data$tx)="TX"
label(data$tx_mv)="Missing value verified for: TX"
label(data$sacaton_no)="Sacaton number"
label(data$sacaton_no_mv)="Missing value verified for: Sacaton number"
label(data$phx_number)="PIMC number"
label(data$phx_number_mv)="Missing value verified for: PIMC number"
label(data$dob)="DOB"
label(data$dob_mv)="Missing value verified for: DOB"
label(data$gender)="Gender"
label(data$gender_mv)="Missing value verified for: Gender"
label(data$height)="Height (cm)"
label(data$height_mv)="Missing value verified for: Height (cm)"
label(data$allergies)="Allergies?"
label(data$allergies_mv)="Missing value verified for: Allergies?"
label(data$allerglist)="Allergies (specify)"
label(data$allerglist_mv)="Missing value verified for: Allergies (specify)"
label(data$district)="District"
label(data$district_mv)="Missing value verified for: District"
label(data$phonenum)="Phone number"
label(data$phonenum_mv)="Missing value verified for: Phone number"
label(data$bldprsr)="Blood pressure arm"
label(data$bldprsr_mv)="Missing value verified for: Blood pressure arm"
label(data$dtdiabonset)="Date of Diabetes Onset?"
label(data$dtdiabonset_mv)="Missing value verified for: Date of Diabetes Onset?"
label(data$closedate)="Close date"
label(data$closedate_mv)="Missing value verified for: Close date"
label(data$washoutdate)="Washout date"
label(data$washoutdate_mv)="Missing value verified for: Washout date"
label(data$demoform_complete)="Complete?"
label(data$f2visitdate)="Date of Visit"
label(data$f2visitdate_mv)="Missing value verified for: Date of Visit"
label(data$f2visittype)="Type of Visit"
label(data$f2visittype_mv)="Missing value verified for: Type of Visit"
label(data$f2testintvrl)="Test Interval"
label(data$f2testintvrl_mv)="Missing value verified for: Test Interval"
label(data$f2bldprsrarm)="Blood pressure arm"
label(data$f2bldprsrarm_mv)="Missing value verified for: Blood pressure arm"
label(data$f2frstmsrsys)="1st measure - Systolic"
label(data$f2frstmsrsys_mv)="Missing value verified for: 1st measure - Systolic"
label(data$f2frstmsrdia)="1st measure - Diastolic"
label(data$f2frstmsrdia_mv)="Missing value verified for: 1st measure - Diastolic"
label(data$f2scndmsrsys)="2nd measure - Systolic"
label(data$f2scndmsrsys_mv)="Missing value verified for: 2nd measure - Systolic"
label(data$f2scndmsrdia)="2nd measure - Diastolic"
label(data$f2scndmsrdia_mv)="Missing value verified for: 2nd measure - Diastolic"
label(data$f2heartrate)="Heart rate (beats/min)"
label(data$f2heartrate_mv)="Missing value verified for: Heart rate (beats/min)"
label(data$subject_screening_3_and_2_form_02_complete)="Complete?"
label(data$f3visitdate)="Date of Visit"
label(data$f3visitdate_mv)="Missing value verified for: Date of Visit"
label(data$f3visittype)="Type of Visit"
label(data$f3visittype_mv)="Missing value verified for: Type of Visit"
label(data$f3testintvrl)="Test Interval"
label(data$f3testintvrl_mv)="Missing value verified for: Test Interval"
label(data$f3spturncltd)="Urine collected (spot) for albumin/creatinine?"
label(data$f3spturncltd_mv)="Missing value verified for: Urine collected (spot) for albumin/creatinine?"
label(data$f3venipfork)="Venipuncture for K+ and creatinine?"
label(data$f3venipfork_mv)="Missing value verified for: Venipuncture for K+ and creatinine?"
label(data$f3bldprsrarm)="Blood pressure arm"
label(data$f3bldprsrarm_mv)="Missing value verified for: Blood pressure arm"
label(data$f3frstmsrsys)="1st measure - Systolic"
label(data$f3frstmsrsys_mv)="Missing value verified for: 1st measure - Systolic"
label(data$f3frstmsrdia)="1st measure - Diastolic"
label(data$f3frstmsrdia_mv)="Missing value verified for: 1st measure - Diastolic"
label(data$f3scndmsrsys)="2nd measure - Systolic"
label(data$f3scndmsrsys_mv)="Missing value verified for: 2nd measure - Systolic"
label(data$f3scndmsrdia)="2nd measure - Diastolic"
label(data$f3scndmsrdia_mv)="Missing value verified for: 2nd measure - Diastolic"
label(data$f3heartrate)="Heart rate (beats/min)"
label(data$f3heartrate_mv)="Missing value verified for: Heart rate (beats/min)"
label(data$f3prohbmeds)="Any prohibiting study medicines?"
label(data$f3prohbmeds_mv)="Missing value verified for: Any prohibiting study medicines?"
label(data$f3elgstatus)="Elgibility status"
label(data$f3elgstatus_mv)="Missing value verified for: Elgibility status"
label(data$f3inelgibleexplain)="If inelgible, explain"
label(data$f3inelgibleexplain_mv)="Missing value verified for: If inelgible, explain"
label(data$f3medsurghis)="Medical/surgical history"
label(data$f3medsurghis_mv)="Missing value verified for: Medical/surgical history"
label(data$f3_acratio)="A/C ratio"
label(data$f3_acratio_mv)="Missing value verified for: A/C ratio"
label(data$f3srmcreatyn)="Serum creatinine"
label(data$f3srmcreatyn_mv)="Missing value verified for: Serum creatinine"
label(data$f3_pregnancy)="Pregnancy"
label(data$f3_pregnancy_mv)="Missing value verified for: Pregnancy"
label(data$f3_staffpref)="Staff preference"
label(data$f3_staffpref_mv)="Missing value verified for: Staff preference"
label(data$f3_subjpref)="Subject preference"
label(data$f3_subjpref_mv)="Missing value verified for: Subject preference"
label(data$f3_other)="Other"
label(data$f3_other_mv)="Missing value verified for: Other"
label(data$subject_screeningelgibility_form_03_complete)="Complete?"
label(data$f3meds_visitdate)="Visit Date"
label(data$f3meds_visitdate_mv)="Missing value verified for: Visit Date"
label(data$f3drugcode1)="Drug Code"
label(data$f3drugcode1_mv)="Missing value verified for: Drug Code"
label(data$f3dosage1)="Dosage"
label(data$f3dosage1_mv)="Missing value verified for: Dosage"
label(data$f3timesday1)="Number of times taken"
label(data$f3timesday1_mv)="Missing value verified for: Number of times taken"
label(data$f3dwm1)="Times/day, week or month"
label(data$f3dwm1_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode2)="Drug Code"
label(data$f3drugcode2_mv)="Missing value verified for: Drug Code"
label(data$f3dosage2)="Dosage"
label(data$f3dosage2_mv)="Missing value verified for: Dosage"
label(data$f3timesday2)="Number of time taken"
label(data$f3timesday2_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm2)="Times/day, week or month"
label(data$f3dwm2_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode3)="Drug Code"
label(data$f3drugcode3_mv)="Missing value verified for: Drug Code"
label(data$f3dosage3)="Dosage"
label(data$f3dosage3_mv)="Missing value verified for: Dosage"
label(data$f3timesday3)="Number of time taken"
label(data$f3timesday3_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm3)="Times/day, week or month"
label(data$f3dwm3_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode4)="Drug Code"
label(data$f3drugcode4_mv)="Missing value verified for: Drug Code"
label(data$f3dosage4)="Dosage"
label(data$f3dosage4_mv)="Missing value verified for: Dosage"
label(data$f3timesday4)="Number of time taken"
label(data$f3timesday4_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm4)="Times/day, week or month"
label(data$f3dwm4_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode5)="Drug Code"
label(data$f3drugcode5_mv)="Missing value verified for: Drug Code"
label(data$f3dosage5)="Dosage"
label(data$f3dosage5_mv)="Missing value verified for: Dosage"
label(data$f3timesday5)="Number of time taken"
label(data$f3timesday5_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm5)="Times/day, week or month"
label(data$f3dwm5_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode6)="Drug Code"
label(data$f3drugcode6_mv)="Missing value verified for: Drug Code"
label(data$f3dosage6)="Dosage"
label(data$f3dosage6_mv)="Missing value verified for: Dosage"
label(data$f3timesday6)="Number of time taken"
label(data$f3timesday6_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm6)="Times/day, week or month"
label(data$f3dwm6_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode7)="Drug Code"
label(data$f3drugcode7_mv)="Missing value verified for: Drug Code"
label(data$f3dosage7)="Dosage"
label(data$f3dosage7_mv)="Missing value verified for: Dosage"
label(data$f3timesday7)="Number of time taken"
label(data$f3timesday7_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm7)="Times/day, week or month"
label(data$f3dwm7_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode8)="Drug Code"
label(data$f3drugcode8_mv)="Missing value verified for: Drug Code"
label(data$f3dosage8)="Dosage"
label(data$f3dosage8_mv)="Missing value verified for: Dosage"
label(data$f3timesday8)="Number of time taken"
label(data$f3timesday8_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm8)="Times/day, week or month"
label(data$f3dwm8_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode9)="Drug Code"
label(data$f3drugcode9_mv)="Missing value verified for: Drug Code"
label(data$f3dosage9)="Dosage"
label(data$f3dosage9_mv)="Missing value verified for: Dosage"
label(data$f3timesday9)="Number of time taken"
label(data$f3timesday9_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm9)="Times/day, week or month"
label(data$f3dwm9_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode10)="Drug Code"
label(data$f3drugcode10_mv)="Missing value verified for: Drug Code"
label(data$f3dosage10)="Dosage"
label(data$f3dosage10_mv)="Missing value verified for: Dosage"
label(data$f3timesday10)="Number of time taken"
label(data$f3timesday10_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm10)="Times/day, week or month"
label(data$f3dwm10_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode11)="Drug Code"
label(data$f3drugcode11_mv)="Missing value verified for: Drug Code"
label(data$f3dosage11)="Dosage"
label(data$f3dosage11_mv)="Missing value verified for: Dosage"
label(data$f3timesday11)="Number of time taken"
label(data$f3timesday11_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm11)="Times/day, week or month"
label(data$f3dwm11_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode12)="Drug Code"
label(data$f3drugcode12_mv)="Missing value verified for: Drug Code"
label(data$f3dosage12)="Dosage"
label(data$f3dosage12_mv)="Missing value verified for: Dosage"
label(data$f3timesday12)="Number of time taken"
label(data$f3timesday12_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm12)="Times/day, week or month"
label(data$f3dwm12_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode13)="Drug Code"
label(data$f3drugcode13_mv)="Missing value verified for: Drug Code"
label(data$f3dosage13)="Dosage"
label(data$f3dosage13_mv)="Missing value verified for: Dosage"
label(data$f3timesday13)="Number of time taken"
label(data$f3timesday13_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm13)="Times/day, week or month"
label(data$f3dwm13_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode14)="Drug Code"
label(data$f3drugcode14_mv)="Missing value verified for: Drug Code"
label(data$f3dosage14)="Dosage"
label(data$f3dosage14_mv)="Missing value verified for: Dosage"
label(data$f3timesday14)="Number of time taken"
label(data$f3timesday14_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm14)="Times/day, week or month"
label(data$f3dwm14_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode15)="Drug Code"
label(data$f3drugcode15_mv)="Missing value verified for: Drug Code"
label(data$f3dosage15)="Dosage"
label(data$f3dosage15_mv)="Missing value verified for: Dosage"
label(data$f3timesday15)="Number of time taken"
label(data$f3timesday15_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm15)="Times/day, week or month"
label(data$f3dwm15_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode16)="Drug Code"
label(data$f3drugcode16_mv)="Missing value verified for: Drug Code"
label(data$f3dosage16)="Dosage"
label(data$f3dosage16_mv)="Missing value verified for: Dosage"
label(data$f3timesday16)="Number of time taken"
label(data$f3timesday16_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm16)="Times/day, week or month"
label(data$f3dwm16_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode17)="Drug Code"
label(data$f3drugcode17_mv)="Missing value verified for: Drug Code"
label(data$f3dosage17)="Dosage"
label(data$f3dosage17_mv)="Missing value verified for: Dosage"
label(data$f3timesday17)="Number of time taken"
label(data$f3timesday17_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm17)="Times/day, week or month"
label(data$f3dwm17_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode18)="Drug Code"
label(data$f3drugcode18_mv)="Missing value verified for: Drug Code"
label(data$f3dosage18)="Dosage"
label(data$f3dosage18_mv)="Missing value verified for: Dosage"
label(data$f3timesday18)="Number of time taken"
label(data$f3timesday18_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm18)="Times/day, week or month"
label(data$f3dwm18_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode19)="Drug Code"
label(data$f3drugcode19_mv)="Missing value verified for: Drug Code"
label(data$f3dosage19)="Dosage"
label(data$f3dosage19_mv)="Missing value verified for: Dosage"
label(data$f3timesday19)="Number of time taken"
label(data$f3timesday19_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm19)="Times/day, week or month"
label(data$f3dwm19_mv)="Missing value verified for: Times/day, week or month"
label(data$f3drugcode20)="Drug Code"
label(data$f3drugcode20_mv)="Missing value verified for: Drug Code"
label(data$f3dosage20)="Dosage"
label(data$f3dosage20_mv)="Missing value verified for: Dosage"
label(data$f3timesday20)="Number of time taken"
label(data$f3timesday20_mv)="Missing value verified for: Number of time taken"
label(data$f3dwm20)="Times/day, week or month"
label(data$f3dwm20_mv)="Missing value verified for: Times/day, week or month"
label(data$subject_screeningelgibility_meds_form_03_complete)="Complete?"
label(data$clrncvisitdate)="Date of Visit"
label(data$clrncvisitdate_mv)="Missing value verified for: Date of Visit"
label(data$clrncvisittype)="Type of Visit"
label(data$clrncvisittype_mv)="Missing value verified for: Type of Visit"
label(data$clrnctestintvrl)="Test Interval"
label(data$clrnctestintvrl_mv)="Missing value verified for: Test Interval"
label(data$renal_clearance_visit_complete)="Complete?"
label(data$intfu_visitdate)="Date of Visit"
label(data$intfu_visitdate_mv)="Missing value verified for: Date of Visit"
label(data$intfu_visittype)="Type of Visit"
label(data$intfu_visittype_mv)="Missing value verified for: Type of Visit"
label(data$intfu_testintrvl)="Test Interval"
label(data$intfu_testintrvl_mv)="Missing value verified for: Test Interval"
label(data$interim_followup_visit_complete)="Complete?"
label(data$f4visitdate)="Visit Date"
label(data$f4visitdate_mv)="Missing value verified for: Visit Date"
label(data$f4visittype)="Visit type"
label(data$f4visittype_mv)="Missing value verified for: Visit type"
label(data$f4testintvl)="Test Interval"
label(data$f4testintvl_mv)="Missing value verified for: Test Interval"
label(data$f4bldprsr)="Blood pressure arm"
label(data$f4bldprsr_mv)="Missing value verified for: Blood pressure arm"
label(data$f4frstmsrsys)="1st measure - Systolic"
label(data$f4frstmsrsys_mv)="Missing value verified for: 1st measure - Systolic"
label(data$f4frstmsrdia)="1st measure - Diastolic"
label(data$f4frstmsrdia_mv)="Missing value verified for: 1st measure - Diastolic"
label(data$f4scndmsrsys)="2nd measure - Systolic"
label(data$f4scndmsrsys_mv)="Missing value verified for: 2nd measure - Systolic"
label(data$f4scndmsrdia)="2nd measure - Diastolic"
label(data$f4scndmsrdia_mv)="Missing value verified for: 2nd measure - Diastolic"
label(data$f4heartrate)="Heart rate (beats/min)"
label(data$f4heartrate_mv)="Missing value verified for: Heart rate (beats/min)"
label(data$f4pregnant)="Pregnancy Test?"
label(data$f4pregnant_mv)="Missing value verified for: Pregnancy Test?"
label(data$f4medsnow)="Is the subect currently taking any other medicine?"
label(data$f4medsnow_mv)="Missing value verified for: Is the subect currently taking any other medicine?"
label(data$f4srmpotasm)="Results of serum potassium (mmol/L)"
label(data$f4srmpotasm_mv)="Missing value verified for: Results of serum potassium (mmol/L)"
label(data$f4srmglu)="Results of serum glucose (mg/dL)"
label(data$f4srmglu_mv)="Missing value verified for: Results of serum glucose (mg/dL)"
label(data$f4forminit)="Initials of person who completed form"
label(data$f4forminit_mv)="Missing value verified for: Initials of person who completed form"
label(data$blood_pressure_form_04_complete)="Complete?"
label(data$f4_visitdate)="Visit date"
label(data$f4_visitdate_mv)="Missing value verified for: Visit date"
label(data$f4_pillbotldt)="Pill bottle date"
label(data$f4_pillbotldt_mv)="Missing value verified for: Pill bottle date"
label(data$f4_dsplstpill)="Date of last visit at which pills were dispensed"
label(data$f4_dsplstpill_mv)="Missing value verified for: Date of last visit at which pills were dispensed"
label(data$f4_dysbtwndsp)="Number of days between pill-dispensing visits"
label(data$f4_dysbtwndsp_mv)="Missing value verified for: Number of days between pill-dispensing visits"
label(data$f4_nodlydose)="Number of daily doses dispensed at last visit"
label(data$f4_nodlydose_mv)="Missing value verified for: Number of daily doses dispensed at last visit"
label(data$f4_dlydosrtrn)="Number of daily doses returned"
label(data$f4_dlydosrtrn_mv)="Missing value verified for: Number of daily doses returned"
label(data$f4_dlydosetkn)="Number of daily doses actually taken"
label(data$f4_dlydosetkn_mv)="Missing value verified for: Number of daily doses actually taken"
label(data$f4_compliance)="Compliance (%)"
label(data$f4_compliance_mv)="Missing value verified for: Compliance (%)"
label(data$f4_totdlydosn)="Total number of daily doses dispensed for next visit"
label(data$f4_totdlydosn_mv)="Missing value verified for: Total number of daily doses dispensed for next visit"
label(data$f4_curdos)="Current dose of medicine (mg)"
label(data$f4_curdos_mv)="Missing value verified for: Current dose of medicine (mg)"
label(data$f4_dateadjstd)="Date of medicine adjustment"
label(data$f4_dateadjstd_mv)="Missing value verified for: Date of medicine adjustment"
label(data$f4_ajstddose)="Adjusted dose (mg)"
label(data$f4_ajstddose_mv)="Missing value verified for: Adjusted dose (mg)"
label(data$f4_medseffects)="Side effects"
label(data$f4_medseffects_mv)="Missing value verified for: Side effects"
label(data$f4_expeffects)="If any side effects, please explain"
label(data$f4_expeffects_mv)="Missing value verified for: If any side effects, please explain"
label(data$drug_compliance_form_04_complete)="Complete?"
label(data$f4meds_visitdate)="Visit Date"
label(data$f4meds_visitdate_mv)="Missing value verified for: Visit Date"
label(data$f4drugcode1)="Drug Code"
label(data$f4drugcode1_mv)="Missing value verified for: Drug Code"
label(data$f4dosage1)="Dosage"
label(data$f4dosage1_mv)="Missing value verified for: Dosage"
label(data$f4timesday1)="Number of times taken"
label(data$f4timesday1_mv)="Missing value verified for: Number of times taken"
label(data$f4dwm1)="Times/day, week or month"
label(data$f4dwm1_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode2)="Drug Code"
label(data$f4drugcode2_mv)="Missing value verified for: Drug Code"
label(data$f4dosage2)="Dosage"
label(data$f4dosage2_mv)="Missing value verified for: Dosage"
label(data$f4timesday2)="Number of time taken"
label(data$f4timesday2_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm2)="Times/day, week or month"
label(data$f4dwm2_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode3)="Drug Code"
label(data$f4drugcode3_mv)="Missing value verified for: Drug Code"
label(data$f4dosage3)="Dosage"
label(data$f4dosage3_mv)="Missing value verified for: Dosage"
label(data$f4timesday3)="Number of time taken"
label(data$f4timesday3_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm3)="Times/day, week or month"
label(data$f4dwm3_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode4)="Drug Code"
label(data$f4drugcode4_mv)="Missing value verified for: Drug Code"
label(data$f4dosage4)="Dosage"
label(data$f4dosage4_mv)="Missing value verified for: Dosage"
label(data$f4timesday4)="Number of time taken"
label(data$f4timesday4_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm4)="Times/day, week or month"
label(data$f4dwm4_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode5)="Drug Code"
label(data$f4drugcode5_mv)="Missing value verified for: Drug Code"
label(data$f4dosage5)="Dosage"
label(data$f4dosage5_mv)="Missing value verified for: Dosage"
label(data$f4timesday5)="Number of time taken"
label(data$f4timesday5_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm5)="Times/day, week or month"
label(data$f4dwm5_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode6)="Drug Code"
label(data$f4drugcode6_mv)="Missing value verified for: Drug Code"
label(data$f4dosage6)="Dosage"
label(data$f4dosage6_mv)="Missing value verified for: Dosage"
label(data$f4timesday6)="Number of time taken"
label(data$f4timesday6_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm6)="Times/day, week or month"
label(data$f4dwm6_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode7)="Drug Code"
label(data$f4drugcode7_mv)="Missing value verified for: Drug Code"
label(data$f4dosage7)="Dosage"
label(data$f4dosage7_mv)="Missing value verified for: Dosage"
label(data$f4timesday7)="Number of time taken"
label(data$f4timesday7_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm7)="Times/day, week or month"
label(data$f4dwm7_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode8)="Drug Code"
label(data$f4drugcode8_mv)="Missing value verified for: Drug Code"
label(data$f4dosage8)="Dosage"
label(data$f4dosage8_mv)="Missing value verified for: Dosage"
label(data$f4timesday8)="Number of time taken"
label(data$f4timesday8_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm8)="Times/day, week or month"
label(data$f4dwm8_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode9)="Drug Code"
label(data$f4drugcode9_mv)="Missing value verified for: Drug Code"
label(data$f4dosage9)="Dosage"
label(data$f4dosage9_mv)="Missing value verified for: Dosage"
label(data$f4timesday9)="Number of time taken"
label(data$f4timesday9_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm9)="Times/day, week or month"
label(data$f4dwm9_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode10)="Drug Code"
label(data$f4drugcode10_mv)="Missing value verified for: Drug Code"
label(data$f4dosage10)="Dosage"
label(data$f4dosage10_mv)="Missing value verified for: Dosage"
label(data$f4timesday10)="Number of time taken"
label(data$f4timesday10_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm10)="Times/day, week or month"
label(data$f4dwm10_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode11)="Drug Code"
label(data$f4drugcode11_mv)="Missing value verified for: Drug Code"
label(data$f4dosage11)="Dosage"
label(data$f4dosage11_mv)="Missing value verified for: Dosage"
label(data$f4timesday11)="Number of time taken"
label(data$f4timesday11_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm11)="Times/day, week or month"
label(data$f4dwm11_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode12)="Drug Code"
label(data$f4drugcode12_mv)="Missing value verified for: Drug Code"
label(data$f4dosage12)="Dosage"
label(data$f4dosage12_mv)="Missing value verified for: Dosage"
label(data$f4timesday12)="Number of time taken"
label(data$f4timesday12_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm12)="Times/day, week or month"
label(data$f4dwm12_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode13)="Drug Code"
label(data$f4drugcode13_mv)="Missing value verified for: Drug Code"
label(data$f4dosage13)="Dosage"
label(data$f4dosage13_mv)="Missing value verified for: Dosage"
label(data$f4timesday13)="Number of time taken"
label(data$f4timesday13_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm13)="Times/day, week or month"
label(data$f4dwm13_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode14)="Drug Code"
label(data$f4drugcode14_mv)="Missing value verified for: Drug Code"
label(data$f4dosage14)="Dosage"
label(data$f4dosage14_mv)="Missing value verified for: Dosage"
label(data$f4timesday14)="Number of time taken"
label(data$f4timesday14_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm14)="Times/day, week or month"
label(data$f4dwm14_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode15)="Drug Code"
label(data$f4drugcode15_mv)="Missing value verified for: Drug Code"
label(data$f4dosage15)="Dosage"
label(data$f4dosage15_mv)="Missing value verified for: Dosage"
label(data$f4timesday15)="Number of time taken"
label(data$f4timesday15_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm15)="Times/day, week or month"
label(data$f4dwm15_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode16)="Drug Code"
label(data$f4drugcode16_mv)="Missing value verified for: Drug Code"
label(data$f4dosage16)="Dosage"
label(data$f4dosage16_mv)="Missing value verified for: Dosage"
label(data$f4timesday16)="Number of time taken"
label(data$f4timesday16_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm16)="Times/day, week or month"
label(data$f4dwm16_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode17)="Drug Code"
label(data$f4drugcode17_mv)="Missing value verified for: Drug Code"
label(data$f4dosage17)="Dosage"
label(data$f4dosage17_mv)="Missing value verified for: Dosage"
label(data$f4timesday17)="Number of time taken"
label(data$f4timesday17_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm17)="Times/day, week or month"
label(data$f4dwm17_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode18)="Drug Code"
label(data$f4drugcode18_mv)="Missing value verified for: Drug Code"
label(data$f4dosage18)="Dosage"
label(data$f4dosage18_mv)="Missing value verified for: Dosage"
label(data$f4timesday18)="Number of time taken"
label(data$f4timesday18_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm18)="Times/day, week or month"
label(data$f4dwm18_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode19)="Drug Code"
label(data$f4drugcode19_mv)="Missing value verified for: Drug Code"
label(data$f4dosage19)="Dosage"
label(data$f4dosage19_mv)="Missing value verified for: Dosage"
label(data$f4timesday19)="Number of time taken"
label(data$f4timesday19_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm19)="Times/day, week or month"
label(data$f4dwm19_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode20)="Drug Code"
label(data$f4drugcode20_mv)="Missing value verified for: Drug Code"
label(data$f4dosage20)="Dosage"
label(data$f4dosage20_mv)="Missing value verified for: Dosage"
label(data$f4timesday20)="Number of time taken"
label(data$f4timesday20_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm20)="Times/day, week or month"
label(data$f4dwm20_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode21)="Drug Code"
label(data$f4drugcode21_mv)="Missing value verified for: Drug Code"
label(data$f4dosage21)="Dosage"
label(data$f4dosage21_mv)="Missing value verified for: Dosage"
label(data$f4timesday21)="Number of time taken"
label(data$f4timesday21_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm21)="Times/day, week or month"
label(data$f4dwm21_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode22)="Drug Code"
label(data$f4drugcode22_mv)="Missing value verified for: Drug Code"
label(data$f4dosage22)="Dosage"
label(data$f4dosage22_mv)="Missing value verified for: Dosage"
label(data$f4timesday22)="Number of time taken"
label(data$f4timesday22_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm22)="Times/day, week or month"
label(data$f4dwm22_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode23)="Drug Code"
label(data$f4drugcode23_mv)="Missing value verified for: Drug Code"
label(data$f4dosage23)="Dosage"
label(data$f4dosage23_mv)="Missing value verified for: Dosage"
label(data$f4timesday23)="Number of time taken"
label(data$f4timesday23_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm23)="Times/day, week or month"
label(data$f4dwm23_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode24)="Drug Code"
label(data$f4drugcode24_mv)="Missing value verified for: Drug Code"
label(data$f4dosage24)="Dosage"
label(data$f4dosage24_mv)="Missing value verified for: Dosage"
label(data$f4timesday24)="Number of time taken"
label(data$f4timesday24_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm24)="Times/day, week or month"
label(data$f4dwm24_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode25)="Drug Code"
label(data$f4drugcode25_mv)="Missing value verified for: Drug Code"
label(data$f4dosage25)="Dosage"
label(data$f4dosage25_mv)="Missing value verified for: Dosage"
label(data$f4timesday25)="Number of time taken"
label(data$f4timesday25_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm25)="Times/day, week or month"
label(data$f4dwm25_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode26)="Drug Code"
label(data$f4drugcode26_mv)="Missing value verified for: Drug Code"
label(data$f4dosage26)="Dosage"
label(data$f4dosage26_mv)="Missing value verified for: Dosage"
label(data$f4timesday26)="Number of time taken"
label(data$f4timesday26_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm26)="Times/day, week or month"
label(data$f4dwm26_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode27)="Drug Code"
label(data$f4drugcode27_mv)="Missing value verified for: Drug Code"
label(data$f4dosage27)="Dosage"
label(data$f4dosage27_mv)="Missing value verified for: Dosage"
label(data$f4timesday27)="Number of time taken"
label(data$f4timesday27_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm27)="Times/day, week or month"
label(data$f4dwm27_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode28)="Drug Code"
label(data$f4drugcode28_mv)="Missing value verified for: Drug Code"
label(data$f4dosage28)="Dosage"
label(data$f4dosage28_mv)="Missing value verified for: Dosage"
label(data$f4timesday28)="Number of time taken"
label(data$f4timesday28_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm28)="Times/day, week or month"
label(data$f4dwm28_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode29)="Drug Code"
label(data$f4drugcode29_mv)="Missing value verified for: Drug Code"
label(data$f4dosage29)="Dosage"
label(data$f4dosage29_mv)="Missing value verified for: Dosage"
label(data$f4timesday29)="Number of time taken"
label(data$f4timesday29_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm29)="Times/day, week or month"
label(data$f4dwm29_mv)="Missing value verified for: Times/day, week or month"
label(data$f4drugcode30)="Drug Code"
label(data$f4drugcode30_mv)="Missing value verified for: Drug Code"
label(data$f4dosage30)="Dosage"
label(data$f4dosage30_mv)="Missing value verified for: Dosage"
label(data$f4timesday30)="Number of time taken"
label(data$f4timesday30_mv)="Missing value verified for: Number of time taken"
label(data$f4dwm30)="Times/day, week or month"
label(data$f4dwm30_mv)="Missing value verified for: Times/day, week or month"
label(data$meds_form_04_complete)="Complete?"
label(data$f5visitdate)="Visit Date"
label(data$f5visitdate_mv)="Missing value verified for: Visit Date"
label(data$f5rcntillns)="Recent Illnesses Rqring Med Att."
label(data$f5rcntillns_mv)="Missing value verified for: Recent Illnesses Rqring Med Att."
label(data$f5outpat)="Out Patient?"
label(data$f5outpat_mv)="Missing value verified for: Out Patient?"
label(data$f5outpatspec)="Out Patient - Specify"
label(data$f5outpatspec_mv)="Missing value verified for: Out Patient - Specify"
label(data$f5hospwosrg)="Hosp w/o Surgery"
label(data$f5hospwosrg_mv)="Missing value verified for: Hosp w/o Surgery"
label(data$f5hospwoexp)="Hosp w/o Surgery - Exp"
label(data$f5hospwoexp_mv)="Missing value verified for: Hosp w/o Surgery - Exp"
label(data$f5hospwsrg)="Hosp w/Surgery"
label(data$f5hospwsrg_mv)="Missing value verified for: Hosp w/Surgery"
label(data$f5hospwexp)="Hosp w/Surgery - Exp"
label(data$f5hospwexp_mv)="Missing value verified for: Hosp w/Surgery - Exp"
label(data$f5medprbsdoc)="Any Docmntd Meds Prob since Last Phys?"
label(data$f5medprbsdoc_mv)="Missing value verified for: Any Docmntd Meds Prob since Last Phys?"
label(data$f5crnryartds)="Coronary artery disease?"
label(data$f5crnryartds_mv)="Missing value verified for: Coronary artery disease?"
label(data$f5cancer)="Cancer?"
label(data$f5cancer_mv)="Missing value verified for: Cancer?"
label(data$f5cancerexp)="Cancer - Explain"
label(data$f5cancerexp_mv)="Missing value verified for: Cancer - Explain"
label(data$f5crbrlvasds)="Cerebral Vascular disease?"
label(data$f5crbrlvasds_mv)="Missing value verified for: Cerebral Vascular disease?"
label(data$f5peripvasds)="Peripheral vascular disease?"
label(data$f5peripvasds_mv)="Missing value verified for: Peripheral vascular disease?"
label(data$f5hypertensn)="Hypertension?"
label(data$f5hypertensn_mv)="Missing value verified for: Hypertension?"
label(data$f5seizures)="Seizures?"
label(data$f5seizures_mv)="Missing value verified for: Seizures?"
label(data$f5gnitrnryds)="Genitourinary disease?"
label(data$f5gnitrnryds_mv)="Missing value verified for: Genitourinary disease?"
label(data$f5gnitrnyexp)="GD Explain"
label(data$f5gnitrnyexp_mv)="Missing value verified for: GD Explain"
label(data$f5lungdsz)="Lung Disease?"
label(data$f5lungdsz_mv)="Missing value verified for: Lung Disease?"
label(data$f5majsurg)="Major surgery?"
label(data$f5majsurg_mv)="Missing value verified for: Major surgery?"
label(data$f5majsurgexp)="MS Explain"
label(data$f5majsurgexp_mv)="Missing value verified for: MS Explain"
label(data$f5othrmeddia)="Other medical diagnosis?"
label(data$f5othrmeddia_mv)="Missing value verified for: Other medical diagnosis?"
label(data$f5othrmedexp)="OMD Explain"
label(data$f5othrmedexp_mv)="Missing value verified for: OMD Explain"
label(data$f5lungs)="Lungs"
label(data$f5lungs_mv)="Missing value verified for: Lungs"
label(data$f5lungsexp)="Explain"
label(data$f5lungsexp_mv)="Missing value verified for: Explain"
label(data$f5heart)="Heart"
label(data$f5heart_mv)="Missing value verified for: Heart"
label(data$f5heartexp)="Explain"
label(data$f5heartexp_mv)="Missing value verified for: Explain"
label(data$f5skin)="Skin"
label(data$f5skin_mv)="Missing value verified for: Skin"
label(data$f5skinexp)="Explain"
label(data$f5skinexp_mv)="Missing value verified for: Explain"
label(data$f5edema)="Edema"
label(data$f5edema_mv)="Missing value verified for: Edema"
label(data$f5edemaexp)="Explain"
label(data$f5edemaexp_mv)="Missing value verified for: Explain"
label(data$history_and_physical_examination_form_05_complete)="Complete?"
label(data$f6visitdate)="Visit Date"
label(data$f6visitdate_mv)="Missing value verified for: Visit Date"
label(data$f6ivsites)="IV Sites"
label(data$f6ivsites_mv)="Missing value verified for: IV Sites"
label(data$f6weight)="Weight (kg)"
label(data$f6weight_mv)="Missing value verified for: Weight (kg)"
label(data$f6bsa)="BSA"
label(data$f6bsa_mv)="Missing value verified for: BSA"
label(data$f6u0strttime)="U0 start time (24 hr. clock)"
label(data$f6u0strttime_mv)="Missing value verified for: U0 start time (24 hr. clock)"
label(data$f6u0volume)="U0 volume (ml)"
label(data$f6u0volume_mv)="Missing value verified for: U0 volume (ml)"
label(data$f6bldprsrarm)="Bld Prsr (arm used)"
label(data$f6bldprsrarm_mv)="Missing value verified for: Bld Prsr (arm used)"
label(data$f6frstmsrsty)="1st measure (systolic)"
label(data$f6frstmsrsty_mv)="Missing value verified for: 1st measure (systolic)"
label(data$f6frstmsrdia)="1st measure (diastolic)"
label(data$f6frstmsrdia_mv)="Missing value verified for: 1st measure (diastolic)"
label(data$f6scndmsrsys)="2nd measure (systolic)"
label(data$f6scndmsrsys_mv)="Missing value verified for: 2nd measure (systolic)"
label(data$f6scndmsrdia)="2nd measure (diastolic)"
label(data$f6scndmsrdia_mv)="Missing value verified for: 2nd measure (diastolic)"
label(data$f6heartrate)="Heart rate (beats/min)"
label(data$f6heartrate_mv)="Missing value verified for: Heart rate (beats/min)"
label(data$f6u0h20)="U0 H20 + 10cc/kg load (ml)"
label(data$f6u0h20_mv)="Missing value verified for: U0 H20 + 10cc/kg load (ml)"
label(data$f6fnleqlurtm)="Final equilibration urine time (min)"
label(data$f6fnleqlurtm_mv)="Missing value verified for: Final equilibration urine time (min)"
label(data$f6p1time)="P1 time (min)"
label(data$f6p1time_mv)="Missing value verified for: P1 time (min)"
label(data$f6urnvol)="Total equilibration period urine volume (ml)"
label(data$f6urnvol_mv)="Missing value verified for: Total equilibration period urine volume (ml)"
label(data$f6h20vol)="Total equilibration H20 volume (ml)"
label(data$f6h20vol_mv)="Missing value verified for: Total equilibration H20 volume (ml)"
label(data$f6u1time)="U1 time (min)"
label(data$f6u1time_mv)="Missing value verified for: U1 time (min)"
label(data$f6p2time)="P2 time (min)"
label(data$f6p2time_mv)="Missing value verified for: P2 time (min)"
label(data$f6u1vol)="U1 volume (ml)"
label(data$f6u1vol_mv)="Missing value verified for: U1 volume (ml)"
label(data$f6u1h20vol)="U1 H20 volume (ml)"
label(data$f6u1h20vol_mv)="Missing value verified for: U1 H20 volume (ml)"
label(data$f6u2time)="U2 time (min)"
label(data$f6u2time_mv)="Missing value verified for: U2 time (min)"
label(data$f6p3time)="P3 time (min)"
label(data$f6p3time_mv)="Missing value verified for: P3 time (min)"
label(data$f6u2vol)="U2 volume (ml)"
label(data$f6u2vol_mv)="Missing value verified for: U2 volume (ml)"
label(data$f6u2h20vol)="U2 H20 volume (ml)"
label(data$f6u2h20vol_mv)="Missing value verified for: U2 H20 volume (ml)"
label(data$f6u3time)="U3 time (min)"
label(data$f6u3time_mv)="Missing value verified for: U3 time (min)"
label(data$f6p4time)="P4 time (min)"
label(data$f6p4time_mv)="Missing value verified for: P4 time (min)"
label(data$f6u3vol)="U3 volume (ml)"
label(data$f6u3vol_mv)="Missing value verified for: U3 volume (ml)"
label(data$f6u3h20vol)="U3 H20 volume (ml)"
label(data$f6u3h20vol_mv)="Missing value verified for: U3 H20 volume (ml)"
label(data$f6u4time)="U4 time (min)"
label(data$f6u4time_mv)="Missing value verified for: U4 time (min)"
label(data$f6p5time)="P5 time (min)"
label(data$f6p5time_mv)="Missing value verified for: P5 time (min)"
label(data$f6u4vol)="U4 volume (ml)"
label(data$f6u4vol_mv)="Missing value verified for: U4 volume (ml)"
label(data$f6u4h20vol)="U4 H20 volume (ml)"
label(data$f6u4h20vol_mv)="Missing value verified for: U4 H20 volume (ml)"
label(data$f6u5time)="U5 time (min)"
label(data$f6u5time_mv)="Missing value verified for: U5 time (min)"
label(data$f6p6time)="P6 time (min)"
label(data$f6p6time_mv)="Missing value verified for: P6 time (min)"
label(data$f6u5vol)="U5 volume (ml)"
label(data$f6u5vol_mv)="Missing value verified for: U5 volume (ml)"
label(data$f6u5h20vol)="U5 H2O volume (ml)"
label(data$f6u5h20vol_mv)="Missing value verified for: U5 H2O volume (ml)"
label(data$f6urndate)="Urinalysis date"
label(data$f6urndate_mv)="Missing value verified for: Urinalysis date"
label(data$f6glucose)="Glucose"
label(data$f6glucose_mv)="Missing value verified for: Glucose"
label(data$f6bilirubin)="Bilirubin"
label(data$f6bilirubin_mv)="Missing value verified for: Bilirubin"
label(data$f6ketones)="Ketones"
label(data$f6ketones_mv)="Missing value verified for: Ketones"
label(data$f6specgrvty)="Specific gravity"
label(data$f6specgrvty_mv)="Missing value verified for: Specific gravity"
label(data$f6bloodocult)="Blood, occult"
label(data$f6bloodocult_mv)="Missing value verified for: Blood, occult"
label(data$f6ph)="pH"
label(data$f6ph_mv)="Missing value verified for: pH"
label(data$f6protein)="Protein"
label(data$f6protein_mv)="Missing value verified for: Protein"
label(data$f6casts)="Casts"
label(data$f6casts_mv)="Missing value verified for: Casts"
label(data$f6rbccast)="RBC"
label(data$f6rbccast_mv)="Missing value verified for: RBC"
label(data$f6wbccast)="WBC"
label(data$f6wbccast_mv)="Missing value verified for: WBC"
label(data$f6hyalinecast)="Hyaline"
label(data$f6hyalinecast_mv)="Missing value verified for: Hyaline"
label(data$f6granularcast)="Granular"
label(data$f6granularcast_mv)="Missing value verified for: Granular"
label(data$f6othercast)="Other "
label(data$f6othercast_mv)="Missing value verified for: Other "
label(data$f6othercast_exp)="Other (specify)"
label(data$f6othercast_exp_mv)="Missing value verified for: Other (specify)"
label(data$f6wbc)="WBC"
label(data$f6wbc_mv)="Missing value verified for: WBC"
label(data$f6rbc)="RBC"
label(data$f6rbc_mv)="Missing value verified for: RBC"
label(data$f6epithcells)="Epithelial cells"
label(data$f6epithcells_mv)="Missing value verified for: Epithelial cells"
label(data$f6bacteria)="Bacteria"
label(data$f6bacteria_mv)="Missing value verified for: Bacteria"
label(data$renal_clearance_data_form_06_complete)="Complete?"
label(data$f7qcnumber)="Quality control number"
label(data$f7qcnumber_mv)="Missing value verified for: Quality control number"
label(data$f7visitdate)="Date of sample collection"
label(data$f7visitdate_mv)="Missing value verified for: Date of sample collection"
label(data$f7tstintrvl)="Test interval"
label(data$f7tstintrvl_mv)="Missing value verified for: Test interval"
label(data$f7typesample)="Type of sample"
label(data$f7typesample_mv)="Missing value verified for: Type of sample"
label(data$laboratory_quality_control_form_07_complete)="Complete?"
label(data$f8visitdate)="Date sample collected"
label(data$f8visitdate_mv)="Missing value verified for: Date sample collected"
label(data$f8analdate)="Date sample analyzed"
label(data$f8analdate_mv)="Missing value verified for: Date sample analyzed"
label(data$f8alti)="ALTI (U/L)"
label(data$f8alti_mv)="Missing value verified for: ALTI (U/L)"
label(data$f8glucose)="Glucose (mg/dL)"
label(data$f8glucose_mv)="Missing value verified for: Glucose (mg/dL)"
label(data$f8bun)="BUN (mg/dL)"
label(data$f8bun_mv)="Missing value verified for: BUN (mg/dL)"
label(data$f8creatinine)="Creatinine (mg/dL)"
label(data$f8creatinine_mv)="Missing value verified for: Creatinine (mg/dL)"
label(data$f8sodium)="Sodium (mmol/L)"
label(data$f8sodium_mv)="Missing value verified for: Sodium (mmol/L)"
label(data$f8potassium)="Potassium (mmol/L)"
label(data$f8potassium_mv)="Missing value verified for: Potassium (mmol/L)"
label(data$f8chloride)="Chloride (mmol/L)"
label(data$f8chloride_mv)="Missing value verified for: Chloride (mmol/L)"
label(data$f8calcium)="Calcium (mg/dL)"
label(data$f8calcium_mv)="Missing value verified for: Calcium (mg/dL)"
label(data$f8phosphorus)="Phosphorus (mg/dL)"
label(data$f8phosphorus_mv)="Missing value verified for: Phosphorus (mg/dL)"
label(data$f8totlprotn)="Total protein (gm/dL)"
label(data$f8totlprotn_mv)="Missing value verified for: Total protein (gm/dL)"
label(data$f8albumin)="Albumin (gm/dL)"
label(data$f8albumin_mv)="Missing value verified for: Albumin (gm/dL)"
label(data$f8cholestero)="Cholesterol (mg/dL)"
label(data$f8cholestero_mv)="Missing value verified for: Cholesterol (mg/dL)"
label(data$f8triglyceri)="Triglycerides (mg/dL)"
label(data$f8triglyceri_mv)="Missing value verified for: Triglycerides (mg/dL)"
label(data$f8hdl)="HDL Cholesterol (mg/dl)"
label(data$f8hdl_mv)="Missing value verified for: HDL Cholesterol (mg/dl)"
label(data$f8ldl)="LDL Cholesterol (mg/dl)"
label(data$f8ldl_mv)="Missing value verified for: LDL Cholesterol (mg/dl)"
label(data$f8totlbilrbn)="Total bilirubin (mg/dL)"
label(data$f8totlbilrbn_mv)="Missing value verified for: Total bilirubin (mg/dL)"
label(data$f8dbilirubin)="Direct bilirubin (mg/dL)"
label(data$f8dbilirubin_mv)="Missing value verified for: Direct bilirubin (mg/dL)"
label(data$f8alkp04)="ALK P04 (alkaline phosphatase) (U/L)"
label(data$f8alkp04_mv)="Missing value verified for: ALK P04 (alkaline phosphatase) (U/L)"
label(data$f8ast)="AST (SGOT) (U/L)"
label(data$f8ast_mv)="Missing value verified for: AST (SGOT) (U/L)"
label(data$smacchemistry_panel_form_08_complete)="Complete?"
label(data$f9visitdate)="Date of sample collection"
label(data$f9visitdate_mv)="Missing value verified for: Date of sample collection"
label(data$f9dateanal)="Date specimen analyzed"
label(data$f9dateanal_mv)="Missing value verified for: Date specimen analyzed"
label(data$f9wbc)="WBC (x 10^3)"
label(data$f9wbc_mv)="Missing value verified for: WBC (x 10^3)"
label(data$f9rbc)="RBC (x 10^6)"
label(data$f9rbc_mv)="Missing value verified for: RBC (x 10^6)"
label(data$f9hemoglob)="Hemoglobin (g/dL)"
label(data$f9hemoglob_mv)="Missing value verified for: Hemoglobin (g/dL)"
label(data$f9hematocrit)="Hematocrit (%)"
label(data$f9hematocrit_mv)="Missing value verified for: Hematocrit (%)"
label(data$f9mcv)="MCV (fL)"
label(data$f9mcv_mv)="Missing value verified for: MCV (fL)"
label(data$f9mch)="MCH (pg)"
label(data$f9mch_mv)="Missing value verified for: MCH (pg)"
label(data$f9mchc)="MCHC (%)"
label(data$f9mchc_mv)="Missing value verified for: MCHC (%)"
label(data$f9platelets)="Platelets (x 10^3)"
label(data$f9platelets_mv)="Missing value verified for: Platelets (x 10^3)"
label(data$f9lymphinst)="Lymph, Inst%"
label(data$f9lymphinst_mv)="Missing value verified for: Lymph, Inst%"
label(data$f9monoinst)="Mono, Inst%"
label(data$f9monoinst_mv)="Missing value verified for: Mono, Inst%"
label(data$f9neutinst)="Neut, Inst%"
label(data$f9neutinst_mv)="Missing value verified for: Neut, Inst%"
label(data$f9eosinst)="Eos, Inst%"
label(data$f9eosinst_mv)="Missing value verified for: Eos, Inst%"
label(data$f9basoinst)="Baso, Inst%"
label(data$f9basoinst_mv)="Missing value verified for: Baso, Inst%"
label(data$f9gran)="Gran, %"
label(data$f9gran_mv)="Missing value verified for: Gran, %"
label(data$cbc_form_09_complete)="Complete?"
label(data$f10visitdate)="Visit Date"
label(data$f10visitdate_mv)="Missing value verified for: Visit Date"
label(data$f10spturalbu)="Spot urine albumin (mg/L)"
label(data$f10spturalbu_mv)="Missing value verified for: Spot urine albumin (mg/L)"
label(data$f10spturncrea)="Spot urine creatinine (g/L)"
label(data$f10spturncrea_mv)="Missing value verified for: Spot urine creatinine (g/L)"
label(data$f10serumcreat)="Serum creatinine (mg/dL)"
label(data$f10serumcreat_mv)="Missing value verified for: Serum creatinine (mg/dL)"
label(data$f10strgbxnum)="Storage box number"
label(data$f10strgbxnum_mv)="Missing value verified for: Storage box number"
label(data$daes_routine_labs_form_10_complete)="Complete?"
label(data$f11visitdate)="Visit date"
label(data$f11visitdate_mv)="Missing value verified for: Visit date"
label(data$f11hba1c)="HbA1c (%)"
label(data$f11hba1c_mv)="Missing value verified for: HbA1c (%)"
label(data$f11p0glu)="P0 Glucose"
label(data$f11p0glu_mv)="Missing value verified for: P0 Glucose"
label(data$f11p0creatn)="P0 serum creatinine (mg/dL)"
label(data$f11p0creatn_mv)="Missing value verified for: P0 serum creatinine (mg/dL)"
label(data$f11p0albumn)="P0 albumin (mg/L)"
label(data$f11p0albumn_mv)="Missing value verified for: P0 albumin (mg/L)"
label(data$f11p3albumn)="P3 albumin (mg/L)"
label(data$f11p3albumn_mv)="Missing value verified for: P3 albumin (mg/L)"
label(data$f11p0iggsamp)="P0 IgG (mg/L)"
label(data$f11p0iggsamp_mv)="Missing value verified for: P0 IgG (mg/L)"
label(data$f11p3iggsamp)="P3 IgG (mg/L)"
label(data$f11p3iggsamp_mv)="Missing value verified for: P3 IgG (mg/L)"
label(data$f11u0creatn)="U0 creatinine (g/L)"
label(data$f11u0creatn_mv)="Missing value verified for: U0 creatinine (g/L)"
label(data$f11u0albumn)="U0 albumin (mg/L)"
label(data$f11u0albumn_mv)="Missing value verified for: U0 albumin (mg/L)"
label(data$f11u3albumn)="U3 albumin (mg/L)"
label(data$f11u3albumn_mv)="Missing value verified for: U3 albumin (mg/L)"
label(data$f11u0igg)="U0 IgG (mg/L)"
label(data$f11u0igg_mv)="Missing value verified for: U0 IgG (mg/L)"
label(data$f11u3igg)="U3 IgG (mg/L)"
label(data$f11u3igg_mv)="Missing value verified for: U3 IgG (mg/L)"
label(data$daes_lab_clearance_results_form_11_complete)="Complete?"
label(data$f12expdatevis)="Expected date of visit"
label(data$f12expdatevis_mv)="Missing value verified for: Expected date of visit"
label(data$f12typevisit)="Type of visit"
label(data$f12typevisit_mv)="Missing value verified for: Type of visit"
label(data$f12testintrvl)="Test interval"
label(data$f12testintrvl_mv)="Missing value verified for: Test interval"
label(data$f12reasnmissd)="Reason visit was missed"
label(data$f12reasnmissd_mv)="Missing value verified for: Reason visit was missed"
label(data$f12otherspec)="Reason visit was missed, other (specify)"
label(data$f12otherspec_mv)="Missing value verified for: Reason visit was missed, other (specify)"
label(data$missed_visitintercurrent_event_form_12_complete)="Complete?"
label(data$f16visitdate)="Visit Date"
label(data$f16visitdate_mv)="Missing value verified for: Visit Date"
label(data$f16qcnumber)="Quality control number"
label(data$f16qcnumber_mv)="Missing value verified for: Quality control number"
label(data$f16analdate)="Date sample analyzed"
label(data$f16analdate_mv)="Missing value verified for: Date sample analyzed"
label(data$f16glucose)="Glucose (mg/dL)"
label(data$f16glucose_mv)="Missing value verified for: Glucose (mg/dL)"
label(data$f16bun)="BUN (mg/dL)"
label(data$f16bun_mv)="Missing value verified for: BUN (mg/dL)"
label(data$f16creatinine)="Creatinine (mg/dL)"
label(data$f16creatinine_mv)="Missing value verified for: Creatinine (mg/dL)"
label(data$f16sodium)="Sodium (mmol/L)"
label(data$f16sodium_mv)="Missing value verified for: Sodium (mmol/L)"
label(data$f16potassium)="Potassium (mmol/L)"
label(data$f16potassium_mv)="Missing value verified for: Potassium (mmol/L)"
label(data$f16chloride)="Chloride (mmol/L)"
label(data$f16chloride_mv)="Missing value verified for: Chloride (mmol/L)"
label(data$f16calcium)="Calcium (mg/dL)"
label(data$f16calcium_mv)="Missing value verified for: Calcium (mg/dL)"
label(data$f16phosphorus)="Phosphorus (mg/dL)"
label(data$f16phosphorus_mv)="Missing value verified for: Phosphorus (mg/dL)"
label(data$f16totlprotn)="Total protein (gm/dL)"
label(data$f16totlprotn_mv)="Missing value verified for: Total protein (gm/dL)"
label(data$f16albumin)="Albumin (gm/dL)"
label(data$f16albumin_mv)="Missing value verified for: Albumin (gm/dL)"
label(data$f16cholestero)="Cholesterol (mg/dL)"
label(data$f16cholestero_mv)="Missing value verified for: Cholesterol (mg/dL)"
label(data$f16triglyceri)="Triglycerides (mg/dL)"
label(data$f16triglyceri_mv)="Missing value verified for: Triglycerides (mg/dL)"
label(data$f16totlbilrbn)="Total bilirubin (mg/dL)"
label(data$f16totlbilrbn_mv)="Missing value verified for: Total bilirubin (mg/dL)"
label(data$f16alkp04)="ALK P04 (alkaline phosphatase) (U/L)"
label(data$f16alkp04_mv)="Missing value verified for: ALK P04 (alkaline phosphatase) (U/L)"
label(data$f16ast)="AST (SGOT) (U/L)"
label(data$f16ast_mv)="Missing value verified for: AST (SGOT) (U/L)"
label(data$smac_quality_control_form_16_complete)="Complete?"
label(data$f17visitdate)="Visit Date"
label(data$f17visitdate_mv)="Missing value verified for: Visit Date"
label(data$f17qcnumber)="Quality control number"
label(data$f17qcnumber_mv)="Missing value verified for: Quality control number"
label(data$f17analdate)="Date specimen analyzed"
label(data$f17analdate_mv)="Missing value verified for: Date specimen analyzed"
label(data$f17wbc)="WBC (x 10^3)"
label(data$f17wbc_mv)="Missing value verified for: WBC (x 10^3)"
label(data$f17rbc)="RBC (x 10^6)"
label(data$f17rbc_mv)="Missing value verified for: RBC (x 10^6)"
label(data$f17hemoglob)="Hemoglobin (g/dL)"
label(data$f17hemoglob_mv)="Missing value verified for: Hemoglobin (g/dL)"
label(data$f17hematocrit)="Hematocrit (%)"
label(data$f17hematocrit_mv)="Missing value verified for: Hematocrit (%)"
label(data$f17mcv)="MCV (femtoliter)"
label(data$f17mcv_mv)="Missing value verified for: MCV (femtoliter)"
label(data$f17mch)="MCH (pg)"
label(data$f17mch_mv)="Missing value verified for: MCH (pg)"
label(data$f17mchc)="MCHC (%)"
label(data$f17mchc_mv)="Missing value verified for: MCHC (%)"
label(data$f17platelets)="Platelets (x 10^3)"
label(data$f17platelets_mv)="Missing value verified for: Platelets (x 10^3)"
label(data$cbc_quality_control_form_17_complete)="Complete?"
label(data$f18visitdate)="Visit date"
label(data$f18visitdate_mv)="Missing value verified for: Visit date"
label(data$f18qcnumber)="Quality control number"
label(data$f18qcnumber_mv)="Missing value verified for: Quality control number"
label(data$f18hba1c)="HbA1c (%)"
label(data$f18hba1c_mv)="Missing value verified for: HbA1c (%)"
label(data$f18p0creatn)="P0 serum creatinine (mg/dL)"
label(data$f18p0creatn_mv)="Missing value verified for: P0 serum creatinine (mg/dL)"
label(data$f18p0albumin)="P0 albumin (mg/L)"
label(data$f18p0albumin_mv)="Missing value verified for: P0 albumin (mg/L)"
label(data$f18p0igg)="P0 IgG (mg/L)"
label(data$f18p0igg_mv)="Missing value verified for: P0 IgG (mg/L)"
label(data$f18u0creatn)="U0 creatinine (g/L)"
label(data$f18u0creatn_mv)="Missing value verified for: U0 creatinine (g/L)"
label(data$f18u0albumin)="U0 albumin (mg/L)"
label(data$f18u0albumin_mv)="Missing value verified for: U0 albumin (mg/L)"
label(data$f18u0igg)="U0 IgG (mg/L)"
label(data$f18u0igg_mv)="Missing value verified for: U0 IgG (mg/L)"
label(data$daes_clearance_quality_control_results_form_18_complete)="Complete?"
label(data$f19qcnumber)="Quality control number"
label(data$f19qcnumber_mv)="Missing value verified for: Quality control number"
label(data$f19urnflowrate)="Urine flow rate >= 10 ml/min "
label(data$f19urnflowrate_mv)="Missing value verified for: Urine flow rate >= 10 ml/min "
label(data$f19assaydate)="Date of assay"
label(data$f19assaydate_mv)="Missing value verified for: Date of assay"
label(data$f19urnioth3qc)="Urine Ioth. 3 (mg/dl) QC Result"
label(data$f19urnioth3qc_mv)="Missing value verified for: Urine Ioth. 3 (mg/dl) QC Result"
label(data$f19urnioth3df)="Urine Ioth. 3 (mg/dl) Dilution Factor"
label(data$f19urnioth3df_mv)="Missing value verified for: Urine Ioth. 3 (mg/dl) Dilution Factor"
label(data$f19serumioth3qc)="Serum Ioth. 3 (mg/dl) QC Result"
label(data$f19serumioth3qc_mv)="Missing value verified for: Serum Ioth. 3 (mg/dl) QC Result"
label(data$f19serumioth3df)="Serum Ioth. 3 (mg/dl) Dilution Factor"
label(data$f19serumioth3df_mv)="Missing value verified for: Serum Ioth. 3 (mg/dl) Dilution Factor"
label(data$f19serumioth4qc)="Serum Ioth. 4 (mg/dl) QC Result"
label(data$f19serumioth4qc_mv)="Missing value verified for: Serum Ioth. 4 (mg/dl) QC Result"
label(data$f19serumioth4df)="Serum Ioth. 4 (mg/dl) Dilution Factor"
label(data$f19serumioth4df_mv)="Missing value verified for: Serum Ioth. 4 (mg/dl) Dilution Factor"
label(data$f19urnpah3qc)="Urine PAH 3 (mg/dl) QC Result"
label(data$f19urnpah3qc_mv)="Missing value verified for: Urine PAH 3 (mg/dl) QC Result"
label(data$f19urnpah3df)="Urine PAH 3 (mg/dl) Dilution Factor"
label(data$f19urnpah3df_mv)="Missing value verified for: Urine PAH 3 (mg/dl) Dilution Factor"
label(data$f19serumpah3qc)="Serum PAH 3 (mg/dl) QC Result"
label(data$f19serumpah3qc_mv)="Missing value verified for: Serum PAH 3 (mg/dl) QC Result"
label(data$f19serumpah3df)="Serum PAH 3 (mg/dl) Dilution Factor"
label(data$f19serumpah3df_mv)="Missing value verified for: Serum PAH 3 (mg/dl) Dilution Factor"
label(data$f19serumpah4qc)="Serum PAH 4 (mg/dl) QC Result"
label(data$f19serumpah4qc_mv)="Missing value verified for: Serum PAH 4 (mg/dl) QC Result"
label(data$f19serumpah4df)="Serum PAH 4 (mg/dl) Dilution Factor"
label(data$f19serumpah4df_mv)="Missing value verified for: Serum PAH 4 (mg/dl) Dilution Factor"
label(data$renal_function_quality_control_results_form_19_complete)="Complete?"
label(data$f13dtlastvst)="Date of Last Visit"
label(data$f13dtlastvst_mv)="Missing value verified for: Date of Last Visit"
label(data$f13typevist)="Type of last visit?"
label(data$f13typevist_mv)="Missing value verified for: Type of last visit?"
label(data$f13tstintrvl)="Test interval of last visit?"
label(data$f13tstintrvl_mv)="Missing value verified for: Test interval of last visit?"
label(data$f13endstgrnl)="End stage renal disease"
label(data$f13endstgrnl_mv)="Missing value verified for: End stage renal disease"
label(data$f13esrddate)="Date of diagnosis"
label(data$f13esrddate_mv)="Missing value verified for: Date of diagnosis"
label(data$f13nondiabkid)="Non-diabetic kidney disease"
label(data$f13nondiabkid_mv)="Missing value verified for: Non-diabetic kidney disease"
label(data$f13ndkddate)="Date of diagnosis"
label(data$f13ndkddate_mv)="Missing value verified for: Date of diagnosis"
label(data$f13imprdbldr)="Impaired bladder functioning"
label(data$f13imprdbldr_mv)="Missing value verified for: Impaired bladder functioning"
label(data$f13ibfdate)="Date of diagnosis"
label(data$f13ibfdate_mv)="Missing value verified for: Date of diagnosis"
label(data$f13cngstvhrt)="Congestive heart failure"
label(data$f13cngstvhrt_mv)="Missing value verified for: Congestive heart failure"
label(data$f13chfdate)="Date of diagnosis"
label(data$f13chfdate_mv)="Missing value verified for: Date of diagnosis"
label(data$f13ascites)="Ascites"
label(data$f13ascites_mv)="Missing value verified for: Ascites"
label(data$f13ascitesdt)="Date of diagnosis"
label(data$f13ascitesdt_mv)="Missing value verified for: Date of diagnosis"
label(data$f13potreactntomed)="Potential reactions to medication"
label(data$f13potreactntomed_mv)="Missing value verified for: Potential reactions to medication"
label(data$f13potreactntomeddate)="Date of reaction to medication?"
label(data$f13potreactntomeddate_mv)="Missing value verified for: Date of reaction to medication?"
label(data$f13otherdiag)="Other diagnosis"
label(data$f13otherdiag_mv)="Missing value verified for: Other diagnosis"
label(data$f13othrspecfy)="Other (specify)"
label(data$f13othrspecfy_mv)="Missing value verified for: Other (specify)"
label(data$f13othrdiagdt)="Other diagnosis date"
label(data$f13othrdiagdt_mv)="Missing value verified for: Other diagnosis date"
label(data$f13prfmwdclrnc)="Should the withdrawal clearance study be performed?"
label(data$f13prfmwdclrnc_mv)="Missing value verified for: Should the withdrawal clearance study be performed?"
label(data$f13prfmwdclrncdate)="Date withdrawal clearance scheduled?"
label(data$f13prfmwdclrncdate_mv)="Missing value verified for: Date withdrawal clearance scheduled?"
label(data$f13stoppntdt)="Date stop point declared?"
label(data$f13stoppntdt_mv)="Missing value verified for: Date stop point declared?"
label(data$f13dtcompltd)="Date form completed?"
label(data$f13dtcompltd_mv)="Missing value verified for: Date form completed?"
label(data$f13forminit)="Initials of person who completed form"
label(data$f13forminit_mv)="Missing value verified for: Initials of person who completed form"
label(data$stop_point_form_13_complete)="Complete?"
label(data$f14dob)="Date of birth"
label(data$f14dob_mv)="Missing value verified for: Date of birth"
label(data$f14dod)="Date of death"
label(data$f14dod_mv)="Missing value verified for: Date of death"
label(data$f14codprimary)="Cause of death (primary) - ICD code"
label(data$f14codprimary_mv)="Missing value verified for: Cause of death (primary) - ICD code"
label(data$f14codunerlying)="Cause of death (underlying) - ICD code"
label(data$f14codunerlying_mv)="Missing value verified for: Cause of death (underlying) - ICD code"
label(data$f14coddrprimary)="Cause of death (DR review primary) - ICD code"
label(data$f14coddrprimary_mv)="Missing value verified for: Cause of death (DR review primary) - ICD code"
label(data$f14coddrunerlying)="Cause of death (DR review underlying) - ICD code"
label(data$f14coddrunerlying_mv)="Missing value verified for: Cause of death (DR review underlying) - ICD code"
label(data$f14physiciancr)="Physician performing chart review"
label(data$f14physiciancr_mv)="Missing value verified for: Physician performing chart review"
label(data$f14crphysicianothr)="Physician performing chart review (other specify)"
label(data$f14crphysicianothr_mv)="Missing value verified for: Physician performing chart review (other specify)"
label(data$f14autopsyprfrmd)="Has autopsy been performed"
label(data$f14autopsyprfrmd_mv)="Missing value verified for: Has autopsy been performed"
label(data$f14autopsyprimarycause)="Result of autopsy (primary cause of death) - ICD code"
label(data$f14autopsyprimarycause_mv)="Missing value verified for: Result of autopsy (primary cause of death) - ICD code"
label(data$f14autopsyunderlyingcause)="Result of autopsy (underlying cause of death) - ICD code"
label(data$f14autopsyunderlyingcause_mv)="Missing value verified for: Result of autopsy (underlying cause of death) - ICD code"
label(data$f14deathloc)="Location of death"
label(data$f14deathloc_mv)="Missing value verified for: Location of death"
label(data$f4lodother)="Location of death, other (specify)"
label(data$f4lodother_mv)="Missing value verified for: Location of death, other (specify)"
label(data$f14datecompleted)="Date form completed"
label(data$f14datecompleted_mv)="Missing value verified for: Date form completed"
label(data$f14initials)="Initials of person who completed form"
label(data$f14initials_mv)="Missing value verified for: Initials of person who completed form"
label(data$death_notice_form_14_complete)="Complete?"
label(data$phxlabs_visitdate)="Visit Date"
label(data$phxlabs_visitdate_mv)="Missing value verified for: Visit Date"
label(data$scr)="Serum Creatinine"
label(data$scr_mv)="Missing value verified for: Serum Creatinine"
label(data$gfr)="Glomerular Filtration Rate (GFR)"
label(data$gfr_mv)="Missing value verified for: Glomerular Filtration Rate (GFR)"
label(data$hba1a)="Hemoglobin A1a"
label(data$hba1a_mv)="Missing value verified for: Hemoglobin A1a"
label(data$hba1b)="Hemoglobin A1b"
label(data$hba1b_mv)="Missing value verified for: Hemoglobin A1b"
label(data$hba1c)="Hemoglobin A1c"
label(data$hba1c_mv)="Missing value verified for: Hemoglobin A1c"
label(data$hba1o)="Hemoglobin Ao"
label(data$hba1o_mv)="Missing value verified for: Hemoglobin Ao"
label(data$hbf)="Hemoglobin F"
label(data$hbf_mv)="Missing value verified for: Hemoglobin F"
label(data$p0albumin)="P0 Albumin"
label(data$p0albumin_mv)="Missing value verified for: P0 Albumin"
label(data$p0gluc)="P0 Glucose"
label(data$p0gluc_mv)="Missing value verified for: P0 Glucose"
label(data$p0igg)="P0 IgG"
label(data$p0igg_mv)="Missing value verified for: P0 IgG"
label(data$p0oncoticprsr)="P0 Oncotic Pressure"
label(data$p0oncoticprsr_mv)="Missing value verified for: P0 Oncotic Pressure"
label(data$p3albumin)="P3 Albumin"
label(data$p3albumin_mv)="Missing value verified for: P3 Albumin"
label(data$p3igg)="P3 IgG"
label(data$p3igg_mv)="Missing value verified for: P3 IgG"
label(data$p3oncoticprsr)="P3 Oncotic Pressure"
label(data$p3oncoticprsr_mv)="Missing value verified for: P3 Oncotic Pressure"
label(data$pahclrnc)="PAH Clearance"
label(data$pahclrnc_mv)="Missing value verified for: PAH Clearance"
label(data$u0igg)="U0 IgG"
label(data$u0igg_mv)="Missing value verified for: U0 IgG"
label(data$u0igg_is_below_limit)="U0 IgG is below detection limit of assay?"
label(data$u0igg_is_below_limit_mv)="Missing value verified for: U0 IgG is below detection limit of assay?"
label(data$u3flow_gfr)="U3 Flow"
label(data$u3flow_gfr_mv)="Missing value verified for: U3 Flow"
label(data$u3gfr)="U3 GFR"
label(data$u3gfr_mv)="Missing value verified for: U3 GFR"
label(data$u3use_gfr)="U3 Use"
label(data$u3use_gfr_mv)="Missing value verified for: U3 Use"
label(data$u3albumin)="U3 Albumin"
label(data$u3albumin_mv)="Missing value verified for: U3 Albumin"
label(data$u3albumin_is_below_limit)="U3 Albumin is below detection limit of assay?"
label(data$u3albumin_is_below_limit_mv)="Missing value verified for: U3 Albumin is below detection limit of assay?"
label(data$u3igg)="U3 IgG"
label(data$u3igg_mv)="Missing value verified for: U3 IgG"
label(data$u3igg_is_below_limit)="U3 IgG is below detection limit of assay?"
label(data$u3igg_is_below_limit_mv)="Missing value verified for: U3 IgG is below detection limit of assay?"
label(data$ualb)="Urine Albumin"
label(data$ualb_mv)="Missing value verified for: Urine Albumin"
label(data$ualb_is_below_limit)="Urine Albumin is below detection limit of assay?"
label(data$ualb_is_below_limit_mv)="Missing value verified for: Urine Albumin is below detection limit of assay?"
label(data$ucr)="Urine Creatinine"
label(data$ucr_mv)="Missing value verified for: Urine Creatinine"
label(data$phoenix_labs_complete)="Complete?"
label(data$stflabs_visitdate)="Visit date"
label(data$stflabs_visitdate_mv)="Missing value verified for: Visit date"
label(data$s_gfr)="Glomerular Filtration Rate (GFR)"
label(data$s_gfr_mv)="Missing value verified for: Glomerular Filtration Rate (GFR)"
label(data$s_ff)="I/P RATIO GFR"
label(data$s_ff_mv)="Missing value verified for: I/P RATIO GFR"
label(data$s_pahclrnc)="PAH Clearance"
label(data$s_pahclrnc_mv)="Missing value verified for: PAH Clearance"
label(data$s_gfrb)="GFR Adjusted for BSA "
label(data$s_gfrb_mv)="Missing value verified for: GFR Adjusted for BSA "
label(data$s_pahb)="PAH Adjusted for BSA"
label(data$s_pahb_mv)="Missing value verified for: PAH Adjusted for BSA"
label(data$s_rpf)="Renal Plasma Flow"
label(data$s_rpf_mv)="Missing value verified for: Renal Plasma Flow"
label(data$s_rpfb)="RPF Adjusted for BSA"
label(data$s_rpfb_mv)="Missing value verified for: RPF Adjusted for BSA"
label(data$stanford_labs_complete)="Complete?"
label(data$biopsy_dt)="Biopsy date"
label(data$biopsy_dt_mv)="Missing value verified for: Biopsy date"
label(data$biopsy_source)="Source"
label(data$biopsy_source_mv)="Missing value verified for: Source"
label(data$glomerul)="Glomerul"
label(data$glomerul_mv)="Missing value verified for: Glomerul"
label(data$pct_glob)="PCT Glob"
label(data$pct_glob_mv)="Missing value verified for: PCT Glob"
label(data$mean_fia)="Mean Fia"
label(data$mean_fia_mv)="Missing value verified for: Mean Fia"
label(data$mean_nv)="Mean Nv"
label(data$mean_nv_mv)="Missing value verified for: Mean Nv"
label(data$mean_nv0)="Mean Nv0"
label(data$mean_nv0_mv)="Missing value verified for: Mean Nv0"
label(data$mean_nv1)="Mean Nv1"
label(data$mean_nv1_mv)="Missing value verified for: Mean Nv1"
label(data$mean_n_p)="Mean N P"
label(data$mean_n_p_mv)="Missing value verified for: Mean N P"
label(data$mean_n_e)="Mean N E"
label(data$mean_n_e_mv)="Missing value verified for: Mean N E"
label(data$mean_aa)="Mean AA"
label(data$mean_aa_mv)="Missing value verified for: Mean AA"
label(data$mean_sv)="Mean SV"
label(data$mean_sv_mv)="Missing value verified for: Mean SV"
label(data$mean_bmt)="Mean BMT"
label(data$mean_bmt_mv)="Missing value verified for: Mean BMT"
label(data$mean_fsf)="Mean FSF"
label(data$mean_fsf_mv)="Missing value verified for: Mean FSF"
label(data$mean_fpw)="Mean FPW"
label(data$mean_fpw_mv)="Missing value verified for: Mean FPW"
label(data$pct_denu)="Pct Denu"
label(data$pct_denu_mv)="Missing value verified for: Pct Denu"
label(data$pct_fene)="Pct Fene"
label(data$pct_fene_mv)="Missing value verified for: Pct Fene"
label(data$pct_fen)="Pct Fen"
label(data$pct_fen_mv)="Missing value verified for: Pct Fen"
label(data$sa_using)="Sa Using"
label(data$sa_using_mv)="Missing value verified for: Sa Using"
label(data$number_o)="Number O"
label(data$number_o_mv)="Missing value verified for: Number O"
label(data$alb_stat)="Alb Stat"
label(data$alb_stat_mv)="Missing value verified for: Alb Stat"
label(data$comments)="Comments"
label(data$comments_mv)="Missing value verified for: Comments"
label(data$biopsy_complete)="Complete?"
#Setting Units


#Setting Factors(will create new variable for factors)
data$redcap_event_name.factor = factor(data$redcap_event_name,levels=c("ss_interval_3_arm_1","ss_interval_2_arm_1","ss_interval_1_arm_1","rc_interval_0_arm_1","rv_interval_07_arm_1","rc_interval_1_arm_1","rv_interval_3_arm_1","rv_interval_6_arm_1","rv_interval_9_arm_1","rc_interval_12_arm_1","rv_interval_15_arm_1","rv_interval_18_arm_1","rv_interval_21_arm_1","rc_interval_24_arm_1","rv_interval_27_arm_1","rv_interval_30_arm_1","rv_interval_33_arm_1","rc_interval_36_arm_1","rv_interval_39_arm_1","rv_interval_42_arm_1","rv_interval_45_arm_1","rc_interval_48_arm_1","rv_interval_51_arm_1","rv_interval_54_arm_1","rv_interval_57_arm_1","rc_interval_60_arm_1","rv_interval_63_arm_1","rv_interval_66_arm_1","rv_interval_69_arm_1","rc_interval_72_arm_1","rv_interval_75_arm_1","rv_interval_78_arm_1","rv_interval_81_arm_1","rc_interval_84_arm_1","rv_interval_87_arm_1","rv_interval_90_arm_1","rv_interval_93_arm_1","rc_interval_96_arm_1","rv_interval_99_arm_1","rv_interval_102_arm_1","rv_interval_105_arm_1","rc_interval_108_arm_1","rv_interval_111_arm_1","rv_interval_114_arm_1","rv_interval_117_arm_1","rc_interval_120_arm_1","rv_interval_123_arm_1","rv_interval_126_arm_1","rv_interval_129_arm_1","rc_interval_132_arm_1","rv_interval_135_arm_1","rv_interval_138_arm_1","rv_interval_141_arm_1","rc_interval_144_arm_1","rv_interval_147_arm_1","rv_interval_150_arm_1","rv_interval_153_arm_1","rc_interval_156_arm_1","rv_interval_159_arm_1","rv_interval_162_arm_1","rv_interval_165_arm_1","rc_interval_168_arm_1","rv_interval_171_arm_1","rv_interval_174_arm_1","rv_interval_177_arm_1","rc_interval_180_arm_1","rv_interval_183_arm_1","rv_interval_186_arm_1","rv_interval_189_arm_1","rc_interval_192_arm_1","rv_interval_195_arm_1","rv_interval_198_arm_1","rv_interval_201_arm_1","rc_interval_204_arm_1"))
data$redcap_repeat_instrument.factor = factor(data$redcap_repeat_instrument,levels=c(""))
data$lastname_mv.factor = factor(data$lastname_mv,levels=c("1","0"))
data$firstname_mv.factor = factor(data$firstname_mv,levels=c("1","0"))
data$stratum.factor = factor(data$stratum,levels=c("1","2"))
data$stratum_mv.factor = factor(data$stratum_mv,levels=c("1","0"))
data$randomnum_mv.factor = factor(data$randomnum_mv,levels=c("1","0"))
data$tx.factor = factor(data$tx,levels=c("0","1"))
data$tx_mv.factor = factor(data$tx_mv,levels=c("1","0"))
data$sacaton_no_mv.factor = factor(data$sacaton_no_mv,levels=c("1","0"))
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
data$closedate_mv.factor = factor(data$closedate_mv,levels=c("1","0"))
data$washoutdate_mv.factor = factor(data$washoutdate_mv,levels=c("1","0"))
data$demoform_complete.factor = factor(data$demoform_complete,levels=c("0","1","2"))
data$f2visitdate_mv.factor = factor(data$f2visitdate_mv,levels=c("1","0"))
data$f2visittype.factor = factor(data$f2visittype,levels=c("4"))
data$f2visittype_mv.factor = factor(data$f2visittype_mv,levels=c("1","0"))
data$f2testintvrl_mv.factor = factor(data$f2testintvrl_mv,levels=c("1","0"))
data$f2bldprsrarm.factor = factor(data$f2bldprsrarm,levels=c("0","1","2"))
data$f2bldprsrarm_mv.factor = factor(data$f2bldprsrarm_mv,levels=c("1","0"))
data$f2frstmsrsys_mv.factor = factor(data$f2frstmsrsys_mv,levels=c("1","0"))
data$f2frstmsrdia_mv.factor = factor(data$f2frstmsrdia_mv,levels=c("1","0"))
data$f2scndmsrsys_mv.factor = factor(data$f2scndmsrsys_mv,levels=c("1","0"))
data$f2scndmsrdia_mv.factor = factor(data$f2scndmsrdia_mv,levels=c("1","0"))
data$f2heartrate_mv.factor = factor(data$f2heartrate_mv,levels=c("1","0"))
data$subject_screening_3_and_2_form_02_complete.factor = factor(data$subject_screening_3_and_2_form_02_complete,levels=c("0","1","2"))
data$f3visitdate_mv.factor = factor(data$f3visitdate_mv,levels=c("1","0"))
data$f3visittype.factor = factor(data$f3visittype,levels=c("3","4","5"))
data$f3visittype_mv.factor = factor(data$f3visittype_mv,levels=c("1","0"))
data$f3testintvrl_mv.factor = factor(data$f3testintvrl_mv,levels=c("1","0"))
data$f3spturncltd.factor = factor(data$f3spturncltd,levels=c("1","2"))
data$f3spturncltd_mv.factor = factor(data$f3spturncltd_mv,levels=c("1","0"))
data$f3venipfork.factor = factor(data$f3venipfork,levels=c("1","2"))
data$f3venipfork_mv.factor = factor(data$f3venipfork_mv,levels=c("1","0"))
data$f3bldprsrarm.factor = factor(data$f3bldprsrarm,levels=c("1","2"))
data$f3bldprsrarm_mv.factor = factor(data$f3bldprsrarm_mv,levels=c("1","0"))
data$f3frstmsrsys_mv.factor = factor(data$f3frstmsrsys_mv,levels=c("1","0"))
data$f3frstmsrdia_mv.factor = factor(data$f3frstmsrdia_mv,levels=c("1","0"))
data$f3scndmsrsys_mv.factor = factor(data$f3scndmsrsys_mv,levels=c("1","0"))
data$f3scndmsrdia_mv.factor = factor(data$f3scndmsrdia_mv,levels=c("1","0"))
data$f3heartrate_mv.factor = factor(data$f3heartrate_mv,levels=c("1","0"))
data$f3prohbmeds.factor = factor(data$f3prohbmeds,levels=c("1","2"))
data$f3prohbmeds_mv.factor = factor(data$f3prohbmeds_mv,levels=c("1","0"))
data$f3elgstatus.factor = factor(data$f3elgstatus,levels=c("1","2"))
data$f3elgstatus_mv.factor = factor(data$f3elgstatus_mv,levels=c("1","0"))
data$f3inelgibleexplain_mv.factor = factor(data$f3inelgibleexplain_mv,levels=c("1","0"))
data$f3medsurghis.factor = factor(data$f3medsurghis,levels=c("1","2"))
data$f3medsurghis_mv.factor = factor(data$f3medsurghis_mv,levels=c("1","0"))
data$f3_acratio.factor = factor(data$f3_acratio,levels=c("1","2"))
data$f3_acratio_mv.factor = factor(data$f3_acratio_mv,levels=c("1","0"))
data$f3srmcreatyn.factor = factor(data$f3srmcreatyn,levels=c("1","2"))
data$f3srmcreatyn_mv.factor = factor(data$f3srmcreatyn_mv,levels=c("1","0"))
data$f3_pregnancy.factor = factor(data$f3_pregnancy,levels=c("1","2"))
data$f3_pregnancy_mv.factor = factor(data$f3_pregnancy_mv,levels=c("1","0"))
data$f3_staffpref.factor = factor(data$f3_staffpref,levels=c("1","2"))
data$f3_staffpref_mv.factor = factor(data$f3_staffpref_mv,levels=c("1","0"))
data$f3_subjpref.factor = factor(data$f3_subjpref,levels=c("1","2"))
data$f3_subjpref_mv.factor = factor(data$f3_subjpref_mv,levels=c("1","0"))
data$f3_other.factor = factor(data$f3_other,levels=c("1","2"))
data$f3_other_mv.factor = factor(data$f3_other_mv,levels=c("1","0"))
data$subject_screeningelgibility_form_03_complete.factor = factor(data$subject_screeningelgibility_form_03_complete,levels=c("0","1","2"))
data$f3meds_visitdate_mv.factor = factor(data$f3meds_visitdate_mv,levels=c("1","0"))
data$f3drugcode1_mv.factor = factor(data$f3drugcode1_mv,levels=c("1","0"))
data$f3dosage1_mv.factor = factor(data$f3dosage1_mv,levels=c("1","0"))
data$f3timesday1_mv.factor = factor(data$f3timesday1_mv,levels=c("1","0"))
data$f3dwm1.factor = factor(data$f3dwm1,levels=c("D","W","M"))
data$f3dwm1_mv.factor = factor(data$f3dwm1_mv,levels=c("1","0"))
data$f3drugcode2_mv.factor = factor(data$f3drugcode2_mv,levels=c("1","0"))
data$f3dosage2_mv.factor = factor(data$f3dosage2_mv,levels=c("1","0"))
data$f3timesday2_mv.factor = factor(data$f3timesday2_mv,levels=c("1","0"))
data$f3dwm2.factor = factor(data$f3dwm2,levels=c("D","W","M"))
data$f3dwm2_mv.factor = factor(data$f3dwm2_mv,levels=c("1","0"))
data$f3drugcode3_mv.factor = factor(data$f3drugcode3_mv,levels=c("1","0"))
data$f3dosage3_mv.factor = factor(data$f3dosage3_mv,levels=c("1","0"))
data$f3timesday3_mv.factor = factor(data$f3timesday3_mv,levels=c("1","0"))
data$f3dwm3.factor = factor(data$f3dwm3,levels=c("D","W","M"))
data$f3dwm3_mv.factor = factor(data$f3dwm3_mv,levels=c("1","0"))
data$f3drugcode4_mv.factor = factor(data$f3drugcode4_mv,levels=c("1","0"))
data$f3dosage4_mv.factor = factor(data$f3dosage4_mv,levels=c("1","0"))
data$f3timesday4_mv.factor = factor(data$f3timesday4_mv,levels=c("1","0"))
data$f3dwm4.factor = factor(data$f3dwm4,levels=c("D","W","M"))
data$f3dwm4_mv.factor = factor(data$f3dwm4_mv,levels=c("1","0"))
data$f3drugcode5_mv.factor = factor(data$f3drugcode5_mv,levels=c("1","0"))
data$f3dosage5_mv.factor = factor(data$f3dosage5_mv,levels=c("1","0"))
data$f3timesday5_mv.factor = factor(data$f3timesday5_mv,levels=c("1","0"))
data$f3dwm5.factor = factor(data$f3dwm5,levels=c("D","W","M"))
data$f3dwm5_mv.factor = factor(data$f3dwm5_mv,levels=c("1","0"))
data$f3drugcode6_mv.factor = factor(data$f3drugcode6_mv,levels=c("1","0"))
data$f3dosage6_mv.factor = factor(data$f3dosage6_mv,levels=c("1","0"))
data$f3timesday6_mv.factor = factor(data$f3timesday6_mv,levels=c("1","0"))
data$f3dwm6.factor = factor(data$f3dwm6,levels=c("D","W","M"))
data$f3dwm6_mv.factor = factor(data$f3dwm6_mv,levels=c("1","0"))
data$f3drugcode7_mv.factor = factor(data$f3drugcode7_mv,levels=c("1","0"))
data$f3dosage7_mv.factor = factor(data$f3dosage7_mv,levels=c("1","0"))
data$f3timesday7_mv.factor = factor(data$f3timesday7_mv,levels=c("1","0"))
data$f3dwm7.factor = factor(data$f3dwm7,levels=c("D","W","M"))
data$f3dwm7_mv.factor = factor(data$f3dwm7_mv,levels=c("1","0"))
data$f3drugcode8_mv.factor = factor(data$f3drugcode8_mv,levels=c("1","0"))
data$f3dosage8_mv.factor = factor(data$f3dosage8_mv,levels=c("1","0"))
data$f3timesday8_mv.factor = factor(data$f3timesday8_mv,levels=c("1","0"))
data$f3dwm8.factor = factor(data$f3dwm8,levels=c("D","W","M"))
data$f3dwm8_mv.factor = factor(data$f3dwm8_mv,levels=c("1","0"))
data$f3drugcode9_mv.factor = factor(data$f3drugcode9_mv,levels=c("1","0"))
data$f3dosage9_mv.factor = factor(data$f3dosage9_mv,levels=c("1","0"))
data$f3timesday9_mv.factor = factor(data$f3timesday9_mv,levels=c("1","0"))
data$f3dwm9.factor = factor(data$f3dwm9,levels=c("D","W","M"))
data$f3dwm9_mv.factor = factor(data$f3dwm9_mv,levels=c("1","0"))
data$f3drugcode10_mv.factor = factor(data$f3drugcode10_mv,levels=c("1","0"))
data$f3dosage10_mv.factor = factor(data$f3dosage10_mv,levels=c("1","0"))
data$f3timesday10_mv.factor = factor(data$f3timesday10_mv,levels=c("1","0"))
data$f3dwm10.factor = factor(data$f3dwm10,levels=c("D","W","M"))
data$f3dwm10_mv.factor = factor(data$f3dwm10_mv,levels=c("1","0"))
data$f3drugcode11_mv.factor = factor(data$f3drugcode11_mv,levels=c("1","0"))
data$f3dosage11_mv.factor = factor(data$f3dosage11_mv,levels=c("1","0"))
data$f3timesday11_mv.factor = factor(data$f3timesday11_mv,levels=c("1","0"))
data$f3dwm11.factor = factor(data$f3dwm11,levels=c("D","W","M"))
data$f3dwm11_mv.factor = factor(data$f3dwm11_mv,levels=c("1","0"))
data$f3drugcode12_mv.factor = factor(data$f3drugcode12_mv,levels=c("1","0"))
data$f3dosage12_mv.factor = factor(data$f3dosage12_mv,levels=c("1","0"))
data$f3timesday12_mv.factor = factor(data$f3timesday12_mv,levels=c("1","0"))
data$f3dwm12.factor = factor(data$f3dwm12,levels=c("D","W","M"))
data$f3dwm12_mv.factor = factor(data$f3dwm12_mv,levels=c("1","0"))
data$f3drugcode13_mv.factor = factor(data$f3drugcode13_mv,levels=c("1","0"))
data$f3dosage13_mv.factor = factor(data$f3dosage13_mv,levels=c("1","0"))
data$f3timesday13_mv.factor = factor(data$f3timesday13_mv,levels=c("1","0"))
data$f3dwm13.factor = factor(data$f3dwm13,levels=c("D","W","M"))
data$f3dwm13_mv.factor = factor(data$f3dwm13_mv,levels=c("1","0"))
data$f3drugcode14_mv.factor = factor(data$f3drugcode14_mv,levels=c("1","0"))
data$f3dosage14_mv.factor = factor(data$f3dosage14_mv,levels=c("1","0"))
data$f3timesday14_mv.factor = factor(data$f3timesday14_mv,levels=c("1","0"))
data$f3dwm14.factor = factor(data$f3dwm14,levels=c("D","W","M"))
data$f3dwm14_mv.factor = factor(data$f3dwm14_mv,levels=c("1","0"))
data$f3drugcode15_mv.factor = factor(data$f3drugcode15_mv,levels=c("1","0"))
data$f3dosage15_mv.factor = factor(data$f3dosage15_mv,levels=c("1","0"))
data$f3timesday15_mv.factor = factor(data$f3timesday15_mv,levels=c("1","0"))
data$f3dwm15.factor = factor(data$f3dwm15,levels=c("D","W","M"))
data$f3dwm15_mv.factor = factor(data$f3dwm15_mv,levels=c("1","0"))
data$f3drugcode16_mv.factor = factor(data$f3drugcode16_mv,levels=c("1","0"))
data$f3dosage16_mv.factor = factor(data$f3dosage16_mv,levels=c("1","0"))
data$f3timesday16_mv.factor = factor(data$f3timesday16_mv,levels=c("1","0"))
data$f3dwm16.factor = factor(data$f3dwm16,levels=c("D","W","M"))
data$f3dwm16_mv.factor = factor(data$f3dwm16_mv,levels=c("1","0"))
data$f3drugcode17_mv.factor = factor(data$f3drugcode17_mv,levels=c("1","0"))
data$f3dosage17_mv.factor = factor(data$f3dosage17_mv,levels=c("1","0"))
data$f3timesday17_mv.factor = factor(data$f3timesday17_mv,levels=c("1","0"))
data$f3dwm17.factor = factor(data$f3dwm17,levels=c("D","W","M"))
data$f3dwm17_mv.factor = factor(data$f3dwm17_mv,levels=c("1","0"))
data$f3drugcode18_mv.factor = factor(data$f3drugcode18_mv,levels=c("1","0"))
data$f3dosage18_mv.factor = factor(data$f3dosage18_mv,levels=c("1","0"))
data$f3timesday18_mv.factor = factor(data$f3timesday18_mv,levels=c("1","0"))
data$f3dwm18.factor = factor(data$f3dwm18,levels=c("D","W","M"))
data$f3dwm18_mv.factor = factor(data$f3dwm18_mv,levels=c("1","0"))
data$f3drugcode19_mv.factor = factor(data$f3drugcode19_mv,levels=c("1","0"))
data$f3dosage19_mv.factor = factor(data$f3dosage19_mv,levels=c("1","0"))
data$f3timesday19_mv.factor = factor(data$f3timesday19_mv,levels=c("1","0"))
data$f3dwm19.factor = factor(data$f3dwm19,levels=c("D","W","M"))
data$f3dwm19_mv.factor = factor(data$f3dwm19_mv,levels=c("1","0"))
data$f3drugcode20_mv.factor = factor(data$f3drugcode20_mv,levels=c("1","0"))
data$f3dosage20_mv.factor = factor(data$f3dosage20_mv,levels=c("1","0"))
data$f3timesday20_mv.factor = factor(data$f3timesday20_mv,levels=c("1","0"))
data$f3dwm20.factor = factor(data$f3dwm20,levels=c("D","W","M"))
data$f3dwm20_mv.factor = factor(data$f3dwm20_mv,levels=c("1","0"))
data$subject_screeningelgibility_meds_form_03_complete.factor = factor(data$subject_screeningelgibility_meds_form_03_complete,levels=c("0","1","2"))
data$clrncvisitdate_mv.factor = factor(data$clrncvisitdate_mv,levels=c("1","0"))
data$clrncvisittype.factor = factor(data$clrncvisittype,levels=c("3","5"))
data$clrncvisittype_mv.factor = factor(data$clrncvisittype_mv,levels=c("1","0"))
data$clrnctestintvrl_mv.factor = factor(data$clrnctestintvrl_mv,levels=c("1","0"))
data$renal_clearance_visit_complete.factor = factor(data$renal_clearance_visit_complete,levels=c("0","1","2"))
data$intfu_visitdate_mv.factor = factor(data$intfu_visitdate_mv,levels=c("1","0"))
data$intfu_visittype.factor = factor(data$intfu_visittype,levels=c("3","5"))
data$intfu_visittype_mv.factor = factor(data$intfu_visittype_mv,levels=c("1","0"))
data$intfu_testintrvl_mv.factor = factor(data$intfu_testintrvl_mv,levels=c("1","0"))
data$interim_followup_visit_complete.factor = factor(data$interim_followup_visit_complete,levels=c("0","1","2"))
data$f4visitdate_mv.factor = factor(data$f4visitdate_mv,levels=c("1","0"))
data$f4visittype_mv.factor = factor(data$f4visittype_mv,levels=c("1","0"))
data$f4testintvl_mv.factor = factor(data$f4testintvl_mv,levels=c("1","0"))
data$f4bldprsr.factor = factor(data$f4bldprsr,levels=c("0","1","2"))
data$f4bldprsr_mv.factor = factor(data$f4bldprsr_mv,levels=c("1","0"))
data$f4frstmsrsys_mv.factor = factor(data$f4frstmsrsys_mv,levels=c("1","0"))
data$f4frstmsrdia_mv.factor = factor(data$f4frstmsrdia_mv,levels=c("1","0"))
data$f4scndmsrsys_mv.factor = factor(data$f4scndmsrsys_mv,levels=c("1","0"))
data$f4scndmsrdia_mv.factor = factor(data$f4scndmsrdia_mv,levels=c("1","0"))
data$f4heartrate_mv.factor = factor(data$f4heartrate_mv,levels=c("1","0"))
data$f4pregnant.factor = factor(data$f4pregnant,levels=c("0","1","2"))
data$f4pregnant_mv.factor = factor(data$f4pregnant_mv,levels=c("1","0"))
data$f4medsnow.factor = factor(data$f4medsnow,levels=c("1","2"))
data$f4medsnow_mv.factor = factor(data$f4medsnow_mv,levels=c("1","0"))
data$f4srmpotasm_mv.factor = factor(data$f4srmpotasm_mv,levels=c("1","0"))
data$f4srmglu_mv.factor = factor(data$f4srmglu_mv,levels=c("1","0"))
data$f4forminit_mv.factor = factor(data$f4forminit_mv,levels=c("1","0"))
data$blood_pressure_form_04_complete.factor = factor(data$blood_pressure_form_04_complete,levels=c("0","1","2"))
data$f4_visitdate_mv.factor = factor(data$f4_visitdate_mv,levels=c("1","0"))
data$f4_pillbotldt_mv.factor = factor(data$f4_pillbotldt_mv,levels=c("1","0"))
data$f4_dsplstpill_mv.factor = factor(data$f4_dsplstpill_mv,levels=c("1","0"))
data$f4_dysbtwndsp_mv.factor = factor(data$f4_dysbtwndsp_mv,levels=c("1","0"))
data$f4_nodlydose_mv.factor = factor(data$f4_nodlydose_mv,levels=c("1","0"))
data$f4_dlydosrtrn_mv.factor = factor(data$f4_dlydosrtrn_mv,levels=c("1","0"))
data$f4_dlydosetkn_mv.factor = factor(data$f4_dlydosetkn_mv,levels=c("1","0"))
data$f4_compliance_mv.factor = factor(data$f4_compliance_mv,levels=c("1","0"))
data$f4_totdlydosn_mv.factor = factor(data$f4_totdlydosn_mv,levels=c("1","0"))
data$f4_curdos_mv.factor = factor(data$f4_curdos_mv,levels=c("1","0"))
data$f4_dateadjstd_mv.factor = factor(data$f4_dateadjstd_mv,levels=c("1","0"))
data$f4_ajstddose_mv.factor = factor(data$f4_ajstddose_mv,levels=c("1","0"))
data$f4_medseffects.factor = factor(data$f4_medseffects,levels=c("0","1","2","3","4","5","6","7","8","9","10","11"))
data$f4_medseffects_mv.factor = factor(data$f4_medseffects_mv,levels=c("1","0"))
data$f4_expeffects_mv.factor = factor(data$f4_expeffects_mv,levels=c("1","0"))
data$drug_compliance_form_04_complete.factor = factor(data$drug_compliance_form_04_complete,levels=c("0","1","2"))
data$f4meds_visitdate_mv.factor = factor(data$f4meds_visitdate_mv,levels=c("1","0"))
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
data$meds_form_04_complete.factor = factor(data$meds_form_04_complete,levels=c("0","1","2"))
data$f5visitdate_mv.factor = factor(data$f5visitdate_mv,levels=c("1","0"))
data$f5rcntillns.factor = factor(data$f5rcntillns,levels=c("1","2"))
data$f5rcntillns_mv.factor = factor(data$f5rcntillns_mv,levels=c("1","0"))
data$f5outpat.factor = factor(data$f5outpat,levels=c("0","1","2"))
data$f5outpat_mv.factor = factor(data$f5outpat_mv,levels=c("1","0"))
data$f5outpatspec_mv.factor = factor(data$f5outpatspec_mv,levels=c("1","0"))
data$f5hospwosrg.factor = factor(data$f5hospwosrg,levels=c("0","1","2"))
data$f5hospwosrg_mv.factor = factor(data$f5hospwosrg_mv,levels=c("1","0"))
data$f5hospwoexp_mv.factor = factor(data$f5hospwoexp_mv,levels=c("1","0"))
data$f5hospwsrg.factor = factor(data$f5hospwsrg,levels=c("0","1","2"))
data$f5hospwsrg_mv.factor = factor(data$f5hospwsrg_mv,levels=c("1","0"))
data$f5hospwexp_mv.factor = factor(data$f5hospwexp_mv,levels=c("1","0"))
data$f5medprbsdoc.factor = factor(data$f5medprbsdoc,levels=c("0","1","2"))
data$f5medprbsdoc_mv.factor = factor(data$f5medprbsdoc_mv,levels=c("1","0"))
data$f5crnryartds.factor = factor(data$f5crnryartds,levels=c("0","1","2"))
data$f5crnryartds_mv.factor = factor(data$f5crnryartds_mv,levels=c("1","0"))
data$f5cancer.factor = factor(data$f5cancer,levels=c("0","1","2"))
data$f5cancer_mv.factor = factor(data$f5cancer_mv,levels=c("1","0"))
data$f5cancerexp_mv.factor = factor(data$f5cancerexp_mv,levels=c("1","0"))
data$f5crbrlvasds.factor = factor(data$f5crbrlvasds,levels=c("0","1","2"))
data$f5crbrlvasds_mv.factor = factor(data$f5crbrlvasds_mv,levels=c("1","0"))
data$f5peripvasds.factor = factor(data$f5peripvasds,levels=c("0","1","2"))
data$f5peripvasds_mv.factor = factor(data$f5peripvasds_mv,levels=c("1","0"))
data$f5hypertensn.factor = factor(data$f5hypertensn,levels=c("0","1","2"))
data$f5hypertensn_mv.factor = factor(data$f5hypertensn_mv,levels=c("1","0"))
data$f5seizures.factor = factor(data$f5seizures,levels=c("0","1","2"))
data$f5seizures_mv.factor = factor(data$f5seizures_mv,levels=c("1","0"))
data$f5gnitrnryds.factor = factor(data$f5gnitrnryds,levels=c("0","1","2"))
data$f5gnitrnryds_mv.factor = factor(data$f5gnitrnryds_mv,levels=c("1","0"))
data$f5gnitrnyexp_mv.factor = factor(data$f5gnitrnyexp_mv,levels=c("1","0"))
data$f5lungdsz.factor = factor(data$f5lungdsz,levels=c("0","1","2"))
data$f5lungdsz_mv.factor = factor(data$f5lungdsz_mv,levels=c("1","0"))
data$f5majsurg.factor = factor(data$f5majsurg,levels=c("0","1","2"))
data$f5majsurg_mv.factor = factor(data$f5majsurg_mv,levels=c("1","0"))
data$f5majsurgexp_mv.factor = factor(data$f5majsurgexp_mv,levels=c("1","0"))
data$f5othrmeddia.factor = factor(data$f5othrmeddia,levels=c("0","1","2"))
data$f5othrmeddia_mv.factor = factor(data$f5othrmeddia_mv,levels=c("1","0"))
data$f5othrmedexp_mv.factor = factor(data$f5othrmedexp_mv,levels=c("1","0"))
data$f5lungs.factor = factor(data$f5lungs,levels=c("0","1","2"))
data$f5lungs_mv.factor = factor(data$f5lungs_mv,levels=c("1","0"))
data$f5lungsexp_mv.factor = factor(data$f5lungsexp_mv,levels=c("1","0"))
data$f5heart.factor = factor(data$f5heart,levels=c("0","1","2"))
data$f5heart_mv.factor = factor(data$f5heart_mv,levels=c("1","0"))
data$f5heartexp_mv.factor = factor(data$f5heartexp_mv,levels=c("1","0"))
data$f5skin.factor = factor(data$f5skin,levels=c("0","1","2"))
data$f5skin_mv.factor = factor(data$f5skin_mv,levels=c("1","0"))
data$f5skinexp_mv.factor = factor(data$f5skinexp_mv,levels=c("1","0"))
data$f5edema.factor = factor(data$f5edema,levels=c("0","1","2"))
data$f5edema_mv.factor = factor(data$f5edema_mv,levels=c("1","0"))
data$f5edemaexp_mv.factor = factor(data$f5edemaexp_mv,levels=c("1","0"))
data$history_and_physical_examination_form_05_complete.factor = factor(data$history_and_physical_examination_form_05_complete,levels=c("0","1","2"))
data$f6visitdate_mv.factor = factor(data$f6visitdate_mv,levels=c("1","0"))
data$f6ivsites_mv.factor = factor(data$f6ivsites_mv,levels=c("1","0"))
data$f6weight_mv.factor = factor(data$f6weight_mv,levels=c("1","0"))
data$f6bsa_mv.factor = factor(data$f6bsa_mv,levels=c("1","0"))
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
data$f6glucose.factor = factor(data$f6glucose,levels=c("0","T","t","1","2","3","4"))
data$f6glucose_mv.factor = factor(data$f6glucose_mv,levels=c("1","0"))
data$f6bilirubin.factor = factor(data$f6bilirubin,levels=c("0","1","2","3"))
data$f6bilirubin_mv.factor = factor(data$f6bilirubin_mv,levels=c("1","0"))
data$f6ketones.factor = factor(data$f6ketones,levels=c("0","T","t","1","2","3","4"))
data$f6ketones_mv.factor = factor(data$f6ketones_mv,levels=c("1","0"))
data$f6specgrvty_mv.factor = factor(data$f6specgrvty_mv,levels=c("1","0"))
data$f6bloodocult.factor = factor(data$f6bloodocult,levels=c("1","2"))
data$f6bloodocult_mv.factor = factor(data$f6bloodocult_mv,levels=c("1","0"))
data$f6ph_mv.factor = factor(data$f6ph_mv,levels=c("1","0"))
data$f6protein.factor = factor(data$f6protein,levels=c("0","T","t","1","2","3","4"))
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
data$renal_clearance_data_form_06_complete.factor = factor(data$renal_clearance_data_form_06_complete,levels=c("0","1","2"))
data$f7qcnumber_mv.factor = factor(data$f7qcnumber_mv,levels=c("1","0"))
data$f7visitdate_mv.factor = factor(data$f7visitdate_mv,levels=c("1","0"))
data$f7tstintrvl_mv.factor = factor(data$f7tstintrvl_mv,levels=c("1","0"))
data$f7typesample.factor = factor(data$f7typesample,levels=c("1","2","3","4","5","6"))
data$f7typesample_mv.factor = factor(data$f7typesample_mv,levels=c("1","0"))
data$laboratory_quality_control_form_07_complete.factor = factor(data$laboratory_quality_control_form_07_complete,levels=c("0","1","2"))
data$f8visitdate_mv.factor = factor(data$f8visitdate_mv,levels=c("1","0"))
data$f8analdate_mv.factor = factor(data$f8analdate_mv,levels=c("1","0"))
data$f8alti_mv.factor = factor(data$f8alti_mv,levels=c("1","0"))
data$f8glucose_mv.factor = factor(data$f8glucose_mv,levels=c("1","0"))
data$f8bun_mv.factor = factor(data$f8bun_mv,levels=c("1","0"))
data$f8creatinine_mv.factor = factor(data$f8creatinine_mv,levels=c("1","0"))
data$f8sodium_mv.factor = factor(data$f8sodium_mv,levels=c("1","0"))
data$f8potassium_mv.factor = factor(data$f8potassium_mv,levels=c("1","0"))
data$f8chloride_mv.factor = factor(data$f8chloride_mv,levels=c("1","0"))
data$f8calcium_mv.factor = factor(data$f8calcium_mv,levels=c("1","0"))
data$f8phosphorus_mv.factor = factor(data$f8phosphorus_mv,levels=c("1","0"))
data$f8totlprotn_mv.factor = factor(data$f8totlprotn_mv,levels=c("1","0"))
data$f8albumin_mv.factor = factor(data$f8albumin_mv,levels=c("1","0"))
data$f8cholestero_mv.factor = factor(data$f8cholestero_mv,levels=c("1","0"))
data$f8triglyceri_mv.factor = factor(data$f8triglyceri_mv,levels=c("1","0"))
data$f8hdl_mv.factor = factor(data$f8hdl_mv,levels=c("1","0"))
data$f8ldl_mv.factor = factor(data$f8ldl_mv,levels=c("1","0"))
data$f8totlbilrbn_mv.factor = factor(data$f8totlbilrbn_mv,levels=c("1","0"))
data$f8dbilirubin_mv.factor = factor(data$f8dbilirubin_mv,levels=c("1","0"))
data$f8alkp04_mv.factor = factor(data$f8alkp04_mv,levels=c("1","0"))
data$f8ast_mv.factor = factor(data$f8ast_mv,levels=c("1","0"))
data$smacchemistry_panel_form_08_complete.factor = factor(data$smacchemistry_panel_form_08_complete,levels=c("0","1","2"))
data$f9visitdate_mv.factor = factor(data$f9visitdate_mv,levels=c("1","0"))
data$f9dateanal_mv.factor = factor(data$f9dateanal_mv,levels=c("1","0"))
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
data$cbc_form_09_complete.factor = factor(data$cbc_form_09_complete,levels=c("0","1","2"))
data$f10visitdate_mv.factor = factor(data$f10visitdate_mv,levels=c("1","0"))
data$f10spturalbu_mv.factor = factor(data$f10spturalbu_mv,levels=c("1","0"))
data$f10spturncrea_mv.factor = factor(data$f10spturncrea_mv,levels=c("1","0"))
data$f10serumcreat_mv.factor = factor(data$f10serumcreat_mv,levels=c("1","0"))
data$f10strgbxnum_mv.factor = factor(data$f10strgbxnum_mv,levels=c("1","0"))
data$daes_routine_labs_form_10_complete.factor = factor(data$daes_routine_labs_form_10_complete,levels=c("0","1","2"))
data$f11visitdate_mv.factor = factor(data$f11visitdate_mv,levels=c("1","0"))
data$f11hba1c_mv.factor = factor(data$f11hba1c_mv,levels=c("1","0"))
data$f11p0glu_mv.factor = factor(data$f11p0glu_mv,levels=c("1","0"))
data$f11p0creatn_mv.factor = factor(data$f11p0creatn_mv,levels=c("1","0"))
data$f11p0albumn_mv.factor = factor(data$f11p0albumn_mv,levels=c("1","0"))
data$f11p3albumn_mv.factor = factor(data$f11p3albumn_mv,levels=c("1","0"))
data$f11p0iggsamp_mv.factor = factor(data$f11p0iggsamp_mv,levels=c("1","0"))
data$f11p3iggsamp_mv.factor = factor(data$f11p3iggsamp_mv,levels=c("1","0"))
data$f11u0creatn_mv.factor = factor(data$f11u0creatn_mv,levels=c("1","0"))
data$f11u0albumn_mv.factor = factor(data$f11u0albumn_mv,levels=c("1","0"))
data$f11u3albumn_mv.factor = factor(data$f11u3albumn_mv,levels=c("1","0"))
data$f11u0igg_mv.factor = factor(data$f11u0igg_mv,levels=c("1","0"))
data$f11u3igg_mv.factor = factor(data$f11u3igg_mv,levels=c("1","0"))
data$daes_lab_clearance_results_form_11_complete.factor = factor(data$daes_lab_clearance_results_form_11_complete,levels=c("0","1","2"))
data$f12expdatevis_mv.factor = factor(data$f12expdatevis_mv,levels=c("1","0"))
data$f12typevisit.factor = factor(data$f12typevisit,levels=c("3","5"))
data$f12typevisit_mv.factor = factor(data$f12typevisit_mv,levels=c("1","0"))
data$f12testintrvl_mv.factor = factor(data$f12testintrvl_mv,levels=c("1","0"))
data$f12reasnmissd.factor = factor(data$f12reasnmissd,levels=c("1","2","3","4","5","6","7","8","9"))
data$f12reasnmissd_mv.factor = factor(data$f12reasnmissd_mv,levels=c("1","0"))
data$f12otherspec_mv.factor = factor(data$f12otherspec_mv,levels=c("1","0"))
data$missed_visitintercurrent_event_form_12_complete.factor = factor(data$missed_visitintercurrent_event_form_12_complete,levels=c("0","1","2"))
data$f16visitdate_mv.factor = factor(data$f16visitdate_mv,levels=c("1","0"))
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
data$smac_quality_control_form_16_complete.factor = factor(data$smac_quality_control_form_16_complete,levels=c("0","1","2"))
data$f17visitdate_mv.factor = factor(data$f17visitdate_mv,levels=c("1","0"))
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
data$cbc_quality_control_form_17_complete.factor = factor(data$cbc_quality_control_form_17_complete,levels=c("0","1","2"))
data$f18visitdate_mv.factor = factor(data$f18visitdate_mv,levels=c("1","0"))
data$f18qcnumber_mv.factor = factor(data$f18qcnumber_mv,levels=c("1","0"))
data$f18hba1c_mv.factor = factor(data$f18hba1c_mv,levels=c("1","0"))
data$f18p0creatn_mv.factor = factor(data$f18p0creatn_mv,levels=c("1","0"))
data$f18p0albumin_mv.factor = factor(data$f18p0albumin_mv,levels=c("1","0"))
data$f18p0igg_mv.factor = factor(data$f18p0igg_mv,levels=c("1","0"))
data$f18u0creatn_mv.factor = factor(data$f18u0creatn_mv,levels=c("1","0"))
data$f18u0albumin_mv.factor = factor(data$f18u0albumin_mv,levels=c("1","0"))
data$f18u0igg_mv.factor = factor(data$f18u0igg_mv,levels=c("1","0"))
data$daes_clearance_quality_control_results_form_18_complete.factor = factor(data$daes_clearance_quality_control_results_form_18_complete,levels=c("0","1","2"))
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
data$renal_function_quality_control_results_form_19_complete.factor = factor(data$renal_function_quality_control_results_form_19_complete,levels=c("0","1","2"))
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
data$f13potreactntomed.factor = factor(data$f13potreactntomed,levels=c("1","2","3","4","5","6","7"))
data$f13potreactntomed_mv.factor = factor(data$f13potreactntomed_mv,levels=c("1","0"))
data$f13potreactntomeddate_mv.factor = factor(data$f13potreactntomeddate_mv,levels=c("1","0"))
data$f13otherdiag.factor = factor(data$f13otherdiag,levels=c("1","2"))
data$f13otherdiag_mv.factor = factor(data$f13otherdiag_mv,levels=c("1","0"))
data$f13othrspecfy_mv.factor = factor(data$f13othrspecfy_mv,levels=c("1","0"))
data$f13othrdiagdt_mv.factor = factor(data$f13othrdiagdt_mv,levels=c("1","0"))
data$f13prfmwdclrnc.factor = factor(data$f13prfmwdclrnc,levels=c("1","2"))
data$f13prfmwdclrnc_mv.factor = factor(data$f13prfmwdclrnc_mv,levels=c("1","0"))
data$f13prfmwdclrncdate_mv.factor = factor(data$f13prfmwdclrncdate_mv,levels=c("1","0"))
data$f13stoppntdt_mv.factor = factor(data$f13stoppntdt_mv,levels=c("1","0"))
data$f13dtcompltd_mv.factor = factor(data$f13dtcompltd_mv,levels=c("1","0"))
data$f13forminit_mv.factor = factor(data$f13forminit_mv,levels=c("1","0"))
data$stop_point_form_13_complete.factor = factor(data$stop_point_form_13_complete,levels=c("0","1","2"))
data$f14dob_mv.factor = factor(data$f14dob_mv,levels=c("1","0"))
data$f14dod_mv.factor = factor(data$f14dod_mv,levels=c("1","0"))
data$f14codprimary_mv.factor = factor(data$f14codprimary_mv,levels=c("1","0"))
data$f14codunerlying_mv.factor = factor(data$f14codunerlying_mv,levels=c("1","0"))
data$f14coddrprimary_mv.factor = factor(data$f14coddrprimary_mv,levels=c("1","0"))
data$f14coddrunerlying_mv.factor = factor(data$f14coddrunerlying_mv,levels=c("1","0"))
data$f14physiciancr.factor = factor(data$f14physiciancr,levels=c("1","2"))
data$f14physiciancr_mv.factor = factor(data$f14physiciancr_mv,levels=c("1","0"))
data$f14crphysicianothr_mv.factor = factor(data$f14crphysicianothr_mv,levels=c("1","0"))
data$f14autopsyprfrmd.factor = factor(data$f14autopsyprfrmd,levels=c("1","2","9"))
data$f14autopsyprfrmd_mv.factor = factor(data$f14autopsyprfrmd_mv,levels=c("1","0"))
data$f14autopsyprimarycause_mv.factor = factor(data$f14autopsyprimarycause_mv,levels=c("1","0"))
data$f14autopsyunderlyingcause_mv.factor = factor(data$f14autopsyunderlyingcause_mv,levels=c("1","0"))
data$f14deathloc.factor = factor(data$f14deathloc,levels=c("1","2","3","4","5","9"))
data$f14deathloc_mv.factor = factor(data$f14deathloc_mv,levels=c("1","0"))
data$f4lodother_mv.factor = factor(data$f4lodother_mv,levels=c("1","0"))
data$f14datecompleted_mv.factor = factor(data$f14datecompleted_mv,levels=c("1","0"))
data$f14initials_mv.factor = factor(data$f14initials_mv,levels=c("1","0"))
data$death_notice_form_14_complete.factor = factor(data$death_notice_form_14_complete,levels=c("0","1","2"))
data$phxlabs_visitdate_mv.factor = factor(data$phxlabs_visitdate_mv,levels=c("1","0"))
data$scr_mv.factor = factor(data$scr_mv,levels=c("1","0"))
data$gfr_mv.factor = factor(data$gfr_mv,levels=c("1","0"))
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
data$u0igg_mv.factor = factor(data$u0igg_mv,levels=c("1","0"))
data$u0igg_is_below_limit.factor = factor(data$u0igg_is_below_limit,levels=c("1","0"))
data$u0igg_is_below_limit_mv.factor = factor(data$u0igg_is_below_limit_mv,levels=c("1","0"))
data$u3flow_gfr_mv.factor = factor(data$u3flow_gfr_mv,levels=c("1","0"))
data$u3gfr_mv.factor = factor(data$u3gfr_mv,levels=c("1","0"))
data$u3use_gfr_mv.factor = factor(data$u3use_gfr_mv,levels=c("1","0"))
data$u3albumin_mv.factor = factor(data$u3albumin_mv,levels=c("1","0"))
data$u3albumin_is_below_limit.factor = factor(data$u3albumin_is_below_limit,levels=c("1","0"))
data$u3albumin_is_below_limit_mv.factor = factor(data$u3albumin_is_below_limit_mv,levels=c("1","0"))
data$u3igg_mv.factor = factor(data$u3igg_mv,levels=c("1","0"))
data$u3igg_is_below_limit.factor = factor(data$u3igg_is_below_limit,levels=c("1","0"))
data$u3igg_is_below_limit_mv.factor = factor(data$u3igg_is_below_limit_mv,levels=c("1","0"))
data$ualb_mv.factor = factor(data$ualb_mv,levels=c("1","0"))
data$ualb_is_below_limit.factor = factor(data$ualb_is_below_limit,levels=c("1","0"))
data$ualb_is_below_limit_mv.factor = factor(data$ualb_is_below_limit_mv,levels=c("1","0"))
data$ucr_mv.factor = factor(data$ucr_mv,levels=c("1","0"))
data$phoenix_labs_complete.factor = factor(data$phoenix_labs_complete,levels=c("0","1","2"))
data$stflabs_visitdate_mv.factor = factor(data$stflabs_visitdate_mv,levels=c("1","0"))
data$s_gfr_mv.factor = factor(data$s_gfr_mv,levels=c("1","0"))
data$s_ff_mv.factor = factor(data$s_ff_mv,levels=c("1","0"))
data$s_pahclrnc_mv.factor = factor(data$s_pahclrnc_mv,levels=c("1","0"))
data$s_gfrb_mv.factor = factor(data$s_gfrb_mv,levels=c("1","0"))
data$s_pahb_mv.factor = factor(data$s_pahb_mv,levels=c("1","0"))
data$s_rpf_mv.factor = factor(data$s_rpf_mv,levels=c("1","0"))
data$s_rpfb_mv.factor = factor(data$s_rpfb_mv,levels=c("1","0"))
data$stanford_labs_complete.factor = factor(data$stanford_labs_complete,levels=c("0","1","2"))
data$biopsy_dt_mv.factor = factor(data$biopsy_dt_mv,levels=c("1","0"))
data$biopsy_source.factor = factor(data$biopsy_source,levels=c("P","C"))
data$biopsy_source_mv.factor = factor(data$biopsy_source_mv,levels=c("1","0"))
data$glomerul_mv.factor = factor(data$glomerul_mv,levels=c("1","0"))
data$pct_glob_mv.factor = factor(data$pct_glob_mv,levels=c("1","0"))
data$mean_fia_mv.factor = factor(data$mean_fia_mv,levels=c("1","0"))
data$mean_nv_mv.factor = factor(data$mean_nv_mv,levels=c("1","0"))
data$mean_nv0_mv.factor = factor(data$mean_nv0_mv,levels=c("1","0"))
data$mean_nv1_mv.factor = factor(data$mean_nv1_mv,levels=c("1","0"))
data$mean_n_p_mv.factor = factor(data$mean_n_p_mv,levels=c("1","0"))
data$mean_n_e_mv.factor = factor(data$mean_n_e_mv,levels=c("1","0"))
data$mean_aa_mv.factor = factor(data$mean_aa_mv,levels=c("1","0"))
data$mean_sv_mv.factor = factor(data$mean_sv_mv,levels=c("1","0"))
data$mean_bmt_mv.factor = factor(data$mean_bmt_mv,levels=c("1","0"))
data$mean_fsf_mv.factor = factor(data$mean_fsf_mv,levels=c("1","0"))
data$mean_fpw_mv.factor = factor(data$mean_fpw_mv,levels=c("1","0"))
data$pct_denu_mv.factor = factor(data$pct_denu_mv,levels=c("1","0"))
data$pct_fene_mv.factor = factor(data$pct_fene_mv,levels=c("1","0"))
data$pct_fen_mv.factor = factor(data$pct_fen_mv,levels=c("1","0"))
data$sa_using_mv.factor = factor(data$sa_using_mv,levels=c("1","0"))
data$number_o_mv.factor = factor(data$number_o_mv,levels=c("1","0"))
data$alb_stat_mv.factor = factor(data$alb_stat_mv,levels=c("1","0"))
data$comments_mv.factor = factor(data$comments_mv,levels=c("1","0"))
data$biopsy_complete.factor = factor(data$biopsy_complete,levels=c("0","1","2"))

levels(data$redcap_event_name.factor)=c("SS Interval -3","SS Interval -2","SS Interval -1","RC Interval 0","RV Interval 0.7","RC Interval 1","RV Interval 3","RV Interval 6","RV Interval 9","RC Interval 12","RV Interval 15","RV Interval 18","RV Interval 21","RC Interval 24","RV Interval 27","RV Interval 30","RV Interval 33","RC Interval 36","RV Interval 39","RV Interval 42","RV Interval 45","RC Interval 48","RV Interval 51","RV Interval 54","RV Interval 57","RC Interval 60","RV Interval 63","RV Interval 66","RV Interval 69","RC Interval 72","RV Interval 75","RV Interval 78","RV Interval 81","RC Interval 84","RV Interval 87","RV Interval 90","RV Interval 93","RC Interval 96","RV Interval 99","RV Interval 102","RV Interval 105","RC Interval 108","RV Interval 111","RV Interval 114","RV Interval 117","RC Interval 120","RV Interval 123","RV Interval 126","RV Interval 129","RC Interval 132","RV Interval 135","RV Interval 138","RV Interval 141","RC Interval 144","RV Interval 147","RV Interval 150","RV Interval 153","RC Interval 156","RV Interval 159","RV Interval 162","RV Interval 165","RC Interval 168","RV Interval 171","RV Interval 174","RV Interval 177","RC Interval 180","RV Interval 183","RV Interval 186","RV Interval 189","RC Interval 192","RV Interval 195","RV Interval 198","RV Interval 201","RC Interval 204")
levels(data$redcap_repeat_instrument.factor)=c("")
levels(data$lastname_mv.factor)=c("Yes","No")
levels(data$firstname_mv.factor)=c("Yes","No")
levels(data$stratum.factor)=c("Normal UAE","Micro UAE")
levels(data$stratum_mv.factor)=c("Yes","No")
levels(data$randomnum_mv.factor)=c("Yes","No")
levels(data$tx.factor)=c("Placebo","Losartan")
levels(data$tx_mv.factor)=c("Yes","No")
levels(data$sacaton_no_mv.factor)=c("Yes","No")
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
levels(data$closedate_mv.factor)=c("Yes","No")
levels(data$washoutdate_mv.factor)=c("Yes","No")
levels(data$demoform_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f2visitdate_mv.factor)=c("Yes","No")
levels(data$f2visittype.factor)=c("Screening")
levels(data$f2visittype_mv.factor)=c("Yes","No")
levels(data$f2testintvrl_mv.factor)=c("Yes","No")
levels(data$f2bldprsrarm.factor)=c("Unknown","Right","Left")
levels(data$f2bldprsrarm_mv.factor)=c("Yes","No")
levels(data$f2frstmsrsys_mv.factor)=c("Yes","No")
levels(data$f2frstmsrdia_mv.factor)=c("Yes","No")
levels(data$f2scndmsrsys_mv.factor)=c("Yes","No")
levels(data$f2scndmsrdia_mv.factor)=c("Yes","No")
levels(data$f2heartrate_mv.factor)=c("Yes","No")
levels(data$subject_screening_3_and_2_form_02_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f3visitdate_mv.factor)=c("Yes","No")
levels(data$f3visittype.factor)=c("Clearance","Screening","Follow-up")
levels(data$f3visittype_mv.factor)=c("Yes","No")
levels(data$f3testintvrl_mv.factor)=c("Yes","No")
levels(data$f3spturncltd.factor)=c("No","Yes")
levels(data$f3spturncltd_mv.factor)=c("Yes","No")
levels(data$f3venipfork.factor)=c("No","Yes")
levels(data$f3venipfork_mv.factor)=c("Yes","No")
levels(data$f3bldprsrarm.factor)=c("Right","Left")
levels(data$f3bldprsrarm_mv.factor)=c("Yes","No")
levels(data$f3frstmsrsys_mv.factor)=c("Yes","No")
levels(data$f3frstmsrdia_mv.factor)=c("Yes","No")
levels(data$f3scndmsrsys_mv.factor)=c("Yes","No")
levels(data$f3scndmsrdia_mv.factor)=c("Yes","No")
levels(data$f3heartrate_mv.factor)=c("Yes","No")
levels(data$f3prohbmeds.factor)=c("No","Yes")
levels(data$f3prohbmeds_mv.factor)=c("Yes","No")
levels(data$f3elgstatus.factor)=c("not elgible","elgible")
levels(data$f3elgstatus_mv.factor)=c("Yes","No")
levels(data$f3inelgibleexplain_mv.factor)=c("Yes","No")
levels(data$f3medsurghis.factor)=c("No","Yes")
levels(data$f3medsurghis_mv.factor)=c("Yes","No")
levels(data$f3_acratio.factor)=c("No","Yes")
levels(data$f3_acratio_mv.factor)=c("Yes","No")
levels(data$f3srmcreatyn.factor)=c("No","Yes")
levels(data$f3srmcreatyn_mv.factor)=c("Yes","No")
levels(data$f3_pregnancy.factor)=c("No","Yes")
levels(data$f3_pregnancy_mv.factor)=c("Yes","No")
levels(data$f3_staffpref.factor)=c("No","Yes")
levels(data$f3_staffpref_mv.factor)=c("Yes","No")
levels(data$f3_subjpref.factor)=c("No","Yes")
levels(data$f3_subjpref_mv.factor)=c("Yes","No")
levels(data$f3_other.factor)=c("No","Yes")
levels(data$f3_other_mv.factor)=c("Yes","No")
levels(data$subject_screeningelgibility_form_03_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f3meds_visitdate_mv.factor)=c("Yes","No")
levels(data$f3drugcode1_mv.factor)=c("Yes","No")
levels(data$f3dosage1_mv.factor)=c("Yes","No")
levels(data$f3timesday1_mv.factor)=c("Yes","No")
levels(data$f3dwm1.factor)=c("Day","Week","Month")
levels(data$f3dwm1_mv.factor)=c("Yes","No")
levels(data$f3drugcode2_mv.factor)=c("Yes","No")
levels(data$f3dosage2_mv.factor)=c("Yes","No")
levels(data$f3timesday2_mv.factor)=c("Yes","No")
levels(data$f3dwm2.factor)=c("Day","Week","Month")
levels(data$f3dwm2_mv.factor)=c("Yes","No")
levels(data$f3drugcode3_mv.factor)=c("Yes","No")
levels(data$f3dosage3_mv.factor)=c("Yes","No")
levels(data$f3timesday3_mv.factor)=c("Yes","No")
levels(data$f3dwm3.factor)=c("Day","Week","Month")
levels(data$f3dwm3_mv.factor)=c("Yes","No")
levels(data$f3drugcode4_mv.factor)=c("Yes","No")
levels(data$f3dosage4_mv.factor)=c("Yes","No")
levels(data$f3timesday4_mv.factor)=c("Yes","No")
levels(data$f3dwm4.factor)=c("Day","Week","Month")
levels(data$f3dwm4_mv.factor)=c("Yes","No")
levels(data$f3drugcode5_mv.factor)=c("Yes","No")
levels(data$f3dosage5_mv.factor)=c("Yes","No")
levels(data$f3timesday5_mv.factor)=c("Yes","No")
levels(data$f3dwm5.factor)=c("Day","Week","Month")
levels(data$f3dwm5_mv.factor)=c("Yes","No")
levels(data$f3drugcode6_mv.factor)=c("Yes","No")
levels(data$f3dosage6_mv.factor)=c("Yes","No")
levels(data$f3timesday6_mv.factor)=c("Yes","No")
levels(data$f3dwm6.factor)=c("Day","Week","Month")
levels(data$f3dwm6_mv.factor)=c("Yes","No")
levels(data$f3drugcode7_mv.factor)=c("Yes","No")
levels(data$f3dosage7_mv.factor)=c("Yes","No")
levels(data$f3timesday7_mv.factor)=c("Yes","No")
levels(data$f3dwm7.factor)=c("Day","Week","Month")
levels(data$f3dwm7_mv.factor)=c("Yes","No")
levels(data$f3drugcode8_mv.factor)=c("Yes","No")
levels(data$f3dosage8_mv.factor)=c("Yes","No")
levels(data$f3timesday8_mv.factor)=c("Yes","No")
levels(data$f3dwm8.factor)=c("Day","Week","Month")
levels(data$f3dwm8_mv.factor)=c("Yes","No")
levels(data$f3drugcode9_mv.factor)=c("Yes","No")
levels(data$f3dosage9_mv.factor)=c("Yes","No")
levels(data$f3timesday9_mv.factor)=c("Yes","No")
levels(data$f3dwm9.factor)=c("Day","Week","Month")
levels(data$f3dwm9_mv.factor)=c("Yes","No")
levels(data$f3drugcode10_mv.factor)=c("Yes","No")
levels(data$f3dosage10_mv.factor)=c("Yes","No")
levels(data$f3timesday10_mv.factor)=c("Yes","No")
levels(data$f3dwm10.factor)=c("Day","Week","Month")
levels(data$f3dwm10_mv.factor)=c("Yes","No")
levels(data$f3drugcode11_mv.factor)=c("Yes","No")
levels(data$f3dosage11_mv.factor)=c("Yes","No")
levels(data$f3timesday11_mv.factor)=c("Yes","No")
levels(data$f3dwm11.factor)=c("Day","Week","Month")
levels(data$f3dwm11_mv.factor)=c("Yes","No")
levels(data$f3drugcode12_mv.factor)=c("Yes","No")
levels(data$f3dosage12_mv.factor)=c("Yes","No")
levels(data$f3timesday12_mv.factor)=c("Yes","No")
levels(data$f3dwm12.factor)=c("Day","Week","Month")
levels(data$f3dwm12_mv.factor)=c("Yes","No")
levels(data$f3drugcode13_mv.factor)=c("Yes","No")
levels(data$f3dosage13_mv.factor)=c("Yes","No")
levels(data$f3timesday13_mv.factor)=c("Yes","No")
levels(data$f3dwm13.factor)=c("Day","Week","Month")
levels(data$f3dwm13_mv.factor)=c("Yes","No")
levels(data$f3drugcode14_mv.factor)=c("Yes","No")
levels(data$f3dosage14_mv.factor)=c("Yes","No")
levels(data$f3timesday14_mv.factor)=c("Yes","No")
levels(data$f3dwm14.factor)=c("Day","Week","Month")
levels(data$f3dwm14_mv.factor)=c("Yes","No")
levels(data$f3drugcode15_mv.factor)=c("Yes","No")
levels(data$f3dosage15_mv.factor)=c("Yes","No")
levels(data$f3timesday15_mv.factor)=c("Yes","No")
levels(data$f3dwm15.factor)=c("Day","Week","Month")
levels(data$f3dwm15_mv.factor)=c("Yes","No")
levels(data$f3drugcode16_mv.factor)=c("Yes","No")
levels(data$f3dosage16_mv.factor)=c("Yes","No")
levels(data$f3timesday16_mv.factor)=c("Yes","No")
levels(data$f3dwm16.factor)=c("Day","Week","Month")
levels(data$f3dwm16_mv.factor)=c("Yes","No")
levels(data$f3drugcode17_mv.factor)=c("Yes","No")
levels(data$f3dosage17_mv.factor)=c("Yes","No")
levels(data$f3timesday17_mv.factor)=c("Yes","No")
levels(data$f3dwm17.factor)=c("Day","Week","Month")
levels(data$f3dwm17_mv.factor)=c("Yes","No")
levels(data$f3drugcode18_mv.factor)=c("Yes","No")
levels(data$f3dosage18_mv.factor)=c("Yes","No")
levels(data$f3timesday18_mv.factor)=c("Yes","No")
levels(data$f3dwm18.factor)=c("Day","Week","Month")
levels(data$f3dwm18_mv.factor)=c("Yes","No")
levels(data$f3drugcode19_mv.factor)=c("Yes","No")
levels(data$f3dosage19_mv.factor)=c("Yes","No")
levels(data$f3timesday19_mv.factor)=c("Yes","No")
levels(data$f3dwm19.factor)=c("Day","Week","Month")
levels(data$f3dwm19_mv.factor)=c("Yes","No")
levels(data$f3drugcode20_mv.factor)=c("Yes","No")
levels(data$f3dosage20_mv.factor)=c("Yes","No")
levels(data$f3timesday20_mv.factor)=c("Yes","No")
levels(data$f3dwm20.factor)=c("Day","Week","Month")
levels(data$f3dwm20_mv.factor)=c("Yes","No")
levels(data$subject_screeningelgibility_meds_form_03_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$clrncvisitdate_mv.factor)=c("Yes","No")
levels(data$clrncvisittype.factor)=c("Clearance","Follow-up")
levels(data$clrncvisittype_mv.factor)=c("Yes","No")
levels(data$clrnctestintvrl_mv.factor)=c("Yes","No")
levels(data$renal_clearance_visit_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$intfu_visitdate_mv.factor)=c("Yes","No")
levels(data$intfu_visittype.factor)=c("Clearance","Follow-up")
levels(data$intfu_visittype_mv.factor)=c("Yes","No")
levels(data$intfu_testintrvl_mv.factor)=c("Yes","No")
levels(data$interim_followup_visit_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f4visitdate_mv.factor)=c("Yes","No")
levels(data$f4visittype_mv.factor)=c("Yes","No")
levels(data$f4testintvl_mv.factor)=c("Yes","No")
levels(data$f4bldprsr.factor)=c("Unknown","Right","Left")
levels(data$f4bldprsr_mv.factor)=c("Yes","No")
levels(data$f4frstmsrsys_mv.factor)=c("Yes","No")
levels(data$f4frstmsrdia_mv.factor)=c("Yes","No")
levels(data$f4scndmsrsys_mv.factor)=c("Yes","No")
levels(data$f4scndmsrdia_mv.factor)=c("Yes","No")
levels(data$f4heartrate_mv.factor)=c("Yes","No")
levels(data$f4pregnant.factor)=c("Unknown","No","Yes")
levels(data$f4pregnant_mv.factor)=c("Yes","No")
levels(data$f4medsnow.factor)=c("No","Yes")
levels(data$f4medsnow_mv.factor)=c("Yes","No")
levels(data$f4srmpotasm_mv.factor)=c("Yes","No")
levels(data$f4srmglu_mv.factor)=c("Yes","No")
levels(data$f4forminit_mv.factor)=c("Yes","No")
levels(data$blood_pressure_form_04_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f4_visitdate_mv.factor)=c("Yes","No")
levels(data$f4_pillbotldt_mv.factor)=c("Yes","No")
levels(data$f4_dsplstpill_mv.factor)=c("Yes","No")
levels(data$f4_dysbtwndsp_mv.factor)=c("Yes","No")
levels(data$f4_nodlydose_mv.factor)=c("Yes","No")
levels(data$f4_dlydosrtrn_mv.factor)=c("Yes","No")
levels(data$f4_dlydosetkn_mv.factor)=c("Yes","No")
levels(data$f4_compliance_mv.factor)=c("Yes","No")
levels(data$f4_totdlydosn_mv.factor)=c("Yes","No")
levels(data$f4_curdos_mv.factor)=c("Yes","No")
levels(data$f4_dateadjstd_mv.factor)=c("Yes","No")
levels(data$f4_ajstddose_mv.factor)=c("Yes","No")
levels(data$f4_medseffects.factor)=c("None","nausea","diarrhea","headache","vertigo","myalgias","cough","skin rash","symptoms of postural hypotension","anorexia","fatigue","other")
levels(data$f4_medseffects_mv.factor)=c("Yes","No")
levels(data$f4_expeffects_mv.factor)=c("Yes","No")
levels(data$drug_compliance_form_04_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f4meds_visitdate_mv.factor)=c("Yes","No")
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
levels(data$meds_form_04_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f5visitdate_mv.factor)=c("Yes","No")
levels(data$f5rcntillns.factor)=c("No","Yes")
levels(data$f5rcntillns_mv.factor)=c("Yes","No")
levels(data$f5outpat.factor)=c("Unknown","No","Yes")
levels(data$f5outpat_mv.factor)=c("Yes","No")
levels(data$f5outpatspec_mv.factor)=c("Yes","No")
levels(data$f5hospwosrg.factor)=c("Unknown","No","Yes")
levels(data$f5hospwosrg_mv.factor)=c("Yes","No")
levels(data$f5hospwoexp_mv.factor)=c("Yes","No")
levels(data$f5hospwsrg.factor)=c("Unknown","No","Yes")
levels(data$f5hospwsrg_mv.factor)=c("Yes","No")
levels(data$f5hospwexp_mv.factor)=c("Yes","No")
levels(data$f5medprbsdoc.factor)=c("Unknown","No","Yes")
levels(data$f5medprbsdoc_mv.factor)=c("Yes","No")
levels(data$f5crnryartds.factor)=c("Unknown","No","Yes")
levels(data$f5crnryartds_mv.factor)=c("Yes","No")
levels(data$f5cancer.factor)=c("Unknown","No","Yes")
levels(data$f5cancer_mv.factor)=c("Yes","No")
levels(data$f5cancerexp_mv.factor)=c("Yes","No")
levels(data$f5crbrlvasds.factor)=c("Unknown","No","Yes")
levels(data$f5crbrlvasds_mv.factor)=c("Yes","No")
levels(data$f5peripvasds.factor)=c("Unknown","No","Yes")
levels(data$f5peripvasds_mv.factor)=c("Yes","No")
levels(data$f5hypertensn.factor)=c("Unknown","No","Yes")
levels(data$f5hypertensn_mv.factor)=c("Yes","No")
levels(data$f5seizures.factor)=c("Unknown","No","Yes")
levels(data$f5seizures_mv.factor)=c("Yes","No")
levels(data$f5gnitrnryds.factor)=c("Unknown","No","Yes")
levels(data$f5gnitrnryds_mv.factor)=c("Yes","No")
levels(data$f5gnitrnyexp_mv.factor)=c("Yes","No")
levels(data$f5lungdsz.factor)=c("Unknown","No","Yes")
levels(data$f5lungdsz_mv.factor)=c("Yes","No")
levels(data$f5majsurg.factor)=c("Unknown","No","Yes")
levels(data$f5majsurg_mv.factor)=c("Yes","No")
levels(data$f5majsurgexp_mv.factor)=c("Yes","No")
levels(data$f5othrmeddia.factor)=c("Unknown","No","Yes")
levels(data$f5othrmeddia_mv.factor)=c("Yes","No")
levels(data$f5othrmedexp_mv.factor)=c("Yes","No")
levels(data$f5lungs.factor)=c("Unknown","Normal","Abnormal")
levels(data$f5lungs_mv.factor)=c("Yes","No")
levels(data$f5lungsexp_mv.factor)=c("Yes","No")
levels(data$f5heart.factor)=c("Unknown","Normal","Abnormal")
levels(data$f5heart_mv.factor)=c("Yes","No")
levels(data$f5heartexp_mv.factor)=c("Yes","No")
levels(data$f5skin.factor)=c("Unknown","Normal","Abnormal")
levels(data$f5skin_mv.factor)=c("Yes","No")
levels(data$f5skinexp_mv.factor)=c("Yes","No")
levels(data$f5edema.factor)=c("Unknown","Normal","Abnormal")
levels(data$f5edema_mv.factor)=c("Yes","No")
levels(data$f5edemaexp_mv.factor)=c("Yes","No")
levels(data$history_and_physical_examination_form_05_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f6visitdate_mv.factor)=c("Yes","No")
levels(data$f6ivsites_mv.factor)=c("Yes","No")
levels(data$f6weight_mv.factor)=c("Yes","No")
levels(data$f6bsa_mv.factor)=c("Yes","No")
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
levels(data$f6glucose.factor)=c("Neg","Trace (100 mg/dl)","Trace (100 mg/dl)","250 mg/dl","500 mg/dl","1000 mg/dl",">2000 mg/dl")
levels(data$f6glucose_mv.factor)=c("Yes","No")
levels(data$f6bilirubin.factor)=c("Neg","1+ (small)","2+ (moderate)","3+ (large)")
levels(data$f6bilirubin_mv.factor)=c("Yes","No")
levels(data$f6ketones.factor)=c("Neg","Trace (5)","Trace (5)","Small (15)","Moderate (40)","80",">160")
levels(data$f6ketones_mv.factor)=c("Yes","No")
levels(data$f6specgrvty_mv.factor)=c("Yes","No")
levels(data$f6bloodocult.factor)=c("No","Yes")
levels(data$f6bloodocult_mv.factor)=c("Yes","No")
levels(data$f6ph_mv.factor)=c("Yes","No")
levels(data$f6protein.factor)=c("Neg","Trace","Trace","1+ 30 mg/dl","2+ 100 mg/dl","3+ 300 mg/dl","4+ > 2000 mg/dl")
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
levels(data$renal_clearance_data_form_06_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f7qcnumber_mv.factor)=c("Yes","No")
levels(data$f7visitdate_mv.factor)=c("Yes","No")
levels(data$f7tstintrvl_mv.factor)=c("Yes","No")
levels(data$f7typesample.factor)=c("Stanford Lab GFR sample (iothalamate and PAH)","NIH Lab P0 creatinine, albumin, IgG - DAES clearance results","SMAC","CBC","U0 albumin, IgG, and creatinine (spot urine) - DAES clearance result","NIH Lab glycosylated hemoglobin (HbA1c)")
levels(data$f7typesample_mv.factor)=c("Yes","No")
levels(data$laboratory_quality_control_form_07_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f8visitdate_mv.factor)=c("Yes","No")
levels(data$f8analdate_mv.factor)=c("Yes","No")
levels(data$f8alti_mv.factor)=c("Yes","No")
levels(data$f8glucose_mv.factor)=c("Yes","No")
levels(data$f8bun_mv.factor)=c("Yes","No")
levels(data$f8creatinine_mv.factor)=c("Yes","No")
levels(data$f8sodium_mv.factor)=c("Yes","No")
levels(data$f8potassium_mv.factor)=c("Yes","No")
levels(data$f8chloride_mv.factor)=c("Yes","No")
levels(data$f8calcium_mv.factor)=c("Yes","No")
levels(data$f8phosphorus_mv.factor)=c("Yes","No")
levels(data$f8totlprotn_mv.factor)=c("Yes","No")
levels(data$f8albumin_mv.factor)=c("Yes","No")
levels(data$f8cholestero_mv.factor)=c("Yes","No")
levels(data$f8triglyceri_mv.factor)=c("Yes","No")
levels(data$f8hdl_mv.factor)=c("Yes","No")
levels(data$f8ldl_mv.factor)=c("Yes","No")
levels(data$f8totlbilrbn_mv.factor)=c("Yes","No")
levels(data$f8dbilirubin_mv.factor)=c("Yes","No")
levels(data$f8alkp04_mv.factor)=c("Yes","No")
levels(data$f8ast_mv.factor)=c("Yes","No")
levels(data$smacchemistry_panel_form_08_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f9visitdate_mv.factor)=c("Yes","No")
levels(data$f9dateanal_mv.factor)=c("Yes","No")
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
levels(data$cbc_form_09_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f10visitdate_mv.factor)=c("Yes","No")
levels(data$f10spturalbu_mv.factor)=c("Yes","No")
levels(data$f10spturncrea_mv.factor)=c("Yes","No")
levels(data$f10serumcreat_mv.factor)=c("Yes","No")
levels(data$f10strgbxnum_mv.factor)=c("Yes","No")
levels(data$daes_routine_labs_form_10_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f11visitdate_mv.factor)=c("Yes","No")
levels(data$f11hba1c_mv.factor)=c("Yes","No")
levels(data$f11p0glu_mv.factor)=c("Yes","No")
levels(data$f11p0creatn_mv.factor)=c("Yes","No")
levels(data$f11p0albumn_mv.factor)=c("Yes","No")
levels(data$f11p3albumn_mv.factor)=c("Yes","No")
levels(data$f11p0iggsamp_mv.factor)=c("Yes","No")
levels(data$f11p3iggsamp_mv.factor)=c("Yes","No")
levels(data$f11u0creatn_mv.factor)=c("Yes","No")
levels(data$f11u0albumn_mv.factor)=c("Yes","No")
levels(data$f11u3albumn_mv.factor)=c("Yes","No")
levels(data$f11u0igg_mv.factor)=c("Yes","No")
levels(data$f11u3igg_mv.factor)=c("Yes","No")
levels(data$daes_lab_clearance_results_form_11_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f12expdatevis_mv.factor)=c("Yes","No")
levels(data$f12typevisit.factor)=c("Clearance","Follow-up")
levels(data$f12typevisit_mv.factor)=c("Yes","No")
levels(data$f12testintrvl_mv.factor)=c("Yes","No")
levels(data$f12reasnmissd.factor)=c("illness","hospitalization","personal/family business","employment","weather","unknown","pregnancy","convalescence","Other")
levels(data$f12reasnmissd_mv.factor)=c("Yes","No")
levels(data$f12otherspec_mv.factor)=c("Yes","No")
levels(data$missed_visitintercurrent_event_form_12_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f16visitdate_mv.factor)=c("Yes","No")
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
levels(data$smac_quality_control_form_16_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f17visitdate_mv.factor)=c("Yes","No")
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
levels(data$cbc_quality_control_form_17_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f18visitdate_mv.factor)=c("Yes","No")
levels(data$f18qcnumber_mv.factor)=c("Yes","No")
levels(data$f18hba1c_mv.factor)=c("Yes","No")
levels(data$f18p0creatn_mv.factor)=c("Yes","No")
levels(data$f18p0albumin_mv.factor)=c("Yes","No")
levels(data$f18p0igg_mv.factor)=c("Yes","No")
levels(data$f18u0creatn_mv.factor)=c("Yes","No")
levels(data$f18u0albumin_mv.factor)=c("Yes","No")
levels(data$f18u0igg_mv.factor)=c("Yes","No")
levels(data$daes_clearance_quality_control_results_form_18_complete.factor)=c("Incomplete","Unverified","Complete")
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
levels(data$renal_function_quality_control_results_form_19_complete.factor)=c("Incomplete","Unverified","Complete")
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
levels(data$f13potreactntomed.factor)=c("Hyperkalemia","Angioedema","Neutropenia","Thrombocytopenia","Symptomatic hypotension","Intractable cough","Other")
levels(data$f13potreactntomed_mv.factor)=c("Yes","No")
levels(data$f13potreactntomeddate_mv.factor)=c("Yes","No")
levels(data$f13otherdiag.factor)=c("No","Yes")
levels(data$f13otherdiag_mv.factor)=c("Yes","No")
levels(data$f13othrspecfy_mv.factor)=c("Yes","No")
levels(data$f13othrdiagdt_mv.factor)=c("Yes","No")
levels(data$f13prfmwdclrnc.factor)=c("No","Yes")
levels(data$f13prfmwdclrnc_mv.factor)=c("Yes","No")
levels(data$f13prfmwdclrncdate_mv.factor)=c("Yes","No")
levels(data$f13stoppntdt_mv.factor)=c("Yes","No")
levels(data$f13dtcompltd_mv.factor)=c("Yes","No")
levels(data$f13forminit_mv.factor)=c("Yes","No")
levels(data$stop_point_form_13_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f14dob_mv.factor)=c("Yes","No")
levels(data$f14dod_mv.factor)=c("Yes","No")
levels(data$f14codprimary_mv.factor)=c("Yes","No")
levels(data$f14codunerlying_mv.factor)=c("Yes","No")
levels(data$f14coddrprimary_mv.factor)=c("Yes","No")
levels(data$f14coddrunerlying_mv.factor)=c("Yes","No")
levels(data$f14physiciancr.factor)=c("Nelson","Sievers")
levels(data$f14physiciancr_mv.factor)=c("Yes","No")
levels(data$f14crphysicianothr_mv.factor)=c("Yes","No")
levels(data$f14autopsyprfrmd.factor)=c("No","Yes","Unknown")
levels(data$f14autopsyprfrmd_mv.factor)=c("Yes","No")
levels(data$f14autopsyprimarycause_mv.factor)=c("Yes","No")
levels(data$f14autopsyunderlyingcause_mv.factor)=c("Yes","No")
levels(data$f14deathloc.factor)=c("in hospital","at work","at home","en route to hospital","Other (specify)","Unknown")
levels(data$f14deathloc_mv.factor)=c("Yes","No")
levels(data$f4lodother_mv.factor)=c("Yes","No")
levels(data$f14datecompleted_mv.factor)=c("Yes","No")
levels(data$f14initials_mv.factor)=c("Yes","No")
levels(data$death_notice_form_14_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$phxlabs_visitdate_mv.factor)=c("Yes","No")
levels(data$scr_mv.factor)=c("Yes","No")
levels(data$gfr_mv.factor)=c("Yes","No")
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
levels(data$u0igg_mv.factor)=c("Yes","No")
levels(data$u0igg_is_below_limit.factor)=c("Yes","No")
levels(data$u0igg_is_below_limit_mv.factor)=c("Yes","No")
levels(data$u3flow_gfr_mv.factor)=c("Yes","No")
levels(data$u3gfr_mv.factor)=c("Yes","No")
levels(data$u3use_gfr_mv.factor)=c("Yes","No")
levels(data$u3albumin_mv.factor)=c("Yes","No")
levels(data$u3albumin_is_below_limit.factor)=c("Yes","No")
levels(data$u3albumin_is_below_limit_mv.factor)=c("Yes","No")
levels(data$u3igg_mv.factor)=c("Yes","No")
levels(data$u3igg_is_below_limit.factor)=c("Yes","No")
levels(data$u3igg_is_below_limit_mv.factor)=c("Yes","No")
levels(data$ualb_mv.factor)=c("Yes","No")
levels(data$ualb_is_below_limit.factor)=c("Yes","No")
levels(data$ualb_is_below_limit_mv.factor)=c("Yes","No")
levels(data$ucr_mv.factor)=c("Yes","No")
levels(data$phoenix_labs_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$stflabs_visitdate_mv.factor)=c("Yes","No")
levels(data$s_gfr_mv.factor)=c("Yes","No")
levels(data$s_ff_mv.factor)=c("Yes","No")
levels(data$s_pahclrnc_mv.factor)=c("Yes","No")
levels(data$s_gfrb_mv.factor)=c("Yes","No")
levels(data$s_pahb_mv.factor)=c("Yes","No")
levels(data$s_rpf_mv.factor)=c("Yes","No")
levels(data$s_rpfb_mv.factor)=c("Yes","No")
levels(data$stanford_labs_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$biopsy_dt_mv.factor)=c("Yes","No")
levels(data$biopsy_source.factor)=c("University of Minnesota","Stanford")
levels(data$biopsy_source_mv.factor)=c("Yes","No")
levels(data$glomerul_mv.factor)=c("Yes","No")
levels(data$pct_glob_mv.factor)=c("Yes","No")
levels(data$mean_fia_mv.factor)=c("Yes","No")
levels(data$mean_nv_mv.factor)=c("Yes","No")
levels(data$mean_nv0_mv.factor)=c("Yes","No")
levels(data$mean_nv1_mv.factor)=c("Yes","No")
levels(data$mean_n_p_mv.factor)=c("Yes","No")
levels(data$mean_n_e_mv.factor)=c("Yes","No")
levels(data$mean_aa_mv.factor)=c("Yes","No")
levels(data$mean_sv_mv.factor)=c("Yes","No")
levels(data$mean_bmt_mv.factor)=c("Yes","No")
levels(data$mean_fsf_mv.factor)=c("Yes","No")
levels(data$mean_fpw_mv.factor)=c("Yes","No")
levels(data$pct_denu_mv.factor)=c("Yes","No")
levels(data$pct_fene_mv.factor)=c("Yes","No")
levels(data$pct_fen_mv.factor)=c("Yes","No")
levels(data$sa_using_mv.factor)=c("Yes","No")
levels(data$number_o_mv.factor)=c("Yes","No")
levels(data$alb_stat_mv.factor)=c("Yes","No")
levels(data$comments_mv.factor)=c("Yes","No")
levels(data$biopsy_complete.factor)=c("Incomplete","Unverified","Complete")
