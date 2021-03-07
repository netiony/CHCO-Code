# Read in
data=read.csv('./FicollUniversityOfMi_DATA_2021-02-04_1618.csv',na.strings = "")
#Setting Labels
# label(data$record_id)="Patient ID"
# label(data$redcap_event_name)="Event Name"
# label(data$redcap_repeat_instrument)="Repeat Instrument"
# label(data$redcap_repeat_instance)="Repeat Instance"
# label(data$lastname)="Lastname"
# label(data$lastname_mv)="Missing value verified for: Lastname"
# label(data$firstname)="Firstname"
# label(data$firstname_mv)="Missing value verified for: Firstname"
# label(data$sacaton_no)="Sacaton number"
# label(data$sacaton_no_mv)="Missing value verified for: Sacaton number"
# label(data$phx_number)="PIMC number"
# label(data$phx_number_mv)="Missing value verified for: PIMC number"
# label(data$dob)="DOB"
# label(data$dob_mv)="Missing value verified for: DOB"
# label(data$gender)="Gender"
# label(data$gender_mv)="Missing value verified for: Gender"
# label(data$height)="Height (cm)"
# label(data$height_mv)="Missing value verified for: Height (cm)"
# label(data$allergies)="Allergies?"
# label(data$allergies_mv)="Missing value verified for: Allergies?"
# label(data$allerglist)="Allergies (specify)"
# label(data$allerglist_mv)="Missing value verified for: Allergies (specify)"
# label(data$district)="District"
# label(data$district_mv)="Missing value verified for: District"
# label(data$phonenum)="Phone number"
# label(data$phonenum_mv)="Missing value verified for: Phone number"
# label(data$dtdiabonset)="Date of Diabetes Onset?"
# label(data$dtdiabonset_mv)="Missing value verified for: Date of Diabetes Onset?"
# label(data$demoform_complete)="Complete?"
# label(data$clrncvisitdate)="Date of Visit"
# label(data$clrncvisitdate_mv)="Missing value verified for: Date of Visit"
# label(data$clrncvisittype)="Type of Visit"
# label(data$clrncvisittype_mv)="Missing value verified for: Type of Visit"
# label(data$clrnctestintvrl)="Test Interval"
# label(data$clrnctestintvrl_mv)="Missing value verified for: Test Interval"
# label(data$clearancevisit_complete)="Complete?"
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
# label(data$f4medseffcts)="Med effects?"
# label(data$f4medseffcts_mv)="Missing value verified for: Med effects?"
# label(data$f4pregnant)="Pregnancy Test?"
# label(data$f4pregnant_mv)="Missing value verified for: Pregnancy Test?"
# label(data$f4medsnow)="Is the subect currently taking any other medicine?"
# label(data$f4medsnow_mv)="Missing value verified for: Is the subect currently taking any other medicine?"
# label(data$f4srmpotasm)="Serum potassium"
# label(data$f4srmpotasm_mv)="Missing value verified for: Serum potassium"
# label(data$f4forminit)="Initials of person who completed form"
# label(data$f4forminit_mv)="Missing value verified for: Initials of person who completed form"
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
# label(data$f4meds_complete)="Complete?"
# label(data$f5rcntillns)="Recent Illnesses Rqring Med Att."
# label(data$f5rcntillns_mv)="Missing value verified for: Recent Illnesses Rqring Med Att."
# label(data$f5outpat)="Out Patient?"
# label(data$f5outpat_mv)="Missing value verified for: Out Patient?"
# label(data$f5outpatspec)="Out Patient - Specify"
# label(data$f5outpatspec_mv)="Missing value verified for: Out Patient - Specify"
# label(data$f5hospwosrg)="Hosp w/o Surgery"
# label(data$f5hospwosrg_mv)="Missing value verified for: Hosp w/o Surgery"
# label(data$f5hospwoexp)="Hosp w/o Surgery - Exp"
# label(data$f5hospwoexp_mv)="Missing value verified for: Hosp w/o Surgery - Exp"
# label(data$f5hospwsrg)="Hosp w/Surgery"
# label(data$f5hospwsrg_mv)="Missing value verified for: Hosp w/Surgery"
# label(data$f5hospwexp)="Hosp w/Surgery - Exp"
# label(data$f5hospwexp_mv)="Missing value verified for: Hosp w/Surgery - Exp"
# label(data$f5medprbsdoc)="Any Docmntd Meds Prob since Last Phys?"
# label(data$f5medprbsdoc_mv)="Missing value verified for: Any Docmntd Meds Prob since Last Phys?"
# label(data$f5crnryartds)="Coronary artery disease?"
# label(data$f5crnryartds_mv)="Missing value verified for: Coronary artery disease?"
# label(data$f5cancer)="Cancer?"
# label(data$f5cancer_mv)="Missing value verified for: Cancer?"
# label(data$f5cancerexp)="Cancer - Explain"
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
# label(data$f5gnitrnyexp)="GD Explain"
# label(data$f5gnitrnyexp_mv)="Missing value verified for: GD Explain"
# label(data$f5lungdsz)="Lung Disease?"
# label(data$f5lungdsz_mv)="Missing value verified for: Lung Disease?"
# label(data$f5majsurg)="Major surgery?"
# label(data$f5majsurg_mv)="Missing value verified for: Major surgery?"
# label(data$f5majsurgexp)="MS Explain"
# label(data$f5majsurgexp_mv)="Missing value verified for: MS Explain"
# label(data$f5othrmeddia)="Other medical diagnosis?"
# label(data$f5othrmeddia_mv)="Missing value verified for: Other medical diagnosis?"
# label(data$f5othrmedexp)="OMD Explain"
# label(data$f5othrmedexp_mv)="Missing value verified for: OMD Explain"
# label(data$f5lungs)="Lungs"
# label(data$f5lungs_mv)="Missing value verified for: Lungs"
# label(data$f5lungsexp)="Explain"
# label(data$f5lungsexp_mv)="Missing value verified for: Explain"
# label(data$f5heart)="Heart"
# label(data$f5heart_mv)="Missing value verified for: Heart"
# label(data$f5heartexp)="Explain"
# label(data$f5heartexp_mv)="Missing value verified for: Explain"
# label(data$f5skin)="Skin"
# label(data$f5skin_mv)="Missing value verified for: Skin"
# label(data$f5skinexp)="Explain"
# label(data$f5skinexp_mv)="Missing value verified for: Explain"
# label(data$f5edema)="Edema"
# label(data$f5edema_mv)="Missing value verified for: Edema"
# label(data$f5edemaexp)="Explain"
# label(data$f5edemaexp_mv)="Missing value verified for: Explain"
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
# label(data$f6othercast_exp)="Other (specify)"
# label(data$f6othercast_exp_mv)="Missing value verified for: Other (specify)"
# label(data$f6wbc)="WBC"
# label(data$f6wbc_mv)="Missing value verified for: WBC"
# label(data$f6rbc)="RBC"
# label(data$f6rbc_mv)="Missing value verified for: RBC"
# label(data$f6epithcells)="Epithelial cells"
# label(data$f6epithcells_mv)="Missing value verified for: Epithelial cells"
# label(data$f6bacteria)="Bacteria"
# label(data$f6bacteria_mv)="Missing value verified for: Bacteria"
# label(data$f6renalclrnc_complete)="Complete?"
# label(data$f8visitdate)="Date sample collected"
# label(data$f8visitdate_mv)="Missing value verified for: Date sample collected"
# label(data$f8analdate)="Date sample analyzed"
# label(data$f8analdate_mv)="Missing value verified for: Date sample analyzed"
# label(data$f8glucose)="Glucose (mg/dL)"
# label(data$f8glucose_mv)="Missing value verified for: Glucose (mg/dL)"
# label(data$f8bun)="BUN (mg/dL)"
# label(data$f8bun_mv)="Missing value verified for: BUN (mg/dL)"
# label(data$f8creatinine)="Creatinine (mg/dL)"
# label(data$f8creatinine_mv)="Missing value verified for: Creatinine (mg/dL)"
# label(data$f8sodium)="Sodium (mmol/L)"
# label(data$f8sodium_mv)="Missing value verified for: Sodium (mmol/L)"
# label(data$f8potassium)="Potassium (mmol/L)"
# label(data$f8potassium_mv)="Missing value verified for: Potassium (mmol/L)"
# label(data$f8chloride)="Chloride (mmol/L)"
# label(data$f8chloride_mv)="Missing value verified for: Chloride (mmol/L)"
# label(data$f8calcium)="Calcium (mg/dL)"
# label(data$f8calcium_mv)="Missing value verified for: Calcium (mg/dL)"
# label(data$f8phosphorus)="Phosphorus (mg/dL)"
# label(data$f8phosphorus_mv)="Missing value verified for: Phosphorus (mg/dL)"
# label(data$f8totlprotn)="Total protein (gm/dL)"
# label(data$f8totlprotn_mv)="Missing value verified for: Total protein (gm/dL)"
# label(data$f8albumin)="Albumin (gm/dL)"
# label(data$f8albumin_mv)="Missing value verified for: Albumin (gm/dL)"
# label(data$f8cholestero)="Cholesterol (mg/dL)"
# label(data$f8cholestero_mv)="Missing value verified for: Cholesterol (mg/dL)"
# label(data$f8triglyceri)="Triglycerides (mg/dL)"
# label(data$f8triglyceri_mv)="Missing value verified for: Triglycerides (mg/dL)"
# label(data$f8totlbilrbn)="Total bilirubin (mg/dL)"
# label(data$f8totlbilrbn_mv)="Missing value verified for: Total bilirubin (mg/dL)"
# label(data$f8alkp04)="ALK P04 (alkaline phosphatase) (U/L)"
# label(data$f8alkp04_mv)="Missing value verified for: ALK P04 (alkaline phosphatase) (U/L)"
# label(data$f8ast)="AST (SGOT) (U/L)"
# label(data$f8ast_mv)="Missing value verified for: AST (SGOT) (U/L)"
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
# label(data$f9lymphinst)="Lymph, Inst%"
# label(data$f9lymphinst_mv)="Missing value verified for: Lymph, Inst%"
# label(data$f9monoinst)="Mono, Inst%"
# label(data$f9monoinst_mv)="Missing value verified for: Mono, Inst%"
# label(data$f9neutinst)="Neut, Inst%"
# label(data$f9neutinst_mv)="Missing value verified for: Neut, Inst%"
# label(data$f9eosinst)="Eos, Inst%"
# label(data$f9eosinst_mv)="Missing value verified for: Eos, Inst%"
# label(data$f9basoinst)="Baso, Inst%"
# label(data$f9basoinst_mv)="Missing value verified for: Baso, Inst%"
# label(data$f9gran)="Gran, %"
# label(data$f9gran_mv)="Missing value verified for: Gran, %"
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
# label(data$f13reactionspecify)="Other (specify)"
# label(data$f13reactionspecify_mv)="Missing value verified for: Other (specify)"
# label(data$f13datereactn)="Date of reaction"
# label(data$f13datereactn_mv)="Missing value verified for: Date of reaction"
# label(data$f13stoppnt_complete)="Complete?"
# label(data$f7qcnumber)="Quality control number"
# label(data$f7qcnumber_mv)="Missing value verified for: Quality control number"
# label(data$f7visitdate)="Date of sample collection"
# label(data$f7visitdate_mv)="Missing value verified for: Date of sample collection"
# label(data$f7tstintrvl)="Test interval"
# label(data$f7tstintrvl_mv)="Missing value verified for: Test interval"
# label(data$f7typesample)="Type of sample"
# label(data$f7typesample_mv)="Missing value verified for: Type of sample"
# label(data$f7labqc_complete)="Complete?"
# label(data$f12expdatevis)="Expected date of visit"
# label(data$f12expdatevis_mv)="Missing value verified for: Expected date of visit"
# label(data$f12typevisit)="Type of visit"
# label(data$f12typevisit_mv)="Missing value verified for: Type of visit"
# label(data$f12testintrvl)="Test interval"
# label(data$f12testintrvl_mv)="Missing value verified for: Test interval"
# label(data$f12reasnmissd)="Reason visit was missed"
# label(data$f12reasnmissd_mv)="Missing value verified for: Reason visit was missed"
# label(data$f12otherspec)="Reason visit was missed, other (specify)"
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
# label(data$f19renlfuncqcrslt_complete)="Complete?"
# label(data$scr)="Serum Creatinine"
# label(data$scr_mv)="Missing value verified for: Serum Creatinine"
# label(data$gfr)="Glomerular Filtration Rate (GFR)"
# label(data$gfr_mv)="Missing value verified for: Glomerular Filtration Rate (GFR)"
# label(data$hba1a)="Hemoglobin A1a"
# label(data$hba1a_mv)="Missing value verified for: Hemoglobin A1a"
# label(data$hba1b)="Hemoglobin A1b"
# label(data$hba1b_mv)="Missing value verified for: Hemoglobin A1b"
# label(data$hba1c)="Hemoglobin A1c"
# label(data$hba1c_mv)="Missing value verified for: Hemoglobin A1c"
# label(data$hba1o)="Hemoglobin Ao"
# label(data$hba1o_mv)="Missing value verified for: Hemoglobin Ao"
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
# label(data$u0igg_mv)="Missing value verified for: U0 IgG"
# label(data$u0igg_is_below_limit)="U0 IgG is below detection limit of assay?"
# label(data$u0igg_is_below_limit_mv)="Missing value verified for: U0 IgG is below detection limit of assay?"
# label(data$u3flow_gfr)="U3 Flow"
# label(data$u3flow_gfr_mv)="Missing value verified for: U3 Flow"
# label(data$u3gfr)="U3 GFR"
# label(data$u3gfr_mv)="Missing value verified for: U3 GFR"
# label(data$u3use_gfr)="U3 Use"
# label(data$u3use_gfr_mv)="Missing value verified for: U3 Use"
# label(data$u3albumin)="U3 Albumin"
# label(data$u3albumin_mv)="Missing value verified for: U3 Albumin"
# label(data$u3albumin_is_below_limit)="U3 Albumin is below detection limit of assay?"
# label(data$u3albumin_is_below_limit_mv)="Missing value verified for: U3 Albumin is below detection limit of assay?"
# label(data$u3igg)="U3 IgG"
# label(data$u3igg_mv)="Missing value verified for: U3 IgG"
# label(data$u3igg_is_below_limit)="U3 IgG is below detection limit of assay?"
# label(data$u3igg_is_below_limit_mv)="Missing value verified for: U3 IgG is below detection limit of assay?"
# label(data$ualb)="Urine Albumin"
# label(data$ualb_mv)="Missing value verified for: Urine Albumin"
# label(data$ualb_is_below_limit)="Urine Albumin is below detection limit of assay?"
# label(data$ualb_is_below_limit_mv)="Missing value verified for: Urine Albumin is below detection limit of assay?"
# label(data$ucr)="Urine Creatinine"
# label(data$ucr_mv)="Missing value verified for: Urine Creatinine"
# label(data$phoenix_labs_complete)="Complete?"
# label(data$s_gfr)="Glomerular Filtration Rate (GFR)"
# label(data$s_gfr_mv)="Missing value verified for: Glomerular Filtration Rate (GFR)"
# label(data$s_ff)="I/P RATIO GFR"
# label(data$s_ff_mv)="Missing value verified for: I/P RATIO GFR"
# label(data$s_pahclrnc)="PAH Clearance"
# label(data$s_pahclrnc_mv)="Missing value verified for: PAH Clearance"
# label(data$s_gfrb)="GFR Adjusted for BSA "
# label(data$s_gfrb_mv)="Missing value verified for: GFR Adjusted for BSA "
# label(data$s_pahb)="PAH Adjusted for BSA"
# label(data$s_pahb_mv)="Missing value verified for: PAH Adjusted for BSA"
# label(data$s_rpf)="Renal Plasma Flow"
# label(data$s_rpf_mv)="Missing value verified for: Renal Plasma Flow"
# label(data$s_rpfb)="RPF Adjusted for BSA"
# label(data$s_rpfb_mv)="Missing value verified for: RPF Adjusted for BSA"
# label(data$s_p0oncoticprsr)="P0 Oncotic Pressure"
# label(data$s_p0oncoticprsr_mv)="Missing value verified for: P0 Oncotic Pressure"
# label(data$s_p3oncoticprsr)="P3 Oncotic Pressure"
# label(data$s_p3oncoticprsr_mv)="Missing value verified for: P3 Oncotic Pressure"
# label(data$s_pie)="EfferentArterial Oncotic Pressure"
# label(data$s_pie_mv)="Missing value verified for: EfferentArterial Oncotic Pressure"
# label(data$s_pigc)="Glomerular Capillary Oncotic Pressure"
# label(data$s_pigc_mv)="Missing value verified for: Glomerular Capillary Oncotic Pressure"
# label(data$stanford_labs_complete)="Complete?"
#Setting Factors(will create new variable for factors)
data$redcap_event_name.factor = factor(data$redcap_event_name,levels=c("rc_interval_1_arm_1","rc_interval_2_arm_1","rc_interval_3_arm_1","rc_interval_4_arm_1","rc_interval_5_arm_1","rc_interval_6_arm_1"))
data$redcap_repeat_instrument.factor = factor(data$redcap_repeat_instrument,levels=c(""))
data$lastname_mv.factor = factor(data$lastname_mv,levels=c("1","0"))
data$firstname_mv.factor = factor(data$firstname_mv,levels=c("1","0"))
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
data$dtdiabonset_mv.factor = factor(data$dtdiabonset_mv,levels=c("1","0"))
data$demoform_complete.factor = factor(data$demoform_complete,levels=c("0","1","2"))
data$clrncvisitdate_mv.factor = factor(data$clrncvisitdate_mv,levels=c("1","0"))
data$clrncvisittype.factor = factor(data$clrncvisittype,levels=c("3","5"))
data$clrncvisittype_mv.factor = factor(data$clrncvisittype_mv,levels=c("1","0"))
data$clrnctestintvrl_mv.factor = factor(data$clrnctestintvrl_mv,levels=c("1","0"))
data$clearancevisit_complete.factor = factor(data$clearancevisit_complete,levels=c("0","1","2"))
data$f4bldprsr.factor = factor(data$f4bldprsr,levels=c("1","2"))
data$f4bldprsr_mv.factor = factor(data$f4bldprsr_mv,levels=c("1","0"))
data$f4frstmsrsys_mv.factor = factor(data$f4frstmsrsys_mv,levels=c("1","0"))
data$f4frstmsrdia_mv.factor = factor(data$f4frstmsrdia_mv,levels=c("1","0"))
data$f4scndmsrsys_mv.factor = factor(data$f4scndmsrsys_mv,levels=c("1","0"))
data$f4scndmsrdia_mv.factor = factor(data$f4scndmsrdia_mv,levels=c("1","0"))
data$f4heartrate_mv.factor = factor(data$f4heartrate_mv,levels=c("1","0"))
data$f4medseffcts.factor = factor(data$f4medseffcts,levels=c("0","1","2"))
data$f4medseffcts_mv.factor = factor(data$f4medseffcts_mv,levels=c("1","0"))
data$f4pregnant.factor = factor(data$f4pregnant,levels=c("1","2"))
data$f4pregnant_mv.factor = factor(data$f4pregnant_mv,levels=c("1","0"))
data$f4medsnow.factor = factor(data$f4medsnow,levels=c("1","2"))
data$f4medsnow_mv.factor = factor(data$f4medsnow_mv,levels=c("1","0"))
data$f4srmpotasm_mv.factor = factor(data$f4srmpotasm_mv,levels=c("1","0"))
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
data$f8totlbilrbn_mv.factor = factor(data$f8totlbilrbn_mv,levels=c("1","0"))
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
data$f13reactionspecify_mv.factor = factor(data$f13reactionspecify_mv,levels=c("1","0"))
data$f13datereactn_mv.factor = factor(data$f13datereactn_mv,levels=c("1","0"))
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
data$s_gfr_mv.factor = factor(data$s_gfr_mv,levels=c("1","0"))
data$s_ff_mv.factor = factor(data$s_ff_mv,levels=c("1","0"))
data$s_pahclrnc_mv.factor = factor(data$s_pahclrnc_mv,levels=c("1","0"))
data$s_gfrb_mv.factor = factor(data$s_gfrb_mv,levels=c("1","0"))
data$s_pahb_mv.factor = factor(data$s_pahb_mv,levels=c("1","0"))
data$s_rpf_mv.factor = factor(data$s_rpf_mv,levels=c("1","0"))
data$s_rpfb_mv.factor = factor(data$s_rpfb_mv,levels=c("1","0"))
data$s_p0oncoticprsr_mv.factor = factor(data$s_p0oncoticprsr_mv,levels=c("1","0"))
data$s_p3oncoticprsr_mv.factor = factor(data$s_p3oncoticprsr_mv,levels=c("1","0"))
data$s_pie_mv.factor = factor(data$s_pie_mv,levels=c("1","0"))
data$s_pigc_mv.factor = factor(data$s_pigc_mv,levels=c("1","0"))
data$stanford_labs_complete.factor = factor(data$stanford_labs_complete,levels=c("0","1","2"))

levels(data$redcap_event_name.factor)=c("RC Interval 1","RC Interval 2","RC Interval 3","RC Interval 4","RC Interval 5","RC Interval 6")
levels(data$redcap_repeat_instrument.factor)=c("")
levels(data$lastname_mv.factor)=c("Yes","No")
levels(data$firstname_mv.factor)=c("Yes","No")
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
levels(data$dtdiabonset_mv.factor)=c("Yes","No")
levels(data$demoform_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$clrncvisitdate_mv.factor)=c("Yes","No")
levels(data$clrncvisittype.factor)=c("Clearance","Follow-up")
levels(data$clrncvisittype_mv.factor)=c("Yes","No")
levels(data$clrnctestintvrl_mv.factor)=c("Yes","No")
levels(data$clearancevisit_complete.factor)=c("Incomplete","Unverified","Complete")
levels(data$f4bldprsr.factor)=c("Right","Left")
levels(data$f4bldprsr_mv.factor)=c("Yes","No")
levels(data$f4frstmsrsys_mv.factor)=c("Yes","No")
levels(data$f4frstmsrdia_mv.factor)=c("Yes","No")
levels(data$f4scndmsrsys_mv.factor)=c("Yes","No")
levels(data$f4scndmsrdia_mv.factor)=c("Yes","No")
levels(data$f4heartrate_mv.factor)=c("Yes","No")
levels(data$f4medseffcts.factor)=c("Unknown","No","Yes")
levels(data$f4medseffcts_mv.factor)=c("Yes","No")
levels(data$f4pregnant.factor)=c("No","Yes")
levels(data$f4pregnant_mv.factor)=c("Yes","No")
levels(data$f4medsnow.factor)=c("No","Yes")
levels(data$f4medsnow_mv.factor)=c("Yes","No")
levels(data$f4srmpotasm_mv.factor)=c("Yes","No")
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
levels(data$f8totlbilrbn_mv.factor)=c("Yes","No")
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
levels(data$f13reactionspecify_mv.factor)=c("Yes","No")
levels(data$f13datereactn_mv.factor)=c("Yes","No")
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
levels(data$s_gfr_mv.factor)=c("Yes","No")
levels(data$s_ff_mv.factor)=c("Yes","No")
levels(data$s_pahclrnc_mv.factor)=c("Yes","No")
levels(data$s_gfrb_mv.factor)=c("Yes","No")
levels(data$s_pahb_mv.factor)=c("Yes","No")
levels(data$s_rpf_mv.factor)=c("Yes","No")
levels(data$s_rpfb_mv.factor)=c("Yes","No")
levels(data$s_p0oncoticprsr_mv.factor)=c("Yes","No")
levels(data$s_p3oncoticprsr_mv.factor)=c("Yes","No")
levels(data$s_pie_mv.factor)=c("Yes","No")
levels(data$s_pigc_mv.factor)=c("Yes","No")
levels(data$stanford_labs_complete.factor)=c("Incomplete","Unverified","Complete")

ficoll = data
rm(data)