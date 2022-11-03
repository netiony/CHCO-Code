harmonize_data = function(new_export = T,data_dir = "~/Documents/Work/CHCO/Petter Bjornstad/Data Harmonization/Data Exports"){
  library(knitr)
  library(childsds)
  library(tidyverse)
  library(parsedate)
  source("~/GitHub/shared-resources/Data Cleaning/Calculated Variables/eGFR.R")
  source("~/GitHub/shared-resources/Data Cleaning/Calculated Variables/hemodynamics.R")
  # Re-import from REDCap if necessary
  if(new_export){
    library(redcapAPI)
    # API import
    tokens <- read.csv("/Users/timvigers/Documents/Work/CHCO/Petter Bjornstad/Data Harmonization/api_tokens.csv")
    uri <- "https://redcap.ucdenver.edu/api/"
    renalheir <- exportRecords(
      redcapConnection(
        url = uri, token = tokens$Token[tokens$Study == "Renal-HEIR"]
      ),
      labels = F
    )
    renalheir$study <- "RENAL-HEIR"
    renalheir$visit <- "Baseline"
    renalheir[renalheir == -99] <- NA
    renalheir[renalheir == -999] <- NA
    renalheir[renalheir == -9999] <- NA
    write.csv(renalheir,file = paste0(data_dir,"/renalheir.csv"),
              row.names = F,na = "")
    
    penguin <- exportRecords(
      redcapConnection(
        url = uri, token = tokens$Token[tokens$Study == "PENGUIN"]
      ),
      labels = F
    )
    penguin$study <- "PENGUIN"
    penguin$visit <- "Baseline"
    penguin[,c("co_enroll","co_enroll_id")] = NA
    penguin[penguin == -99] <- NA
    penguin[penguin == -999] <- NA
    penguin[penguin == -9999] <- NA
    write.csv(penguin,file = paste0(data_dir,"/penguin.csv"),
              row.names = F,na = "")
    
    crocodile <- exportRecords(
      redcapConnection(
        url = uri, token = tokens$Token[tokens$Study == "CROCODILE"]
      ),
      labels = F
    )
    crocodile$study <- "CROCODILE"
    crocodile$visit <- "Baseline"
    crocodile[,c("co_enroll","co_enroll_id")] = NA
    crocodile[crocodile == -99] <- NA
    crocodile[crocodile == -999] <- NA
    crocodile[crocodile == -9999] <- NA
    write.csv(crocodile,file = paste0(data_dir,"/crocodile.csv"),
              row.names = F,na = "")
    
    coffee <- exportRecords(
      redcapConnection(
        url = uri, token = tokens$Token[tokens$Study == "COFFEE"]
      ),
      labels = F
    )
    coffee$study <- "COFFEE"
    coffee$visit <- "Baseline"
    coffee[,c("co_enroll","co_enroll_id")] = NA
    coffee[coffee == -99] <- NA
    coffee[coffee == -999] <- NA
    coffee[coffee == -9999] <- NA
    write.csv(coffee,file = paste0(data_dir,"/coffee.csv"),
              row.names = F,na = "")
    
    casper <- exportRecords(
      redcapConnection(
        url = uri, token = tokens$Token[tokens$Study == "CASPER"]
      ),
      labels = F
    )
    casper$study <- "CASPER"
    casper$visit <- "Baseline"
    casper[,c("co_enroll","co_enroll_id")] = NA
    casper[casper == -99] <- NA
    casper[casper == -999] <- NA
    casper[casper == -9999] <- NA
    write.csv(casper,file = paste0(data_dir,"/casper.csv"),
              row.names = F,na = "")
    
    improve <- exportRecords(
      redcapConnection(
        url = uri, token = tokens$Token[tokens$Study == "IMPROVE"]
      ),
      labels = F
    )
    
    improve$study <- "IMPROVE"
    improve$visit=as.character(improve$study_visit)
    improve$visit[is.na(improve$visit)] = "Baseline"
    improve[improve == -99] <- NA
    improve[improve == -999] <- NA
    improve[improve == -9999] <- NA
    write.csv(improve,file = paste0(data_dir,"/improve.csv"),
              row.names = F,na = "")
  }
  # Read in with all columns as character
  renalheir = read.csv(paste0(data_dir,"/renalheir.csv"),
                       colClasses = "character",na.strings = "")
  penguin = read.csv(paste0(data_dir,"/penguin.csv"),
                       colClasses = "character",na.strings = "")
  improve = read.csv(paste0(data_dir,"/improve.csv"),
                       colClasses = "character",na.strings = "")
  crocodile = read.csv(paste0(data_dir,"/crocodile.csv"),
                       colClasses = "character",na.strings = "")
  coffee = read.csv(paste0(data_dir,"/coffee.csv"),
                       colClasses = "character",na.strings = "")
  casper = read.csv(paste0(data_dir,"/casper.csv"),
                       colClasses = "character",na.strings = "")
  # Fix ID column name
  improve$record_id = improve$subject_id
  renalheir$record_id = renalheir$subject_id
  coffee$record_id = coffee$subject_id
  casper$record_id = casper$subject_id
  # Base dataframe
  base_vars = c("record_id","study","co_enroll","co_enroll_id","visit")
  df = do.call(rbind,
               list(renalheir[,base_vars],penguin[,base_vars],improve[,base_vars],
               crocodile[,base_vars],coffee[,base_vars],casper[,base_vars]))
  # Add outcomes one at a time
  renalheir$hba1c
  crocodile 
  

  
  # Fill
  fill_vars <- c("group", "dob", "gender", "race", "ethnicity", "diagnosis_date","sglt2i")
  df <- df %>%
    group_by(subject_id) %>%
    fill(all_of(fill_vars), .direction = "downup")
  # Only T2D are on SGLT2i
  df$sglt2i[df$group != "T2D"] = "No"
  
  ###############################################################################
  # Final formatting and calculated fields
  ###############################################################################
  
  # HbA1c
  df$hba1c = coalesce(df$hba1c,df$mmtt_hba1c_base)
  
  # Age
  df$age_clamp = round(as.numeric(difftime(df$clamp_date,df$dob,units = "days"))/365.25)
  df$age_mri = round(as.numeric(difftime(df$mri_date,df$dob,units = "days"))/365.25)
  df$age = coalesce(df$age_consent,df$age_biopsy,df$age_clamp,df$age_mri)
  
  # BMI
  df$bmi = coalesce(df$screen_bmi,df$vitals_bmi,df$clamp_bmi,df$mri_bmi)
  # BMI percentile
  ## Excluding adults
  df$bmi_z <- sds(
    value = df$bmi,
    age = df$age,
    sex = df$gender, male = "Male", female = "Female",
    item = "bmi", type = "SDS",
    ref = cdc.ref
  )
  df$bmi_percentile <- sds(
    value = df$bmi,
    age = df$age,
    sex = df$gender, male = "Male", female = "Female",
    item = "bmi", type = "perc",
    ref = cdc.ref
  )
  
  # Diabetes duration
  df$disease_duration = round(as.numeric(difftime(df$diagnosis_date,df$dob,units = "days"))/365.25)
  
  # Various eGFRs
  df <- data.frame(cbind(df, egfr_calc(
    age = df$age, serum_creatinine = df$serum_creatinine,
    cystatin_c = df$cystatin_c, height = df$clamp_height, sex = df$gender
  )))
  
  # Hemodynamics
  hemodynamics <- c("Pglo", "Ra", "Re", "RVR", "FF", "RBF")
  
  hemo <- hemodynamics_calc(
    total_protein = df$total_protein,
    gfr = df$gfr, rpf = df$rpf, map = df$map,
    hematocrit_minus_10 = df$hematocrit_minus_10,
    hematocrit_minus_5 = df$hematocrit_minus_5,
    hematocrit_90 = df$hematocrit_90,
    hematocrit_120 = df$hematocrit_120,
    group = df$group
  )
  
  hemo <- hemo %>% rename(
    Pglo = Pglo_abs, Ra = Ra_abs, Re = Re_abs, RVR = RVRAbs,
    FF = FFabs, RBF = RBFabs
  )
  
  df <- data.frame(cbind(df, hemo[, hemodynamics]))
  
  # Co-enroll "No" values
  df$co_enroll[is.na(df$co_enroll)] <- "No"
  
  # Visit names
  df$visit <- factor(df$visit,
                     levels = c(
                       "Baseline", "Pre-Surgery",
                       "3 Months Post-Surgery", "12 Months Post-Surgery"
                     )
  )
  
  # Fix CROCODILE IDs
  df$subject_id[df$study=="CROCODILE"] = paste0("CRC-",str_pad(df$subject_id[df$study=="CROCODILE"],2,"left","0"))
  
  # Sort and return!
  df <- df %>% arrange(study, subject_id, visit)
  return(df)
}
