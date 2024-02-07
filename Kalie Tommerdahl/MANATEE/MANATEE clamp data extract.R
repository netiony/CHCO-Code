library(readxl)
library(reshape2)

# Path to Excel file
file_path <- "/Volumes/Peds Endo/Kalie Tommerdahl/MANATEE/Data raw/MANATEE-T1D_Clamp Spreadsheet 1.xlsx"
sheet_names <- excel_sheets(file_path)

# Function to extract specific cells
excel_extract <- function(file_path, sheet_names, range) {
  # Create an empty list to store values
  values_list <- list()
  
  # Loop through each sheet and read data
  for (i in seq_along(sheet_names)) {
    # Read the Excel file
    sheet_data <- read_excel(file_path, sheet = sheet_names[i], range = range, col_names = FALSE)
    
    # Check if any data is returned
    if (length(sheet_data) > 0) {
      # If data is not empty, assign it to the corresponding sheet name in the list
      values_list[[sheet_names[i]]] <- sheet_data[[1]]
    } else {
      # If the range is empty or contains no data, assign NA to the corresponding sheet name
      values_list[[sheet_names[i]]] <- NA
    }
  }
  
  return(values_list)
}  

# Function to extract numeric columns
excel_extract_numeric_column <- function(file_path, range) {
  # Get sheet names
  sheet_names <- excel_sheets(file_path)
  
  # Initialize an empty list to store vectors of numeric entries from column Q
  values_list <- list()
  
  # Loop through each sheet and extract numeric entries from column Q
  for (sheet_name in sheet_names) {
    options(scipen = 999)
    # Read data from column Q of the current sheet
    sheet_data <- readxl::read_excel(file_path, sheet = sheet_name, range = range, col_names = FALSE, col_types = "numeric")
    
    # Subset sheet_data to numeric entries using grepl
    numeric_entries <- sheet_data[grepl("^\\d+\\.?\\d*$", sheet_data$...1), ]
    
    # Convert numeric entries to a vector
    numeric_vector <- as.vector(numeric_entries$...1)
    
    # Store numeric entries in the values_list
    values_list[[sheet_name]] <- numeric_vector
  }
  
  return(values_list)
}

# Function to extract raw ranges for clamp
excel_extract_ranges <- function(file_path, range) {
  # Get sheet names
  sheet_names <- excel_sheets(file_path)
  
  # Initialize an empty list to store data frames
  combined_data <- NULL
  
  # Loop through each sheet and extract data
  for (sheet_name in sheet_names) {
    options(scipen = 999)
    # Read data from current sheet
    sheet_data <- readxl::read_excel(file_path, sheet = sheet_name, range = range, col_names = TRUE)
    
    # Modify column names: lowercase, replace spaces with underscores, remove parentheses
    colnames(sheet_data) <- tolower(gsub("\\s|/", "_", gsub("\\(|\\)", "", colnames(sheet_data))))
    
    # Filter non-NA study time
    sheet_data <- sheet_data[!is.na(sheet_data$study_time), ]
    
    record_id <- as.numeric(gsub("\\D", "", sheet_name))
    visit_id <- substr(sheet_name, nchar(sheet_name) - 1, nchar(sheet_name))
    
    # Add record_id column
    sheet_data$record_id <- record_id
    sheet_data$visit_id <- visit_id
    
    # Combine data
    if (is.null(combined_data)) {
      combined_data <- sheet_data
    } else {
      combined_data <- rbind(combined_data, sheet_data)
    }
  }
  
  return(combined_data)
}

# Date of clamp
clamp_visit_date <- excel_extract(file_path, sheet_names, range = "C2")

# Extract bg and insulin
bg_insulin <- excel_extract_ranges(file_path, range = "G7:M100")
bg_clean <- bg_insulin %>%
  dplyr::select(record_id, visit_id, patient_glucose_mg_dl, study_time) %>%
  rename("bg" = "patient_glucose_mg_dl") %>%
  mutate(study_time = case_when(study_time == 202 ~ 200,
                                study_time == 232 ~ 230,
                                study_time == 262 ~ 265,
                                T ~ study_time))

bg_clean_wide <- dcast(bg_clean, 
                       record_id + visit_id ~ paste0("bg_", study_time), 
                       value.var = "bg")

# Numeric values from column Q
clamp_numeric <- excel_extract_numeric_column(file_path, range = "Q1:Q100")

clamp_df <- as.data.frame(clamp_numeric,
                          row.names = c("clamp_wt",
                                        "clamp_d20",
                                        "clamp_ht",
                                        "clamp_bsa",
                                        "p1_raw_m",
                                        "p1_raw_leanm",
                                        "p1_gc_m",
                                        "p1_gc_leanm",
                                        "p2_raw_m",
                                        "p2_raw_leanm",
                                        "p2_gc_m",
                                        "p2_gc_leanm"))

# Participant arms
control <- c(6, 16, 20, 21, 22)
treatment <- c(2, 3, 5, 8, 9, 10, 11, 12, 14, 15, 17, 19, 23, 24, 25, 26, 27, 28, 29)

# Combine all extracted data
clean_df <- as.data.frame(t(clamp_df)) %>%
  mutate(record_id = as.numeric(gsub("\\D", "", colnames(clamp_df))),
         visit_id = substr(colnames(clamp_df), nchar(colnames(clamp_df)) - 1, nchar(colnames(clamp_df))),
         clamp_visit_date = clamp_visit_date, 
         clamp_yn = 1,
         clamp_bg = 1,
         redcap_event_name = 
           case_when(visit_id == "BL" & record_id %in% treatment ~ "study_visit_baseli_arm_1",
                     visit_id == "FL" & record_id %in% treatment ~ "study_visit_follow_arm_1",
                     visit_id == "BL" & record_id %in% control ~ "study_visit_baseli_arm_2",
                     visit_id == "FL" & record_id %in% control ~ "study_visit_follow_arm_2")) %>%
  unnest(cols = clamp_visit_date) %>%
  left_join(bg_clean_wide)

clean_df <- clean_df[, 
                     c("record_id", "redcap_event_name", "clamp_visit_date", setdiff(names(clean_df), 
                                                                              c("record_id", "redcap_event_name", "clamp_visit_date")))] %>%
  dplyr::select(-visit_id, -bg_91) %>%
  filter(!is.na(redcap_event_name)) %>%
  as.data.frame()

write.csv(clean_df, "/Volumes/Peds Endo/Kalie Tommerdahl/MANATEE/Data clean/MANATEE-T1D_Clamp Spreadsheet.csv", row.names = F)
