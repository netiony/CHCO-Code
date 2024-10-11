# Define the main folder path
main_folder <- "/Volumes/Peds Endo/Petter Bjornstad/CROCODILE/Data_Raw/CROCODILE_RawData"

# Get a list of all subfolders within the main folder
subfolders <- list.dirs(main_folder, recursive = FALSE)

# Create an empty data frame to store the information
file_info <- data.frame(
  Subfolder = character(),
  Barcodes = character(),
  Features = character(),
  Matrix = character(),
  stringsAsFactors = FALSE
)

# Loop through each subfolder
for (subfolder in subfolders) {
  # Extract the subfolder name and ensure underscores instead of hyphens
  subfolder_name <- gsub("-", "_", basename(subfolder))
  
  # List the .gz files in the subfolder
  files <- list.files(subfolder, pattern = "*.gz", full.names = TRUE)
  
  # Create empty variables for barcodes, features, and matrix
  barcodes_name <- NA
  features_name <- NA
  matrix_name <- NA
  
  # Loop through each file and rename it
  for (file in files) {
    # Get the file name without extension
    file_name <- tools::file_path_sans_ext(basename(file))
    
    # Replace repetitive occurrences of the sample ID in the file name
    new_file_name <- paste0(subfolder_name, "_", gsub(paste0(subfolder_name, "_"), "", file_name), ".gz")
    
    # Define the new file path
    new_file_path <- file.path(subfolder, new_file_name)
    
    # Rename the file
    file.rename(file, new_file_path)
    
    # Assign the new file name to the appropriate variable
    if (grepl("barcodes", file_name)) {
      barcodes_name <- new_file_name
    } else if (grepl("features", file_name)) {
      features_name <- new_file_name
    } else if (grepl("matrix", file_name)) {
      matrix_name <- new_file_name
    }
  }
  
  # Add the information to the data frame
  file_info <- rbind(file_info, data.frame(
    Subfolder = subfolder_name,
    Barcodes = barcodes_name,
    Features = features_name,
    Matrix = matrix_name,
    stringsAsFactors = FALSE
  ))
}

# Save the data frame to a CSV file
write.csv(file_info, file = "/Volumes/Peds Endo/Petter Bjornstad/CROCODILE/Data_Raw/CROCODILE_RawData/GEO_individual_file_info.csv", row.names = FALSE)

# Optional: Print a message when done
cat("Files have been renamed and information saved to CSV.")

########################################################################
########################### MD5 checksums ##############################
########################################################################

# Load the required package
library(digest)

# Get a list of all subfolders within the main folder
subfolders <- list.dirs(main_folder, recursive = FALSE)

# Create an empty data frame to store the information
file_info <- data.frame(
  Subfolder = character(),
  FileName = character(),
  MD5_Checksum = character(),
  stringsAsFactors = FALSE
)

# Loop through each subfolder
for (subfolder in subfolders) {
  # List the .gz files in the subfolder
  files <- list.files(subfolder, pattern = "*.gz", full.names = TRUE)
  
  # Loop through each file and calculate its MD5 checksum
  for (file in files) {
    # Get the MD5 checksum of the file
    checksum <- digest(file, algo = "md5", file = T)
    print(file)
    print(checksum)
    # Add the information to the data frame
    file_info <- rbind(file_info, data.frame(
      Subfolder = basename(subfolder),
      FileName = basename(file),
      MD5_Checksum = checksum,
      stringsAsFactors = FALSE
    ))
  }
}

# Save the data frame with checksums to a CSV file
write.csv(file_info, file = "/Volumes/Peds Endo/Petter Bjornstad/CROCODILE/Data_Raw/CROCODILE_RawData/md5_checksums.csv", row.names = FALSE)

# Optional: Print a message when done
cat("MD5 checksums calculated and saved to CSV.")

########################################################################
############################## RDS -> TSV ##############################
########################################################################

# Load the combined RDS file
rds_file <- "/Volumes/Peds Endo/Petter Bjornstad/scRNA/data_clean/seurat_data_CRC.RDS"
seurat_object <- readRDS(rds_file)  # Assuming your RDS file contains a list of participant data

# Load necessary libraries
library(Seurat)
library(Matrix)
library(digest)

# Ensure the output directory exists
output_dir <- "/Volumes/Peds Endo/Petter Bjornstad/CROCODILE/Data_Raw/CROCODILE_ProcessedData"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Extract participant names (assuming participant IDs are stored in the 'orig.ident' metadata column)
participants <- unique(seurat_object@meta.data$orig.ident)

# Initialize a data frame to store participant IDs, file names, and checksums
file_info <- data.frame(
  Participant_ID = character(),
  MTX_File = character(),
  MTX_MD5 = character(),
  Barcodes_File = character(),
  Barcodes_MD5 = character(),
  Features_File = character(),
  Features_MD5 = character(),
  stringsAsFactors = FALSE
)

# Function to save data in MTX format and barcodes/features in TSV format, and compute MD5 checksums
save_data <- function(counts_matrix, barcodes, features, participant_id, output_dir) {
  # Replace hyphens with underscores in participant ID
  participant_id_fixed <- gsub("-", "_", participant_id)
  
  # Create a new folder for this participant
  participant_folder <- file.path(output_dir, participant_id_fixed)
  dir.create(participant_folder, showWarnings = FALSE)
  
  # Set up file paths with "_processed" suffix
  mtx_path <- file.path(participant_folder, paste0(participant_id_fixed, "_matrix_processed.mtx"))
  barcodes_path <- file.path(participant_folder, paste0(participant_id_fixed, "_barcodes_processed.tsv"))
  features_path <- file.path(participant_folder, paste0(participant_id_fixed, "_features_processed.tsv"))
  
  # Save counts matrix as MTX (Matrix Market format)
  Matrix::writeMM(counts_matrix, file = mtx_path)
  
  # Save barcodes (cell names)
  write.table(barcodes, file = barcodes_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Save features (gene names)
  write.table(features, file = features_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Calculate MD5 checksums
  mtx_md5 <- digest(mtx_path, algo = "md5", file = TRUE)
  barcodes_md5 <- digest(barcodes_path, algo = "md5", file = TRUE)
  features_md5 <- digest(features_path, algo = "md5", file = TRUE)
  
  # Return the file names (without paths) and MD5 checksums for recording
  return(list(
    mtx_file = basename(mtx_path), mtx_md5 = mtx_md5,
    barcodes_file = basename(barcodes_path), barcodes_md5 = barcodes_md5,
    features_file = basename(features_path), features_md5 = features_md5
  ))
}

# Loop through each participant and save their data
for (participant_id in participants) {
  # Subset the Seurat object for the specific participant
  participant_subset <- subset(seurat_object, subset = orig.ident == participant_id)
  
  # Extract the counts matrix for the participant
  counts_matrix <- GetAssayData(participant_subset, slot = "counts")
  
  # Extract barcodes (cell names) and features (gene names)
  barcodes <- colnames(counts_matrix)
  features <- rownames(counts_matrix)
  
  # Save the counts matrix, barcodes, and features, and calculate MD5 checksums
  file_paths <- save_data(counts_matrix, barcodes, features, participant_id, output_dir)
  
  # Add participant ID, file names, and MD5 checksums to the data frame
  file_info <- rbind(file_info, data.frame(
    Participant_ID = gsub("-", "_", participant_id),
    MTX_File = file_paths$mtx_file,
    MTX_MD5 = file_paths$mtx_md5,
    Barcodes_File = file_paths$barcodes_file,
    Barcodes_MD5 = file_paths$barcodes_md5,
    Features_File = file_paths$features_file,
    Features_MD5 = file_paths$features_md5,
    stringsAsFactors = FALSE
  ))
}

# Save the file info data frame to a CSV file, including the MD5 checksums
write.csv(file_info, file = file.path(output_dir, "processed_participant_file_info_with_md5.csv"), row.names = FALSE)

# Optional: Print a message when done
cat("All participant data, barcodes, and features saved to MTX and TSV formats, and MD5 checksums compiled into CSV.")