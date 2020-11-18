library(readxl)
library(stringi)

# read in separate sheets
pft <- read_xlsx("H:\\Endocrinology\\Chan\\Functional data analysis\\Prospective data\\Reconsent Preliminary Data.xlsx", 
                sheet="PFTs")
admit <- read_xlsx("H:\\Endocrinology\\Chan\\Functional data analysis\\Prospective data\\Reconsent Preliminary Data.xlsx", 
                   sheet="AdmitsAndPEX")
glyc <- read_xlsx("H:\\Endocrinology\\Chan\\Functional data analysis\\Prospective data\\Reconsent Preliminary Data.xlsx", 
          sheet="OGTTandA1C", col_types = c("text","date","date","text","numeric",
                                            "numeric","numeric","numeric"))
mod <- read_xlsx("H:\\Endocrinology\\Chan\\Functional data analysis\\Prospective data\\Reconsent Preliminary Data.xlsx", 
          sheet="Modulators")

# deal with duplicate patients
# first need to fix typo
fix_ids <- function(df) {
  temp <- df
  temp$ID[temp$ID=="017CA and 097"] <- "017CC and 097"
  dups <- temp[grepl("and",temp$ID,fixed=TRUE),]
  temp <- temp[!grepl("and",temp$ID,fixed=TRUE),]
  dups1 <- dups
  dups2 <- dups
  # in dup1, keep first ID
  dups1$ID <- gsub( " .*$", "", dups1$ID)
  # in dup2, keep second ID
  dups2$ID <- sub(".*d", "", dups2$ID)
  dups2$ID <- sub(".* ", "", dups2$ID)
  df <- rbind(temp,dups1,dups2)
  return(df)
}

pft <- fix_ids(pft)
admit <- fix_ids(admit)
glyc <- fix_ids(glyc)
mod <- fix_ids(mod)

# now do the other data cleaning steps in fpca fdapace file, fixing dates, etc.
# make sure IDs match across files
pft$subjectid <- pft$ID
admit$subjectid <- admit$ID
glyc$subjectid <- glyc$ID
mod$subjectid <- mod$ID
pft$ID <- NULL
admit$ID <- NULL
glyc$ID <- NULL
mod$ID <- NULL

for (i in 1:nrow(glyc)) {
  if (!grepl("/",glyc$`Date of Test`[i],fixed=TRUE)) {
    glyc$testdate[i] <- as.Date(as.numeric(glyc$`Date of Test`[i]), origin = "1899-12-30")
  }
  else{
    glyc$testdate[i] <- as.Date(glyc$`Date of Test`[i],format="%m/%d/%Y")
  }
}
glyc$testdate <- as.Date(glyc$testdate, origin = "1970-01-01")
glyc$`Date of Test` <- NULL

mod$subjectid <- stri_pad(mod$subjectid,width = 3,pad="0")

# write.csv(pft,"H:\\Endocrinology\\Chan\\Functional data analysis\\Prospective data\\pft_clean.csv")
# write.csv(admit,"H:\\Endocrinology\\Chan\\Functional data analysis\\Prospective data\\admit_clean.csv")
# write.csv(glyc,"H:\\Endocrinology\\Chan\\Functional data analysis\\Prospective data\\glyc_clean.csv")
# write.csv(mod,"H:\\Endocrinology\\Chan\\Functional data analysis\\Prospective data\\mod_clean.csv")
