library(dplyr)
library(tidyr)
library(table1)
library(purrr)
library(Hmisc)

dat <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))

# withdrew/LTFU
exclude <- read.csv("/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/PANTHER/Data_Cleaned/panther_withdrew_ltfu_list.csv")$record_id

ier <- dat %>%
  filter(study == "PANTHER") %>%
  dplyr::summarise(dplyr::across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   dplyr::across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  fill(hba1c) %>% fill(acr_u) %>% ungroup() %>%
  dplyr::mutate(race_ethnicity_condensed = case_when(race == "White" & 
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                                     race == "Black or African American" & 
                                                       ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                                     ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                                     T ~ "Not Hispanic or Latino Other"),
                race = case_when(race == "Black/African American & White" ~ "More than one", T ~ race),
                status = case_when(record_id %in% exclude ~ "Removed",
                                   T ~ "Currently enrolled"),
                age = case_when(record_id == "PAN-18-O" ~ 10, T ~ age)) %>% # 10 at 2025, manually editing because no date of consent?
  filter(visit == "baseline") %>%
  filter(record_id %nin% exclude) %>%
  dplyr::select(race, ethnicity, sex, age, status) %>%
  dplyr::rename(Race = race,
                Ethnicity = ethnicity,
                Gender = sex,
                Age = age,
                Status = status) %>%
  dplyr::mutate(Age = floor(Age), 
                "Age Unit" = "Years") %>%
  arrange(Status)

write.csv(ier, "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/PANTHER/Data_Cleaned/2025_panther_ier.csv", row.names = F)

summary_table <- ier %>%
  group_by(Race, Ethnicity, Gender) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = c(Ethnicity, Gender),
    values_from = count,
    values_fill = 0
  )
