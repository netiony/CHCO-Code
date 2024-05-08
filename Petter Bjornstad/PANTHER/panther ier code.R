library(dplyr)
library(tidyr)
library(table1)
library(purrr)
library(Hmisc)

dat <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", na.strings = c(" ", "", "-9999",-9999))

ier <- dat %>%
  filter(study == "PANTHER") %>%
  dplyr::summarise(dplyr::across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   dplyr::across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit)) %>%
  fill(hba1c) %>% fill(acr_u) %>% ungroup() %>%
  filter(participation_status != "Removed"|is.na(participation_status)) %>%
  mutate(race_ethnicity_condensed = case_when(race == "White" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
                                              race == "Black or African American" & 
                                                ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
                                              ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
                                              T ~ "Not Hispanic or Latino Other"),
         race = case_when(race == "Black/African American & White" ~ "More than one", T ~ race)) %>%
  filter(visit == "baseline") %>%
  dplyr::select(race, ethnicity, sex, age) %>%
  mutate(age = floor(age), 
         "Age Unit" = "Years") %>%
  dplyr::rename("Race" = "race",
                "Ethnicity" = "ethnicity",
                "Gender" = "sex",
                "Age" = "age")

write.csv(ier, "/Volumes/Peds Endo/Petter Bjornstad/PANTHER/2024_panther_ier.csv", row.names = F)
