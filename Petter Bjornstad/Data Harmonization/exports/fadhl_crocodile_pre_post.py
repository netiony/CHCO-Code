"""
This code is designed pulls our harmonized data and selects the variables Fadhl 
is using in his pre-/post-surgery analysis. 
"""
__author__ = "Tim Vigers"
__credits__ = ["Tim Vigers"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"

import os
os.chdir(
    "/Users/timvigers/Documents/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
from data_harmonization import harmonize_data
# Get dataset
df = harmonize_data()
df.to_csv("~/df.csv", index=False)
# Variable list
vars = ["record_id", "co_enroll_id", "study", "dob", "diabetes_dx_date",
        "race", "ethnicity", "race_ethnicity", "visit", "procedure", "date",
        "sglt2i_ever", "sglt2i_timepoint", "age", "sex", "bmi", "creatinine_s",
        "cystatin_c_s", "gfr", "gfr_bsa", "acr_u", "hba1c",
        'dexa_body_fat', 'dexa_fat_kg', 'dexa_trunk_kg',
        'bod_pod_body_fat', 'bod_pod_fat_kg', 'bod_pod_trunk_kg']

fadhl = df[vars]
fadhl.to_csv("~/fadhl.csv", index=False)
