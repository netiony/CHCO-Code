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
    "C:/Users/timbv/Documents/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
from data_harmonization import harmonize_data
# Get dataset
df = harmonize_data()
# df.to_csv("~/df.csv", index=False)
# Variable list
vars = ["record_id", "co_enroll_id", "study", "dob", "diabetes_dx_date",
        "race", "ethnicity", "race_ethnicity", "visit", "procedure", "date",
        "age", "sex", "bmi", "creatinine_s",
        "cystatin_c_s", "gfr", "gfr_bsa", "acr_u", "hba1c",
        'dexa_body_fat', 'dexa_fat_kg', 'dexa_trunk_kg',
        'bod_pod_body_fat', 'bod_pod_fat_kg',
        "bold_r_bl_cortex", "bold_r_bl_medulla", "bold_r_bl_kidney",
        "bold_r_pf_cortex", "bold_r_pf_medulla", "bold_r_pf_kidney",
        "bold_l_bl_cortex", "bold_l_bl_medulla", "bold_l_bl_kidney",
        "bold_l_pf_cortex", "bold_l_pf_medulla", "bold_l_pf_kidney",
        "gloms", "gloms_gs", "ifta", "vessels_other", "fia", "glom_tuft_area", "glom_volume_weibel", "glom_volume_wiggins", "glom_volume_con", "mes_matrix_area", "mes_index", "mes_volume_weibel", "mes_volume_wiggins", "mes_volume_con", "glom_nuc_count", "mes_nuc_count", "art_intima", "art_media", "pod_nuc_density", "pod_cell_volume"]

fadhl = df[vars]
fadhl.to_csv("~/fadhl.csv", index=False)
