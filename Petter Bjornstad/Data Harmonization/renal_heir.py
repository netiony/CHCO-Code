"""
This code is designed to pull data from the CROCODILE REDCap project and output data in a "semi-long" format with one row per study procedure, and a visit column for longitudinal clustering when combined with other studies.
"""
__author__ = "Tim Vigers"
__credits__ = ["Tim Vigers"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"

# Libraries
import os
import redcap
import pandas as pd
import numpy as np
os.chdir("/Users/timvigers/Documents/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
from harmonization_functions import combine_checkboxes
from harmonization_functions import find_duplicate_columns

# REDCap project variables
tokens = pd.read_csv(
    "~/Dropbox/Work/CHCO/Petter Bjornstad/Data Harmonization/api_tokens.csv")
uri = "https://redcap.ucdenver.edu/api/"
token = tokens.loc[tokens["Study"] == "Renal-HEIR", "Token"].iloc[0]
proj = redcap.Project(url=uri, token=token)
# Get project metadata
meta = pd.DataFrame(proj.metadata)

# ------------------------------------------------------------------------------
# Demographics
# ------------------------------------------------------------------------------

dem_cols = ["subject_id", "co_enroll_id", "dob", "diagnosis",
            "group", "gender", "race", "ethnicity"]
# Export
demo = pd.DataFrame(proj.export_records(fields=dem_cols))
demo.drop(redcap_cols, axis=1, inplace=True)
demo.rename({"gender": "sex", "diagnosis": "diabetes_dx_date"},
            inplace=True, axis=1)
dem_cols[3] = "diabetes_dx_date"
dem_cols[4] = "sex"
# Race columns combined into one
demo = combine_checkboxes(demo, base_name="race", levels=[
                          "American Indian or Alaskan Native", "Asian", "Hawaiian or Pacific Islander", "Black or African American", "White", "Unknown", "Other"])
# Same for ethnicity
demo = combine_checkboxes(demo,
                          base_name="ethnicity",
                          levels=["Hispanic or Latino", "Not Hispanic or Latino", "Unknown/Not Reported"])
# Relevel sex and group
demo["sex"].replace({1: "Male", 0: "Female", 2: "Other",
                     "1": "Male", "0": "Female", "2": "Other"}, inplace=True)
demo["group"].replace({2: "Type 2 Diabetes", 3: "Obese Control",
                       4: "Lean Control",
                       "2": "Type 2 Diabetes", "3": "Obese Control",
                       "4": "Lean Control"}, inplace=True)

# ------------------------------------------------------------------------------
# Medications
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Physical exam
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "physical_exam", "field_name"]]
phys = pd.DataFrame(proj.export_records(fields=var))
phys["procedure"] = "physical_exam"


phys.drop(["phys_normal", "phys_abnormal"], axis=1, inplace=True)
phys.columns = phys.columns.str.replace(r"phys_", "", regex=True)
phys.rename({"sysbp": "sbp", "diasbp": "dbp"}, inplace=True, axis=1)

# ------------------------------------------------------------------------------
# Screening labs
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "screening_labs", "field_name"]]
screen = pd.DataFrame(proj.export_records(fields=var))
screen.drop(["prescreen_a1c", "prescreen_a1c_date",
            "screen_menstrual", "screen_upt"], axis=1, inplace=True)
screen.columns = screen.columns.str.replace(
    r"labs_|screen_", "", regex=True)
screen.rename({"creat_s": "creatinine_s", "uacr": "acr_u",
              "creat_u": "creatinine_u", "hg": "hemoglobin"}, axis=1, inplace=True)
screen["procedure"] = "screening"

# ------------------------------------------------------------------------------
# Labs
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "study_visit_baseline_vitalslabs", "field_name"]]
labs = pd.DataFrame(proj.export_records(fields=var))
labs.drop(["baseline_vitals", "visit_upt",
          "visit_uptresult", "baseline_labs", "pilabs_yn", "pi_copeptin", "pi_renin", "pi_angiotensin2", "pi_osmo_s", "pi_osmo_u", "pi_lithium_s", "pi_lithium_u", "metabolomics_yn", "kim_yn", "pi_kim_ykl40", "pi_kim_ngal", "pi_kim_kim1", "pi_kim_il18", "pi_kim_tnfr1", "pi_kim_tnfr2"], axis=1, inplace=True)
labs.columns = labs.columns.str.replace(
    r"visit_|bl_", "", regex=True)
labs.rename({"uacr": "acr_u"}, axis=1, inplace=True)
labs["procedure"] = "labs"

# ------------------------------------------------------------------------------
# BOLD/ASL MRI
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "study_visit_boldasl_mri", "field_name"]]
mri = pd.DataFrame(proj.export_records(fields=var))
mri.columns = mri.columns.str.replace(
    r"mri_", "", regex=True)
mri["procedure"] = "bold_mri"

# ------------------------------------------------------------------------------
# DXA Scan
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "study_visit_dxa_scan", "field_name"]]
dxa = pd.DataFrame(proj.export_records(fields=var))
dxa.columns = dxa.columns.str.replace(
    r"dxa_", "", regex=True)
dxa["procedure"] = "dxa"

# ------------------------------------------------------------------------------
# Clamp
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "study_visit_he_clamp", "field_name"]]
clamp = pd.DataFrame(proj.export_records(fields=var))
clamp.drop(["clamp_yn", "clamp_d20", "clamp_ffa",
           "clamp_insulin", "hct_yn", "clamp_bg"], axis=1, inplace=True)
clamp.rename({"clamp_wt": "weight", "clamp_ht": "height"},
             inplace=True, axis=1)
clamp.columns = clamp.columns.str.replace(r"clamp_", "", regex=True)
clamp["procedure"] = "clamp"

# ------------------------------------------------------------------------------
# Renal Clearance Testing
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "study_visit_renal_clearance_testing", "field_name"]]
rct = pd.DataFrame(proj.export_records(fields=var))
rct.drop(["iohexol_yn", "pah_yn"], axis=1, inplace=True)
rct["procedure"] = "renal_clearance_testing"

# ------------------------------------------------------------------------------
# Kidney Biopsy
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "optional_kidney_biopsy_56ba", "field_name"]]
var = var + ["gloms", "gloms_gs", "ifta", "vessels_other", "fia",
             "glom_tuft_area", "glom_volume_weibel", "glom_volume_wiggins",
             "glom_volume_con", "mes_matrix_area",
             "mes_index", "mes_volume_weibel", "mes_volume_wiggins",
             "mes_volume_con", "glom_nuc_count", "mes_nuc_count", "art_intima",
             "art_media", "pod_nuc_density", "pod_cell_volume"]
biopsy = pd.DataFrame(proj.export_records(fields=var))
biopsy.drop([col for col in biopsy.columns if '_yn' in col] +
            [col for col in biopsy.columns if 'procedure_' in col],
            axis=1, inplace=True)
biopsy.columns = biopsy.columns.str.replace(r"bx_", "", regex=True)
biopsy.columns = biopsy.columns.str.replace(r"labs_", "", regex=True)
biopsy.columns = biopsy.columns.str.replace(r"vitals_", "", regex=True)
biopsy["procedure"] = "kidney_biopsy"

# ------------------------------------------------------------------------------
# PET scan
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "optional_pet_scan", "field_name"]]
pet = pd.DataFrame(proj.export_records(fields=var))
pet.drop(["petcon_yn"], axis=1, inplace=True)
pet.columns = pet.columns.str.replace(r"pet_", "", regex=True)
pet["procedure"] = "pet_scan"

# MERGE
crocodile = pd.merge(phys, screen, how="outer")
crocodile = pd.merge(crocodile, labs, how="outer")
crocodile = pd.merge(crocodile, mri, how="outer")
crocodile = pd.merge(crocodile, dxa, how="outer")
crocodile = pd.merge(crocodile, clamp, how="outer")
crocodile = pd.merge(crocodile, rct, how="outer")
crocodile = pd.merge(crocodile, biopsy, how="outer")
crocodile = pd.merge(crocodile, pet, how="outer")
crocodile = pd.merge(crocodile, demo, how="outer")
# REORGANIZE
crocodile["visit"] = "baseline"
crocodile["study"] = "CROCODILE"
id_cols = ["subject_id", "co_enroll_id", "study"] + \
    dem_cols[1:] + ["visit", "procedure", "date"]
other_cols = crocodile.columns.difference(id_cols, sort=False).tolist()
other_cols.sort()
crocodile = crocodile[id_cols + other_cols]
# SORT
crocodile.sort_values(["subject_id", "date", "procedure"], inplace=True)
# Check for duplicated column names
dups = find_duplicate_columns(crocodile)
dups.to_csv("~/croc_duplicate_columns.csv", index=False)
# Print final data
crocodile.to_csv("~/crocodile.csv", index=False)
