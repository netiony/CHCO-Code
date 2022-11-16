"""
This code is designed to pull data from the IMPROVE REDCap project and output data in a long format with one row per study procedure per visit.
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
token = tokens.loc[tokens["Study"] == "IMPROVE", "Token"].iloc[0]
proj = redcap.Project(url=uri, token=token)
# Get project metadata
meta = pd.DataFrame(proj.metadata)
# Columns to drop
redcap_cols = ["redcap_event_name",
               "redcap_repeat_instrument", "redcap_repeat_instance"]

# ------------------------------------------------------------------------------
# Demographics
# ------------------------------------------------------------------------------

dem_cols = ["subject_id", "co_enroll_id", "dob", "diagnosis",
            "gender", "race", "ethnicity"]
# Export
demo = pd.DataFrame(proj.export_records(fields=dem_cols,
                                        events=["screening_arm_1"]))
demo["group"] = "Type 2 Diabetes"
demo.drop(redcap_cols, axis=1, inplace=True)
demo.rename({"gender": "sex", "diagnosis": "diabetes_dx_date"},
            inplace=True, axis=1)
# Race columns combined into one
demo = combine_checkboxes(demo, base_name="race", levels=[
                          "American Indian or Alaskan Native", "Asian", "Hawaiian or Pacific Islander", "Black or African American", "White", "Unknown", "Other"])
# Same for ethnicity
demo = combine_checkboxes(demo,
                          base_name="ethnicity",
                          levels=["Hispanic or Latino", "Not Hispanic or Latino", "Unknown/Not Reported"])
# Relevel sex and group
demo["sex"].replace({1: "Male", 2: "Female", 3: "Other",
                     "1": "Male", "2": "Female", "3": "Other"}, inplace=True)

# ------------------------------------------------------------------------------
# Medications
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Physical exam
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "physical_exam", "field_name"]]
phys = pd.DataFrame(proj.export_records(
    fields=var, events=["screening_arm_1"]))
phys["visit"] = "baseline"
phys["procedure"] = "physical_exam"
phys.drop(redcap_cols + ["phys_norm", "phys_no", "breast_tanner",
                         "testicular_volume", "lmp", "screen_bmi_percentile", "activity_factor_male", "activity_factor_female", "schofield_male", "schofield_female"], axis=1, inplace=True)
phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
phys.rename({"sys_bp": "sbp", "dys_bp": "dbp", "waist_circumference": "waistcm",
            "hip_circumference": "hipcm"}, inplace=True, axis=1)

# ------------------------------------------------------------------------------
# Screening labs
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "screening_labs", "field_name"]]
screen = pd.DataFrame(proj.export_records(fields=var,
                                          events=["screening_arm_1"]))
screen.drop(redcap_cols + ['a1c_pre', 'a1c_pre_date', "screen_pregnant"],
            axis=1, inplace=True)
screen.columns = screen.columns.str.replace(
    r"screen_|_of_screen", "", regex=True)
screen.rename({"serum_creatinine": "creatinine_s", "urine_acr": "acr_u",
              "urine_cre": "creatinine_u", "urine_mab": "mab_u"},
              axis=1, inplace=True)
screen["procedure"] = "screening"
screen["visit"] = "baseline"

# ------------------------------------------------------------------------------
# Accelerometry
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "accelerometry", "field_name"]]
var = var + ["study_visit"]
accel = pd.DataFrame(proj.export_records(fields=var))
accel = accel.loc[accel["acc_wear_percent"] != ""]
accel.drop(redcap_cols + ["study_visit_accel"], axis=1, inplace=True)
accel.columns = accel.columns.str.replace(
    r"acc_|accel_", "", regex=True)
accel["procedure"] = "accelerometry"

# ------------------------------------------------------------------------------
# Cardio/Abdominal MRI
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "cardioabdominal_mri", "field_name"]]
var = var + ["study_visit"]
mri = pd.DataFrame(proj.export_records(fields=var))
mri.drop(redcap_cols + ["mri_cardio", "mri_abdo",
                        "mri_aortic", "study_visit_mri"],
         axis=1, inplace=True)
mri = mri.loc[mri["mri_visit_date"] != ""]
mri.columns = mri.columns.str.replace(
    r"mri_|visit_", "", regex=True)
mri["procedure"] = "cardio_abdominal_mri"

# ------------------------------------------------------------------------------
# MMTT + Metabolic Cart
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "mmtt_metabolic_cart", "field_name"]]
mmtt = pd.DataFrame(proj.export_records(fields=var))
mmtt.drop(redcap_cols + ["study_visit_mttt", "mmtt_vitals", "mmtt_pregnant",
                         "mmtt_lmp", "mmtt_brmr", "mmtt_60rmr", "mmtt_base_labs", "mmtt_ffa_labs", "mmtt_insulin",
                         "mmtt_glp1", "mmtt_cpep", "mmtt_yy", "mmtt_glucagon",
                         "mmtt_gluc"],
          axis=1, inplace=True)
mmtt = mmtt.loc[mmtt["mmtt_date"] != ""]
mmtt.columns = mmtt.columns.str.replace(
    r"mmtt_", "", regex=True)
mmtt.rename({"wt": "weight", "ht": "height", "waist": "waistcm",
             "hip": "hipcm", "hr": "pulse", "sys_bp": "sbp", "dia_bp": "dbp"},
            inplace=True, axis=1)
mmtt["procedure"] = "mmtt"

# ------------------------------------------------------------------------------
# DXA
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "body_composition_dxa_bod_pod", "field_name"]]
dxa = pd.DataFrame(proj.export_records(fields=var))
dxa.drop(redcap_cols + ["study_visit_bodycomp", "dxa_complete",
                        "bodpod_complete"],
         axis=1, inplace=True)
dxa = dxa.loc[dxa["bodcomp_date"] != ""]
dxa.columns = dxa.columns.str.replace(
    r"dxa_|bp_|bodcomp_", "", regex=True)

# ------------------------------------------------------------------------------
# Clamp
# ------------------------------------------------------------------------------

var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                            == "clamp", "field_name"]]
clamp = pd.DataFrame(proj.export_records(fields=var))
clamp.drop(redcap_cols + ["study_visit_clamp", "baseline", "fasting_labs",
                          "bg_labs", "ns_bolus", "urine_labs"],
           axis=1, inplace=True)
clamp = clamp.loc[clamp["clamp_date"] != ""]
clamp.columns = clamp.columns.str.replace(
    r"clamp_", "", regex=True)
clamp.to_csv("~/clamp.csv")


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
id_cols = ["subject_id", "study"] + \
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
crocodile.to_csv("crocodile.csv", index=False)
