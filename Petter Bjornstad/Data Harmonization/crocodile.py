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
import redcap
import pandas as pd
import numpy as np
from combine_checkboxes import combine_checkboxes

# For testing
import os
os.chdir("/Users/timvigers/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")


# REDCap project variables
tokens = pd.read_csv(
    "~/Dropbox/Work/CHCO/Petter Bjornstad/Data Harmonization/api_tokens.csv")
uri = "https://redcap.ucdenver.edu/api/"
token = tokens.loc[tokens["Study"] == "CROCODILE", "Token"].iloc[0]
proj = redcap.Project(url=uri, token=token)
# Get project metadata
meta = pd.DataFrame(proj.metadata)
df = pd.DataFrame(proj.export_records())

# ------------------------------------------------------------------------------
# Demographics
# ------------------------------------------------------------------------------

dem_cols = ["record_id", "dob", "diabetes_dx_date",
            "group", "sex", "race", "ethnicity"]
# Export
demo = proj.export_records(format_type="df", fields=dem_cols)
# Race columns combined into one
demo = combine_checkboxes(demo, base_name="race",
                          levels=[
                              'American Indian/Alaska Native', 'Asian', 'Hawaiian/Pacific Islander', 'Black/African American', 'White', 'Other', 'Unknown'])
# Same for ethnicity
demo = combine_checkboxes(demo,
                          base_name="ethnicity",
                          levels=["Hispanic or Latino", "Not Hispanic or Latino", "Unknown/Not Reported"])
# Relevel sex and group
demo["sex"].replace({1: "Male", 2: "Female", 3: "Other"}, inplace=True)
demo["group"].replace({1: "Type 1 Diabetes", 2: "Control"}, inplace=True)

# ------------------------------------------------------------------------------
# Medications
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Physical exam
# ------------------------------------------------------------------------------

var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                           == "physical_exam", "field_name"]]
phys = df.filter(var, axis=1)
phys["procedure"] = "physical_exam"
phys.drop(["phys_normal", "phys_abnormal"], axis=1, inplace=True)
phys.columns = phys.columns.str.replace(r"phys_", "")

# ------------------------------------------------------------------------------
# Screening labs
# ------------------------------------------------------------------------------

var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                           == "screening_labs", "field_name"]]
screen = df.filter(var, axis=1)
screen.drop(["prescreen_a1c", "prescreen_a1c_date",
            "screen_menstrual", "screen_upt"], axis=1, inplace=True)
screen.columns = screen.columns.str.replace(
    r"labs_|screen_", "", regex=True)
screen.rename({"creat_s": "creatinine_s",
              "creat_u": "creatinine_u"}, axis=1, inplace=True)
screen["procedure"] = "screening"

# ------------------------------------------------------------------------------
# Labs
# ------------------------------------------------------------------------------

var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                           == "study_visit_baseline_vitalslabs", "field_name"]]
labs = df.filter(var, axis=1)
labs.drop(["baseline_vitals", "visit_upt",
          "visit_uptresult", "baseline_labs", "pilabs_yn", "pi_copeptin", "pi_renin", "pi_angiotensin2", "pi_osmo_s", "pi_osmo_u", "pi_lithium_s", "pi_lithium_u", "metabolomics_yn", "kim_yn", "pi_kim_ykl40", "pi_kim_ngal", "pi_kim_kim1", "pi_kim_il18", "pi_kim_tnfr1", "pi_kim_tnfr2"], axis=1, inplace=True)
labs.columns = labs.columns.str.replace(
    r"visit_|bl_", "", regex=True)
labs["procedure"] = "labs"


# ------------------------------------------------------------------------------
# BOLD/ASL MRI
# ------------------------------------------------------------------------------

var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                           == "study_visit_boldasl_mri", "field_name"]]
mri = df.filter(var, axis=1)
mri.columns = mri.columns.str.replace(
    r"mri_", "", regex=True)
mri["procedure"] = "bold_mri"

# ------------------------------------------------------------------------------
# DXA Scan
# ------------------------------------------------------------------------------
var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                           == "study_visit_dxa_scan", "field_name"]]
dxa = df.filter(var, axis=1)
dxa.columns = dxa.columns.str.replace(
    r"dxa_", "", regex=True)
dxa["procedure"] = "dxa"

# MERGE
crocodile = pd.merge(phys, screen, how="outer")
crocodile = pd.merge(crocodile, labs, how="outer")
crocodile = pd.merge(crocodile, mri, how="outer")
crocodile = pd.merge(crocodile, dxa, how="outer")
# REORGANIZE
crocodile["visit"] = "baseline"
crocodile["study"] = "CROCODILE"
id_cols = ["record_id", "study", "visit", "procedure", "date"]
other_cols = crocodile.columns.difference(id_cols, sort=False).tolist()
crocodile = crocodile[id_cols + other_cols]
# SORT
crocodile.sort_values(["record_id", "date", "procedure"], inplace=True)
# Print final data
crocodile.to_csv("~/crocodile.csv", index=False)
