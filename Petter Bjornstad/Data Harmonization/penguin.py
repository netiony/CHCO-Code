"""
This code is designed to pull data from the PENGUIN REDCap project and output data in a "semi-long" format with one row per study procedure, and a visit column for longitudinal clustering when combined with other studies.
"""
__author__ = "Tim Vigers"
__credits__ = ["Tim Vigers"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_penguin():
    # Libraries
    import os
    import redcap
    import pandas as pd
    import numpy as np
    os.chdir(
        "/Users/timvigers/Documents/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
    from harmonization_functions import combine_checkboxes
    from harmonization_functions import find_duplicate_columns
    # REDCap project variables
    tokens = pd.read_csv(
        "~/Dropbox/Work/CHCO/Petter Bjornstad/Data Harmonization/api_tokens.csv")
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "PENGUIN", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)

# ------------------------------------------------------------------------------
# Demographics
# ------------------------------------------------------------------------------

    dem_cols = ["record_id", "dob", "group", "sex", "race", "ethnicity"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    # Race columns combined into one
    demo = combine_checkboxes(demo, base_name="race", levels=[
        "American Indian or Alaskan Native",
        "Asian",
        "Hawaiian or Pacific Islander",
        "Black or African American",
        "White",
        "Unknown",
        "Other"])
    # Same for ethnicity
    demo = combine_checkboxes(demo,
                              base_name="ethnicity",
                              levels=["Hispanic or Latino", "Not Hispanic or Latino", "Unknown/Not Reported"])
    # Relevel sex and group
    demo["sex"].replace({2: "Male", 1: "Female", 3: "Other",
                        "2": "Male", "1": "Female", "3": "Other"}, inplace=True)
    demo["diabetes_dx_date"] = np.nan
    demo["co_enroll_id"] = np.nan

# ------------------------------------------------------------------------------
# Medications
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Physical exam
# ------------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "physical_exam", "field_name"]]
    phys = pd.DataFrame(proj.export_records(fields=var))
    phys.drop(["phys_age", "phys_normal", "phys_abnormal"],
              axis=1, inplace=True)
    phys["procedure"] = "physical_exam"
    phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
    phys.rename({"sysbp": "sbp", "diasbp": "dbp"}, inplace=True, axis=1)

# ------------------------------------------------------------------------------
# Screening labs
# ------------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "screening_labs", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var))
    screen.drop(["screen_creat_lab", "screen_upt", "screen_menstrual"],
                axis=1, inplace=True)
    screen.columns = screen.columns.str.replace(
        r"screen_|labs_", "", regex=True)
    screen.rename({"creat_s": "creatinine_s"}, axis=1, inplace=True)
    screen["procedure"] = "screening"

# ------------------------------------------------------------------------------
# Baseline labs
# ------------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_baseline_vitalslabs", "field_name"]]
    labs = pd.DataFrame(proj.export_records(fields=var))
    labs.drop(["baseline_vitals", "visit_upt", "visit_uptresult",
               "baseline_labs", "u24_labs", "pilabs_yn", "metabolomics_yn", "kim_yn"], axis=1, inplace=True)
    labs.columns = labs.columns.str.replace(
        r"visit_|bl_", "", regex=True)
    labs.rename({"a1c": "hba1c", "uacr": "acr_u"}, axis=1, inplace=True)
    labs["procedure"] = "labs"

# ------------------------------------------------------------------------------
# DXA Scan
# ------------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_dxa_scan", "field_name"]]
    dxa = pd.DataFrame(proj.export_records(fields=var))
    dxa.columns = dxa.columns.str.replace(
        r"dxa_", "", regex=True)
    dxa["procedure"] = "dxa"

# ------------------------------------------------------------------------------
# Clamp
# ------------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
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

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_renal_clearance_testing", "field_name"]]
    rct = pd.DataFrame(proj.export_records(fields=var))
    rct.drop(["iohexol_yn", "pah_yn", "egfr"], axis=1, inplace=True)
    rct["procedure"] = "renal_clearance_testing"

# ------------------------------------------------------------------------------
# PET scan
# ------------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "pet_scan", "field_name"]]
    pet = pd.DataFrame(proj.export_records(fields=var))
    pet.drop(["petcom_yn"], axis=1, inplace=True)
    pet.columns = pet.columns.str.replace(r"pet_", "", regex=True)
    pet["procedure"] = "pet_scan"

# ------------------------------------------------------------------------------
# fMRI
# ------------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "fmri", "field_name"]]
    mri = pd.DataFrame(proj.export_records(fields=var))
    mri["procedure"] = "fmri"

    # MERGE
    df = pd.merge(phys, screen, how="outer")
    df = pd.merge(df, labs, how="outer")
    df = pd.merge(df, dxa, how="outer")
    df = pd.merge(df, clamp, how="outer")
    df = pd.merge(df, rct, how="outer")
    df = pd.merge(df, pet, how="outer")
    df = pd.merge(df, mri, how="outer")
    df = pd.merge(df, demo, how="outer")
    # REORGANIZE
    df["visit"] = "baseline"
    df["study"] = "PENGUIN"
    id_cols = ["record_id", "co_enroll_id", "study"] + \
        dem_cols[1:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols).tolist()
    df = df[id_cols + other_cols]
    # SORT
    df.sort_values(["record_id", "date", "procedure"], inplace=True)
    # Check for duplicated column names
    # dups = find_duplicate_columns(df)
    # dups.to_csv("~/df_duplicate_columns.csv", index=False)
    # Return final data
    return df
