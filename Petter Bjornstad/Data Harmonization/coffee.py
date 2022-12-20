"""
This code is designed to pull data from the COFFEE REDCap project and output data in a "semi-long" format with one row per study procedure, and a visit column for longitudinal clustering when combined with other studies.
"""
__author__ = ["Tim Vigers","Ye Ji Choi"]
__credits__ = ["Tim Vigers","Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_coffee():
    # Libraries
    import redcap
    import pandas as pd
    from natsort import natsorted, ns
    from harmonization_functions import combine_checkboxes
    # REDCap project variables
    tokens = pd.read_csv(
        "/home/timvigers/Work/CHCO/Petter Bjornstad/Data Harmonization/api_tokens.csv")
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "COFFEE", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)

    # --------------------------------------------------------------------------
    # Demographics
    # --------------------------------------------------------------------------

    dem_cols = ["subject_id", "dob", "diagnosis",
                "gender", "race", "ethnicity"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    demo["co_enroll_id"] = ""
    demo.rename({"gender": "sex", "diagnosis": "diabetes_dx_date"},
                inplace=True, axis=1)
    dem_cols[2] = "diabetes_dx_date"
    dem_cols[3] = "sex"
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
    demo["sex"].replace({1: "Male", 0: "Female", 3: "Other",
                        "1": "Male", "0": "Female", "3": "Other"}, inplace=True)
    demo["group"] = "Type 1 Diabetes"

    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "medical_history_casper", "field_name"]]
    med = pd.DataFrame(proj.export_records(fields=var))
    # Just SGLT2i for now
    med = med[["subject_id", "diabetes_med_other___4"]]
    med["diabetes_med_other___4"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"diabetes_med_other___4": "sglti_timepoint"},
               axis=1, inplace=True)

    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "physical_exam_casper", "field_name"]]
    phys = pd.DataFrame(proj.export_records(fields=var))
    phys["procedure"] = "physical_exam"
    phys.drop(["phys_norm", "phys_no", "breast_tanner",
               "testicular_volume", "lmp", "screen_bmi_percentile", "male_activity_factor", "fem_activity_factor", "schofield_male", "schofield_female"], axis=1, inplace=True)
    phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
    phys.rename({"sys_bp": "sbp", "dys_bp": "dbp", "waist_circumference": "waistcm",
                "hip_circumference": "hipcm"}, inplace=True, axis=1)

    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "screening_labs_casper", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var))
    screen.drop(['a1c_pre', 'a1c_pre_date', "screen_pregnant"],
                axis=1, inplace=True)
    screen.columns = screen.columns.str.replace(
        r"screen_|_of_screen", "", regex=True)
    screen.rename({"serum_creatinine": "creatinine_s", "urine_acr": "acr_u",
                   "urine_cre": "creatinine_u", "urine_mab": "microalbumin_u"},
                  axis=1, inplace=True)
    screen["procedure"] = "screening"

    # --------------------------------------------------------------------------
    # Clamp
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "clamp", "field_name"]]
    clamp = pd.DataFrame(proj.export_records(fields=var))
    clamp.drop(["baseline", "fasting_labs", "bg_labs", "urine_labs", "hct_lab",
                "cs_clamp_date"],
               axis=1, inplace=True)
    clamp.columns = clamp.columns.str.replace(
        r"clamp_|cf_", "", regex=True)
    clamp.rename({"cystatin_c": "cystatin_c_s",
                  "serum_creatinine": "creatinine_s"
                  }, inplace=True, axis=1)
    clamp["procedure"] = "clamp"

    # --------------------------------------------------------------------------
    # Outcomes
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "outcomes", "field_name"]]
    out = pd.DataFrame(proj.export_records(fields=var))
    out.drop(["kidney_outcomes", "egfr", "metab_outcomes",
              "asl_outcomes", "bold_outcomes"],
             axis=1, inplace=True)
    out.columns = out.columns.str.replace(
        r"mri_", "", regex=True)
    out["procedure"] = "kidney_outcomes"
    # MERGE
    df = pd.merge(phys, screen, how="outer")
    df = pd.merge(df, med, how="outer")
    df = pd.merge(df, clamp, how="outer")
    df = pd.merge(df, out, how="outer")
    df = pd.merge(df, demo, how="outer")
    # REORGANIZE
    df["visit"] = "baseline"
    df["study"] = "COFFEE"
    id_cols = ["subject_id", "co_enroll_id", "study"] + \
        dem_cols[1:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    df = df[id_cols + other_cols]
    # SORT
    df.sort_values(["subject_id", "date", "procedure"], inplace=True)
    # Rename subject identifier
    df.rename({"subject_id": "record_id"}, axis=1, inplace=True)
    # Return final data
    return df
