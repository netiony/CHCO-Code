"""
This code is designed to pull data from the PENGUIN REDCap project and output data in a "semi-long" format with one row per study procedure, and a visit column for longitudinal clustering when combined with other studies.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_penguin():
    # Libraries
    import os
    home_dir = os.path.expanduser("~")
    os.chdir(home_dir + "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
    import redcap
    import pandas as pd
    import numpy as np
    from natsort import natsorted, ns
    from harmonization_functions import combine_checkboxes
    # REDCap project variables
    try:
        tokens = pd.read_csv(
            "/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/api_tokens.csv")
    except FileNotFoundError:
        tokens = pd.read_csv(
            "/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/api_tokens.csv")
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "PENGUIN", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)
    # Replace missing values
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999]
    rep = rep + [str(r) for r in rep]

    # --------------------------------------------------------------------------
    # Demographics
    # --------------------------------------------------------------------------

    dem_cols = ["record_id", "dob", "group", "sex", "race", "ethnicity"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    demo.replace(rep, "", inplace=True)  # Replace missing values
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
    demo["diabetes_dx_date"] = ""
    demo["co_enroll_id"] = ""
    demo["group"] = "PKD"

    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------
    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "medical_history", "field_name"]]
    med = pd.DataFrame(proj.export_records(fields=var))
    med.replace(rep, "", inplace=True)  # Replace missing values
    med_list = {"htn_med___1": "ace_inhibitor",
                "htn_med___2": "angiotensin_receptor_blocker",
                "htn_med___3": "beta_blocker",
                "htn_med___4": "ca_channel_blocker",
                "htn_med___5": "diuretic",
                "htn_med___6": "statin"}
    og_names = list(med_list.keys())
    med = med[["record_id"] + og_names]
    med.rename(med_list, axis=1, inplace=True)
    # RAASi
    med = med.assign(raasi_timepoint=np.maximum(pd.to_numeric(
        med["ace_inhibitor"]), pd.to_numeric(med["angiotensin_receptor_blocker"])))
    # Replace 0/1 values with yes/no
    med.iloc[:, 1:] = med.iloc[:, 1:].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"})
    med["procedure"] = "medications"
    med["visit"] = "screening"

    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "physical_exam", "field_name"]]
    phys = pd.DataFrame(proj.export_records(fields=var))
    phys.replace(rep, "", inplace=True)  # Replace missing values
    phys.drop(["phys_age", "phys_normal", "phys_abnormal"],
              axis=1, inplace=True)
    phys["procedure"] = "physical_exam"
    phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
    phys.rename({"sysbp": "sbp", "diasbp": "dbp"}, inplace=True, axis=1)
    phys["visit"] = "screening"

    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "screening_labs", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var))
    screen.replace(rep, "", inplace=True)  # Replace missing values
    screen.drop(["screen_creat_lab", "screen_upt", "screen_menstrual"],
                axis=1, inplace=True)
    screen.columns = screen.columns.str.replace(
        r"screen_|labs_", "", regex=True)
    screen.rename({"creat_s": "creatinine_s"}, axis=1, inplace=True)
    screen["procedure"] = "screening"
    screen["visit"] = "screening"

    # --------------------------------------------------------------------------
    # Baseline labs
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_baseline_vitalslabs", "field_name"]]
    labs = pd.DataFrame(proj.export_records(fields=var))
    labs.replace(rep, "", inplace=True)  # Replace missing values
    labs.drop(["baseline_vitals", "visit_upt", "visit_uptresult",
               "baseline_labs", "u24_labs", "pilabs_yn", "metabolomics_yn", "kim_yn"], axis=1, inplace=True)
    labs.columns = labs.columns.str.replace(
        r"visit_|bl_", "", regex=True)
    labs.rename({"a1c": "hba1c", "uacr": "acr_u"}, axis=1, inplace=True)
    labs["procedure"] = "labs"
    labs["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # DXA Scan
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_dxa_scan", "field_name"]]
    dxa = pd.DataFrame(proj.export_records(fields=var))
    dxa.replace(rep, "", inplace=True)  # Replace missing values
    dxa.rename({"bodyfat_percent": "body_fat",
                "leanmass_percent": "lean_mass", "fatmass_kg": "fat_kg",
                "leanmass_kg": "lean_kg", "trunkmass_kg": "trunk_kg",
                "trunkmass_percent": "trunk_mass",
                "bmd": "bone_mineral_density"}, axis=True, inplace=True)
    dxa.columns = dxa.columns.str.replace(
        r"dxa_", "", regex=True)
    dxa_cols = dxa.columns[2:].to_list()
    dxa.rename(dict(zip(dxa_cols, ["dexa_" + d for d in dxa_cols])),
               axis=1, inplace=True)
    dxa["procedure"] = "dxa"
    dxa["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # Clamp
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_he_clamp", "field_name"]]
    clamp = pd.DataFrame(proj.export_records(fields=var))
    clamp.replace(rep, "", inplace=True)  # Replace missing values
    # Format
    clamp.drop(["clamp_yn", "clamp_d20", "clamp_ffa",
                "clamp_insulin", "hct_yn", "clamp_bg"], axis=1, inplace=True)
    clamp.rename({"clamp_wt": "weight", "clamp_ht": "height"},
                 inplace=True, axis=1)
    clamp.columns = clamp.columns.str.replace(r"clamp_", "", regex=True)
    clamp.columns = clamp.columns.str.replace(
        r"insulin_minus", "insulin_minus_", regex=True)
    clamp.columns = clamp.columns.str.replace(
        r"ffa_minus", "ffa_minus_", regex=True)
    clamp.columns = clamp.columns.str.replace(r"bg_", "glucose_", regex=True)
    clamp.columns = clamp.columns.str.replace(
        r"glucose_minus", "glucose_minus_", regex=True)
    clamp["procedure"] = "clamp"
    clamp["visit"] = "baseline"
    clamp["insulin_sensitivity_method"] = "hyperinsulinemic_euglycemic_clamp"
    # FFA
    # See /home/timvigers/Work/CHCO/Petter Bjornstad/IHD/Background/Renal Heir Equations.docx
    ffa = [c for c in clamp.columns if "ffa_" in c]
    clamp[ffa] = clamp[ffa].apply(
        pd.to_numeric, errors='coerce')
    clamp["baseline_ffa"] = \
        clamp[['ffa_minus_20', 'ffa_minus_10', 'ffa_0']].mean(axis=1)
    clamp["p1_steady_state_ffa"] = \
        clamp[['ffa_70', 'ffa_80', 'ffa_90']].mean(axis=1)
    clamp["p2_steady_state_ffa"] = \
        clamp[['ffa_250', 'ffa_260', 'ffa_270']].mean(axis=1)
    clamp["p1_ffa_suppression"] = (
        (clamp["baseline_ffa"] - clamp["p1_steady_state_ffa"]) / clamp["baseline_ffa"]) * 100
    clamp["p2_ffa_suppression"] = (
        (clamp["baseline_ffa"] - clamp["p2_steady_state_ffa"]) / clamp["baseline_ffa"]) * 100
    clamp["ffa_method"] = "hyperinsulinemic_euglycemic_clamp"

    # --------------------------------------------------------------------------
    # Renal Clearance Testing
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_renal_clearance_testing", "field_name"]]
    rct = pd.DataFrame(proj.export_records(fields=var))
    rename = {"gfr": "gfr_raw_plasma_urine", "gfr_bsa": "gfr_bsa_plasma_urine",
              "erpf": "erpf_raw_plasma_urine", "erpf_bsa": "erpf_bsa_plasma_urine",
              "gfr_15mgmin": "gfr_raw_plasma", "gfrbsa": "gfr_bsa_plasma",
              "erpf_pah_85": "erpf_raw_plasma", "erpfbsa": "erpf_bsa_plasma"}
    rct.rename(rename, axis=1, inplace=True)
    rct = rct[["record_id"] + list(rename.values())]
    rct["procedure"] = "renal_clearance_testing"
    rct["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # PET scan
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "pet_scan", "field_name"]]
    pet = pd.DataFrame(proj.export_records(fields=var))
    pet.replace(rep, "", inplace=True)  # Replace missing values
    pet.drop(["petcom_yn"], axis=1, inplace=True)
    pet.columns = pet.columns.str.replace(r"pet_", "", regex=True)
    pet["procedure"] = "pet_scan"
    pet["visit"] = "baseline"

    # --------------------------------------------------------------------------
    # fMRI
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "fmri", "field_name"]]
    mri = pd.DataFrame(proj.export_records(fields=var))
    mri.replace(rep, "", inplace=True)  # Replace missing values
    mri["procedure"] = "fmri"
    mri["visit"] = "baseline"

    # MERGE
    df = pd.concat([phys, screen], join='outer', ignore_index=True)
    df = pd.concat([df, med], join='outer', ignore_index=True)
    df = pd.concat([df, labs], join='outer', ignore_index=True)
    df = pd.concat([df, dxa], join='outer', ignore_index=True)
    df = pd.concat([df, clamp], join='outer', ignore_index=True)
    df = pd.concat([df, rct], join='outer', ignore_index=True)
    df = pd.concat([df, pet], join='outer', ignore_index=True)
    df = pd.concat([df, mri], join='outer', ignore_index=True)
    df = pd.merge(df, demo, how="outer")
    df = df.copy()
    # REORGANIZE
    df["study"] = "PENGUIN"
    id_cols = ["record_id", "co_enroll_id", "study"] + \
        dem_cols[1:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    df = df[id_cols + other_cols]
    # SORT
    df.sort_values(["record_id", "date", "procedure"], inplace=True)
    # Drop empty columns
    df.replace("", np.nan, inplace=True)
    df.dropna(how='all', axis=1, inplace=True)
    # Return final data
    return df
