"""
This code is designed to pull data from the CROCODILE REDCap project and output data in a "semi-long" format with one row per study procedure, and a visit column for longitudinal clustering when combined with other studies.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_crocodile():
    # Libraries
    import os
    home_dir = os.path.expanduser("~")
    os.chdir(home_dir + "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
    import sys
    import redcap
    import pandas as pd
    from natsort import natsorted, ns
    from harmonization_functions import combine_checkboxes
    # REDCap project variables
    try:
      tokens = pd.read_csv("/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/api_tokens.csv")
    except FileNotFoundError:
      tokens = pd.read_csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/api_tokens.csv")
    uri = "https://redcap.ucdenver.edu/api/"
    token = tokens.loc[tokens["Study"] == "CROCODILE", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)

    # --------------------------------------------------------------------------
    # Demographics
    # --------------------------------------------------------------------------

    dem_cols = ["record_id", "dob", "diabetes_dx_date",
                "group", "sex", "race", "ethnicity"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    demo["co_enroll_id"] = ""
    # Race columns combined into one
    demo = combine_checkboxes(demo, base_name="race", levels=[
        'American Indian/Alaska Native', 'Asian', 'Hawaiian/Pacific Islander', 'Black/African American', 'White', 'Other', 'Unknown'])
    # Same for ethnicity
    demo = combine_checkboxes(demo,
                              base_name="ethnicity",
                              levels=["Hispanic or Latino", "Not Hispanic or Latino", "Unknown/Not Reported"])
    # Relevel sex and group
    demo["sex"].replace({1: "Male", 2: "Female", 3: "Other",
                        "1": "Male", "2": "Female", "3": "Other"}, inplace=True)
    demo["group"].replace({1: "Type 1 Diabetes", 2: "Lean Control",
                           "1": "Type 1 Diabetes", "2": "Lean Control"}, inplace=True)

    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "medical_history", "field_name"]]
    med = pd.DataFrame(proj.export_records(fields=var))
    med_list = {"diabetes_meds_other___1": "metformin_timepoint",
                "diabetes_meds_other___2": "tzd_timepoint",
                "diabetes_meds_other___3": "glp1_agonist_timepoint",
                "diabetes_meds_other___4": "sglti_timepoint",
                "diabetes_meds_other___5": "other_diabetes_med_timepoint",
                "htn_med___1": "ace_inhibitor",
                "htn_med___2": "angiotensin_receptor_blocker",
                "htn_med___3": "beta_blocker",
                "htn_med___4": "ca_channel_blocker",
                "htn_med___5": "diuretic",
                "htn_med___6": "statin"}
    og_names = list(med_list.keys())
    med = med[["record_id"] + og_names]
    med[og_names] = med[og_names].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"})
    med.rename(med_list, axis=1, inplace=True)

    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "physical_exam", "field_name"]]
    phys = pd.DataFrame(proj.export_records(fields=var))
    phys["procedure"] = "physical_exam"
    phys.drop(["phys_normal", "phys_abnormal"], axis=1, inplace=True)
    phys.columns = phys.columns.str.replace(r"phys_", "", regex=True)
    phys.rename({"sysbp": "sbp", "diasbp": "dbp"}, inplace=True, axis=1)

    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "screening_labs", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var))
    screen.drop(["prescreen_a1c", "prescreen_a1c_date",
                "screen_menstrual", "screen_upt"], axis=1, inplace=True)
    screen.columns = screen.columns.str.replace(
        r"labs_|screen_", "", regex=True)
    screen.rename({"creat_s": "creatinine_s", "uacr": "acr_u", "a1c": "hba1c",
                   "creat_u": "creatinine_u", "hg": "hemoglobin"}, axis=1, inplace=True)
    screen["procedure"] = "screening"

    # --------------------------------------------------------------------------
    # Labs
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_baseline_vitalslabs", "field_name"]]
    labs = pd.DataFrame(proj.export_records(fields=var))
    labs.drop(["baseline_vitals", "visit_upt",
               "visit_uptresult", "baseline_labs", "pilabs_yn", "pi_copeptin", "pi_renin", "pi_angiotensin2", "pi_osmo_s", "pi_osmo_u", "pi_lithium_s", "pi_lithium_u", "metabolomics_yn", "kim_yn", "pi_kim_ykl40", "pi_kim_ngal", "pi_kim_kim1", "pi_kim_il18", "pi_kim_tnfr1", "pi_kim_tnfr2"], axis=1, inplace=True)
    labs.columns = labs.columns.str.replace(
        r"visit_|bl_", "", regex=True)
    labs.rename({"uacr": "acr_u", "a1c": "hba1c"}, axis=1, inplace=True)
    labs["procedure"] = "labs"

    # --------------------------------------------------------------------------
    # BOLD/ASL MRI
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_boldasl_mri", "field_name"]]
    mri = pd.DataFrame(proj.export_records(fields=var))
    mri.columns = mri.columns.str.replace(
        r"mri_", "", regex=True)
    mri.rename({"volume_right": "right_kidney_volume_ml", "volume_left": "left_kidney_volume_ml"},
               axis=1, inplace=True)
    mri["procedure"] = "bold_mri"

    # --------------------------------------------------------------------------
    # DXA Scan
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_dxa_scan", "field_name"]]
    dxa = pd.DataFrame(proj.export_records(fields=var))
    dxa.columns = dxa.columns.str.replace(
        r"dxa_|_percent", "", regex=True)
    dxa.rename({"bodyfat": "body_fat", "leanmass": "lean_mass",
                "trunkmass": "trunk_mass", "fatmass_kg": "fat_kg",
                "leanmass_kg": "lean_kg", "trunkmass_kg": "trunk_kg",
                "bmd": "bone_mineral_density"}, axis=1, inplace=True)
    dxa_cols = dxa.columns[2:].to_list()
    dxa.rename(dict(zip(dxa_cols, ["dexa_" + d for d in dxa_cols])),
               axis=1, inplace=True)
    dxa["procedure"] = "dxa"

    # --------------------------------------------------------------------------
    # Clamp
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "study_visit_he_clamp", "field_name"]]
    clamp = pd.DataFrame(proj.export_records(fields=var))
    # Format
    clamp.drop(["clamp_yn", "clamp_d20", "clamp_ffa",
                "clamp_insulin", "hct_yn", "clamp_bg"], axis=1, inplace=True)
    clamp.rename({"clamp_wt": "weight", "clamp_ht": "height",
                  "cystatin_c": "cystatin_c_s"},
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
    clamp["insulin_sensitivity_method"] = "hyperinsulinemic_euglycemic_clamp"
    # FFA
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
    rename = {"gfr_raw": "gfr_raw_plasma_urine", "gfr_bsa": "gfr_bsa_plasma_urine",
              "erpf_raw": "erpf_raw_plasma_urine", "erpf": "erpf_bsa_plasma_urine",
              "gfr_15mgmin": "gfr_raw_plasma", "gfrbsa": "gfr_bsa_plasma",
              "erpf_pah_85": "erpf_raw_plasma", "erpfbsa": "erpf_bsa_plasma"}
    rct.rename(rename, axis=1, inplace=True)
    rct = rct[["record_id"] + list(rename.values())]
    rct["procedure"] = "renal_clearance_testing"

    # --------------------------------------------------------------------------
    # Kidney Biopsy
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
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

    # --------------------------------------------------------------------------
    # PET scan
    # --------------------------------------------------------------------------

    var = ["record_id"] + [v for v in meta.loc[meta["form_name"]
                                               == "optional_pet_scan", "field_name"]]
    pet = pd.DataFrame(proj.export_records(fields=var))
    pet.drop(["petcon_yn"], axis=1, inplace=True)
    pet.columns = pet.columns.str.replace(r"pet_", "", regex=True)
    pet["procedure"] = "pet_scan"

    # MERGE
    df = pd.concat([phys, screen], join='outer', ignore_index=True)
    df = pd.concat([df, med], join='outer', ignore_index=True)
    df = pd.concat([df, labs], join='outer', ignore_index=True)
    df = pd.concat([df, mri], join='outer', ignore_index=True)
    df = pd.concat([df, dxa], join='outer', ignore_index=True)
    df = pd.concat([df, clamp], join='outer', ignore_index=True)
    df = pd.concat([df, rct], join='outer', ignore_index=True)
    df = pd.concat([df, biopsy], join='outer', ignore_index=True)
    df = pd.concat([df, pet], join='outer', ignore_index=True)
    df = pd.merge(df, demo, how="outer")
    df = df.copy()
    # REORGANIZE
    df["visit"] = "baseline"
    df["study"] = "CROCODILE"
    id_cols = ["record_id", "co_enroll_id", "study"] + \
        dem_cols[1:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    df = df[id_cols + other_cols]
    # SORT
    df.sort_values(["record_id", "date", "procedure"], inplace=True)
    # Rename IDs
    df["record_id"] = ["CRC-" + str(i).zfill(2) for i in df["record_id"]]
    # Replace missing values
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999]
    rep = rep + [str(r) for r in rep]
    df.replace(rep, "", inplace=True)
    # Print final data
    return df

