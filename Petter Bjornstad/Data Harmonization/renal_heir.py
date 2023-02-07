"""
This code is designed to pull data from the RENAL HEIR REDCap project and output data in a "semi-long" format with one row per study procedure, and a visit column for longitudinal clustering when combined with other studies.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_renal_heir():
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
    token = tokens.loc[tokens["Study"] == "Renal-HEIR", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)

    # --------------------------------------------------------------------------
    # Demographics
    # --------------------------------------------------------------------------

    dem_cols = ["subject_id", "co_enroll_id", "dob", "diagnosis",
                "group", "gender", "race", "ethnicity"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols))
    demo.rename({"gender": "sex", "diagnosis": "diabetes_dx_date"},
                inplace=True, axis=1)
    dem_cols[3] = "diabetes_dx_date"
    dem_cols[5] = "sex"
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

    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------

    var = ["subject_id"] + ["sglt2i"] + [v for v in meta.loc[meta["form_name"]
                                                == "medical_history", "field_name"]]
    med = pd.DataFrame(proj.export_records(fields=var))
    # Just SGLT2i for now
    med = med[["subject_id", "sglt2i", "diabetes_med_other___3", "diabetes_med___3"]]
    med.replace({0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"diabetes_med_other___3": "sglti_timepoint",
                "sglt2i": "sglt2i_ever"}, axis=1, inplace=True)
    # Insulin
    med["diabetes_med___3"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"diabetes_med___3": "insulin_med_timepoint"},
               axis=1, inplace=True)

    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "physical_exam", "field_name"]]
    phys = pd.DataFrame(proj.export_records(fields=var))
    phys["procedure"] = "physical_exam"
    phys.drop(["male_activity_factor", "fem_activity_factor", "schofield_male",
               "schofield_female", "phys_norm", "phys_no", "breast_tanner", "testicular_volume", "lmp", "screen_bmi_percentile"], axis=1, inplace=True)
    phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
    phys.rename({"sys_bp": "sbp", "dys_bp": "dbp",
                "waist_circumference": "waistcm", "hip_circumference": "hipcm"}, inplace=True, axis=1)

    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "screening_labs", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var))
    screen.drop(["a1c_pre", "a1c_pre_date", "screen_pregnant"],
                axis=1, inplace=True)
    screen.columns = screen.columns.str.replace(
        r"_of_screen|screen_", "", regex=True)
    screen.rename({"serum_creatinine": "creatinine_s", "urine_acr": "acr_u",
                   "urine_mab": "microalbumin_u", "urine_cre": "creatinine_u"},
                  axis=1, inplace=True)
    screen["procedure"] = "screening"

    # --------------------------------------------------------------------------
    # DXA Scan
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "body_composition_dxa", "field_name"]]
    dxa = pd.DataFrame(proj.export_records(fields=var))
    dxa.columns = dxa.columns.str.replace(
        r"dexa_", "", regex=True)
    dxa_cols = dxa.columns[2:].to_list()
    dxa.rename(dict(zip(dxa_cols, ["dexa_" + d for d in dxa_cols])),
               axis=1, inplace=True)
    dxa["procedure"] = "dxa"

    # --------------------------------------------------------------------------
    # Clamp
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "clamp", "field_name"]]
    clamp = pd.DataFrame(proj.export_records(fields=var))
    # Format
    clamp.drop(["baseline", "fasting_labs", "urine_labs", "hct_lab",
                "bg_labs", "ffa_lab", "cpep_lab", "insulin_labs"], axis=1, inplace=True)
    clamp.columns = clamp.columns.str.replace(
        r"clamp_", "", regex=True)
    clamp.rename({"serum_creatinine": "creatinine_s",
                  "serum_sodium": "sodium_s",
                  "cystatin_c": "cystatin_c_s",
                  "urine_mab_baseline": "microalbumin_u",
                  "urine_cre_baseline": "creatinine_u"},
                 inplace=True, axis=1)
    clamp.columns = clamp.columns.str.replace(r"clamp_", "", regex=True)
    clamp["procedure"] = "clamp"
    clamp["insulin_sensitivity_method"] = "hyperglycemic_clamp"
    # M
    num_vars = ["d20_infusion", "weight"]
    clamp[num_vars] = clamp[num_vars].apply(
        pd.to_numeric, errors='coerce')
    clamp["raw_m"] = (clamp["d20_infusion"] * 190 / 60) / clamp["weight"]
    # FFA
    # See /home/timvigers/Work/CHCO/Petter Bjornstad/IHD/Background/Renal Heir Equations.docx
    ffa = [c for c in clamp.columns if "ffa_" in c]
    clamp[ffa] = clamp[ffa].apply(
        pd.to_numeric, errors='coerce')
    clamp["baseline_ffa"] = \
        clamp[['ffa_minus_10', 'ffa_minus_5']].mean(axis=1)
    clamp["steady_state_ffa"] = \
        clamp[['ffa_220', 'ffa_230', 'ffa_240', 'ffa_250']].mean(axis=1)
    clamp["ffa_suppression"] = (
        (clamp["baseline_ffa"] - clamp["steady_state_ffa"]) / clamp["baseline_ffa"]) * 100
    clamp["ffa_method"] = "hyperglycemic_clamp"
    # Insulin
    ins = ['insulin_minus_10', 'insulin_minus_5', 'insulin_2', 'insulin_4',
           'insulin_6', 'insulin_8', 'insulin_10', 'insulin_120',
           'insulin_220', 'insulin_230', 'insulin_240', 'insulin_250']
    clamp[ins] = clamp[ins].apply(
        pd.to_numeric, errors='coerce')
    clamp["steady_state_insulin"] = clamp[['insulin_220',
                                           'insulin_230', 'insulin_240', 'insulin_250']].mean(axis=1) * 6
    # C peptide
    cpep = [c for c in clamp.columns if "cpeptide" in c]
    clamp[cpep] = clamp[cpep].apply(
        pd.to_numeric, errors='coerce')
    clamp["steady_state_cpeptide"] = clamp[['cpeptide_220',
                                            'cpeptide_230', 'cpeptide_240', 'cpeptide_250']].mean(axis=1)
    # ACPRg
    clamp["acprg"] = clamp[['cpeptide_2', 'cpeptide_4',
                            'cpeptide_6', 'cpeptide_8', 'cpeptide_10']].mean(axis=1) - clamp[['cpeptide_minus_10', 'cpeptide_minus_5']].mean(axis=1)
    # AIRg
    clamp["airg"] = clamp[['insulin_2', 'insulin_4',
                           'insulin_6', 'insulin_8', 'insulin_10']].mean(axis=1) * 6 - clamp[['insulin_minus_10', 'insulin_minus_5']].mean(axis=1) * 6
    # DI
    clamp["di"] = \
        (clamp["raw_m"] / clamp["steady_state_insulin"]) * clamp["airg"]

    # --------------------------------------------------------------------------
    # Outcomes
    # --------------------------------------------------------------------------

    var = ["subject_id"] + [v for v in meta.loc[meta["form_name"]
                                                == "outcomes", "field_name"]]
    out = pd.DataFrame(proj.export_records(fields=var))
    out.drop(["kidney_outcomes", "egfr", "metab_outcomes", "asl_outcomes"],
             axis=1, inplace=True)
    out.columns = out.columns.str.replace(
        r"mri_", "", regex=True)
    rename = {"gfr": "gfr_raw_plasma", "gfr_bsa": "gfr_bsa_plasma",
              "rpf": "erpf_raw_plasma", "rpf_bsa": "erpf_bsa_plasma"}
    out.rename(rename, axis=1, inplace=True)
    out["procedure"] = "kidney_outcomes"

    # --------------------------------------------------------------------------
    # Kidney Biopsy
    # --------------------------------------------------------------------------

    var = ["subject_id", ] + [v for v in meta.loc[meta["form_name"]
                                                  == "kidney_biopsy", "field_name"]]
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

    # MERGE
    df = pd.concat([phys, screen], join='outer', ignore_index=True)
    df = pd.concat([df, med], join='outer', ignore_index=True)
    df = pd.concat([df, dxa], join='outer', ignore_index=True)
    df = pd.concat([df, clamp], join='outer', ignore_index=True)
    df = pd.concat([df, out], join='outer', ignore_index=True)
    df = pd.concat([df, biopsy], join='outer', ignore_index=True)
    df = pd.merge(df, demo, how="outer")
    df = df.copy()
    # REORGANIZE
    df["visit"] = "baseline"
    df["study"] = "RENAL-HEIR"
    id_cols = ["subject_id", "co_enroll_id", "study"] + \
        dem_cols[2:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    df = df[id_cols + other_cols]
    # SORT
    df.sort_values(["subject_id", "date", "procedure"], inplace=True)
    # Rename subject identifier
    df.rename({"subject_id": "record_id"}, axis=1, inplace=True)
    # Replace missing values
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999]
    rep = rep + [str(r) for r in rep]
    df.replace(rep, "", inplace=True)
    # Return final data
    return df
