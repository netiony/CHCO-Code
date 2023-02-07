"""
This code is designed to pull data from the IMPROVE REDCap project and output data in a long format with one row per study procedure per visit.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def clean_improve():
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
    token = tokens.loc[tokens["Study"] == "IMPROVE", "Token"].iloc[0]
    proj = redcap.Project(url=uri, token=token)
    # Get project metadata
    meta = pd.DataFrame(proj.metadata)
    # Columns to drop
    redcap_cols = ["redcap_event_name",
                   "redcap_repeat_instrument", "redcap_repeat_instance"]

    # --------------------------------------------------------------------------
    # Demographics
    # --------------------------------------------------------------------------

    dem_cols = ["subject_id", "co_enroll_id", "dob", "diagnosis",
                "gender", "race", "ethnicity", "sglt2i"]
    # Export
    demo = pd.DataFrame(proj.export_records(fields=dem_cols,
                                            events=["screening_arm_1"]))
    demo["group"] = "Type 2 Diabetes"
    demo.drop(redcap_cols, axis=1, inplace=True)
    demo.rename({"gender": "sex", "diagnosis": "diabetes_dx_date"},
                inplace=True, axis=1)
    dem_cols[3] = "diabetes_dx_date"
    dem_cols[4] = "sex"
    dem_cols[7] = "sglt2i_ever"
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
    demo["sglt2i"].replace({1: "Yes", 0: "No", "1": "Yes", "0": "No"},
                           inplace=True)
    demo.rename({"sglt2i": "sglt2i_ever"}, axis=1, inplace=True)

    # --------------------------------------------------------------------------
    # Medications
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit", "diabetes_med_other"]
    med = pd.DataFrame(proj.export_records(fields=var))
    # SGLT2i (diabetes_med_other___4), RAASi (htn_med_type___1, htn_med_type___2), Metformin (diabetes_med_other___1)
    med = med[["subject_id", "diabetes_med_other___4", "htn_med_type___1", "htn_med_type___2", "diabetes_med_other___1", "diabetes_med___1", "diabetes_med___2"]]
    # SGLT2i
    med["diabetes_med_other___4"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"diabetes_med_other___4": "sglti_timepoint"},
               axis=1, inplace=True)
    # RAASi
    med = med.assign(raasi = np.maximum(pd.to_numeric(med["htn_med_type___1"]), pd.to_numeric(med["htn_med_type___2"])))
    med.drop(med[['htn_med_type___1', 'htn_med_type___2']], axis=1, inplace=True)
    med["raasi_timepoint"].replace(
    {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    # Metformin
    med["diabetes_med___1"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"diabetes_med___1": "metformin_timepoint"},
               axis=1, inplace=True)
    # Insulin
    med["diabetes_med_2"].replace(
        {0: "No", "0": "No", 1: "Yes", "1": "Yes"}, inplace=True)
    med.rename({"diabetes_med_2": "insulin_med_timepoint"},
               axis=1, inplace=True)


    # --------------------------------------------------------------------------
    # Physical exam
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "physical_exam", "field_name"]]
    phys = pd.DataFrame(proj.export_records(
        fields=var, events=["screening_arm_1"]))
    phys["procedure"] = "physical_exam"
    phys.drop(redcap_cols + ["phys_norm", "phys_no", "breast_tanner",
                             "testicular_volume", "lmp", "screen_bmi_percentile", "activity_factor_male", "activity_factor_female", "schofield_male", "schofield_female"], axis=1, inplace=True)
    phys.columns = phys.columns.str.replace(r"phys_|screen_", "", regex=True)
    phys.rename({"sys_bp": "sbp", "dys_bp": "dbp",
                 "waist_circumference": "waistcm",
                 "hip_circumference": "hipcm"}, inplace=True, axis=1)

    # --------------------------------------------------------------------------
    # Screening labs
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "screening_labs", "field_name"]]
    screen = pd.DataFrame(proj.export_records(fields=var,
                                              events=["screening_arm_1"]))
    screen.drop(redcap_cols + ['a1c_pre', 'a1c_pre_date', "screen_pregnant"],
                axis=1, inplace=True)
    screen.columns = screen.columns.str.replace(
        r"screen_|_of_screen", "", regex=True)
    screen.rename({"serum_creatinine": "creatinine_s", "urine_acr": "acr_u",
                   "urine_cre": "creatinine_u", "urine_mab": "microalbumin_u"},
                  axis=1, inplace=True)
    screen["procedure"] = "screening"

    # --------------------------------------------------------------------------
    # Accelerometry
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "accelerometry", "field_name"]]
    accel = pd.DataFrame(proj.export_records(fields=var))
    accel = accel.loc[accel["acc_wear_percent"] != ""]
    accel.drop(redcap_cols + ["study_visit_accel"], axis=1, inplace=True)
    accel.columns = accel.columns.str.replace(
        r"acc_|accel_", "", regex=True)
    accel["procedure"] = "accelerometry"

    # --------------------------------------------------------------------------
    # Cardio/Abdominal MRI
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "cardioabdominal_mri", "field_name"]]
    mri = pd.DataFrame(proj.export_records(fields=var))
    mri.drop(redcap_cols + ["mri_cardio", "mri_abdo",
                            "mri_aortic", "study_visit_mri"],
             axis=1, inplace=True)
    mri.columns = mri.columns.str.replace(
        r"mri_|visit_", "", regex=True)
    mri["procedure"] = "cardio_abdominal_mri"

    # --------------------------------------------------------------------------
    # MMTT + Metabolic Cart
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "mmtt_metabolic_cart", "field_name"]]
    mmtt = pd.DataFrame(proj.export_records(fields=var))
    # Drop unnecessary columns
    mmtt.drop(redcap_cols + ["study_visit_mttt", "mmtt_vitals", "mmtt_pregnant",
                             "mmtt_lmp", "mmtt_brmr", "mmtt_60rmr", "mmtt_base_labs", "mmtt_ffa_labs", "mmtt_insulin",
                             "mmtt_glp1", "mmtt_cpep", "mmtt_yy", "mmtt_glucagon",
                             "mmtt_gluc"],
              axis=1, inplace=True)
    mmtt.columns = mmtt.columns.str.replace(
        r"mmtt_", "", regex=True)
    mmtt.columns = mmtt.columns.str.replace(
        r"_neg_", "_minus_", regex=True)
    mmtt.rename({"wt": "weight", "ht": "height", "waist": "waistcm",
                "hip": "hipcm", "hr": "pulse", "sys_bp": "sbp",
                 "dia_bp": "dbp", "hba1c_base": "hba1c"},
                inplace=True, axis=1)
    mmtt["procedure"] = "mmtt"
    # FFA
    # See /home/timvigers/Work/CHCO/Petter Bjornstad/IHD/Background/Renal Heir Equations.docx
    ffa = [c for c in mmtt.columns if "ffa_" in c]
    mmtt[ffa] = mmtt[ffa].apply(
        pd.to_numeric, errors='coerce')
    mmtt["baseline_ffa"] = mmtt[['ffa_minus_10', 'ffa_0']].mean(axis=1)
    mmtt["steady_state_ffa"] = mmtt['ffa_120']
    mmtt["ffa_suppression"] = (
        (mmtt["baseline_ffa"] - mmtt["steady_state_ffa"]) / mmtt["baseline_ffa"]) * 100
    mmtt["ffa_method"] = "mmtt"

    # --------------------------------------------------------------------------
    # DXA
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "body_composition_dxa_bod_pod", "field_name"]]
    dxa = pd.DataFrame(proj.export_records(fields=var))
    dxa.drop(redcap_cols + ["study_visit_bodycomp", "dxa_complete",
                            "bodpod_complete"],
             axis=1, inplace=True)
    dxa.columns = dxa.columns.str.replace(
        r"dxa_", "dexa_", regex=True)
    dxa.columns = dxa.columns.str.replace(
        r"bp_", "bod_pod_", regex=True)
    dxa.rename({"dexa_bmd": "dexa_bone_mineral_density",
                "bod_pod_fat_mass": "bod_pod_fat_kg",
                "bodcomp_date": "date"}, axis=1, inplace=True)
    dxa["procedure"] = "dxa"

    # --------------------------------------------------------------------------
    # Clamp
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "clamp", "field_name"]]
    clamp = pd.DataFrame(proj.export_records(fields=var))
    # Format
    clamp.drop(redcap_cols + ["study_visit_clamp", "baseline", "fasting_labs",
                              "bg_labs", "ns_bolus", "urine_labs"],
               axis=1, inplace=True)
    clamp.columns = clamp.columns.str.replace(
        r"clamp_", "", regex=True)
    clamp.rename({"cystatin_c": "cystatin_c_s", "urine_mab": "microalbumin_u",
                  "serum_creatinine": "creatinine_s", "acr_baseline": "acr_u",
                  "urine_mab_baseline": "microalbumin_u",
                  "urine_cre_baseline": "creatinine_u"
                  }, inplace=True, axis=1)
    clamp["procedure"] = "clamp"
    clamp["insulin_sensitivity_method"] = "hyperglycemic_clamp"
    # M
    num_vars = ["d20_infusion", "weight"]
    clamp[num_vars] = clamp[num_vars].apply(
        pd.to_numeric, errors='coerce')
    clamp["raw_m"] = (clamp["d20_infusion"] * 190 / 60) / clamp["weight"]
    # No FFA
    # Insulin
    ins = ['insulin_minus_10', 'insulin_minus_5', 'insulin_2', 'insulin_4',
           'insulin_6', 'insulin_8', 'insulin_10', 'insulin_120', 'insulin_220',
           'insulin_230', 'insulin_240', 'insulin_245', 'insulin_249',
           'insulin_250', 'insulin_252', 'insulin_253', 'insulin_254',
           'insulin_255']
    clamp[ins] = clamp[ins].apply(
        pd.to_numeric, errors='coerce')
    clamp["steady_state_insulin"] = clamp[['insulin_220',
                                           'insulin_230', 'insulin_240', 'insulin_245', 'insulin_249',
                                           'insulin_250', 'insulin_252', 'insulin_253', 'insulin_254',
                                           'insulin_255']].mean(axis=1) * 6
    # C peptide
    cpep = [c for c in clamp.columns if "cpeptide" in c]
    clamp[cpep] = clamp[cpep].apply(
        pd.to_numeric, errors='coerce')
    clamp["steady_state_cpeptide"] = clamp[['cpeptide_220', 'cpeptide_230',
                                            'cpeptide_240', 'cpeptide_245',
                                            'cpeptide_249', 'cpeptide_250', 'cpeptide_252', 'cpeptide_253', 'cpeptide_254', 'cpeptide_255']].mean(axis=1)
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

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "outcomes", "field_name"]]
    out = pd.DataFrame(proj.export_records(fields=var))
    out.drop(redcap_cols + ["kidney_outcomes", "egfr", "metab_outcomes",
                            "asl_outcomes", "bold_outcomes", "adc_outcomes"],
             axis=1, inplace=True)
    out.columns = out.columns.str.replace(
        r"mri_", "", regex=True)
    rename = {"gfr": "gfr_raw_plasma", "gfr_bsa": "gfr_bsa_plasma",
              "rpf": "erpf_raw_plasma", "erpf_bsa": "erpf_bsa_plasma"}
    out.rename(rename, axis=1, inplace=True)
    out["procedure"] = "kidney_outcomes"

    # --------------------------------------------------------------------------
    # Kidney Biopsy
    # --------------------------------------------------------------------------

    var = ["subject_id", "study_visit"] + [v for v in meta.loc[meta["form_name"]
                                                               == "kidney_biopsy", "field_name"]]
    var = var + ["gloms", "gloms_gs", "ifta", "vessels_other", "fia",
                 "glom_tuft_area", "glom_volume_weibel", "glom_volume_wiggins",
                 "glom_volume_con", "mes_matrix_area",
                 "mes_index", "mes_volume_weibel", "mes_volume_wiggins",
                 "mes_volume_con", "glom_nuc_count", "mes_nuc_count", "art_intima",
                 "art_media", "pod_nuc_density", "pod_cell_volume"]
    biopsy = pd.DataFrame(proj.export_records(fields=var))
    biopsy.drop(redcap_cols + [col for col in biopsy.columns if '_yn' in col] +
                [col for col in biopsy.columns if 'procedure_' in col],
                axis=1, inplace=True)
    biopsy.columns = biopsy.columns.str.replace(r"bx_", "", regex=True)
    biopsy.columns = biopsy.columns.str.replace(r"labs_", "", regex=True)
    biopsy.columns = biopsy.columns.str.replace(r"vitals_", "", regex=True)
    biopsy["procedure"] = "kidney_biopsy"

    # MERGE
    df = pd.concat([phys, screen], join='outer', ignore_index=True)
    df = pd.concat([df, med], join='outer', ignore_index=True)
    df = pd.concat([df, accel], join='outer', ignore_index=True)
    df = pd.concat([df, mri], join='outer', ignore_index=True)
    df = pd.concat([df, mmtt], join='outer', ignore_index=True)
    df = pd.concat([df, dxa], join='outer', ignore_index=True)
    df = pd.concat([df, clamp], join='outer', ignore_index=True)
    df = pd.concat([df, out], join='outer', ignore_index=True)
    df = pd.concat([df, biopsy], join='outer', ignore_index=True)
    df = pd.merge(df, demo, how="outer")
    df = df.copy()
    # REORGANIZE
    df.rename({"study_visit": "visit"}, axis=1, inplace=True)
    df["study"] = "IMPROVE"
    id_cols = ["subject_id", "study"] + \
        dem_cols[1:] + ["visit", "procedure", "date"]
    other_cols = df.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    df = df[id_cols + other_cols]
    # SORT
    df.sort_values(
        ["subject_id", "visit", "date", "procedure"], inplace=True)
    # Change study visit names
    df["visit"].replace({'': "baseline", '1': "pre_surgery",
                         '2': "3_months_post_surgery",
                              '3': "12_months_post_surgery"}, inplace=True)
    # Rename subject identifier
    df.rename({"subject_id": "record_id"}, axis=1, inplace=True)
    # Replace missing values
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999, -99999]
    rep = rep + [str(r) for r in rep]
    df.replace(rep, "", inplace=True)
    # Return final data
    return df
