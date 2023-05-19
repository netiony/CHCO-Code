"""
This code is designed to pull data from multiple REDCap projects and harmonize
the data in a single dataset. Some studies are cross-sectional but include
measures at multiple visits, and some studies are longitudinal. So, this code
outputs data in a "semi-long" format with one row per study procedure, and a
visit column for longitudinal clustering. The data cleaning process for each 
individual dataset is a separate function, and this function puts them together 
and performs some formatting tweaks.
"""
__author__ = ["Tim Vigers", "Ye Ji Choi"]
__credits__ = ["Tim Vigers", "Ye Ji Choi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def harmonize_data():
    # Libraries
    import os
    import sys
    sys.path.insert(0, os.path.expanduser('~') +
                    "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
    import pandas as pd
    import numpy as np
    from natsort import natsorted, ns
    from casper import clean_casper
    from coffee import clean_coffee
    from crocodile import clean_crocodile
    from improve import clean_improve
    from penguin import clean_penguin
    from renal_heir import clean_renal_heir
    from renal_heiritage import clean_renal_heiritage
    from panther import clean_panther
    from harmonization_functions import calc_egfr
    # Use individual data functions to import cleaned DFs
    casper = clean_casper()
    coffee = clean_coffee()
    crocodile = clean_crocodile()
    improve = clean_improve()
    penguin = clean_penguin()
    renal_heir = clean_renal_heir()
    renal_heiritage = clean_renal_heiritage()
    # panther = clean_panther()
    # Merge
    harmonized = pd.concat([casper, coffee], join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, crocodile],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, improve],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, penguin],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, renal_heir],
                           join='outer', ignore_index=True)
    harmonized = pd.concat([harmonized, renal_heiritage],
                           join='outer', ignore_index=True)
    # harmonized = pd.concat([harmonized, panther],
    #                        join='outer', ignore_index=True)
    # harmonized = pd.concat([harmonized, panda],
    #                        join='outer', ignore_index=True)

    # Fix levels of categorical variables
    harmonized["visit"] = \
        pd.Categorical(harmonized["visit"],
                       categories=['screening', 'baseline', 'pre_surgery',
                                   '3_months_post_surgery', '12_months_post_surgery'],
                       ordered=True)
    harmonized["race"].replace(
        ["American Indian or Alaskan Native & White",
         "Black or African American & White",
         'American Indian or Alaskan Native & Black or African American',
         'Asian & White'], "More Than One", inplace=True)
    harmonized["race"].replace(
        {'Black/African American': 'Black or African American',
         "": "Unknown"}, inplace=True)
    harmonized["ethnicity"].replace({"": "Unknown"}, inplace=True)
    race_ethnicity = harmonized["race"] + \
        ", " + harmonized["ethnicity"]
    harmonized = pd.concat([harmonized, race_ethnicity], axis=1)
    harmonized.rename({0: "race_ethnicity"}, axis=1, inplace=True)
    # Date variables
    dates = ["dob", "date", "diabetes_dx_date"]
    harmonized[dates] = \
        harmonized[dates].apply(pd.to_datetime, errors='coerce')
    # Replace blanks with missing
    harmonized.replace("", np.nan, inplace=True)
    # Convert to numeric
    num_vars = ["height", "total_kidney_volume_ml", "left_kidney_volume_ml",
                "right_kidney_volume_ml"]
    harmonized[num_vars] = harmonized[num_vars].apply(
        pd.to_numeric, errors='coerce')
    # Calculated variables
    # Age
    age = round((harmonized["date"] - harmonized["dob"]).dt.days / 365.25, 2)
    harmonized = pd.concat([harmonized, age], axis=1)
    harmonized.rename({0: "age"}, axis=1, inplace=True)
    # BMI
    harmonized["bmi"] = pd.to_numeric(harmonized["weight"])/((pd.to_numeric(harmonized["height"])/100)**2)
    # Diabetes duration
    disease_duration = \
        round((harmonized["date"] -
              harmonized["diabetes_dx_date"]).dt.days / 365.25, 2)
    harmonized = pd.concat([harmonized, disease_duration], axis=1)
    harmonized.rename({0: "diabetes_duration"}, axis=1, inplace=True)
    # eGFR
    harmonized = calc_egfr(harmonized, age="age",
                           serum_creatinine="creatinine_s", cystatin_c="cystatin_c_s",
                           bun="bun", height="height", sex="sex", male="Male", female="Female", alpha=0.5)
    # Kidney volume
    harmonized["total_kidney_volume_ml"] = \
        harmonized["left_kidney_volume_ml"] + \
        harmonized["right_kidney_volume_ml"]
    # Height adjusted kidney volume
    harmonized["ht_adj_tkv"] = harmonized["total_kidney_volume_ml"] / \
        (harmonized["height"] / 100)
    # Calculate FSOC = bl_bold - pf_bold
    cols = [c for c in harmonized.columns if "_bl_" in c] + \
        [c for c in harmonized.columns if "_pf_" in c]
    harmonized[cols] = harmonized[cols].apply(
        pd.to_numeric, errors='coerce', axis=1)
    harmonized = harmonized.assign(
        fsoc_r_cortex=harmonized["bold_r_bl_cortex"] -
        harmonized["bold_r_pf_cortex"],
        fsoc_r_medulla=harmonized["bold_r_bl_medulla"] -
        harmonized["bold_r_pf_medulla"],
        fsoc_r_kidney=harmonized["bold_r_bl_kidney"] -
        harmonized["bold_r_pf_kidney"],
        fsoc_l_cortex=harmonized["bold_l_bl_cortex"] -
        harmonized["bold_l_pf_cortex"],
        fsoc_l_medulla=harmonized["bold_l_bl_medulla"] -
        harmonized["bold_l_pf_medulla"],
        fsoc_l_kidney=harmonized["bold_l_bl_kidney"] -
        harmonized["bold_l_pf_kidney"])
    # UACR
    harmonized["acr_u"] = \
        pd.to_numeric(harmonized["microalbumin_u"], errors="coerce") * 100 / \
        pd.to_numeric(harmonized["creatinine_u"], errors="coerce")
    # Albuminuria
    alb = []
    for a in harmonized["acr_u"]:
        if a < 30:
            alb.append("A1")
        elif a >= 30 and a <= 300:
            alb.append("A2")
        elif a > 300:
            alb.append("A3")
        else:
            alb.append(np.nan)
    harmonized["albuminuria_cat"] = alb
    harmonized["elevated_albuminuria"] = pd.cut(
        harmonized["acr_u"], [-float("inf"), 30, float("inf")], right=False, labels=["No", "Yes"])
    # FFA suppression negative to be 0
    harmonized["ffa_suppression"] = np.where(
        harmonized["ffa_suppression"] < 0, 0, harmonized["ffa_suppression"])
    harmonized["p1_ffa_suppression"] = np.where(
        harmonized["p1_ffa_suppression"] < 0, 0, harmonized["p1_ffa_suppression"])
    harmonized["p2_ffa_suppression"] = np.where(
        harmonized["p2_ffa_suppression"] < 0, 0, harmonized["p2_ffa_suppression"])
    # FFA suppression combined
    harmonized = \
        harmonized.assign(ffa_suppression_combined=harmonized["ffa_suppression"].where(
            harmonized["ffa_suppression"].notnull(), harmonized["p2_ffa_suppression"]))
    # Fasting Insulin
    harmonized[["insulin_minus_20", "insulin_minus_10", "insulin_minus_5", "insulin_0"]] = harmonized[[
        "insulin_minus_20", "insulin_minus_10", "insulin_minus_5", "insulin_0"]].apply(pd.to_numeric)
    harmonized["fasting_insulin"] = \
        harmonized[["insulin_minus_20", "insulin_minus_10",
                    "insulin_minus_5", "insulin_0"]].apply(lambda x: x.mean(), axis=1)
    # Fasting FFA
    harmonized["fasting_ffa"] = \
        harmonized[["ffa_minus_20", "ffa_minus_10", "ffa_minus_5", "ffa_0"]].apply(
            lambda x: x.mean(), axis=1)
    # Adipose IR (fasting_insulin * fasting_ffa)
    harmonized = harmonized.assign(
        adipose_ir=harmonized["fasting_ffa"] * harmonized["fasting_insulin"])

    # Sort columns
    id_cols = ["record_id", "co_enroll_id", "study", "dob", "diabetes_dx_date",
               "sex", "race", "ethnicity", "visit", "procedure", "date", "group"]
    other_cols = harmonized.columns.difference(id_cols, sort=False).tolist()
    other_cols = natsorted(other_cols, alg=ns.IGNORECASE)
    harmonized = harmonized[id_cols + other_cols]
    # Sort rows
    harmonized.sort_values(
        ["study", "record_id", "visit", "procedure", "date"], inplace=True)
    # Format dates nicely
    harmonized[dates] = harmonized[dates].apply(
        lambda x: x.dt.strftime('%Y-%m-%d'))
    # Return
    harmonized = harmonized.astype(object)
    return harmonized
