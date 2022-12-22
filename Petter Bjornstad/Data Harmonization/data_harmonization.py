"""
This code is designed to pull data from multiple REDCap projects and harmonize
the data in a single dataset. Some studies are cross-sectional but include
measures at multiple visits, and some studies are longitudinal. So, this code
outputs data in a "semi-long" format with one row per study procedure, and a
visit column for longitudinal clustering. The data cleaning process for each individual dataset is a separate function, and this function puts them together and performs some formatting tweaks.
"""
__author__ = "Tim Vigers"
__credits__ = ["Tim Vigers"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"


def harmonize_data():
    # Libraries
    import os
    os.chdir(
        "C:/Users/timbv/Documents/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
    import pandas as pd
    import numpy as np
    from casper import clean_casper
    from coffee import clean_coffee
    from crocodile import clean_crocodile
    from improve import clean_improve
    from penguin import clean_penguin
    from renal_heir import clean_renal_heir
    from harmonization_functions import calc_egfr
    from natsort import natsorted, ns
    # Use individual data functions to import cleaned DFs
    casper = clean_casper()
    coffee = clean_coffee()
    crocodile = clean_crocodile()
    improve = clean_improve()
    penguin = clean_penguin()
    renal_heir = clean_renal_heir()
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
    # Replace missing values
    rep = [-97, -98, -99, -997, -998, -999, -9997, -9998, -9999]
    rep = rep + [str(r) for r in rep]
    harmonized.replace(rep, "", inplace=True)
    # Fix levels of categorical variables
    harmonized["visit"] = \
        pd.Categorical(harmonized["visit"],
                       categories=['baseline', 'pre_surgery',
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
    harmonized["race_ethnicity"] = harmonized["race"] + \
        ", " + harmonized["ethnicity"]
    # Date variables
    harmonized[["dob", "date"]] = harmonized[["dob", "date"]].apply(
        pd.to_datetime, errors='coerce')
    # Replace blanks with missing
    harmonized.replace("", np.nan, inplace=True)
    # Convert to numeric
    num_vars = ["height", "total_kidney_volume_ml", "left_kidney_volume_ml",
                "right_kidney_volume_ml"]
    harmonized[num_vars] = harmonized[num_vars].apply(
        pd.to_numeric, errors='coerce')
    # Calculated variables
    age = round((harmonized["date"] -
                 harmonized["dob"]).dt.days / 365.25, 2)
    harmonized = pd.concat([harmonized, age], axis=1)
    harmonized.rename({0: "age"}, axis=1, inplace=True)
    harmonized = calc_egfr(harmonized, age="age",
                           serum_creatinine="creatinine_s", cystatin_c="cystatin_c_s",
                           bun="bun", height="height", sex="sex", male="Male", female="Female", alpha=0.5)
    # Kidney volume
    harmonized["total_kidney_volume_ml"] = \
        harmonized["left_kidney_volume_ml"] + \
        harmonized["right_kidney_volume_ml"]
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
    # Sort
    harmonized.sort_values(
        ["study", "record_id", "visit", "procedure", "date"], inplace=True)
    # Return
    harmonized = harmonized.astype(object)
    return harmonized
