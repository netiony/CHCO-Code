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

    import os
    os.chdir(
        "/Users/timvigers/Documents/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
    # Libraries
    import pandas as pd
    import numpy as np
    from casper import clean_casper
    from coffee import clean_coffee
    from crocodile import clean_crocodile
    from improve import clean_improve
    from penguin import clean_penguin
    from renal_heir import clean_renal_heir
    from harmonization_functions import calc_egfr
    # Use individual data functions to import cleaned DFs
    casper = clean_casper()
    coffee = clean_coffee()
    crocodile = clean_crocodile()
    improve = clean_improve()
    penguin = clean_penguin()
    renal_heir = clean_renal_heir()
    # Merge
    harmonized = pd.merge(casper, coffee, how="outer")
    harmonized = pd.merge(harmonized, crocodile, how="outer")
    harmonized = pd.merge(harmonized, improve, how="outer")
    harmonized = pd.merge(harmonized, penguin, how="outer")
    harmonized = pd.merge(harmonized, renal_heir, how="outer")
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
    # Calculated variables
    age = round((harmonized["date"] -
                 harmonized["dob"]).dt.days / 365.25, 2)
    harmonized = pd.concat([harmonized, age], axis=1)
    harmonized.rename({0: "age"}, axis=1, inplace=True)
    harmonized = calc_egfr(harmonized, age="age",
                           serum_creatinine="creatinine_s", cystatin_c="cystatin_c_s",
                           bun="bun", height="height", sex="sex", male="Male", female="Female", alpha=0.5)
    # Calculate FSOC = bl_bold - pf_bold
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
    # Sort
    harmonized.sort_values(
        ["study", "record_id", "visit", "procedure", "date"], inplace=True)
    # Return
    return harmonized
