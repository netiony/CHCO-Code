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
    os.chdir(os.path.expanduser('~'))
    os.chdir("/home/timvigers/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
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
    # FFA
    # See /home/timvigers/Work/CHCO/Petter Bjornstad/IHD/Background/Renal Heir Equations.docx
    ffa = [c for c in harmonized.columns if "ffa_" in c]
    ffa = natsorted(ffa, alg=ns.IGNORECASE)
    harmonized[ffa] = harmonized[ffa].apply(
        pd.to_numeric, errors='coerce')
    harmonized["baseline_ffa"] = \
        harmonized[['ffa_minus_10', 'ffa_minus_20', 'ffa_minus_5']].mean(axis=1)
    harmonized["steady_state_ffa"] = \
        harmonized[['ffa_220', 'ffa_230', 'ffa_240', 'ffa_250', 'ffa_260', 'ffa_270']].mean(axis=1)
    harmonized["ffa_supression"] = ((harmonized["baseline_ffa"]-harmonized["steady_state_ffa"])/harmonized["baseline_ffa"])*100
    # Insulin
    ins = ['insulin_220', 'insulin_230', 'insulin_240', 'insulin_245',
        'insulin_249', 'insulin_252', 'insulin_253', 'insulin_254', 'insulin_255',
        'insulin_250','insulin_260','insulin_270']
    harmonized[ins] = harmonized[ins].apply(
        pd.to_numeric, errors='coerce')    
    harmonized["steady_state_insulin"] = harmonized[ins].mean(axis=1) * 6
    # C peptide
    cpep = ['cpeptide_220', 'cpeptide_230', 'cpeptide_240', 'cpeptide_245', 'cpeptide_249', 'cpeptide_250', 
    'cpeptide_252', 'cpeptide_253', 'cpeptide_254', 'cpeptide_255']
    harmonized[ins] = harmonized[ins].apply(
        pd.to_numeric, errors='coerce')    
    harmonized["steady_state_insulin"] = harmonized[ins].mean(axis=1) * 6
    # Sort
    harmonized.sort_values(
        ["study", "record_id", "visit", "procedure", "date"], inplace=True)
    # Return
    return harmonized
