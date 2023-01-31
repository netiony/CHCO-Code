"""
Helper functions for data harmonization.
"""


def combine_checkboxes(df, base_name="", levels=[], sep=" & ",
                       drop=True):
    # Libraries
    import pandas as pd
    import numpy as np
    # Make a copy of the dataframe
    data = df.copy()
    # A list of integer strings (one for each level)
    num = [str(s) for s in range(1, len(levels) + 1)]
    # Get column names
    cols = [base_name + "___" + s for s in num]
    # Which columns are checked
    values = data[cols].apply(lambda x: " & ".join(
        [str(s + 1) for s in np.where((x == 1) | (x == "1") | (x == "Checked"))[0]]), axis=1)
    # Dictionary for number to level translation
    column_translation = dict(zip(num, levels))
    # Replace numbers with levels
    for word, replacement in column_translation.items():
        values = [s.replace(word, replacement) for s in values]
    # Add to dataframe
    data[base_name] = values
    # Drop old columns
    if drop:
        data.drop(cols, axis=1, inplace=True)
    # Return data
    return data


def calc_egfr(df, age="age", serum_creatinine="creatinine_s",
              cystatin_c="cystatin_c_s", bun="bun", height="height", sex="sex",
              male="Male", female="Female", alpha=0.5):
    import pandas as pd
    import numpy as np
    # Make a copy of the dataframe
    data = df.copy()
    # Format input
    sex = data[sex].copy()
    sex.replace({male: "M", female: "F", "": np.nan}, inplace=True)
    qcr = np.floor(data[age])
    serum_creatinine = pd.to_numeric(data[serum_creatinine], errors="coerce")
    cystatin_c = pd.to_numeric(data[cystatin_c], errors="coerce")
    height = pd.to_numeric(data[height], errors="coerce")
    bun = pd.to_numeric(data[bun], errors="coerce")
    age = pd.to_numeric(data[age], errors="coerce")
    # Younger participants
    qcr.replace({8: 0.46, 9: 0.49, 10: 0.51, 11: 0.53,
                12: 0.57, 13: 0.59, 14: 0.61}, inplace=True)
    # Females
    qcr.loc[(qcr == 15) & (sex == "F")] = 0.64
    qcr.loc[(qcr == 16) & (sex == "F")] = 0.67
    qcr.loc[(qcr == 17) & (sex == "F")] = 0.69
    qcr.loc[(qcr == 18) & (sex == "F")] = 0.69
    qcr.loc[(qcr >= 19) & (sex == "F")] = 0.70
    # Males
    qcr.loc[(qcr == 15) & (sex == "M")] = 0.72
    qcr.loc[(qcr == 16) & (sex == "M")] = 0.78
    qcr.loc[(qcr == 17) & (sex == "M")] = 0.82
    qcr.loc[(qcr == 18) & (sex == "M")] = 0.85
    qcr.loc[(qcr == 19) & (sex == "M")] = 0.88
    qcr.loc[(qcr > 19) & (sex == "M")] = 0.90
    # Calculate final metrics
    eGFR_fas_cr = 107.3 / (serum_creatinine / qcr)
    # eGFR FAS combined creatinine and cystatin-C
    f1 = serum_creatinine / qcr
    f2 = 1 - alpha
    f3 = cystatin_c / 0.82
    eGFR_fas_cr_cysc = 107.3 / ((0.5 * f1) + (f2 * f3))
    # eGFR Zappitelli
    eGFR_Zap = (507.76 * np.exp((0.003 * height)) /
                ((cystatin_c ** 0.635) * ((serum_creatinine * 88.4) ** 0.547)))
    # eGFR Schwartz
    m = sex.replace({"M": 1, "F": 0})
    eGFR_Schwartz = 39.1 * ((height / serum_creatinine) ** 0.516) * ((1.8 / cystatin_c) ** 0.294) * ((30 / bun) ** 0.169) * \
        (1.099 ** m) * ((height / 1.4) ** 0.188)
    # eGFR bedside Schwartz
    eGFR_bedside_Schwartz = (41.3 * (height / 100)) / serum_creatinine
    # CKD-EPI Creatinine 2021 https://www.kidney.org/content/ckd-epi-creatinine-equation-2021
    a1_f = -0.241
    a1_m = -0.302
    a2 = -1.200
    f = sex.replace({"M": 0, "F": 1})
    a = sex.replace({"M": -0.302, "F": -0.241})
    k = sex.replace({"M": 0.9, "F": 0.7})
    # Calculate
    eGFR_CKD_epi = 142 * (np.minimum(serum_creatinine / k, 1)**a) * \
        (np.maximum(serum_creatinine / k, 1)**-1.200) * \
        (0.9938**age) * (1.012 * f + (1 - f))
    # Add results to dataframe
    egfr = pd.concat(
        [eGFR_Schwartz, eGFR_bedside_Schwartz, eGFR_Zap, eGFR_fas_cr, eGFR_fas_cr_cysc, eGFR_CKD_epi], axis=1)
    egfr.columns = ["eGFR_Schwartz", "eGFR_bedside_Schwartz", "eGFR_Zap",
                    "eGFR_fas_cr", "eGFR_fas_cr_cysc", "eGFR_CKD_epi"]
    data = pd.concat([data, egfr], axis=1)
    # Return
    return data
