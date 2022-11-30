"""
This code is designed pulls our harmonized data and selects the variables Allie Shapiro is using in her neuroimaging analysis.
"""
__author__ = "Tim Vigers"
__credits__ = ["Tim Vigers"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"

import os
os.chdir(
    "/Users/timvigers/Documents/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
import pandas as pd
from data_harmonization import harmonize_data
# Get dataset
df = harmonize_data()
# Filter and select
mri = df[df["record_id"].isin(
    ["CRC-10", "CRC-12", "CRC-13", "CRC-18", "CRC-20", "CRC-24", "CRC-28", "CRC-31", "CRC-33", "CRC-34", "CRC-35"])]
mri = mri.groupby(by=["record_id"]).agg("last").reset_index()

cols = ["subject_id", "group", "age", "gender", "bmi",
    'pasl2d_left', 'pasl2d_right', 'pcasl3d_left', 'pcasl3d_right',
    "volume_left", "volume_right"] + \
        [c for c in mri.columns if "_bl_" in c] + \
            [c for c in mri.columns if "fsoc" in c]

cols = [c for c in df.columns if "adc" in c]

subject_id % in % paste0("CRC -", c())) % > %
select(,
       contains("adc"), volume_left, volume_right, contains("_bl_"),
       contains("fsoc"))
# Write CSV
write.csv(mri, na = "", row.names = F,
          file = paste0("~/Desktop/bold_mri_",
                      format(Sys.time(), "%Y_%m_%d_%H%M"), ".csv"))
