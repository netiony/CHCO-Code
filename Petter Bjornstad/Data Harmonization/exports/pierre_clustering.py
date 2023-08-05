"""
This code is designed pulls our harmonized data and selects the variables Fadhl
is using in his pre-/post-surgery analysis.
"""
__author__ = "Tim Vigers"
__credits__ = ["Tim Vigers"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Tim Vigers"
__email__ = "timothy.vigers@cuanschutz.edu"
__status__ = "Dev"
import os
import sys

sys.path.insert(
    0, os.path.expanduser("~") + "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization"
)
from data_harmonization import harmonize_data
import pandas as pd

# Get dataset
df = harmonize_data()
# Variable list from Pierre
dict = pd.read_csv(
    "/Volumes/Work/Petter Bjornstad/Data Harmonization/Data Exports/variables_for_Pierre_clustering.csv"
)
dict = dict[dict["Unnamed: 2"] == "x"]
cols = [c for c in dict["variable_name"]]
# Skip non-existent variables
skip = [
    "ht_adj_tkv",
    "avg_pcascl",
    "avg_k_r2",
    "avg_c_r2",
    "avg_m_r2",
    "avg_k_fsoc",
    "avg_c_fsoc",
    "avg_m_fsoc",
    "avg_c_adc",
    "avg_m_k2_wo_cyst_vw",
    "avg_c_k2_wo_cyst_vw",
    "gbm_thick_artmean",
    "gbm_thick_harmmean",
]
cols = [c for c in cols if c not in skip]
# Select columns
df = df[cols]
# Select rows
df = df[df["study"].isin(["CROCODILE", "PANDA", "RENAL-HEIR", "RENAL-HEIRitage"])]
# Write
df.to_csv(
    "/Volumes/Work/Petter Bjornstad/Data Harmonization/Data Exports/pierre_clustering.csv",
    index=False,
)
