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
os.chdir(os.path.expanduser('~'))
os.chdir("GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
from data_harmonization import harmonize_data
from datetime import datetime
import pandas as pd
# Get dataset
df = harmonize_data()
# Write
df.to_csv("~/temp.csv", index=False)
# df.to_csv("/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", index=False)
# # Update for Michigan
# time = datetime.now().strftime('%Y_%m_%d_%I%M%p')
# # Re-save dictionary
# dictionary = pd.read_csv(
#     "/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Data Harmonization/data_dictionary_master.csv")
# dictionary.to_csv("/Users/timvigers/Dropbox/Shared/Michigan/CHCO_Sample_IDs_Clinical/chco_data_dictionary_" +
#                   time + ".csv", index=False)
# # Save cleaned data (de-identified)
# df = df.drop(["dob"], axis=1)
# df.to_csv("/Users/timvigers/Dropbox/Shared/Michigan/CHCO_Sample_IDs_Clinical/chco_harmonized_dataset_" +
#           time + ".csv", index=False)
