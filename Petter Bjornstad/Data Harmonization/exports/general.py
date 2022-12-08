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
os.chdir("Documents/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
import pandas as pd
import numpy as np
from data_harmonization import harmonize_data
# Get dataset
df = harmonize_data()
# Write
df.to_csv("/Volumes/Work/CHCO/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", index=False)
