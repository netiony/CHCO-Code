# Libraries
import os
import sys
sys.path.insert(0, os.path.expanduser('~') +
                  "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
from data_harmonization import harmonize_data

clean = harmonize_data()
clean.to_csv("/Volumes/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv", index=False)
