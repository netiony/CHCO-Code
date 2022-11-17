"""
This function takes REDCap metadata from Petter's studies and tries to match column names across all studies. For example, serum creatinine is sometimes called "serum_creatinine" and sometimes called "creat_s."
"""

# Libraries
import redcap
import pandas as pd
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
# Import REDCap metadata
tokens = pd.read_csv(
    "~/Dropbox/Work/CHCO/Petter Bjornstad/Data Harmonization/api_tokens.csv")
uri = "https://redcap.ucdenver.edu/api/"
# Get project metadata
croc_meta = pd.DataFrame(redcap.Project(url=uri,
                                        token=tokens.loc[tokens["Study"] == "CROCODILE", "Token"].iloc[0]).metadata)
croc_meta = croc_meta[~croc_meta["form_name"].isin(
    ['metabolomics', 'lipidomics', 'proteomics', 'voxelwise'])]
imp_meta = pd.DataFrame(redcap.Project(url=uri,
                                       token=tokens.loc[tokens["Study"] == "IMPROVE", "Token"].iloc[0]).metadata)
# Find each column name's closest match
matches = []
for var in croc_meta["field_label"]:
    closest = process.extractBests(var, imp_meta["field_label"].to_list(),
                                   scorer=fuzz.token_set_ratio, limit=1)
    closest = [(var,) + c for c in closest]
    matches.append(closest)
# To dataframe
matches = [item for sublist in matches for item in sublist]
matches = pd.DataFrame(matches,
                       columns=["original_variable", "closest_match", "score"])
matches.to_csv("~/croc_imp_matches.csv", index=False)
