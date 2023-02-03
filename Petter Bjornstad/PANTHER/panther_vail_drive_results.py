# import OS module
import os
import pandas as pd
import numpy as np
# Get the list of all files
wd = "/Users/timvigers/Dropbox/Work/Panther/"
files = os.listdir(wd + "Data_Raw")
files.sort()
# Combine
res = {"MRN": [], "TestName": [], "Result Numeric": []}
for f in files:
    df = pd.read_csv(wd + "Data_Raw/" + f,
                     usecols=["MRN", "TestName", "Result Numeric"])
    res["MRN"].extend(list(df["MRN"]))
    res["TestName"].extend(list(df["TestName"]))
    res["Result Numeric"].extend(list(df["Result Numeric"]))
# As DF
res = pd.DataFrame(res)
# Format
res = res[res["TestName"] != "R MISC LOB"]
res = res[res["TestName"] != "UrProc"]
# Get IDs
mrns = pd.read_csv(wd + "Data_Cleaned/panther_mrns.csv")
ids = []
for i in res["MRN"]:
    if i in mrns["mrn"].to_list():
        ids.append(mrns.loc[mrns["mrn"] == i, "record_id"].values[0])
    else:
        ids.append("")
res["record_id"] = ids
# Drop rows and columns
res = res[["record_id", 'TestName', 'Result Numeric']]
res = res[res["record_id"] != ""]
# Fix missing
res['Result Numeric'].replace([-999.99, ""], np.nan, inplace=True)
res.dropna(inplace=True)
res.drop_duplicates(subset=['record_id', 'TestName'], inplace=True)
res = res.pivot(index="record_id", columns="TestName",
                values="Result Numeric")
# Write
res.to_csv(wd + "Data_Cleaned/" + "panther_results.csv")
