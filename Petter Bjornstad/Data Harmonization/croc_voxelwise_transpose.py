import os
import pandas as pd
import numpy as np
wd = "/Users/timvigers/Desktop/CRC-VoxelWise/CRC-voxelwiseVectors/"
# List files
files = os.listdir(wd)
files.sort()
# Need to make sure that all voxel lists are the same length (some people have more voxels than others, likely due to larger kidneys)
ns = []
for file in files:
	voxels = pd.read_csv(wd+file,header=None)
	ns.append(voxels.shape[0])
length = max(ns)
# Create dictionary with the correct number of columns
vox_cols = ["voxel_"+str(i).zfill(1) for i in range(length)]
cols = ["record_id","region","parameter"]+vox_cols
res = {k: [] for k in cols}
for file in files:
	voxels = pd.read_csv(wd+file,header=None)
	voxels = voxels.iloc[:,0].tolist()
	voxels += [''] * (length - len(voxels))
	f = file.replace("-","_").split("_")
	res["record_id"].append(f[1])
	res["region"].append(f[2].lower())
	res["parameter"].append(f[3].lower().replace(".csv",""))
	for i in range(len(voxels)):
		res[vox_cols[i]].append(voxels[i])
# Convert to DF
df = pd.DataFrame(res)
df[df == ''] = np.nan
df.to_csv("/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/CROCODILE/Data_Cleaned/voxelwise_rows.csv",index=False)