import pandas as pd
import scanpy as sc
import anndata as ad
t = sc.read_h5ad("/Users/timvigers/Dropbox/Work/CHCO/Petter Bjornstad/scRNA/Data_Raw/"
                 "CU_Anschutz_scRNAseq_data/PB_40datasets_soupxcleaned070821_RPC2_Sean_selected.h5ad")
temp = ad.AnnData.to_df(t)
