library(SeuratDisk)
so = readRDS("/Users/timvigers/Library/CloudStorage/Dropbox/Work/Petter Bjornstad/scRNA/data_clean/seurat_data_no_computations.RDS")
SaveH5Seurat(so, filename = "/Users/timvigers/Library/CloudStorage/Dropbox/Work/Petter Bjornstad/scRNA/data_clean/anndata_no_computations.h5Seurat")
Convert("/Users/timvigers/Library/CloudStorage/Dropbox/Work/Petter Bjornstad/scRNA/data_clean/anndata_no_computations.h5Seurat", dest = "h5ad")
