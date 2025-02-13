#Load in seurat .h5seurat object& convert to .rds
install.packages("Seurat")
install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk",force=T)

library(SeuratDisk)
library(Seurat)

# Load your .h5seurat file
seurat_obj <- LoadH5Seurat("/Users/hhampson/Downloads/KPMP.h5Seurat")

# Save as .RDS
saveRDS(seurat_obj, file = "/Users/hhampson/Downloads/KPMP.rds")
