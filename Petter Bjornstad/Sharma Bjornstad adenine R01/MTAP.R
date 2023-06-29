library(Seurat)
# Read in scRNA object
so <- readRDS("/Volumes/Peds Endo/Petter Bjornstad/scRNA/Data_Clean/seurat_data_no_computations.RDS")
# Limit to RH/IMPROVE baseline visit
so <- so[, !grepl("_12M", so$michigan_id)]
# Combined groups
so$SGLT2i <- factor(so$T2D_HC_Phil,
                    levels = c("HC", "T2D", "T2Di"),
                    labels = c("SGLT2i-", "SGLT2i-", "SGLT2i+")
)
# Normalize and scale
so <- NormalizeData(so)
so <- ScaleData(so, features = rownames(so))
# PCA
so <- RunPCA(so, features = VariableFeatures(object = so))
# Cluster cells
so <- FindNeighbors(so)
so <- FindClusters(so)
# Perform UMAP
so <- RunUMAP(so, dims = 1:20)
# General cell types as identifiers
so$generaltype <- sub("_.*", "", so$LR_clusters)
Idents(so) <- so$generaltype
# Compare MTAP between SLGT2i- and SGLT2i+
markers <- FindMarkers(so, ident.1 = "SGLT2i-", group.by = 'SGLT2i')
