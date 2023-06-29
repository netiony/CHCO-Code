library(Seurat)
library(ggplot2)
# Read in scRNA object
so <- readRDS("/home/tim/UCD/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/scRNA/Data_Clean/seurat_data_no_computations.RDS")
# Combined groups
so$SGLT2i <- factor(so$T2D_HC_Phil,
                    levels = c("HC", "T2D", "T2Di"),
                    labels = c("SGLT2i-", "SGLT2i-", "SGLT2i+")
)
# T2D only
so <- so[, so$T2D_HC_Phil %in% c("T2D", "T2D_post",     "T2Di" )]
# Remove 12 months
so <- so[, !grepl("_12M",so$michigan_id)]
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
# UMAP plot
p <- FeaturePlot(so, features = "MTAP")
p <- LabelClusters(plot = p, id = "ident")
ggsave("~/umap.png",plot=p)
# Ridgeplot
p = VlnPlot(so, features = "MTAP",group.by = "SGLT2i")
ggsave("~/ridge.png",plot=p)
# Violin
p = RidgePlot(so, features = "MTAP")
ggsave("~/violin.png",plot=p)
# Compare MTAP between SLGT2i- and SGLT2i+
markers <- FindMarkers(so,features = "MTAP", ident.1 = "SGLT2i-",
                       group.by = 'SGLT2i',logfc.threshold = 0)
markers
