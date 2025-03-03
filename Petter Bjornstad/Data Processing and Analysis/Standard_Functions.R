#Single Cell RNA Sequencing Quality Control & Pre-Processing Function 
qc_function <- function(so,cut_low,cut_high,dims,mt_pct) {
so <- subset(so, subset = nFeature_RNA > cut_low & nFeature_RNA < cut_high) # implementing this based on the responses from Dylan and Matteo
so <- NormalizeData(so)
so <- ScaleData(so, features = VariableFeatures(so))
so <- RunPCA(so, features = VariableFeatures(so))
ElbowPlot(so)
so <- FindNeighbors(so, dims = 1:dims)
so <- FindClusters(so)
so <- RunUMAP(so, dims = 1:dims, reduction.key = "UMAP_")
}

so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 5000) # implementing this based on the responses from Dylan and Matteo
so <- NormalizeData(so)
so <- ScaleData(so, features = VariableFeatures(so))
so <- RunPCA(so, features = VariableFeatures(so))
ElbowPlot(so)
so <- FindNeighbors(so, dims = 1:10)
so <- FindClusters(so)
so <- RunUMAP(so, dims = 1:30, reduction.key = "UMAP_")