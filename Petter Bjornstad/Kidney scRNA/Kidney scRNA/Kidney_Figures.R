#Figures
#By sex
jpeg(file = fs::path(dir.home,"Results_and_Figures","Umap_by_Sex.jpeg"))
so_kidney_sc <- RunUMAP(so_kidney_sc, dims = 1:15)
DimPlot(so_kidney_sc, reduction = "umap", group.by = "sex") 
dev.off()
#By sglti2
jpeg(file = fs::path(dir.home,"Results_and_Figures","Umap_by_SGLT2i.jpeg"))
DimPlot(so_kidney_sc, reduction = "umap", group.by = "sglt2i_ever") 
dev.off()
#by disease status
jpeg(file = fs::path(dir.home,"Results_and_Figures","Umap_by_disease_status.jpeg"))
DimPlot(so_kidney_sc, reduction = "umap", group.by = "group") 
dev.off()
#General 
jpeg(file = fs::path(dir.home,"Results_and_Figures","Umap.jpeg"))
DimPlot(so_kidney_sc, reduction = "umap")
dev.off()

#Diff exp by sex
jpeg(file = fs::path(dir.home,"Results_and_Figures","Umap.jpeg"))
EnhancedVolcano(m_top,
                lab = rownames(m_top),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Differentially Expressed Genes by Sex',
                subtitle = "Positive Log2 FC = Greater Expression in Females vs. Males",
                pCutoff = 0.05,
                # FCcutoff = 0,
                pointSize = 1.5,
                labSize = 4)


