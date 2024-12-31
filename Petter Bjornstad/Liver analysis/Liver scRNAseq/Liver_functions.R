#Define genes
sens_genes <- c("ACVR1B","ANG","ANGPT1","ANGPTL4","AREG","AXL","BEX3","BMP2","BMP6","C3","CCL1","CCL13",
                "CCL16","CCL2","CCL20","CCL24","CCL26","CCL3","CCL4","CCL5","CCL7","CCL8","CD55",
                "CD9","CSF1","CSF2","CSF2RB","CST4","CTNNB1","CTSB","CXCL1","CXCL10","CXCL12",
                "CXCL16","CXCL2","CXCL3","CXCL8","CXCR2","DKK1","EDN1","EGF","EGFR","EREG","ESM1",
                "ETS2","FAS","FGF1","FGF2","FGF7","GDF15","GEM","GMFG","HGF","HMGB1","ICAM1","ICAM3",
                "IGF1","IGFBP1","IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7","IL10","IL13",
                "IL15","IL18","IL1A","IL1B","IL2","IL32","IL6","IL6ST","IL7","INHA","IQGAP2","ITGA2",
                "ITPKA","JUN","KITLG","LCP1","MIF","MMP1","MMP10","MMP12","MMP13","MMP14","MMP2",
                "MMP3","MMP9","NAP1L4","NRG1","PAPPA","PECAM1","PGF","PIGF","PLAT","PLAU","PLAUR",
                "PTBP1","PTGER2","PTGES","RPS6KA5","SCAMP4","SELPLG","SEMA3F","SERPINB4","SERPINE1",
                "SERPINE2","SPP1","SPX","TIMP2","TNF","TNFRSF10C","TNFRSF11B","TNFRSF1A","TNFRSF1B",
                "TUBGCP2","VEGFA","VEGFC","VGF","WNT16","WNT2")
sens_genes <- c(sens_genes,"CDKN1A")

#DEGS ----
# de.markers(so_liver_sn, genes, "diagnosis_of_diabetes", id2 = "No", id1 = "Yes", NULL, "_top")
# seurat_object <- so_kidney_sc
# group.by <- "epic_glp1ra_1"
# id1 <- "No"
# id2 <- "Yes"
# celltype="PT"
# genes <- sens_genes
# extension="_top"
# # function for de.markers
de.markers <- function(seurat_object, genes, group.by, id1, id2, celltype, extension,logfc.threshold,min.pct){
  m = FindMarkers(seurat_object, features = genes,group.by = group.by,ident.1 = id1, 
                  ident.2 = id2, subset.ident = celltype, verbose = F, logfc.threshold=0.001,
                  min.pct = 0.001)
  m$p_val_adj = p.adjust(m$p_val,method = "fdr")
  m <- m %>% 
    rownames_to_column('gene') %>%
    arrange(p_val) %>%
    column_to_rownames('gene') %>%
    dplyr::select(avg_log2FC,pct.1,pct.2,p_val,p_val_adj) %>%
    filter(!is.na(p_val))
  
  genes_subset <- rownames(m)[m$p_val <= 0.05]
  
  if (length(genes_subset) > 0){
    assign(paste0("genes_subset", extension), genes_subset, envir = .GlobalEnv)
  }
  assign(paste0("m", extension), m, envir = .GlobalEnv)
  return(knitr::kable(m, digits = 3
  ))
}

#DEG & GSEA Function ----
degs_fxn <- function(so,cell,exposure,gene_set,exp_group,ref_group,enrichment,top_gsea) {
  DefaultAssay(so) <- "RNA"
  
  #Conidtion names
  condition <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group," vs. ",ref_group,")")
  
  #Differential Expression by Group
  if (!is.null(cell)) {
  Idents(so) <- so$celltype2
  cell_name <- str_replace_all(cell,"/","_")
  }
  # sens_genes <- c(sens_genes,"CDKN1A")
  # de.markers(so, gene_set, "group", id2 = "neither", id1 = "both", "PT", "")
  de.markers(so, gene_set, exposure, id2 = ref_group, id1 = exp_group, cell, "")
  colnames(m)[2] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group,")")
  colnames(m)[3] <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",ref_group,")")
  
  #Plot DEGs
  m_top <- m
  significant_genes <- m_top %>% filter(p_val_adj < 0.05)
  
  # # Select the top 10 positive and top 10 negative log2FC genes that are significant
  # top_genes <- rbind(
  #   significant_genes %>% arrange(desc(avg_log2FC)) %>% head(40),  # Top 10 positive log2FC
  #   significant_genes %>% arrange(avg_log2FC) %>% head(40)         # Top 10 negative log2FC
  # )
  top_genes <- rbind(
    significant_genes %>% filter(avg_log2FC > 0) %>% arrange(p_val_adj) %>% head(20),  # Top 10 positive significant genes
    significant_genes %>% filter(avg_log2FC < 0) %>% arrange(p_val_adj) %>% head(20)  # Top 10 negative significant genes
  )
  
  if (!is.null(cell)){
    title <- paste0("DEGs in ",cell_name," cells for ",condition)
  } else {
    title <- paste0("Bulk DEGs for ",condition)
  }
  labels <- ifelse(rownames(m_top) %in% rownames(top_genes), rownames(m_top), NA)
  p <- EnhancedVolcano(m_top,
                       lab = labels,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = title,
                       subtitle = paste0("Positive Log2 FC = Greater Expression in ", condition,"\n",
                                         "(Significant at FDR-P<0.05, FC Threshold = 0.5)"),
                       pCutoff = 0.05,
                       FCcutoff = 0.5,
                       labFace = 'bold',
                       pointSize = 4,
                       labSize = 5,
                       drawConnectors = TRUE,
                       widthConnectors = 1.0,
                       colConnectors = 'black',
                       legendPosition=NULL,
                       boxedLabels = TRUE,
                       max.overlaps=60)
  if (!is.null(cell)){
  filename <- paste0("DEGs_in_",cell_name,"_cells_for_",condition,".pdf")
  } else {
    filename <- paste0("Bulk_DEGs_for_",condition,".pdf") 
  }
  pdf(fs::path(dir.results,filename),width=20,height=20)
  plot(p)
  dev.off()
  
  #GSEA
  if (enrichment=="Yes") {
    #Gene set enrichment analysis
    gc()
    sce_sn_hep <- as.SingleCellExperiment(so)
    ## C2 category is according to canonical pathways: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707969/pdf/nihms-743907.pdf
    geneSets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
    ### filter background to only include genes that we assessed
    geneSets$gene_symbol <- toupper(geneSets$gene_symbol)
    geneSets <- geneSets[geneSets$gene_symbol %in% names(sce_sn_hep),]
    m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
    stats <- m$p_val_adj
    names(stats) <- rownames(m)
    eaRes <- fgsea(pathways = m_list, stats = na.omit(stats))
    #ooEA <- order(eaRes$pval, decreasing = FALSE)
    #kable(head(eaRes[ooEA, 1:7], n = 20))
    # Convert the leadingEdge column to comma-separated strings
    eaRes$leadingEdge <- sapply(eaRes$leadingEdge, function(x) paste(x, collapse = ", "))
    gc()
    #Significant pathways
    sig <- eaRes %>% 
      filter(padj<0.05)
    
    #Plot top pathways
    # Subset top pathways for visualization
    top_pathways <- eaRes[order(-abs(eaRes$NES)), ][1:top_gsea, ]  # Top 10 pathways based on adjusted p-value
    
    # Define a significance threshold
    significance_threshold <- 0.05
    
    # Add a significance column for coloring
    top_pathways$significance <- ifelse(
      top_pathways$padj > significance_threshold, "Non-significant",
      ifelse(top_pathways$NES > 0, "Positive Significant", "Negative Significant")
    )
    
    if (!is.null(cell)) {
      title <- paste0("Top Enriched Pathways by NES (GSEA) in ",cell_name," cells among ",condition)
    } else {
      title <- paste0("Top Enriched Pathways by NES (GSEA) for ",condition," (Bulk)")
    }
    # Create a bar plot
    gsea_plot <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = significance)) +
      geom_bar(stat = "identity") +
      coord_flip() +  # Flip coordinates for better readability
      scale_fill_manual(
        values = c(
          "Positive Significant" = "red",
          "Negative Significant" = "blue",
          "Non-significant" = "gray"
        ),
        name = "Significance"
      ) +
      labs(
        title = title,
        x = "Pathway",
        y = "Normalized Enrichment Score (NES)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold")
      )
    
  
    if (!is.null(cell)){
      filename <- paste0("GSEA_top_",top_gsea,"_pathways_",cell_name,"_cells_for",condition,".pdf")
    } else {
      filename <- paste0("Bulk_GSEA_for_",condition,".pdf") 
    }
    pdf(fs::path(dir.results,filename),width=20,height=20)
    plot(gsea_plot)
    dev.off()
    
  }
  #Make sure gene nemaes are in printed file
  deg_results <- m 
  deg_results$Gene <- rownames(deg_results)
    
  #save results to excel file
  write_multiple_sheets <- function(output_file, sheet_data) {
    # Create a new workbook
    wb <- createWorkbook()
    
    # Loop over the list of data frames
    for (sheet_name in names(sheet_data)) {
      # Add a new sheet to the workbook
      addWorksheet(wb, sheet_name)
      
      # Write the data frame to the sheet
      writeData(wb, sheet_name, sheet_data[[sheet_name]])
    }
    # Save the workbook to the specified file
    saveWorkbook(wb, file = output_file, overwrite = TRUE)
    
    # Print a message
    message("Excel file with multiple sheets created: ", output_file)
  } 
  # Example usage
  df1 <- deg_results
  df2 <- eaRes
  
  # Specify the file name and data
  if (!is.null(cell)){
  output_file <- fs::path(dir.results,paste0("Results_",cell_name,"_cells_for_",condition,".xlsx"))  
  } else {
  output_file <- fs::path(dir.results,paste0("Bulk_Results_for_",condition,".xlsx"))
  }

  sheet_data <- list(
    "DEG" = df1,
    "Pathway_Results" = df2
  )
  
  # Call the function
  write_multiple_sheets(output_file, sheet_data)
}

#Visualize Function ----
##a. Bulk ----


##b. Single Cell ----
visualize_function <- function(exposure,cell) {
  
  # List all files in the results directory
  all_files <- dir_ls(path = dir.results, glob = "*.xlsx")  # List only .xlsx files
  
  # Filter the files that contain the current cell type and exposure in the file name
  matching_files <- all_files[str_detect(all_files, cell) & str_detect(all_files, exposure)]
  
  # Check if matching files exist for the current cell type and exposure
  if (length(matching_files) > 0) {
    # Ensure we have at least 3 files
    if (length(matching_files) >= 3) {
      # Read in the first, second, and third files as group1, group2, and group3
      group1 <- read.xlsx(matching_files[1])
      group2 <- read.xlsx(matching_files[2])
      group3 <- read.xlsx(matching_files[3])
      
    } 
  }
  # Create an empty list to store titles
  group_titles <- list()
  
  # Loop through each file in matching_files
  for (i in 1:length(matching_files)) {
    file_path <- matching_files[i]
    
    # Extract the comparison group from the file name
    comparison_group <- sub(".*\\((.*?)\\).*", "\\1", basename(file_path))
    
    # Assign titles to group1_title, group2_title, etc.
    group_titles[[paste0("group", i, "_title")]] <- comparison_group
  }
  
  # Now you can access the group titles like this
  group1_title <- group_titles$group1_title
  group2_title <- group_titles$group2_title
  group3_title <- group_titles$group3_title
  
  # Define significance threshold
  significance_threshold <- 0.05
  
  # Filter significant genes and classify as upregulated or downregulated
  get_gene_sets <- function(data) {
    upregulated <- data %>%
      filter(p_val_adj < significance_threshold & avg_log2FC > 0) %>%
      pull(Gene) # Replace 'gene' with the actual column name for gene IDs
    downregulated <- data %>%
      filter(p_val_adj < significance_threshold & avg_log2FC < 0) %>%
      pull(Gene) # Replace 'gene' with the actual column name for gene IDs
    list(up = upregulated, down = downregulated)
  }
  
  # Calculate gene sets for each group
  group1_sets <- get_gene_sets(group1)
  group2_sets <- get_gene_sets(group2)
  group3_sets <- get_gene_sets(group3)
  
  # Extract upregulated and downregulated sets
  up1 <- group1_sets$up
  down1 <- group1_sets$down
  up2 <- group2_sets$up
  down2 <- group2_sets$down
  up3 <- group3_sets$up
  down3 <- group3_sets$down
  
  # Find mutual upregulated genes between groups
  mutual_up12 <- intersect(up1, up2)  # Upregulated in both Group 1 and Group 2
  mutual_up13 <- intersect(up1, up3)  # Upregulated in both Group 1 and Group 3
  mutual_up23 <- intersect(up2, up3)  # Upregulated in both Group 2 and Group 3
  mutual_up_all <- Reduce(intersect, list(up1, up2, up3))  # Upregulated in all three groups
  
  # Find mutual downregulated genes between groups
  mutual_down12 <- intersect(down1, down2)  # Downregulated in both Group 1 and Group 2
  mutual_down13 <- intersect(down1, down3)  # Downregulated in both Group 1 and Group 3
  mutual_down23 <- intersect(down2, down3)  # Downregulated in both Group 2 and Group 3
  mutual_down_all <- Reduce(intersect, list(down1, down2, down3))  # Downregulated in all three groups
  
  # Find unique upregulated genes in each group
  unique_up1 <- setdiff(up1, union(up2, up3))  # Upregulated only in Group 1
  unique_up2 <- setdiff(up2, union(up1, up3))  # Upregulated only in Group 2
  unique_up3 <- setdiff(up3, union(up1, up2))  # Upregulated only in Group 3
  
  # Find unique downregulated genes in each group
  unique_down1 <- setdiff(down1, union(down2, down3))  # Downregulated only in Group 1
  unique_down2 <- setdiff(down2, union(down1, down3))  # Downregulated only in Group 2
  unique_down3 <- setdiff(down3, union(down1, down2))  # Downregulated only in Group 3
  
  
  
  # Function to draw Venn diagrams for upregulated and downregulated genes
  # Function to draw Venn diagrams for upregulated and downregulated genes
  draw_venn <- function(up1, up2, up3, down1, down2, down3) {
    
    # Plot Venn diagram for upregulated genes
    venn_up <- venn.diagram(
      x = list(
        Group1 = up1,
        Group2 = up2,
        Group3 = up3
      ),
      category.names = c("Group 1", "Group 2", "Group 3"),
      filename = NULL, # Don't save as file, display in plot window
      output = TRUE,
      main = "Upregulated Genes",
      col = "black",
      fill = c("lightblue", "lightgreen", "lightcoral"),
      alpha = 0.5,
      cex = 1.5,
      cat.cex = 1.5,
      cat.pos = c(0, 0, 0),  # Move Group 3 label to the bottom
      cat.dist = 0.05,
      margin = 0.1
    )
    
    # Plot Venn diagram for downregulated genes
    venn_down <- venn.diagram(
      x = list(
        Group1 = down1,
        Group2 = down2,
        Group3 = down3
      ),
      category.names = c("Group 1", "Group 2", "Group 3"),
      filename = NULL, # Don't save as file, display in plot window
      output = TRUE,
      main = "Downregulated Genes",
      col = "black",
      fill = c("lightblue", "lightgreen", "lightcoral"),
      alpha = 0.5,
      cex = 1.5,
      cat.cex = 1.5,
      cat.pos = c(0, 0, 0),  # Move Group 3 label to the bottom
      cat.dist = 0.05,
      margin = 0.1
    )
    
    # Draw the upregulated genes Venn diagram
    grid.draw(venn_up)
    
    # Draw a new page for downregulated Venn diagram
    grid.newpage() 
    
    # Draw the downregulated genes Venn diagram
    grid.draw(venn_down)
  }
  draw_venn(up1, up2, up3, down1, down2, down3)
  
 print(paste0("Group 1 = ",group1_title),paste0("Group 2 = ",group2_title),paste0("Group 3 = ",group3_title))
}

#so=so_liver_sn
#cell=NULL
#exposure=exp
#covariate="diagnosis_of_diabetes"
#gene_set=rownames(so_liver_sn)
#batch_size=100
#exp_group=NULL
#ref_group=NULL

#Mast Function----
mast_fxn <- function(so,cell,exposure,covariate,gene_set,batch_size,exp_group,ref_group) {
  DefaultAssay(so) <- "RNA"

  #Conidtion names
  if (!is.null(exp_group)) {
  condition <- paste0(str_to_title(str_replace_all(exposure,"_"," "))," (",exp_group," vs. ",ref_group,")")
  } else {
    condition <- paste0(str_to_title(exposure))
  }

  #Filter so to specified cell type
  if (!is.null(cell)){
    so <- subset(so,celltype2==cell)
    DefaultAssay(so) <- "RNA"
  } 
  
  #Association by Group
  combined_results <- data.frame()
  num_genes <- length(gene_set)
  for (start_batch in seq(1,num_genes,by=batch_size)) {
    
    end_batch <- min(start_batch + batch_size - 1, num_genes)
    batch <- gene_set[start_batch:end_batch]
    
    # Subset the Seurat object to include only genes in the gene set
    subset <- intersect(batch, rownames(so))
    so_sub <- subset(so_sub, features = subset)
    
    # Extract the expression data matrix from so (e.g., normalized counts)
    expression_matrix <- as.matrix(GetAssayData(so_sub, layer = "data"))
    
    # Extract metadata
    cell_metadata <- so_sub@meta.data
    #cell_metadata[[exposure]] <- factor(cell_metadata[[exposure]])
    
    # Create SingleCellAssay object
    sca <- FromMatrix(exprsArray = expression_matrix, cData = cell_metadata)
    
    # Assuming gene_set is a vector of gene names
    sca_gene_set <- sca[rownames(sca) %in% gene_set, ]
    
    #Define the formula
    model_formula <- as.formula(paste0("~", exposure))
    model_formula <- ~ast
    exp <- paste0(exposure,exp_group)
    
    # Run the linear model with zlm
    zlm_results <- zlm(formula = model_formula, sca = sca_gene_set)
    
    # Summarize results and perform likelihood ratio test (LRT) for outcome
    summary_zlm <- summary(zlm_results, doLRT = "ast")
    
    # Convert the summary to a data table for easy manipulation
    summary_dt <- summary_zlm$datatable
    
    #Overall with Hurdle and logFC
    fcHurdle1 <- merge(summary_dt[contrast==exp & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summary_dt[contrast==exp & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    fcHurdle1[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    fcHurdle1$component <- "FC_H"
    fcHurdle2 <- merge(summary_dt[contrast==exp & component=='C',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                       summary_dt[contrast==exp & component=='C', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    fcHurdle2[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    fcHurdle2$component <- "C"
    fcHurdle3 <- merge(summary_dt[contrast==exp & component=='D',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                       summary_dt[contrast==exp & component=='D', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    fcHurdle3[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    fcHurdle3$component <- "D"
    
    fcHurdle <- rbind(fcHurdle1,fcHurdle2)
    fcHurdle <- rbind(fcHurdle,fcHurdle3)
    
    result1 <- fcHurdle
    result1 <- result1 %>%
      filter(fdr<0.05) %>% 
      dplyr::rename(LogFoldChange=coef) %>% 
      dplyr::rename(Gene=primerid)
    
    # Append the batch results to the combined results using rbindlist
    # combined_results <- rbind(combined_results, batch_results)
    combined_results <- rbindlist(list(combined_results, result1), use.names = TRUE, fill = TRUE)
    
  }

  # Specify the file name and data
  if (!is.null(cell)){
    output_file <- paste0("Results_",cell_name,"_cells_for_",condition,".xlsx")
  } else {
    output_file <- paste0("Results_for_",condition,".xlsx")
  }

  write.csv(combined_results,fs::path(output_file))
  
  # Filter top 10 positive and negative log2FoldChange
  top_pos <- as.data.frame(combined_results) %>%
    filter(component=="FC_H") %>% 
    filter(fdr<0.05) %>%
    filter(LogFoldChange>0) %>% 
    arrange(desc(LogFoldChange)) %>%
    dplyr::slice(1:30) 
    # dplyr::rename(Gene=primerid,
    #               LogFC=coef)
  
  top_neg <- as.data.frame(combined_results) %>%
    filter(component=="FC_H") %>% 
    filter(fdr<0.05) %>%
    filter(LogFoldChange<0)  %>% 
      arrange(LogFoldChange) %>%
      dplyr::slice(1:30) 
    # dplyr::rename(Gene=primerid,
    #               LogFC=coef)
  if(dim(top_pos)[[1]] > 0 | dim(top_neg)[[1]] > 0) {
  top_pos$Direction <- "Positive"
  top_neg$Direction <- "Negative"
  
  # Combine and prepare for plotting
  top_genes <- bind_rows(top_pos, top_neg)
  
  # Set sorting order: all negative genes first, then all positive genes
  top_genes <- top_genes %>%
    mutate(Gene = factor(Gene, levels = Gene[order(Direction, LogFoldChange)]))  # Order by Direction and log2FC
  
  # # Bar chart with flipped x and y axes
  p <- ggplot(top_genes, aes(x = LogFoldChange, y = Gene, fill = Direction)) +
    geom_bar(stat = "identity") +
    labs(x = "logFC", y = "Gene",title=paste0(condition,"vs. ",reference," among", cell," cells in youth with Type 2 Diabetes")) +
    scale_fill_manual(values = c("blue", "red"), labels = c("Lower Expression", "Higher Expression")) +
    theme_minimal()+
    theme(element_text(family="Times"))+
    theme(legend.position="none")
  
  if (!is.null(cell)){
    figure_file<- paste0("MAST_Barchart_",cell_name,"_cells_for_",condition,".pdf")
  } else {
    figure_file <- paste0("MAST_Barchart_for_",condition,".pdf")
  }
  
  pdf(fs::path(dir.results,figure_file),width=20,height=20)
  plot(p)
  dev.off()
  }
}


#Ideas Function ----
ideas_fxn <- function(so, exp_group, ref_group, covariates, var_type, id, method, exposure) {
  
  #Filter so to groups of interest
  so$filt_condition <- ifelse(so$group2==paste0(exp_group) | so$group2==paste0(ref_group),"Yes","No")
  so <- subset(so,filt_condition=="Yes")
  so$group3 <- ifelse(so$group2==paste0(exp_group),paste0(exp_group),paste0(ref_group))
  so$group3 <- factor(so$group3)
  
  # Create individual variable dynamically based on the input 'id'
  so$individual <- so[[id]]  # Assign the column specified by 'id' as the individual variable
  
  # Split genes into batches
  batch_size <- 1000  # Number of genes per batch
  count_matrix <- GetAssayData(so, assay = "RNA", layer = "counts")
  gene_batches <- split(rownames(count_matrix), ceiling(seq_along(rownames(count_matrix)) / batch_size))
  
  # Loop through each gene batch
  all_results <- list()  # To store results from all batches
  
  for (i in seq_along(gene_batches)) {
    cat("Processing batch", i, "of", length(gene_batches), "\n")
    
    # Get the list of genes for the current batch
    batch_genes <- gene_batches[[i]]
    
    # Extract the count matrix for this batch of genes
    batch_matrix <- count_matrix[batch_genes, , drop = FALSE]  # Drop = FALSE ensures it stays a matrix
    batch_matrix <- as.matrix(batch_matrix)
    # Round the count matrix to ensure integer values
    batch_matrix <- round(batch_matrix)  # Round values to ensure integer counts
    
    # Check if the matrix contains non-integer values
    if (any(batch_matrix != floor(batch_matrix))) {
      stop("Batch matrix contains non-integer values!")
    }
    
    # Extract the metadata for the cells corresponding to the batch
    meta_cell <- so@meta.data
    meta_cell$cell_id <- rownames(meta_cell)  # Add cell IDs if not already present
    
    # Ensure the metadata includes the required column (e.g., 'Kit_Lot' for individuals)
    meta_cell$cell_rd <- colSums(batch_matrix)  # Total counts per cell (read depth)
    
    # Add a small constant to ensure no zero or negative values
    meta_cell$cell_rd <- meta_cell$cell_rd + 1  # Adding 1 to avoid log(0) errors
    
    #    # Ensure the metadata includes the required column (e.g., 'Kit_Lot' for individuals)
    # meta_cell$cell_rd[meta_cell$cell_rd == 0] <- 1  # Ensure no zero read depth
    
    # Create individual-level metadata (e.g., group information, exposure)
    meta_ind <- unique(meta_cell[, c("individual", "group3")])  # Use dynamic 'id' for individual variable
    
    # Set variables for differential expression analysis
    var2test <- "group3"  # Group variable to test
    var2adjust <- NULL     # Adjusting variable (e.g., RIN, if applicable)
    var2test_type <- var_type  # "binary" or "continuous"
    var_per_cell <- "cell_rd"  # Cell-level adjustment variable (read depth)
    test_method <- paste0(method)
    
    # Perform IDEAS distance calculation for this batch
    dist_batch <- ideas_dist(
      count_input = batch_matrix,
      meta_cell = meta_cell,
      meta_ind = meta_ind,
      var_per_cell = var_per_cell,
      var2test = var2test,
      var2test_type = var2test_type,
      d_metric = "Was",  # Wasserstein distance
      fit_method = test_method  # Negative binomial fit
    )
    
    # Perform PERMANOVA for this batch
    pval_batch <- permanova(
      dist_array = dist_batch,
      meta_ind = meta_ind,
      var2test = var2test,
      var2adjust = NULL,
      var2test_type = var2test_type,
      n_perm = 999,
      r.seed = 903
    )
    
    # Combine results for this batch
    batch_results <- data.frame(
      gene = batch_genes,
      pvalue = pval_batch
    )
    
    # Store results
    all_results[[i]] <- batch_results
  }
  
  # Combine all batches into one final result
  final_results <- do.call(rbind, all_results)
  
  # Save or further analyze the results
  write.csv(final_results, "DEA_results_ideas.csv", row.names = FALSE)
  
  # Example mock data for log-fold changes (LFC); Replace with actual data
  set.seed(42)
  final_results$logFC <- rnorm(nrow(final_results), mean = 0, sd = 2)  # Mock log fold change values
  
  # Set thresholds for significance
  pval_threshold <- 0.05
  logFC_threshold <- 1
  
  # Classify genes as upregulated, downregulated, or not significant
  final_results <- final_results %>%
    mutate(
      significance = case_when(
        pvalue < pval_threshold & logFC > logFC_threshold ~ "Upregulated",
        pvalue < pval_threshold & logFC < -logFC_threshold ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    )
  
  # Create the volcano plot
  volcano_plot <- ggplot(final_results, aes(x = logFC, y = -log10(pvalue), color = significance)) +
    geom_point(alpha = 0.8, size = 1.5) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
    labs(
      title = "Volcano Plot",
      x = "Log2 Fold Change (logFC)",
      y = "-log10(p-value)"
    ) +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  # Save the plot
  ggsave("volcano_plot_ideas.png", volcano_plot, width = 8, height = 6)
  
  # Optionally, print the plot to the console
  print(volcano_plot)
}


# GeomSplitViolin <- ggproto(
#   "GeomSplitViolin", 
#   GeomViolin, 
#   draw_group = function(self, data, ..., draw_quantiles = NULL) {
#     data <- transform(data, 
#                       xminv = x - violinwidth * (x - xmin), 
#                       xmaxv = x + violinwidth * (xmax - x))
#     grp <- data[1,'group']
#     newdata <- plyr::arrange(
#       transform(data, x = if(grp%%2==1) xminv else xmaxv), 
#       if(grp%%2==1) y else -y
#     )
#     newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
#     newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
#     if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
#       stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
#       quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
#       aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
#       aesthetics$alpha <- rep(1, nrow(quantiles))
#       both <- cbind(quantiles, aesthetics)
#       quantile_grob <- GeomPath$draw_panel(both, ...)
#       ggplot2:::ggname("geom_split_violin", 
#                        grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
#     } else {
#       ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
#     }
#   }
# )
# 
# geom_split_violin <- function (mapping = NULL, 
#                                data = NULL, 
#                                stat = "ydensity", 
#                                position = "identity", ..., 
#                                draw_quantiles = NULL, 
#                                trim = TRUE, 
#                                scale = "area", 
#                                na.rm = FALSE, 
#                                show.legend = NA, 
#                                inherit.aes = TRUE) {
#   layer(data = data, 
#         mapping = mapping, 
#         stat = stat, 
#         geom = GeomSplitViolin, 
#         position = position, 
#         show.legend = show.legend, 
#         inherit.aes = inherit.aes, 
#         params = list(trim = trim, 
#                       scale = scale, 
#                       draw_quantiles = draw_quantiles, 
#                       na.rm = na.rm, ...)
#   )
# }
# 
# split.vp <- function(seurat_object, genes, filepath, color1 = "#6c9a8b", color2 = "#e8998d") {
#   for (i in 1:length(genes)){
#     cat("\n")
#     cat("###", genes[i])
#     cat("\n")  
#     d = VlnPlot(seurat_object, features = genes[i], split.by = "Group", idents = "PT", split.plot = F, pt.size = 0) 
#     d = d$data
#     p = ggplot(d,aes(x=ident, y = !!sym(genes[i]), fill=split))+
#       geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), 
#                   size = 0.2, alpha = 0.3, show.legend = T, aes(color = split)) +
#       geom_split_violin(trim = T) +
#       theme_bw()+
#       theme(legend.title = element_blank(),
#             axis.title.x = element_blank(), 
#             axis.text.x=element_blank(),
#             plot.title = element_text()) +
#       labs(title = (paste0(genes[i]," in  PT Cells")),
#            y = "Expression") +
#       scale_fill_manual(values=c(color1, color2)) + 
#       scale_color_manual(values=c(color1, color2))
#     print(p)
#     cat("\n")
#     # Save
#     ggsave(filename = paste0(filepath,"Violin_",genes[i],".jpeg"),plot = p,scale = 5,
#            width = 800,height = 600,units = "px")
#   }
# }
# 
# split.vp.combined <- function(seurat_object, genes, filepath, color1 = "#6c9a8b", color2 = "#e8998d", idents = NULL) {
#   compiled_d = data.frame()
#   for (i in 1:length(genes)){
#     d = VlnPlot(seurat_object, features = genes[i], split.by = "Group", idents = idents, split.plot = F, pt.size = 0) 
#     d = d$data
#     d$genename = colnames(d)[1]
#     colnames(d)[1] <- "expression"
#     compiled_d = rbind(compiled_d, d)
#   }
#   p = 
#     ggplot(compiled_d,aes(x=genename, y = expression, fill=split)) +
#     geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), 
#                 size = 0.2, alpha = 0.3, show.legend = T, aes(color = split)) +    
#     geom_split_violin(scale = "width", trim = T) +
#     theme_bw() +
#     theme(legend.title = element_blank(),
#           axis.title.x = element_blank(), 
#           plot.title = element_text()) +
#     labs(y = "Expression") +
#     scale_fill_manual(values=c(color1, color2)) + 
#     scale_color_manual(values=c(color1, color2))
#   print(p)
#   cat("\n")
#   # Save
#   ggsave(filename = paste0(filepath,"Violin_combined",".jpeg"),plot = p,scale = 5,
#          width = 1000,height = 600,units = "px")
#   
# }
# 
# # Formatted dot plot function
# # colorlow = "#8ecae6", colormid = "#fcbf49", colorhigh = "#d90429"
# dp.formatted <- function(seurat_object, genes, celltype, group.by, m,
#                          colorlow = "#83c5be", colormid = "#f4f1bb", colorhigh = "#d90429")
# {
#   pt.combined <- DotPlot(seurat_object,
#                          features = genes,idents = celltype, group.by = group.by,
#                          scale = F, cols = "RdYlBu"
#   )$data 
#   
#   pt.plot <- pt.combined %>% 
#     ggplot(aes(x=features.plot, y = id, color = avg.exp.scaled, size = pct.exp)) + 
#     geom_point() +
#     theme_bw() +
#     scale_color_gradient2(low = colorlow, mid = colormid, high = colorhigh, midpoint = 2,
#                           guide = guide_colorbar(label.vjust = 0.8, ticks = F, draw.ulim = T, draw.llim = T),
#                           limits = c(0,4)) +
#     scale_size(range = c(0,4), 
#                limits = c(1,80)) +
#     theme(panel.grid = element_blank(),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#           legend.text = element_text(size = 8),
#           legend.title = element_text(size = 8, vjust = 0.5),
#           legend.spacing.x = unit(.1, "cm"),
#           legend.direction = "horizontal") +
#     guides(size = guide_legend(label.position = "bottom",
#                                title.position = "top"),
#            color = guide_colorbar(label.position = "bottom",
#                                   title.position = "top")) +
#     labs(color = "Scaled average expression",
#          size = "Expression (%) ") + 
#     scale_y_discrete(limits=rev)
#   
#   pt.table <- m %>%
#     filter(rownames(m) %in% genes) %>%
#     dplyr::mutate(p_val_rounded = round(p_val, 4),
#                   p_val = p_format(p_val_rounded, trailing.zero = T, accuracy = 0.001, digits = 3),
#                   p_val_adj_rounded = round(p_val_adj, 4),
#                   p_val_adj = p_format(p_val_adj_rounded, trailing.zero = T, accuracy = 0.001, digits = 3),
#                   pct.1 = sprintf("%.3f", pct.1),
#                   pct.2 = sprintf("%.3f", pct.2),
#                   avg_log2FC = sprintf("%.3f", avg_log2FC)) %>% 
#     dplyr::select(pct.1, pct.2, avg_log2FC, p_val, p_val_adj) 
#   gg.pt.table <- ggtexttable(pt.table,
#                              cols = c("T1D", "HC", "Log2FC", "p-value", "q-value"),
#                              theme = ttheme("blank")) %>%
#     tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1) %>%
#     tab_add_title("% Expressed in PT Cells")
#   
#   pt.plot_table <- ggarrange(pt.plot, NULL, gg.pt.table,
#                              nrow = 1, widths = c(1,-0.1,1), common.legend = F,
#                              legend = "top")
#   
# }
# 
# metabo_sc <- function(data, gene, transcripts, gene_name) {
#   plots <- list()
#   for (i in 1:length(transcripts)) {
#     plot <- subset(data, apply(!is.na(data[, c(transcripts[i], gene)]), 1, all)) %>%
#       ggplot(aes_string(y = gene, x = transcripts[i], color = "group")) +
#       geom_smooth(method = lm, se = F, aes(color=NULL), color = "black",
#                   linetype = "dashed") +
#       geom_point() +
#       geom_smooth(method = lm, se = F) +
#       labs(color = "Group",
#            y = gene_name,
#            x = if(is.null(label(data[[transcripts[i]]]))) (transcripts[i]) else label(data[[transcripts[i]]])) +
#       theme_bw()
#     plots[[i]] <- plot
#   }
#   ggarrange(plotlist = plots, common.legend = T)
# }
# 
# add_direction <- function(df) {
#   df <- df %>%
#     mutate(direction = case_when(
#       (avg_log2FC > 0 & p_val < 0.05) ~ "Upregulated", 
#       (avg_log2FC <= 0 & p_val < 0.05) ~ "Downregulated",
#       TRUE ~ "NS"
#     ))
#   df$direction <- factor(df$direction, levels = c("Downregulated", "Upregulated", "NS") )
#   return(df)
# }
# 
