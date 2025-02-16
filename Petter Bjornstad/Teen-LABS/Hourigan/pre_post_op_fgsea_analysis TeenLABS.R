# Load necessary libraries
library(tidyverse)  # Data manipulation and visualization
library(edgeR)      # Differential expression analysis
library(limma)      # Linear modeling
library(fgsea)      # Gene set enrichment analysis
library(ggplot2)    # Visualization
library(tidyr)      # Data wrangling
library(gridExtra)  # Combine plots

# Set random seed for reproducibility
set.seed(42)

# Load data and filter for Group "BB"
a <- read_tsv("../Step1_Matrix_28-Mar-2023_10_53_54.txt") %>%
  dplyr::filter(Group %in% c("BB"))

# Prepare sample group metadata
b <- a %>%
  select(1:4) %>%  # Select the first four columns
  distinct(SampleID, .keep_all = TRUE) %>%
  mutate(group = paste0(TimePoint, "_", Group, "_", Subject)) %>%
  arrange(group)

# Replace dashes in TimePoint column
b$TimePoint <- gsub("-", "_", b$TimePoint)

# Create design matrix for model
design <- model.matrix(~ 0 + TimePoint, data = b)

# Prepare counts matrix
counts <- a[5:length(a)]
rownames(counts) <- paste0(a$TimePoint, "_", a$Group, "_", a$Subject)

# Transform counts to long format
count2 <- counts %>%
  rownames_to_column() %>%
  gather(var, value, -rowname) %>%
  spread(rowname, value)

# Prepare numeric log-transformed counts
rownames(count2) <- count2$var
count3 <- count2[2:length(count2)] %>%
  mutate_if(is.character, as.numeric) %>%
  log()
rownames(count3) <- count2$var

# Compute duplicate correlation
corfit <- duplicateCorrelation(count3, design, block = b$Subject)

# Define contrasts for differential analysis
cont.matrix <- makeContrasts(
  delta1 = (TimePointPost_Op - TimePointPre_Op),
  levels = design
)

# Fit linear model and apply contrasts
fit <- lmFit(count3, design = design, block = b$Subject, correlation = corfit$consensus)
contrast_fit <- contrasts.fit(fit, contrasts = cont.matrix)

# Apply empirical Bayes smoothing
ebays_fit <- eBayes(contrast_fit, robust = TRUE)

# Extract top table results
res <- as.data.frame(topTable(ebays_fit, coef = 1, number = Inf))

# Load pathways for GSEA
tm <- gmtPathways("./ReactomePathways.gmt")

# Prepare ranked gene list for fgsea
ch <- res %>% arrange(t)
res2 <- cbind(rownames(ch), as.numeric(as.character(ch$t)))
res2 <- as.tibble(res2)
res2$V2 <- as.numeric(res2$V2)
ranks <- deframe(as.data.frame(res2))

# Perform fgsea analysis
fgseaRes <- fgsea(pathways = tm, stats = ranks, maxSize = 500)

# Filter significant pathways (padj < 0.05)
p <- fgseaRes %>%
  arrange(padj) %>%
  filter(padj < 0.05)

# Create dot plot for significant pathways
p1 <- ggplot(p %>%
               mutate(size = case_when(
                 padj > 0.01 ~ "0.01 - 0.05",
                 padj <= 0.01 ~ "<0.01"
               )),
             aes(y = reorder(pathway, desc(padj)), x = "", size = size, fill = NES)) +
  geom_point(shape = 21, alpha = 0.5) +
  labs(size = "padj", fill = "NES") +
  scale_fill_viridis_c() +
  scale_size_manual(values = c(12, 6)) +
  xlab("Significant Enrichments for Pre vs Post-op") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    strip.text = element_text(size = 20),
    strip.background = element_rect(colour = "black", fill = NA),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "right",
    axis.text.x = element_text(color = "black", vjust = 0.5),
    text = element_text(size = 14)
  )

# Define pathways of interest
pathways_interested <- c("IRAK2 mediated activation of TAK1 complex")

# Initialize empty list for plots
plots <- list()

# Create line plots for each pathway of interest
for (path in pathways_interested) {
  filtered_pathway <- p %>%
    filter(pathway == path) %>%
    select(leadingEdge, pathway)
  
  # Extract leading edge genes
  leading_genes <- unlist(filtered_pathway$leadingEdge)
  
  # Subset matrix and reshape for plotting
  subset_matrix <- a[, c("Subject", "TimePoint", leading_genes)]
  subset_long <- subset_matrix %>%
    pivot_longer(cols = -c(Subject, TimePoint), names_to = "Gene", values_to = "RFU")
  subset_long$TimePoint <- factor(subset_long$TimePoint, levels = c("Pre-Op", "Post-Op"))
  
  # Plot data
  plot <- ggplot(subset_long, aes(x = TimePoint, y = RFU, group = Subject, color = Subject)) +
    geom_line() +
    geom_point() +
    facet_wrap(~Gene) +
    theme(text = element_text(size = 14)) +
    labs(title = "Pathways: IRAK2 mediated activation of TAK1 complex")
  
  # Add plot to list
  plots[[path]] <- plot
}

# Combine plots and save as PDF
p1_list <- list(p1)
combined_plots <- c(p1_list, plots)


pdf("pathway_line.pdf", width = 16)
grid.arrange(grobs = combined_plots, ncol = 2)
dev.off()
