library(tidyverse)
library(ggrepel)
library(readxl)
library(clusterProfiler)
library(ReactomePA)
# Import
df = read_excel("/Volumes/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted.xlsx",sheet = "GLYC CPH")
df = df %>% arrange(p.value)
pCutoff = 0.05
# Volcano plot
df$logp <- -log10(df$p.value)
df$sig <- df$adj.p.value <= pCutoff
vp = ggplot(df, aes(x = estimate, y = logp, color = sig)) +
  geom_point(size = 2) +
  geom_label_repel(
    data = df[df$sig, ], aes(label = TargetFullName), color = "black",
    max.overlaps = 10
  ) +
  scale_color_manual(values = c("grey", "#3e6dbf")) +
  xlab("HR") +
  ylab(bquote(~ -Log[10] ~ italic(P))) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("~/Documents/volcano.png",width = 1000,height = 900,units = "px",
       plot = vp,scale = 2.25)
# Forest plot
df %>% filter(adj.p.value <= pCutoff) %>%
  forestplot(labeltext = TargetFullName,
             mean = estimate,
             lower = conf.low,
             upper = conf.high,
             zero = 1,
             xlab = "HR")
# Pathway analysis
gl = log(df$estimate[df$p.value <= 0.05])
names(gl) = df$EntrezGeneID[df$p.value <= 0.05]
gl = sort(gl,decreasing = T)
gl = gl[!duplicated(names(gl))]
gsea = gsePathway(gl,organism = "human")
dp = enrichplot::dotplot(gsea)
ggsave("~/Documents/dot.png",width = 1000,height = 900,units = "px",
       plot = dp,scale = 1.8)
