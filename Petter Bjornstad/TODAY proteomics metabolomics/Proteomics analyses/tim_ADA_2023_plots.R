library(tidyverse)
library(readxl)
library(forestplot)
# Import
df = read_excel("/Volumes/Peds Endo/Petter Bjornstad/TODAY subaward/Results/TODAY proteomics Cox models scaled adjusted 2022-08-31.xlsx",sheet = "GLYC CPH")
pCutoff = 0.001
# Volcano plot
df$logp <- -log10(df$p.value)
df$sig <- df$p.value <= pCutoff
vp = ggplot(df, aes(x = estimate, y = logp, color = sig)) +
  geom_hline(yintercept = -log10(pCutoff), linetype = "dashed") +
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
ggsave("~/Documents/volcano.png",width = 1200,height = 900,units = "px",
       plot = vp,scale = 2.5)
# Forest plot
df %>% filter(row_number() <= 10) %>%
  forestplot(labeltext = TargetFullName,
             mean = estimate,
             lower = conf.low,
             upper = conf.high,
             zero = 1,
             cex  = 2,
             lineheight = "auto",
             xlab = "Lab axis txt")

