library(readxl)
library(ggplot2)
library(dplyr)
library(forcats)
library(ggpubr)

color_table <- tibble(
  up_down = c(-1, 0, 1),
  color = c("indianred", "gray", "steelblue")
)

ipa_plot <- function(data){
  data <- data %>% mutate(up_down = case_when(
    abs(`z-score`) < 0.0001 | is.na(`z-score`) ~ 0,
    `z-score` <= -0.0001 ~ -1,
    `z-score` >= 0.0001 ~ 1,
    .default = 0
  ))
  data$up_down <- factor(data$up_down, levels = unique(data$up_down))
  p <- ggplot(data, aes(reorder(`Ingenuity Canonical Pathways`, `-log(p-value)`), `-log(p-value)`, fill = up_down))+
    geom_col() +
    scale_y_continuous(limits = c(0, 6)) +
    ylab("-log(p-value)") +
    xlab("Pathway") + theme(legend.position="none") +
    scale_fill_manual(values = c("-1" = "steelblue", "0" = "grey", "1" = "indianred")) 
  p <- p + coord_flip()
}

alb <- read_xls("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/IPA/Output from IPA/albuminuria.xls",
                skip = 1)
alb <- alb %>% arrange(desc("-log(p-value)")) 
alb_keep <- alb[1:10,]

hyp <- read_xls("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/IPA/Output from IPA/HYP.xls",
                skip = 1)
hyp <- hyp %>% arrange(desc("-log(p-value)")) 
hyp_keep <- hyp[1:10,]

rapid <- read_xls("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/IPA/Output from IPA/rapid.xls",
                  skip = 1)
rapid <- rapid %>% arrange(desc("-log(p-value)")) 
rapid_keep <- rapid[1:10,]

alb_plot <- ipa_plot(alb_keep)
hyp_plot <- ipa_plot(hyp_keep)
rapid_plot <- ipa_plot(rapid_keep)

p <- ggarrange(alb_plot, hyp_plot, rapid_plot, ncol = 1, nrow = 3, align = "hv", labels = c("A", "B", "C"))

png('/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/IPA/Output from IPA/TODAY_DKD_IPA_pathway.png', 
    res = 600, width = 15, height = 6, units = "in")
p
dev.off()

htn <- read_xls("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/IPA/Output from IPA/HTN.xls",
                  skip = 1)
htn <- htn %>% arrange(desc("-log(p-value)")) 
#htn_keep <- htn[htn$`-log(p-value)` > 1.3,]
htn_keep <- htn[1:20,]
htn_plot <- ipa_plot(htn_keep)
png('/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/IPA/Output from IPA/TODAY_HTN_IPA_pathway.png', 
    res = 600, width = 15, height = 6, units = "in")
htn_plot
dev.off()

glyc <- read_xls("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/IPA/Output from IPA/GLYC.xls",
                skip = 1)
glyc <- glyc %>% arrange(desc("-log(p-value)")) 
glyc_keep <- glyc[1:25,]
glyc_plot <- ipa_plot(glyc_keep)
png('/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/IPA/Output from IPA/TODAY_GLYC_IPA_pathway.png', 
    res = 600, width = 15, height = 6, units = "in")
glyc_plot
dev.off()
file.copy("/Volumes/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/IPA/Output from IPA/TODAY_GLYC_IPA_pathway.png",
          "/Users/pylell/Dropbox/TODAY glycemic manuscript [shared]/Analysis output/TODAY_GLYC_IPA_pathway.png",overwrite = TRUE)