library(readxl)
library(ggplot2)
library(dplyr)
library(forcats)
library(ggpubr)

# make it a function
ipa_plot <- function(data){
  p <- ggplot(data, aes(reorder(`Ingenuity Canonical Pathways`, `-log(p-value)`), `-log(p-value)`))+
    geom_col() +
    scale_y_continuous(limits = c(0, 6)) +
    ylab("-log(p-value)") +
    xlab("Pathway")
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

