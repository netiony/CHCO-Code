library(ggrepel)
library(qpcR)
library(ggpubr)

test_ipa <- readxl::read_xls("/Users/choiyej/Dropbox/PANTHER/IPA results/tanner_stage_12_345_rh_de_table_export.xls", skip = 1)
colnames(test_ipa) <- c("pathway", "neglog_p", "ratio", "zscore", "molecules")



# Volcano plot
sig_cutoff = -log(0.05, base = 10)

test_ipa <- test_ipa %>%
  dplyr::mutate(direction = case_when(zscore < 0 & neglog_p > sig_cutoff ~ "Negative", 
                                      zscore > 0 & neglog_p > sig_cutoff ~ "Positive",
                                      zscore == 0 & neglog_p > sig_cutoff ~ "Unknown",
                                      T ~ "NS")) %>%
  group_by(direction) %>%
  dplyr::mutate(count = row_number(),
                label = case_when(direction != "NS" ~ paste0(str_sub(direction, 1,1), count))) %>%
  ungroup()

test_ipa$direction <- factor(test_ipa$direction, levels = c("Negative", "Positive", "Unknown", "NS"))

ipa.plot <- test_ipa %>% filter(neglog_p > 0) %>%
  ggplot(aes(x = zscore, y = neglog_p)) + 
  geom_point(aes(color = direction), alpha = 0.5, size = 2) +
  geom_hline(aes(yintercept = sig_cutoff), color = "red", linetype = "dashed") +
  labs(x = "Z-Score",
       y = "-log(p-value)",
       color = "Direction") + 
  scale_color_manual(values = c("#bf0603", "#457b9d",
                                "#f6bd60", "#dad7cd")) +
  ggrepel::geom_label_repel(aes(label = label, color = direction),
                            label.size = 0.15,
                            max.overlaps = 100, 
                            min.segment.length = 0.01) +
  theme_bw()

# Table
ipa.table.neg <- test_ipa %>%
  dplyr::filter(direction == "Negative") %>%
  dplyr::select(label, pathway) 
ipa.table.pos <- test_ipa %>%
  dplyr::filter(direction == "Positive") %>%
  dplyr::select(label, pathway) 
ipa.table.unk <- test_ipa %>%
  dplyr::filter(direction == "Unknown") %>%
  dplyr::select(label, pathway) 

cbind.fill <- function(...) {
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function(x) {
    filled_matrix <- ifelse(is.na(x), "", x)
    rbind(filled_matrix, matrix("", n - nrow(filled_matrix), ncol(filled_matrix)))
  }))
}

ipa.table.comb <- cbind.fill(ipa.table.neg,ipa.table.pos,ipa.table.unk)
gg.ipa.table <- ggtexttable(ipa.table.comb, rows = NULL,
                                cols = c("Label", "Pathway", "Label", "Pathway","Label", "Pathway"),
                                theme = ttheme("blank",
                                               tbody.style = tbody_style(fill = "white", 
                                                                         size = 8.5, 
                                                                         hjust = 0,
                                                                         x = 0.1))) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 0.8)



ggarrange(ipa.plot, gg.ipa.table,
          ncol = 2, heights = c(1,1), widths = c(0.5,1),
          legend = "top")

