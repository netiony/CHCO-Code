library(lipidr)
library(ggplot2)

plot_results_volcano_unadj <- function (de.results, show.labels = TRUE) {
  if (!"logFC" %in% colnames(de.results)) {
    message("de.results contains ANOVA-style comparison.", 
            " Average Expression will be plotted instead of logFC.")
    p <- ggplot(de.results, aes(AveExpr, -log10(P.Value), 
                                color = Class, label = Molecule)) + geom_point() + 
      xlab("Average Intensity")
  } else {
    p <- ggplot(de.results, aes(logFC, -log10(P.Value), color = Class, 
                                label = Molecule)) + geom_point() + geom_vline(xintercept = c(1, -1), lty = 2) + theme_bw() 
  }
  p <- p + geom_hline(yintercept = -log10(0.05), lty = 2) + 
    facet_wrap(~contrast)
  
  if (show.labels) {
    sig <- de.results$P.Value < 0.05
    if ("logFC" %in% colnames(de.results)) {
      sig <- sig & abs(de.results$logFC) >= 1
    }
    labels <- de.results[sig, "Molecule"]
        p <- p + geom_label_repel(data = de.results[sig, ], 
                        aes(label = Molecule, 
                            vjust = "inward",  # Center vertically
                            hjust = "inward",  # Center horizontally
                            color = Class), 
                        size = 2.4, 
                        direction = "both", segment.alpha = 0.6, label.padding = 0.15, 
                        force = 0.5, max.overlaps = 1000, show.legend = FALSE)
  }
  
  print(p)
}




significant_molecules_unadj <- function (de.results, p.cutoff = 0.05, logFC.cutoff = 1) 
{
  if (!"logFC" %in% colnames(de.results)) {
    message("de.results contains ANOVA-style comparison.", 
            " LogFC cutoff will be ignored")
    ret <- de.results %>% filter(P.Value < p.cutoff) %>% 
      (function(x) split(x$Molecule, x$contrast))
    return(ret)
  }
  de.results %>% filter(P.Value < p.cutoff, abs(logFC) > logFC.cutoff) %>% 
    (function(x) split(x$Molecule, x$contrast))
}
