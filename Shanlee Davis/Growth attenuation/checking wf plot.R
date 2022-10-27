wfdata <- read.csv("E:/Davis/Growth attenuation analysis/wfdata.csv")


barplot(wfdata$delta_pah,
        col = wfdata$stop_reason2,
        space = 0.5, ylim=c(-50,10),
        main = "Waterfall Plot for Change in Predicted Adult Height",
        ylab = "Change from baseline (cm)",
        cex.axis=1.2, cex.lab=1.2, legend.text = c("Discontinued Early", "Completed Treatment", "On Treatment"))

full_color_list <- palettes_c_names

library(paletteer)
library(ggthemes)
colors <- paletteer_c(palette = "grDevices::Blues", n = 3)
barplot(wfdata$delta_pah,
        col =  colors[ unclass(wfdata$stop_reason2) ]   ,
        space = 0.5, ylim=c(-50,10),
        main = "Waterfall Plot for Change in Predicted Adult Height",
        ylab = "Change from baseline (cm)",
        cex.axis=1.2, cex.lab=1.2, legend.text = c("Discontinued Early", "Completed Treatment", "On Treatment"))


# Scatterplot with categoric color scale
plot(
  x = iris$Petal.Length, 
  y = iris$Petal.Width,
  bg = colors[ unclass(iris$Species) ],
  cex = 3,
  pch=21
)






library(ggplot2)
ggplot(wfdata, aes(y=delta_pah, fill=stop_reason2)) +
  geom_bar(stat="identity")+
  ggtitle("name")


p <- ggplot(wfdata, aes(delta_pah, fill = type)) + geom_rect(aes(x = desc,
                                                             
                                                             xmin = id - 0.45, xmax = id + 0.45, ymin = end,
                                                             
                                                             ymax = start))