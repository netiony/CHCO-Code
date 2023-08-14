plotTrajectory <- function(data, 
                           xvar, 
                           yvar, 
                           grpvar,
                           summary.func,
                           errorbar.func, 
                           length.endcaps,
                           line.colors, 
                           shape.type, 
                           xtitle, 
                           ytitle, 
                           xaxis.breaks=waiver(),
                           legend.position='topright', 
                           legend.text.size=12,
                           legend.title = 'Group',
                           xtitle.size=12, 
                           ytitle.size=12,
                           outputData=FALSE){
  # Load ggplot library
  library(ggplot2)

 
  # Create Plot
  plot.out <- ggplot(data, aes(x=eval(parse(text=paste0('data$',xvar))), y=eval(parse(text=paste0('data$',yvar))), col=eval(parse(text=paste0('data$',grpvar))), shape=eval(parse(text=paste0('data$',grpvar))))) + 
    stat_summary(fun.y = summary.func, geom = "point", size=2, position = position_dodge(0.25)) + 
    stat_summary(fun.data = errorbar.func, geom = "errorbar", width=length.endcaps, linetype='solid', position = position_dodge(0.25)) + 
    stat_summary(fun.y = summary.func, geom = "line", size=1.1, position = position_dodge(0.25)) + 
    scale_x_continuous(breaks = xaxis.breaks) +
    # coord_trans(x = xtrans, y = ytrans) +
    scale_color_manual(name="",labels=levels(data[,grpvar]), values=line.colors) +
    scale_shape_manual(name="",labels=levels(data[,grpvar]), values = shape.type) +
    theme_bw() +
    labs(x=xtitle, y=ytitle, col=legend.title) +
    theme(legend.position = legend.position,
          legend.background = element_rect(color = "white",
                                           fill = "white",
                                           size = 1),
          legend.text = element_text(size = legend.text.size),
          axis.title.x = element_text(size = xtitle.size),
          axis.title.y = element_text(size = ytitle.size),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
          )
  
  # Return plot
  return(plot.out)
}