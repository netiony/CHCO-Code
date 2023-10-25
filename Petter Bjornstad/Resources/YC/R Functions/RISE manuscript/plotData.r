# data - list containing the data for the baseline curve and the data for each
#        timepoint
# 
# groups - Can be used to select treatment groups to plot (e.g. Met or GM)
# 
# legend.names - Specify group names (excluding baseline) for treatment group
#                names in legend
#                
# timepoints - character vector of timepoints to include, NULL if None included.
#              these names must match those used in the input data list (data)
# 
# useBaseline - logical to indicate whether to use baseline data to select x and
#               y axis limits (default is true)
#   
# xpd         - TRUE or FALSE (default FALSE).  If true allows legend to be 
#               printed outside the plot area
#               
# plotcol     - background color of plot area.  Default is white.  Should be 
#               a quoted name for hexcode (i.e. 'white', 'blue', '#f3eff5' etc)
#               
# plotannotcol - Color of axes, plot border, tick marks, axis titles.  
#                Default is black.  Should be a quoted name for hexcode (i.e. 
#                'white', 'blue', '#f3eff5' etc)
plotData <- function(data, timepoints=NULL, timepoint.labels=NULL, groups=NULL, 
                     xaxis.label='X', yaxis.label='Y', xmin=NULL, xmax=NULL, 
                     ymin=NULL, ymax=NULL, useBaseline=T, legend.names=NA, 
                     leg.pos="bottomright", xpd=FALSE, plotbgcol='white',
                     plotannotcol='black'){
# ------------------------------------------------------------------------------
# Transparency parameters for ellipses, points, point text
# ------------------------------------------------------------------------------
# All treatment groups
  # transparency.ellipses <- list(GM=0.15, LM=0.15, Met=0.15, Pla=0.15)
  # transparency.points   <- list(GM=1.00, LM=1.00, Met=1.00, Pla=1.00)
  # transparency.pointtext <- c(GM=1.00, LM=1.00, Met=1.00, Pla=1.00)
  # transparency.arrows    <- c(GM=1.00, LM=1.00, Met=1.00, Pla=1.00)

# GM  
  # transparency.ellipses <- list(GM=0.15, LM=0.00, Met=0.00, Pla=0.00)
  # transparency.points   <- list(GM=1.00, LM=0.00, Met=0.00, Pla=0.00)
  # transparency.pointtext <- c(GM=1.00, LM=0.00, Met=0.00, Pla=0.00)
  # transparency.arrows    <- c(GM=1.00, LM=0.00, Met=0.00, Pla=0.00)

# LM  
  # transparency.ellipses <- list(GM=0.00, LM=0.15, Met=0.00, Pla=0.0)
  # transparency.points   <- list(GM=0.00, LM=1.00, Met=0.00, Pla=0.00)
  # transparency.pointtext <- c(GM=0.00, LM=1.00, Met=0.00, Pla=0.00)
  # transparency.arrows    <- c(GM=0.00, LM=1.00, Met=0.00, Pla=0.00)

# Met  
  transparency.ellipses <- list(GM=0.00, LM=0.00, Met=0.15, Pla=0.00)
  transparency.points   <- list(GM=0.00, LM=0.00, Met=1.00, Pla=0.00)
  transparency.pointtext <- c(GM=0.00, LM=0.00, Met=1.00, Pla=0.00)
  transparency.arrows    <- c(GM=0.00, LM=0.00, Met=1.00, Pla=0.00)

# Pla
  # transparency.ellipses <- list(GM=0.00, LM=0.00, Met=0.00, Pla=0.15)
  # transparency.points   <- list(GM=0.00, LM=0.00, Met=0.00, Pla=1.00)
  # transparency.pointtext <- c(GM=0.00, LM=0.00, Met=0.00, Pla=1.00)
  # transparency.arrows    <- c(GM=0.00, LM=0.00, Met=0.00, Pla=1.00)
  
# Baseline curve and point only
  # transparency.ellipses <- list(GM=0.00, LM=0.00, Met=0.00, Pla=0.00)
  # transparency.points   <- list(GM=0.00, LM=0.00, Met=0.00, Pla=0.00)
  # transparency.pointtext <- c(GM=0.00, LM=0.00, Met=0.00, Pla=0.00)
  # transparency.arrows    <- c(GM=0.00, LM=0.00, Met=0.00, Pla=0.00)

# ------------------------------------------------------------------------------
# Plotting Parameters
# ------------------------------------------------------------------------------

  col.curve.base             <- 'black'
  col.ellipse.base           <- adjustcolor('black', alpha.f = 0.0)
  
  col.ellipses <- list(GM=adjustcolor('darkgreen', alpha.f=transparency.ellipses$GM),
                       LM=adjustcolor('darkorchid4', alpha.f=transparency.ellipses$LM),
                       Met=adjustcolor('saddlebrown', alpha.f=transparency.ellipses$Met),
                       Pla=adjustcolor('royalblue4', alpha.f=transparency.ellipses$Pla))
  
  col.points <-   list(GM=adjustcolor('darkgreen', alpha.f = transparency.points$GM),
                       LM=adjustcolor('darkorchid4', alpha.f = transparency.points$LM),
                       Met=adjustcolor('saddlebrown', alpha.f = transparency.points$Met),
                       Pla=adjustcolor('royalblue4', alpha.f = transparency.points$Pla))

  col.arrows      <- gray(0.1)              
  col.arrow.heads <- gray(0.1, alpha = 0.00)  
  
  cex.basepoint  <- 4       # Controls character expansion for size of baseline plotting point
  cex.timepoints <- 4       # Controls character expansion for size of plotting points for followup
  cex.point.text <- 1     # Controls character expansion for size of text in plotting points
  cex.axis.annotation <- 1.5  # Controls character expansion for size of axis annotation
  cex.labels     <- 1.7     # Controls character expansion for x and y labels
  cex.legend     <- 1      # Controls character expansion for legend labels

  
  font.xylabels <- 2
  font.title    <- 2
  
  lwd.baseline   <- 2
  lwd.arrows     <- 4
  
  pch.points     <- 15    # Plotting symbol for time points
  
  title.label <- ''

# ------------------------------------------------------------------------------
# Data Manipulation in preparation for plotting
# ------------------------------------------------------------------------------
  # Baseline Data
    baseline.data <- data[[1]]
    
  # Follow-up Data
    followup.data <- data[2:length(data)]
    
  # Subset the followup data to keep only specified timepoints
  #     if(!is.null(timepoints)) followup.data <- followup.data[names(followup.data) %in% timepoints]

  # Subset the followup data to include only selected treatment groups
  #  if(!is.null(groups)) followup.data <- lapply(followup.data, function(x){return(x[names(x) %in% groups])}) 

  # If either timepoints or groups is NULL set the followup data to NULL since
  # the indicates that only the baseline curve should be plotted
    if(is.null(timepoints) | is.null(groups)) followup.data <- NULL
    
  # Get timepoint names and treatment group names in the order in which they 
  # appear in the nested list.  These will be used below to properly identify
  # timepoints and treatment groups in the plot
    if(!is.null(followup.data)){
    timepoints.ordered <- names(followup.data)
    groups.ordered     <- lapply(followup.data, function(x){return(names(x))})[[1]]
    } else{
    timepoints.ordered <- NULL
    groups.ordered     <- NULL
    }

  # Determine the x and y plot limits from the full data.  This will force all
  # subplots with various combinations of visit and timepoint to be plotted with
  # identically scaled x and y axes
    # Calculate x and y limits based on the data.  This will only be used if
    # xmin, xmax, ymin, ymax are not passed to the routine
      xy.limits <- xy_limits(data,useBaseline = useBaseline )
      
      xmin.save <- ifelse(is.null(xmin) | is.null(xmax), xy.limits[1], xmin)
      xmax.save <- ifelse(is.null(xmin) | is.null(xmax), xy.limits[2], xmax)
      xinc <- ifelse(is.null(xmin) | is.null(xmax), xy.limits[3], abs(xmax-xmin)/20)

      ymin.save <- ifelse(is.null(ymin) | is.null(ymax), xy.limits[4], ymin)
      ymax.save <- ifelse(is.null(ymin) | is.null(ymax), xy.limits[5], ymax)
      yinc <- ifelse(is.null(ymin) | is.null(ymax), xy.limits[6], abs(ymax-ymin)/20)

      xmin <- xmin.save
      xmax <- xmax.save
      ymin <- ymin.save
      ymax <- ymax.save

# ------------------------------------------------------------------------------
# Plotting the data
# ------------------------------------------------------------------------------
  # Temporary code to make y-label print on two lines.  Couple this par statement
  # with \n in the y-label text to indicate where to break the label across
  # lines
  par(mar=c(5,5,4,2) + 0.1)
  # Initialize plot without data
  plot(baseline.data$curveXY, 
       type='n', 
       xlim = c(xmin,xmax), 
       ylim = c(ymin,ymax),
       main=title.label, 
       font.main=font.title, 
       axes=F, 
       xlab=xaxis.label, 
       ylab=yaxis.label,
       col.lab=plotannotcol,
       col.axis=plotannotcol,
       font.lab=font.xylabels,
       cex.lab=cex.labels,
       mgp=c(2.8,0,0))
  
  # The coordinates of the plot area
    u <- par("usr")
  # Plot background of plot area in desired color
    rect(u[1], u[3], u[2], u[4], col=plotbgcol, border=NA)
  # Overlay plot command to add rest of graph over the current background
    par(new=TRUE)    
    
  # Add axes based on min/max values for x and Y with the required label settings
   axis(1, 
        at=pretty(seq(xmin,xmax,xinc)), 
        cex.axis=cex.axis.annotation,
        cex.lab=cex.labels, 
        tcl=-0.5, 
        xlab=xaxis.label, 
        font.lab=font.xylabels,
        cex.axis=cex.axis.annotation, 
        cex.lab=cex.labels,
        col.axis = plotannotcol,
        col = plotannotcol)
   
   axis(2, 
        at=pretty(seq(ymin,ymax,yinc)), 
        cex.axis=cex.axis.annotation, 
        cex.lab=cex.labels, 
        tcl=-0.5, 
        ylab=yaxis.label, 
        font.lab=font.xylabels,
        cex.axis=cex.axis.annotation, 
        cex.lab=cex.labels,
        col.axis = plotannotcol,
        las=2)
   
  # Add border around the plot
  box(col=plotannotcol)

  if(!is.null(followup.data)){
    # Plot the confidence ellipses, points and arrows between points
      # Plot confidence ellipse for baseline point
        polygon(baseline.data$ellipse,col=col.ellipse.base, border=NA)

    for(i in 1:length(timepoints.ordered)){
      for(j in 1:length(groups.ordered)){
        polygon(followup.data[[c(i,j)]]$ellipse,col=col.ellipses[[match(groups.ordered[j],names(col.ellipses))]], border=NA)
        
        if(i==1) Arrows(baseline.data$meanvec[1], baseline.data$meanvec[2], 
                        followup.data[[c(i,j)]]$meanvec[1],followup.data[[c(i,j)]]$meanvec[2],
           arr.length = 0.2, arr.width = 0.1, arr.adj = 1, arr.type = 'simple', lty='dashed',
           col=adjustcolor(col.arrows, alpha.f = transparency.arrows[groups.ordered[j]]), arr.col = col.arrow.heads, arr.lwd = lwd.arrows)
        if(i > 1) Arrows(followup.data[[c(i-1,j)]]$meanvec[1], followup.data[[c(i-1,j)]]$meanvec[2],
                         followup.data[[c(i,j)]]$meanvec[1],followup.data[[c(i,j)]]$meanvec[2],
           arr.length = 0.2, arr.width = 0.1, arr.adj = 1, arr.type = 'simple', lty='dashed',
           col=adjustcolor(col.arrows, alpha.f = transparency.arrows[groups.ordered[j]]), arr.col = col.arrow.heads, arr.lwd = lwd.arrows)
      }
    }

      # Plot the baseline curve
  lines(baseline.curve, col=col.curve.base, lwd=lwd.baseline)
  # Plot Baseline Point
    points(x=baseline.data$meanvec[1], y=baseline.data$meanvec[2], pch=pch.points, col=col.curve.base, cex=cex.basepoint)
      text(baseline.data$meanvec[1], y=baseline.data$meanvec[2],
           '0',
           col='white', 
           cex=cex.point.text)

    # Plot the timepoints for each treatment group
    for(i in 1:length(timepoints.ordered)){
      for(j in 1:length(groups.ordered)){
         points(x=followup.data[[c(i,j)]]$meanvec[1], 
                y=followup.data[[c(i,j)]]$meanvec[2], 
                pch=pch.points, 
                col=col.points[[match(groups.ordered[j],names(col.points))]], 
                cex=cex.timepoints)  
      }
    }
    
    # Add text to timepoints
    for(i in 1:length(timepoints.ordered)){
      for(j in 1:length(groups.ordered)){
        text(followup.data[[c(i,j)]]$meanvec[1], 
             followup.data[[c(i,j)]]$meanvec[2],
             timepoint.labels[i],
             col=adjustcolor('white',alpha.f = transparency.pointtext[groups.ordered[j]]), 
             cex=cex.point.text)
      }
    }
  }

  # Add legend to plot
  if(is.null(followup.data)){
    legend(x="top", horiz=TRUE, inset=-0.3,legend =  c('Baseline'),pch=c(pch.points), col=c(col.curve.base), bty='n', cex=cex.legend, xjust=1, xpd = TRUE)

  } else
    legend(x="top", horiz=TRUE, inset=-0.3, legend=c('Baseline',ifelse(is.na(legend.names),groups.ordered,legend.names)),pch=c(pch.points), col=unlist(c(col.curve.base,col.points[match(groups,names(col.points))])),bty='n', cex=cex.legend, xjust=1, xpd=TRUE)

}