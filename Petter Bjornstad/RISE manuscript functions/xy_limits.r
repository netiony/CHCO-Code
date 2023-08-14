# data - list containing the followup data for each
#        timepoint
# 
# groups - character vector of treatment groups to include, NULL if None included.
#          these names must match those used in the input data list (data)
# 
# timepoints - character vector of timepoints to include, NULL if None included.
#              these names must match those used in the input data list (data)
# 
xy_limits <- function(data, useBaseline=T){
  # Split the data in baseline and followup visits since the data structure
  # is different at baseline
    base.data <- data[[1]]
    followup.data <- data[2:length(data)]
    
  # Determine the number of timepoints for full data
    ntimepoints <- length(followup.data)
    
  # Determine the number of treatment groups for full data.  Since the number 
  # of treatment groups is the same for all followup visits we only need to take
  # the first element of the returned list
    ngroups <- lapply(followup.data, length)[[1]]
    
  # Loop through and concatentate the baseline curve data and the confidence
  # ellipse data
    if(useBaseline) xylimit.data <- as.matrix(base.data$curveXY)
    for(i in 1:ntimepoints){
      for(j in 1:ngroups){
        if(useBaseline) xylimit.data <- rbind(xylimit.data, followup.data[[c(i,j,2)]])
        if(!useBaseline & i==1 & j==1) xylimit.data <- followup.data[[c(i,j,2)]]
        if(!useBaseline & i>1 | j>1) xylimit.data <- rbind(xylimit.data, followup.data[[c(i,j,2)]])
      }
    }
    
    # Determine the min and max x and y based on the concatenated data
    x.min = min(xylimit.data[,1],na.rm=T)
    x.min = x.min - 0.15*x.min
    x.max = max(xylimit.data[,1],na.rm=T)
    x.max = x.max + 0.15*x.max
    x.inc <- round(0.1*(x.max-x.min),2)

    
    y.min = min(xylimit.data[,2],na.rm=T)
    y.min = y.min - 0.15*y.min
    y.max = max(xylimit.data[,2],na.rm=T)
    y.max = y.max + 0.15*y.max
    y.inc <- round(0.1*(y.max-y.min),2)

    return(c(x.min,x.max,x.inc,y.min,y.max,y.inc))
}