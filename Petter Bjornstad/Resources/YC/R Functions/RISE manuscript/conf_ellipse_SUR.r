# SURmodel  - SUR model for a bivariate outcome
# Cmatx     - Contrast matrix for linear combo of interest
conf_Ellipse_SUR <- function(SURmodel, Cmatx, alpha=0.95){
  library(ellipse)
  mean <- Cmatx %*% coefficients(SURmodel)
  varcov     <- Cmatx %*% summary(SURmodel)$coefCov %*% t(Cmatx)
  confellipse <- ellipse::ellipse(varcov, centre = as.vector(mean))
  returnList <- list(mean,varcov,confellipse)
  names(returnList) <- c('Mean', 'VarCov', 'Ellipse')
  
  return(returnList)
}