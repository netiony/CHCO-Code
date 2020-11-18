library(fdapace)

# load data
data(medfly25)

# Turn the original data into a list of paired amplitude and timing lists
Flies <- MakeFPCAInputs(medfly25$ID, medfly25$Days, medfly25$nEggs)
fpcaObjFlies <- FPCA(Flies$Ly, Flies$Lt, list(plot = TRUE, methodMuCovEst = 'smooth', userBwCov = 2))

A <- FClust(Flies$Ly, Flies$Lt, optnsFPCA = list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.90), k = 2)
# The Neg-Entropy Criterion can be found as: A$clusterObj@bestResult@criterionValue 
CreatePathPlot( fpcaObjFlies, K=2, showObs=FALSE, lty=1, col= A$cluster, xlab = 'Days', ylab = '# of eggs laid')
grid()
