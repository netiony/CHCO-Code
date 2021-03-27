library(FDRsampsize)

# pooled Y-T2D and non-T2D
res=fdr.sampsize(fdr=0.1,
                 ave.pow=0.8,
                 pow.func=power.tcorr,
                 eff.size=rep(c(0.3,0),c(100,6900)),
                 null.effect=0)
res

res=fdr.sampsize(fdr=0.1,
                 ave.pow=0.8,
                 pow.func=power.tcorr,
                 eff.size=rep(c(0.3,0),c(70,6930)),
                 null.effect=0)
res

# Y-T2D only
res=fdr.sampsize(fdr=0.1,
                 ave.pow=0.8,
                 pow.func=power.tcorr,
                 eff.size=rep(c(0.7,0),c(70,6930)),
                 null.effect=0)
res