#Demo of MAST package
BiocManager::install("MAST")
library(MAST)
data(vbeta)
colnames(vbeta)
vbeta <- computeEtFromCt(vbeta)
vbeta.fa <- FromFlatDF(vbeta, idvars=c("Subject.ID", "Chip.Number", "Well"),
                       primerid='Gene', measurement='Et', ncells='Number.of.Cells',
                       geneid="Gene",  cellvars=c('Number.of.Cells', 'Population'),
                       phenovars=c('Stim.Condition','Time'), id='vbeta all', class='FluidigmAssay')
head(rowData(vbeta.fa),3)
sub1 <- vbeta.fa[,1:10]
show(sub1)
sub2 <- subset(vbeta.fa, Well=='A01')
show(sub2)
sub3 <- vbeta.fa[6:10, 1:10]
show(sub3)
sp1 <- split(vbeta.fa, 'Subject.ID')
show(sp1)
vbeta.split<-split(vbeta.fa,"Number.of.Cells")
#see default parameters for plotSCAConcordance
plotSCAConcordance(vbeta.split[[1]],vbeta.split[[2]],
                   filterCriteria=list(nOutlier = 1, sigmaContinuous = 9,
                                       sigmaProportion = 9))
## Split by 'ncells', apply to each component, then recombine
vbeta.filtered <- mast_filter(vbeta.fa, groups='ncells')
## Returned as boolean matrix
was.filtered <- mast_filter(vbeta.fa, apply_filter=FALSE)
## Wells filtered for being discrete outliers
head(subset(was.filtered, pctout))

vbeta.1 <- subset(vbeta.fa, ncells==1)
## Consider the first 20 genes
vbeta.1 <- vbeta.1[1:20,] 
head(colData(vbeta.1))

library(ggplot2)
zlm.output <- zlm(~ Population + Subject.ID, vbeta.1,)
show(zlm.output)

coefAndCI <- summary(zlm.output, logFC=FALSE)$datatable
coefAndCI <- coefAndCI[contrast != '(Intercept)',]
coefAndCI[,contrast:=abbreviate(contrast)]
ggplot(coefAndCI, aes(x=contrast, y=coef, ymin=ci.lo, ymax=ci.hi, col=component))+
  geom_pointrange(position=position_dodge(width=.5)) +facet_wrap(~primerid) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + coord_cartesian(ylim=c(-3, 3))

zlm.lr <- lrTest(zlm.output, 'Population')
pvalue <- ggplot(melt(zlm.lr[,,'Pr(>Chisq)']), aes(x=primerid, y=-log10(value)))+
  geom_bar(stat='identity')+facet_wrap(~test.type) + coord_flip()
print(pvalue)

#With scRNA data
suppressPackageStartupMessages({
  library(ggplot2)
  library(GGally)
  library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(stringr)
  library(NMF)
  library(rsvd)
  library(RColorBrewer)
  library(MAST)
})
#options(mc.cores = detectCores() - 1) #if you have multiple cores to spin
options(mc.cores = 1)
knitr::opts_chunk$set(message = FALSE,error = FALSE,warning = FALSE,cache = FALSE,fig.width=8,fig.height=6)
freq_expressed <- 0.2
FCTHRESHOLD <- log2(1.5)

data(maits, package='MAST')
dim(maits$expressionmat)
head(maits$cdat)
head(maits$fdat)
scaRaw <- FromMatrix(t(maits$expressionmat), maits$cdat, maits$fdat)
