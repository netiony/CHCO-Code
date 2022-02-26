library(sas7bdat)
library(tableone)
library(NormalizeMets)
library(mixOmics)
source("~/GitHub/General-code/foldchange.R")
source("~/GitHub/General-code/editcolnames.R")
setwd(home_dir)
# read in Petter's SAS dataset
alldata <- read.sas7bdat("./Raw data/casperheir_for_laura.sas7bdat")

# read in NG data - those variables are included in the original dataset, but we want to analyze
# them separately
ng_serum <- read.csv("./Raw data/NG serum.csv")
names_ng_serum <- colnames(ng_serum)
ng_urine <- read.csv("./Raw data/ng urine.csv")
names_ng_urine <- colnames(ng_urine)
names_ng <- unique(c(names_ng_serum,names_ng_urine))
names_ng <- names_ng[names_ng != "ID"]

# find variables in alldata that are also in NG
col_alldata <- colnames(alldata)
dups <- col_alldata[col_alldata %in% names_ng]
# take dups out of alldata
alldata <- alldata[,!(colnames(alldata) %in% dups)]

# remove two patients who dropped out and one who could not complete study
alldata <- alldata[!(alldata$ID %in% c("CS-18","CS-45","CS-27")),]

npart <- length(unique(alldata$ID))

targeted <- c("Ala3","Arg3","CSSC3","Glu3","Gln3","His3","Ile3","Leu3","Lys3","Met3","Phe3",
              "Pro3","Ser3","Thr3","Tyr3","Val3","Citrate3","aKG3","succinate3","Fumarate3",
              "Malate3")
nonnorm <- c("Ala3","Ser3")

tab_targeted_group <- CreateTableOne(vars=targeted, data=alldata, strata = "group", test=TRUE)
tab_targeted_group <- print(tab_targeted_group,varLabels=FALSE)

# add FDR corrected p-values
# need to do individal t-tests so I get the actual p-values
p_raw <- NULL
for (i in 1:length(targeted)) {
  temp <- t.test(alldata[alldata$group==1,paste0(targeted[i])],alldata[alldata$group==4,paste0(targeted[i])],data=alldata,var.equal = T)$p.value
  p_raw <- c(p_raw,temp)
}
p_adj <- round(p.adjust(p_raw),3)
p_adj[p_adj<0.001] <- "<0.001"
p_adj <- c("",p_adj)
tab_targeted_group <- cbind(tab_targeted_group,as.matrix(p_adj))
colnames(tab_targeted_group) <- c("1","4","p","test","adj-p")

# comparison of targeted metabolites adjusted for GFR
contrast_adj <- model.matrix(~group + gfr, data=alldata)
ymat <- t(alldata[,targeted])
fit_group <- lmFit(ymat,contrast_adj)
fit_group <- eBayes(fit_group)
# Non-moderated reults
fit_group_reg = lapply(as.data.frame(t(ymat)), function(y){
  mod = lm(y ~ contrast_adj[,-1])
  return(summary(mod)$coefficients[2,4])
})
# format results
results <- as.data.frame(do.call(rbind,fit_group_reg))
colnames(results) = "p"
results$q = p.adjust(results$p,"fdr")
results$p = format.pval(results$p,eps = 0.001,digits = 3)
results$q = format.pval(results$q,eps = 0.001,digits = 3)
write.csv(results,"./Reports/targeted_adjusted_GFR_reg.csv")

# create df for NormalizeMets
metabolite <- as.data.frame(targeted)
row.names(metabolite) <- metabolite$targeted

# featuredata has compound IDs as colnames and sample IDs as rownames
featuredata <- cbind(alldata$ID,alldata[,targeted])
row.names(featuredata) <- featuredata$`alldata$ID`
featuredata$`alldata$ID` <- NULL

sampledata <- alldata[,!(colnames(alldata) %in% targeted)]

# combine the three dataframes
allmetabdata <- list(featuredata=featuredata,sampledata=sampledata,metabolitedata=metabolite)

# Format group
allmetabdata$sampledata$group = factor(allmetabdata$sampledata$group,levels = c(1,4),labels = c("T1D","Control"))

# log transform
logdata <- LogTransform(allmetabdata$featuredata)

# PCA plot
res.pca <- mixOmics::pca(logdata$featuredata,center=T,scale = T)

# create dataset for PLS-DA
plsda1 <- as.data.frame(sampledata[,c("ID","group")])
plsda2 <- as.data.frame(allmetabdata$featuredata)
plsda2 <- as.data.frame(t(plsda2))
#plsda2$COMP.ID <- rownames(plsda2)
#names <- metabolite[,c("COMP.ID","BIOCHEMICAL")]
#plsda2 <- merge(plsda2,names,by="COMP.ID",all.x = TRUE,all.y = TRUE)
#rownames(plsda2) <- plsda2$BIOCHEMICAL
#plsda2$COMP.ID <- NULL
#plsda2$BIOCHEMICAL <- NULL
plsda2 <- as.data.frame(t(plsda2))
plsda2$ID <- rownames(plsda2)
plsda <- merge(plsda1,plsda2,by="ID",all.x = TRUE,all.y = TRUE)
rownames(plsda) <- plsda$ID
plsda$ID <- NULL
plsda[,2:22] <- apply(plsda[2:22],2,as.numeric)

# sPLS-DA analysis
plsda.res = splsda(X = plsda[,c(2:22)], Y=factor(plsda$group,levels = c(1,4),labels = c("T1D","Control")), ncomp = 2)
p = plotIndiv(plsda.res,comp = c(1,2),ellipse = T,legend = T,ind.names = F)
targeted_plsda = ggplot(p$df,aes(x=x,y=y,color = group)) + 
  geom_point(color = p$df$col,aes(shape = group)) + 
  xlab(p$graph$labels$x) + ylab(p$graph$labels$y) +
  ggtitle("Targeted") +
  theme_bw() +
  geom_path(data = p$df.ellipse,aes(x = Col1,y = Col2),inherit.aes = F,color = levels(factor(p$df$col))[1]) +
  geom_path(data = p$df.ellipse,aes(x = Col3,y = Col4),inherit.aes = F,color = levels(factor(p$df$col))[2]) +
  theme(axis.text=element_blank(),axis.ticks = element_blank(),
        legend.position = "none",plot.title = element_text(hjust = 0.5))

plsda.perf = perf(plsda.res, validation = 'Mfold', folds = 5,progressBar = FALSE, 
                  nrepeat = 10, dist = 'max.dist',auc=TRUE)
auc_save <- plsda.perf$auc$comp1[1]
auc_true <- as.numeric(plsda.perf$auc$comp1["AUC.mean"])
# TOP 20 compounds
plsda21 = splsda(X = plsda[,c(2:22)], Y=as.factor(plsda$group), ncomp = 2, keepX = c(21,21))
# get list of top 20 compounds
top21 <- selectVar(plsda21,comp=1)
top21_2 <- selectVar(plsda21,comp=2)

# correlations with continuous variables
tubular = read.csv("./Clean data/tubularinjury.csv",stringsAsFactors = F)
alldata = left_join(alldata,tubular[,c("ID","KIM1","NGAL","MCP1","logYKL40","logIL18")],by = "ID")
#	metabolomics - GIR
#	metabolomics - body composition data 
#	metabolomics - renal measures (GFR, RPF, ACR, renal O2)
corrvars = c("rpf","acr_mean","gfr")
varnames = c("RPF","Mean ACR","GFR")
corrs <- corr.test(alldata[,targeted],alldata[,corrvars],method = "spearman",adjust = "fdr")
corr_sum <- round(corrs$r,3)
corr_p <-   round(corrs$p,3)
corr_p_char <- corr_p
corr_p_char[corr_p>=0.05] <- ""
corr_p_char[corr_p<0.05 & corr_p>=0.01] <- "*"
corr_p_char[!is.na(corr_p) & corr_p<0.01] <- "**"

z<-corr_sum
extra <- corr_p_char
rnames <- colnames(z)
mat_data <- data.matrix(z[,1:ncol(z)])
#rownames(mat_data) <- rnames
#par(cex.main=1.5)
#margins =c(12,20)

c = corr.test(alldata[,targeted],alldata[,corrvars],method = "spearman",adjust = "fdr")
rownames(c$r) <- c("Alanine","Arginine","Cystine","Glutamine","Glycine","Histidine",
                   "Isoleucine","Leucine","Lysine","Methionine","Phenylalanine",
                   "Proline","Serine","Threonine","Tyrosine","Valine","Citrate",
                   "Ketoglutarate","Succinate","Fumarate","Malate")
colnames(c$r) <- varnames

rownames(c$p) <- rownames(c$r)

targeted_heat = pheatmap(c$r,legend = T,main = "Targeted")
colnames(c$p) <- varnames