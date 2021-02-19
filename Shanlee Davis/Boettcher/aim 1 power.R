library(Rfast)

data <- read.csv("E:/Davis/Davis Boettcher/Ross_SCA_GencodeEnsembl_RPKMs1.csv")

keep <- data[data$Gene %in% c("PDZK1","SLC25A1","PLTP","SLC27A1","SLC22A1",
                              "CPT1B","SCARB1"),]
t <- t(keep)
colnames(t) <- t[1,]
final <- as.data.frame(t[4:nrow(t),])

final <- apply(final,2,as.character)
final <- apply(final,2,as.numeric)

cv <- colcvs(final,ln=FALSE)
