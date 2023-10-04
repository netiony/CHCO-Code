###################################
# 
# Filter models by FDR adjustment
#
###################################

dataOut = "/Volumes/RI Biostatistics Core/Shared/Shared Projects/Bothwell/Peds ENDO/Siebert - Multiomics/Data/omes_models_9_22_2023.csv"
#modelsOut = "~/Data/Borengasser/AHA/toAnalyze/ommMbGModels_DXARichApr2020.csv"
modelsOut = "/Volumes/RI Biostatistics Core/Shared/Shared Projects/Bothwell/Peds ENDO/Siebert - Multiomics/Data/omes_models_filtered_9_22_2023.csv"

#dataOut = "~/Data/Borengasser/AHA/toAnalyze/ommMbGP_Apr2020.csv"
#modelsOut = "~/Data/Borengasser/AHA/toAnalyze/models_R01_Apr2020.csv"

d = read.csv(dataOut)
table(d$assays)


# Filter out massive influence points
d = subset(d, maxInfluence < 4) 
table(d$assays)


# Calculate FDR adjusted p-value
d$pAdjFDR = NA

for (i in seq_along(unique(d$assays))){
  pair = unique(d$assays)[i]
  #p.adjust is vectorized
  d$pAdjFDR[d$assays == pair] = p.adjust(d$pCorr[d$assays == pair], method="fdr")
  
}


# organize data
dO = d[order(d$pCorr,decreasing=FALSE),]
dr <- by(dO, dO["assays"], head, n=25)
topNbyAssays = droplevels(Reduce(rbind, dr))


# Filter on FDR 
topN05 = subset(topNbyAssays, pAdjFDR < 0.2)
sort(table(topN05$ana1))
sort(table(topN05$ana2))


# Filter out pa.pg pairs - we already know they're strongly correlated
topN05 <- topN05[!(topN05$assays == "pa.pg"), ]



write.csv(topN05, modelsOut, quote=FALSE, row.names=FALSE)
