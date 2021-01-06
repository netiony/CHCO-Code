library(fishmethods)

data <- read.csv("C:\\Temp\\DATA_2020-04-15 Espinosa TESTO metabolomics jrh - final.csv")
clus.rho(popchar = data$Linoleate, cluster = data$PID)
clus.rho(popchar = data$Octadecatrienoic.acid..linolenic.acid., cluster = data$PID)
clus.rho(popchar = data$Eicosatetraenoic.acid, cluster = data$PID)
clus.rho(popchar = data$Eicosapentaenoic.acid, cluster = data$PID)
clus.rho(popchar = data$Docosahexaenoic.acid, cluster = data$PID)
