#HPC Simulation Script ----
# Load Libraries ----
library(rjags)
library(R2jags)
library(dplyr)
#library(ggplot2)
library(pscl)
library(purrr) 
library(slurmR)

# Set CRAN mirror----
# options(repos = c(CRAN = "https://cloud.r-project.org"))
# install.packages("qgcomp")
library(qgcomp)
# install.packages("pbdMPI",lib="/home1/hhampson/R/x86_64-pc-linux-gnu-library/4.3",type="source")
library(pbdMPI)
library(MASS)


#Set working directory----
setwd("/project/jagoodri_1060/Hailey")

#Load functions----
source("HPC_Simulation_Functions_03_12.R")

#Set Iteration Number (1:15) -----
# for (k in 1:10) {
for (k in 1) {
  # for (scenario in 1:40) {
  for (scenario in 30) {
    
    #Load data----
    Object <- readRDS("FormattedObject.RDS")
    Table <- data.frame(Object$Table)
    Table <- Table[-149]
    P.s.causal.all <- readRDS("Simulation_Scenarios_03_12_25.rds")
    sim.par <- readRDS(paste0("Simulation_Parameter_",scenario,".rds"))  
      # rename(P.s.causal=P.s.scenario)
    # sim.par <- sim.par %>% 
    #   mutate(P.e.causal=ifelse(P.e.causal==6,7,P.e.causal)) 
    
    #Set Parameters----
    P.s <- 206
    P.g <- 92
    P.f <- 33
    P.o <- 22
    P.c <- 11
    P.p <- 7
    P.e <- sim.par$P.e[[1]]
    
    # Load Z matrices----
    # Genus
    # Z.s.g <- readRDS("/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Z.s.g.RDS")
    Z.s.g <- readRDS("Z.s.g.RDS")
    Genus.R <- ncol(Z.s.g)
    #Family
    # Z.g.f <- readRDS("/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Z.g.f.RDS")
    Z.g.f <- readRDS("Z.g.f.RDS")
    Family.R <- ncol(Z.g.f)
    #Order
    # Z.f.o <- readRDS("/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Z.f.o.RDS")
    Z.f.o <- readRDS("Z.f.o.RDS")
    Order.R <- ncol(Z.f.o)
    #Class
    # Z.o.c <- readRDS("/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Z.o.c.RDS")
    Z.o.c <- readRDS("Z.o.c.RDS")
    Class.R <- ncol(Z.o.c)
    #Phylum
    # Z.c.p <- readRDS("/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Z.c.p.RDS")
    Z.c.p <- readRDS("Z.c.p.RDS")
    Phylum.R <- ncol(Z.c.p)
    
    #Create data for simulation
    Y.cont <- Table[, grep("d__", names(Table))]
    Y.names <- colnames(Y.cont)
    Y.obs <- apply(Y.cont, 2, FUN=function(v) { ifelse(v==0, 0, 1)})
    Y.freq <- as.numeric(apply(Y.obs, 2, mean))
    Y.mean <- apply(Y.cont, 2, FUN=function(v) { mean(v[v!=0])})
    Y.mean <- ifelse(is.na(Y.mean), 0.001, Y.mean)
    Y.freq <- Y.freq[Y.freq !=0 & Y.freq !=1]
    Y.freq <- Y.freq[1:P.s] # frequency of species presence vs. absence
    phi <- 5 
    
    # Initiate connections----
    init()
    
    #Function to run Bayesian, Ridge and Zero-Inflated Models ----
    model.fxn <- function(i) {
       # i = 1
      #Set seed for reproducibility
      set.seed(i)
      
      #Simulate Data----
      sim.output <- sim.fxn(i)
      
      # Save the simulated dataset for future use
      # saveRDS(sim.output, file = paste0("Dataset_",scenario,"_scenario_", k, "_", i, ".RDS"))
      
      #Extract outputs from simulated data
      X <- sim.output$X
      Y <- sim.output$Y
      profiles <- sim.output$profiles
      OR_exposure <- sim.par$OR.exposure[i]
      N <- sim.par$N[i]
      Y2 <- data.frame(Y)
      colnames(Y2) <- paste0("species",1:length(colnames(Y2)))
      
      #Format Data for ZING Model 
      Y.g <- as.data.frame(as.matrix(Y2)%*%as.matrix(Z.s.g))
      Y.f <- as.data.frame(as.matrix(Y.g)%*%as.matrix(Z.g.f))
      Y.o <- as.data.frame(as.matrix(Y.f)%*%as.matrix(Z.f.o))
      Y.c <- as.data.frame(as.matrix(Y.o)%*%as.matrix(Z.o.c))
      Y.p <- as.data.frame(as.matrix(Y.c)%*%as.matrix(Z.c.p))
      Y.g$id <- rownames(Y.g)
      Y.f$id <- rownames(Y.f)
      Y.o$id <- rownames(Y.o)
      Y.c$id <- rownames(Y.c)
      Y.p$id <- rownames(Y.p)
      Y2$id <- rownames(Y2)
      Y2 <- full_join(Y2,Y.g,by="id")
      Y2 <- full_join(Y2,Y.f,by="id")
      Y2 <- full_join(Y2,Y.o,by="id")
      Y2 <- full_join(Y2,Y.c,by="id")
      Y2 <- full_join(Y2,Y.p,by="id")
      X2 <- as.data.frame(X)
      X2$id <- rownames(X2)
      colnames(X2) <- c(paste0("X.",1:(length(colnames(X2))-1)),"id")
      Y2 <- full_join(Y2,X2,by="id")
      id_col <- which(colnames(X2)=="id")
      exposure.names <- colnames(X2)[-id_col]
      
      #Names vector
      Taxa.names <- c(rownames(Z.s.g),colnames(Z.s.g),colnames(Z.g.f),colnames(Z.f.o),colnames(Z.o.c),colnames(Z.c.p))
      
      #Run ZING Model ----
      results_list <- map(Taxa.names,dat=Y2,exposure=exposure.names,ZING_Model)
      ZING.results <- bind_rows(results_list) 
      # return(ZING.results)
      
      #Run BaH-ZING Model----
      jdata <- list(N=N, Y=Y, P.s=P.s, X=X, P.e=P.e,
                    Z.s.g=Z.s.g, P.g=P.g,
                    P.f=P.f, Z.g.f=Z.g.f,
                    P.o=P.o, Z.f.o=Z.f.o,
                    P.c=P.c, Z.o.c=Z.o.c,
                    P.p=P.p, Z.c.p=Z.c.p,
                    profiles=profiles)
      var.s <- c("species.beta", "genus.beta", "family.beta", "order.beta", "class.beta", "phylum.beta","species.psi","genus.psi","family.psi","order.psi","class.psi","phylum.psi","species.beta.zero", "genus.beta.zero", "family.beta.zero", "order.beta.zero", "class.beta.zero", "phylum.beta.zero","species.psi.zero","genus.psi.zero","family.psi.zero","order.psi.zero","class.psi.zero","phylum.psi.zero","omega","disp")
      model.fit <- jags.model(file=textConnection(BHRM.microbiome), data=jdata, n.chains=1, n.adapt=100, quiet=F)
      update(model.fit, n.iter=10)
      model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=50, thin=1, progress.bar="none")
      #model.fit <- jags.model(file=textConnection(BHRM.microbiome), data=jdata, n.chains=3, n.adapt=100, quiet=F)
      #update(model.fit, n.iter=1000)
      #model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=5000, thin=1, progress.bar="none")
      # # summarize results
      r <- summary(model.fit)
      results <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))
      results <- results %>%
        mutate(sig = ifelse((X2.5.<0 & X97.5.<0) | (X2.5.>0 & X97.5.>0), "*","N.S."))
      # mutate(Odds.Ratio = exp(Mean),
      #        lcl=exp(X2.5.),
      #        ucl=exp(X97.5.))
      BaHZING.results <- results #5664 - 6348
      rm(model.fit,r,results)
      # return(BaHZING.results)
      # BaHZING.results <- results
      BaHZING.results <- BaHZING.results %>% 
        mutate(Taxa=rownames(BaHZING.results),
               Exposure=rownames(BaHZING.results)) %>%
        mutate(SE=SD/sqrt(N)) %>% 
        rename(P_Value=sig) %>% 
        mutate(Component=ifelse(grepl("zero",Taxa),"Probability","Means")) %>% 
        mutate(Model="BaH-ZING") %>% 
        # rename(OddsRatio=Odds.Ratio) %>% 
        dplyr::select(Taxa,Exposure,Mean,SE,P_Value,Component,Model)
      rownames(BaHZING.results) <- NULL
      
      #Run Ridge BaHZING Model ----
      jdata <- list(N=N, Y=Y, P.s=P.s, X=X, P.e=P.e,
                    P.g=P.g,
                    P.f=P.f,
                    P.o=P.o,
                    P.c=P.c,
                    P.p=P.p,
                    profiles=profiles)
      var.s <- c("species.beta", "genus.beta", "family.beta", "order.beta", "class.beta", "phylum.beta","species.psi","genus.psi","family.psi","order.psi","class.psi","phylum.psi","species.beta.zero", "genus.beta.zero", "family.beta.zero", "order.beta.zero", "class.beta.zero", "phylum.beta.zero","species.psi.zero","genus.psi.zero","family.psi.zero","order.psi.zero","class.psi.zero","phylum.psi.zero","omega","disp")
      model.fit <- jags.model(file=textConnection(RBHRM.microbiome), data=jdata, n.chains=1, n.adapt=100, quiet=F)
      update(model.fit, n.iter=10)
      model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=50, thin=1, progress.bar="none")
      #model.fit <- jags.model(file=textConnection(RBHRM.microbiome), data=jdata, n.chains=3, n.adapt=100, quiet=F)
      #update(model.fit, n.iter=1000)
      #model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=5000, thin=1, progress.bar="none")
      # summarize results
      r <- summary(model.fit)
      results <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))
      results <- results %>%
        filter(grepl("species",rownames(results)))
      results <- results %>%
        mutate(sig = ifelse((X2.5.<0 & X97.5.<0) | (X2.5.>0 & X97.5.>0), "*","N.S."))
      # mutate(Odds.Ratio = exp(Mean),
      #        lcl=exp(X2.5.),
      #        ucl=exp(X97.5.))
      RBaHZING.results <- results #6348
      RBaHZING.results <- RBaHZING.results %>%
        mutate(Taxa=rownames(RBaHZING.results),
               Exposure=rownames(RBaHZING.results)) %>%
        # mutate(SE=SD) %>%
        mutate(SE=SD/sqrt(N)) %>%
        # select(-c(X2.5.,X97.5.,lcl,ucl)) %>%
        rename(P_Value=sig) %>%
        mutate(Component=ifelse(grepl("zero",Taxa),"Probability","Means")) %>%
        mutate(Model="RBaH-ZING") %>%
        # rename(OddsRatio=Odds.Ratio) %>%
        dplyr::select(Taxa,Exposure,Mean,SE,P_Value,Component,Model)
      rownames(RBaHZING.results) <- NULL
      rm(model.fit,r,results)
      # return(RBaHZING.results)
      Results.Output <- rbind(ZING.results,BaHZING.results)
      Results.Output <- rbind(Results.Output,RBaHZING.results)
      return(Results.Output)
    }
    
    # Set number of iterations -------------------------------------------------
   #n_iter <- nrow(sim.par)
    n_iter <- 3
    
    # Run model ---------------------------------------------------------------
    coefs <- pbdLapply(X = 1:n_iter,
                       FUN = model.fxn,
                       pbd.mode = "spmd")
    combined.coefs = do.call(rbind,coefs)
    
    # Save results  -------------------------------------------------------------
    comm.write.csv(combined.coefs,  file = paste0("Small_Test_Full_results_",scenario,"_scenario_",k,"_03_17_25.csv"))
    # Write Success from each connection
    message(paste("SUCCESS from rank", comm.rank()))
    
    # End connections
    finalize(mpi.finalize = TRUE)
    #if (comm.size() > 1) {
    #finalize(mpi.finalize = TRUE)
#}

  }
}

# I have run into an issue where the finalize function does not always 
# close the connection, so the job will run  until it is terminated by the 
# Slurm scheduler, even if it is completed. Therefore, I have included this final 
# line of code. This line throws an error, but will also kill the slurm job without 
# wasting resources (please email me at jagoodri@usc.edu if you find a solution to this)
slurmR::Slurm_clean(coefs)
