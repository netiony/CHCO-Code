library(readr)
library(ggplot2)
library(quantregGrowth)
library(tidyverse)
library(skimr)
#########

# DATA IMPORT & SUBSET 
ds <- read_csv("Desktop/XXY Abstract/Analysis/linked_ks_raw_R.csv")
ds <- subset(ds, select=c(person_id, kline, race, ethnicity, age, 
                          hgt_dup_16, t_height, id_count, id_mw))
ds <- rename(ds, id=person_id, hgt=hgt_dup_16, t=t_height)
ds <- subset(ds, kline=='Case' & age>=2 & age<=20)
ds <- subset(ds, id!='47706' & id!='38877' & id!='34438')
#ds <- subset(ds, id!='47706' & id!='38877' & id!='34438' & id!='13652')

ds <- ds[complete.cases(ds),]
ds <- data.frame(rbind(ds))
fac <- c(2,3,4,7)
ds[,fac] <- lapply(ds[,fac],factor)

str(ds)
View(ds)

# FORMULA - MONOTONE; TESTO & ID_MW AS COVARIATES
tau2=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
mform <-gcrq(data=ds, hgt ~ ps(age, monotone=1) + t + id_mw, tau=tau2) #id_mw as covariate vs "weights" similar 

mform05 <-gcrq(data=ds, hgt ~ ps(age, monotone=1) + t + id_mw, tau=0.05)
mform10 <-gcrq(data=ds, hgt ~ ps(age, monotone=1) + t + id_mw, tau=0.1)
mform25 <-gcrq(data=ds, hgt ~ ps(age, monotone=1) + t + id_mw, tau=0.25)
mform50 <-gcrq(data=ds, hgt ~ ps(age, monotone=1) + t + id_mw, tau=0.5)
mform75 <-gcrq(data=ds, hgt ~ ps(age, monotone=1) + t + id_mw, tau=0.75) 
mform90 <-gcrq(data=ds, hgt ~ ps(age, monotone=1) + t + id_mw, tau=0.9)
mform95 <-gcrq(data=ds, hgt ~ ps(age, monotone=1) + t + id_mw, tau=0.95)

# DATA VISUALIZATION
plot(mform, res=TRUE)
plot(mform, res=FALSE)
#plot(mform, legend=TRUE, conf.level=0.95, shade=TRUE, res=TRUE)
#plot(mform, deriv=TRUE, legend=TRUE, res=TRUE)

########
par(mar=c(3,3,3,2) + 0.1)
mplot95 <- plot(mform95,res=FALSE,xaxt="n",yaxt="n",xlab='',ylab='', 
                ylim=c(70,212), xlim=c(1.5,20.5), lwd=2, tck=0.02)
plot(mform05,res=FALSE,xlab='',ylab='',axes=FALSE,add=TRUE, lwd=2)
plot(mform10,res=FALSE,xlab='',ylab='',axes=FALSE,add=TRUE)
plot(mform25,res=FALSE,xlab='',ylab='',axes=FALSE,add=TRUE)
plot(mform50,res=FALSE,xlab='',ylab='',axes=FALSE,add=TRUE, lwd=2)
plot(mform75,res=FALSE,xlab='',ylab='',axes=FALSE,add=TRUE)
plot(mform90,res=FALSE,xlab='',ylab='',axes=FALSE,add=TRUE)
# Axis labels 
title(xlab='AGE (YEARS)', line=-1.25, cex.lab=0.9, font.lab=2)
title(ylab='STATURE', line=0.25, cex.lab=0.9, font.lab=2)
axis(1, at=2:20, lwd=0.5, tck=1, cex.axis=0.95, padj=-1.7)
axis(2, at=seq(70,210,5), lwd=0.5, tck=1, cex.axis=0.95, las=1, hadj=0.5)
axis(4, at=seq(70,210,5), lwd=0.5, tck=1, cex.axis=0.95, las=1, hadj=0.5)
axis(3, at=3:19, lwd=0.5, tck=1, cex.axis=0.95, padj=3.65)

########
# ADDITION OF CDC DATA 
cdc <- read_csv("Desktop/XXY Abstract/Analysis/cdc_table.csv")
cdc.form05 <- gcrq(data=cdc, five ~ ps(age), tau=0.5)
cdc.form10 <- gcrq(data=cdc, ten ~ ps(age), tau=0.5)
cdc.form25 <- gcrq(data=cdc, two_five ~ ps(age), tau=0.5)
cdc.form50 <- gcrq(data=cdc, five_zero ~ ps(age), tau=0.5)
cdc.form75 <- gcrq(data=cdc, seven_five ~ ps(age), tau=0.5)
cdc.form90 <- gcrq(data=cdc, nine_zero ~ ps(age), tau=0.5)
cdc.form95 <- gcrq(data=cdc, nine_five ~ ps(age), tau=0.5)

plot(cdc.form95,res=FALSE,xlab='',ylab='',axes=FALSE,add=TRUE, lwd=2, col="red")
plot(cdc.form90,res=FALSE,xlab='',ylab='',axes=FALSE,add=TRUE, col="red")
plot(cdc.form75,res=FALSE,xlab='',ylab='',axes=FALSE,add=TRUE, col="red")
plot(cdc.form50,res=FALSE,xlab='',ylab='',axes=FALSE,add=TRUE, lwd=2, col="red")
plot(cdc.form25,res=FALSE,xlab='',ylab='',axes=FALSE,add=TRUE, col="red")
plot(cdc.form10,res=FALSE,xlab='',ylab='',axes=FALSE,add=TRUE, col="red")
plot(cdc.form05,res=FALSE,xlab='',ylab='',axes=FALSE,add=TRUE, lwd=2, col="red")
##########

#FORMULA - MONOTONE & WEIGHTED; NO COVARIATES 
age <- seq(2, 20, by=0.5)
tau2=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
mform <-gcrq(data=ds, hgt ~ ps(age, monotone=1), tau=tau2)
mform.c <- charts(mform, k=age)
mform.c
plot(mform, res=TRUE)
plot(mform, res=FALSE)
####################

####################
library(readr)
library(ggplot2)
library(quantreg.nonpar)
library(tidyverse)
library(skimr)
#########

# DATA IMPORT
ds <- read_csv("Desktop/XXY Abstract/Analysis/linked_ks_raw_R.csv")
ds <- subset(ds, select=c(person_id, kline, race, ethnicity, age, 
                          hgt_dup_16, t_height, id_count, id_mw))
ds <- rename(ds, id=person_id, hgt=hgt_dup_16, t=t_height)
ds <- subset(ds, kline=='Case' & age>=2 & age<=20)
ds <- ds[complete.cases(ds),]

ds <- data.frame(rbind(ds))
fac <- c(2,3,4,7)
ds[,fac] <- lapply(ds[,fac],factor)

facid <- factor(ds$id)
ds$facid <- facid
facid_mw <- factor(ds$id_mw)
ds$facid_mw <- facid_mw
facage <- factor(ds$age)
ds$facage <- facage

str(ds)
View(ds)

# GGPLOT
ggplot(data=ds, aes(x=age, y=hgt, colour=id_mw)) + 
  labs(x="Age (years)",y="Height (cm)",colour="Patient ID") +
  geom_point(size=0.8, alpha=0.7) + 
  geom_smooth(se=FALSE, size=0.5, colour="black")

ggplot(data=ds, aes(x=age, y=hgt, colour=facid_mw)) + 
  labs(x="Age (years)",y="Height (cm)",colour="Patient ID") +
  geom_point(size=0.8, alpha=0.7) + 
  geom_smooth(se=FALSE, size=0.5, colour="black") +
  theme(legend.position="none")

ggplot(data=ds, aes(x=age, y=hgt, colour=t)) + 
  labs(x="Age (years)",y="Height (cm)",colour="Testosterone") +
  geom_point(size=0.8, alpha=0.7) + 
  geom_smooth(se=FALSE, size=0.5, colour="black") + 
  scale_colour_brewer(palette="Paired")

## HISTOGRAM
ggplot(ds, aes(x=age)) + 
  geom_histogram(aes(y=..density..), binwidth=0.5, 
                 colour="black", fill="white") + 
  geom_vline(aes(xintercept=mean(age)), color="black", 
             linetype="dashed", size=1) + 
  geom_density(alpha=.2, fill="#FF6666") 

ggplot(ds, aes(x=hgt)) + 
  geom_histogram(aes(y=..density..), binwidth=2, 
                 colour="black", fill="white") + 
  geom_vline(aes(xintercept=mean(hgt)), color="black", 
             linetype="dashed", size=1) + 
  geom_density(alpha=.2, fill="#FF6666") 

# MODEL STRUCTURE (NOTE id_weight as a factor???)
formula.par <- ds$hgt ~ t + id_mw #formula containing controls in model 
basis.bsp <- create.bspline.basis(breaks = quantile(ds$age, c(0:20)/20)) #quantiles of age for regional inference 

n=length(ds$age)
B=500
alpha=.05
taus=c(1:19)/20 #set of all analyzed quantile indices 
print.taus=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95) #quantile indices that will be print 

# AVERAGE FIRST DERIVITIVE (output is average growth rate at each quantile of height)
piv.bsp <- npqr(formula=formula.par, data=ds, basis=basis.bsp, 
                var="age", taus=taus, print.taus=print.taus, B=B, nderivs=1, 
                average=1, alpha=alpha, process="pivotal", rearrange=FALSE, 
                uniform=TRUE, se="unconditional", printOutput=TRUE, method="fn")
piv.bsp$point.est
piv.bsp$coefficients

# FIRST DERIVITIVE - (growth velocity by quantile)
piv.bsp.firstderiv <- npqr(formula=formula.par, data=ds, basis=basis.bsp, 
                           var="age", taus=taus, B=B, nderivs=1, average=0, 
                           print.taus=print.taus, process="none", printOutput = FALSE)
piv.bsp.firstderiv$point.est
piv.bsp.firstderiv$coefficients

# SECOND DERIVITIVE (growth acceleration by quantile)
piv.bsp.secondderiv <- npqr(formula=formula.par, data=ds, basis=basis.bsp, 
                            var="age", taus=taus, B=B, nderivs=2, average=0, 
                            print.taus=print.taus, process="none", printOutput = FALSE)
piv.bsp.secondderiv$point.est
piv.bsp.secondderiv$coefficients

# DATA VISUALIZATION (VELOCITY & ACCELERATION) only with unique data points??
xsurf1 <- as.vector(piv.bsp.firstderiv$taus)
ysurf1 <- as.vector(piv.bsp.firstderiv$var.unique)
zsurf1 <- t(piv.bsp.firstderiv$point.est)
xsurf2 <- as.vector(piv.bsp.secondderiv$taus)
ysurf2 <- as.vector(piv.bsp.secondderiv$var.unique)
zsurf2 <- t(piv.bsp.secondderiv$point.est)

par(mfrow = c(1, 2))
persp(xsurf1, ysurf1, zsurf1, xlab = "Quantile Index", ylab = "Age (years)",
      zlab = "Growth Rate", ticktype = "detailed", phi=30, theta=120, d=5,
      col = "green", shade = 0.75, main = "Growth Rate (B-splines)")
persp(xsurf2, ysurf2, zsurf2, xlab = "Quantile Index", ylab = "Age (years)",
      zlab = "Growth Acceleration", ticktype = "detailed", phi = 30, theta = 120,
      d = 5, col = "green", shade = 0.75, main = "Growth Acceleration (B-splines)")
#########

# ZERO DERIVITIVE (unclear rearrangement how to plot without forever broken with way too long of output)
piv.fac.fun <- npqr(formula =formula.par, data=ds, basis=facage, var = "age", taus = taus,
                    print.taus = print.taus, B = B, nderivs = 0, average = 0, alpha = alpha,
                    process = "none", rearrange = FALSE, rearrange.vars = "both", se = "conditional",
                    printOutput = TRUE, method = "fn")

piv.fac.fun.re <- npqr(formula =formula.par, data=ds, basis=facage, var = "age", taus = taus,
                       print.taus = print.taus, B = B, nderivs = 0, average = 0, alpha = alpha,
                       process = "none", rearrange = TRUE, rearrange.vars = "both", se = "conditional",
                       printOutput = TRUE, method = "fn")

xsurf <- as.vector(piv.fac.fun$taus)
ysurf <- as.vector(piv.fac.fun$var.unique)
zsurf.fac <- t(piv.fac.fun$point.est)
zsurf.fac.re <- t(piv.fac.fun.re$point.est)
par(mfrow = c(1, 2))
persp(xsurf, ysurf, zsurf.fac, xlab = "Quantile", ylab = "Age (years)",
      zlab = "Height (cm)", ticktype = "detailed", phi = 30, theta = 40, d = 5,
      col = "green", shade = 0.75, main = "XXY Growth Chart")
persp(xsurf, ysurf, zsurf.fac.re, xlab = "Quantile", ylab = "Age (years)",
      zlab = "Height (cm)", ticktype = "detailed", phi = 30, theta = 40, d = 5,
      col = "green", shade = 0.75, main = "XXY Growth Chart (rearranged)")
##########

##########