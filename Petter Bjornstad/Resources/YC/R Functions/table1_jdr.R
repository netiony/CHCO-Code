#JDK
########## FUNCTIONS

#TODO
# Add function reporting for paired data (could be added to TableDesc but paired data comes in many forms, probably better for a seperate function)
# For paired code examples check chirbariatrique (cochran.qtest and mcnemar and friedman)
# Create a package or word doc with all the functions and exemples with description for each parameter
# Add nolevels option to the models ... ?? How would this work, reference is not nolevels
# Add checks in each function if an option is missing 
# Report.models check length of vars and labels 
# Add a bold.pval option in the report models (by default it is fixed at 0.05)
# Add the package infront of all package functions (width from flextable but if BioGenerics loads width no longer works)
# Add in sastor ignore format == TIME
# Add option in each function, word=TRUE to instantly open a word doc with the corresponding table

#Change log:
#2019:
#27/08: Added stop error to check normality function
#27/08: Added in chisq.test the expected<5 then do fisher
#29/08: Changed confint in report.poisson to confint.default (same results but no errors)
#02/09: Added a check for car::Anova() for when Wald statistics cannot be calculated
#04/09: Changed aov without equal variances to Kruskall Wallis
#05/09: Added car package to requiredPackages
#10/09: Created llibrary function which is like library but installs the package if it is not installed
#24/09: Added a progressbar in the TableDesc function and a warning if fisher.test is used with +20 categories
#25/09: Added in the function sastor another reference set for "OUI" and "NON" (capitals)
#26/09: By default if this file is loaded it adds timestamps to the console log
#27/09: Fixed examples and removed SD option for report models
#01/10: Changed sastor() function, no longer set all NA to 0 when only 1 level and new reference to 1 for c(0,1) levels. Also if Oui and Non in levels (+ others) then Oui is set as reference
#02/10: Changed sastor() lapply() for as.numeric or as.factor to make BEST also into numeric
#14/10: Added jdk.desc so now only this file needs to be sourced to load all the functions
#22/10: Added droplevels to almost all functions (at the start)
#22/10: Added a reference TRUE/FALSE option to all the report.model functions
#22/10: Added the reference TRUE/FALSE option to TableDesc and FTdesc to remove the automatic ", n (%)" after the label for categorical variables
#22/10: In the tablemissing function added an option TRT to specificy the number of missing values per TRT level
#23/10: Added (n=..) in the tablemissing function
#23/10: Removed the exact=FALSE from the wilcox.test in comptest function (this changes all pvalues a little if there weren't any ties)
#23/10: Added a check to all the report.model functions to check if the labels length is the same as the number of variables 
#23/10: Removed the totonly parameter from the jdk.desc function (is not an interactable parameter)
#24/10: Added better flextable formatting to the tablemissing function
#29/10: Removed droplevels() from getresult (was returning (%) instead of 0(0%))
#31/10: Replaced an rapply to check for factor or numeric variables in TableDesc and tablemissing (did not work for 1 variable specified)
#08/11: Added line in TableDesc to transform tibble into data.frame objects
#12/11: In report.poisson exponentiated the coefficients and the CI (thus any model before this was incorrect)
#19/11: Changed the tablemissing width specifications to be the same the TableDesc (if TRT is specified, else it's autofitted)
#2020:
#23/01: Added in tablemissing and TableDesc and jdk.desc that if names is missing it is defined as vars 
#23/01: Changed the jdk.desc to have space=TRUE by default
#30/01: jdk.desc now ignores the "group" variable if it was in the varlist
#04/02: Added YES/NO reference set to the sastor function

#Adds time stamps to the console for logs
#Updates the time after the previous line was executed
updatePrompt <- function(...) {options(prompt=paste0(format(Sys.time(), "%H:%M:%S"),"> ")); return(TRUE)}
addTaskCallback(updatePrompt)



'%ni%' <- Negate('%in%')


#Install the required packages
requiredPackages = c("haven","flextable","officer","devEMF","car","dplyr","reshape2")
for(p in requiredPackages){
  if(p %ni% installed.packages()[,"Package"]) {
    install.packages(p)
    library(p,character.only = TRUE)
  } else {library(p,character.only = TRUE)}
}
rm(p)
rm(requiredPackages)



#Function to get table for NA for a variable:
tna <- function(nom) {return(table(is.na(nom)))}
#Function to get table with NA (no need to write useNA="always")
tablena <- function(nom) {return(table(nom,useNA="always"))}


#Function to quickly get the unique length of a variable
lunique <- function(x) {return(length(unique(x)))}




#Function to set as character and then as numeric:
as.numchar <- function(x) {return(suppressWarnings(as.numeric(as.character(x))))} 
#Warning = always introduction of NA if not numeric (<0.05 or >5000 or DM or 0,5)


#Functions to load package or install if not installed
llibrary <- function(p) {
  if (!grepl("\"",deparse(substitute(p)))) {p=deparse(substitute(p))}
  if(p %ni% installed.packages()[,"Package"]) {
    install.packages(p)
    if (class(try(library(p,character.only = TRUE),silent=TRUE))!="try-error") {library(p,character.only = TRUE)}
  } else {library(p,character.only = TRUE)}
} 


#Different View function for viewing large datasets
Vieww <- function(x) {
  utils::View(x)
} 



#Function to get moy+ET (or med(Q1-Q3) if not normal) for numeric and n(%) for factors
##EX: #CTRL+SHIFT+C to comment/uncomment below
# data("mtcars")
# mtcars$am=as.factor(mtcars$am)
# mtcars$mpg=as.numchar(mtcars$mpg)
# getresult(data=mtcars,var="am") #For factors level needs to be specified
# getresult(data=mtcars,var="am",levelj="1")
# getresult(data=mtcars,var="mpg",virg=2) #By default 0 after virg
# #Normality: 1 = normal, 0 = not normal
# getresult(data=mtcars,var="mpg",normal=0) #If normality is not specified it is calculated using checknormality function
getresult <- function(data,var,levelj=NULL,virg=NULL,normal=NULL)
{
  
  x=NULL
  y=NULL
  #By default number of decimals = 0
  if (is.null(virg)) {virg=0}
  
  #If normal is empty then use checknormality function
  if (is.null(normal)) {
    normal=checknormality(data=data,var=var)
  }
  
  #FACTOR
  if (is.factor(data[,var])) {
    #Only NA then return "-":
    if (all(is.na(data[,var]))) {return("-")}
    
    tot=sum(table(data[,var]))
    x=table(data[,var])[which(levels(data[,var])==levelj)]
    y=round(x/tot,2)
    return(paste(x," (",y*100,"%)",sep=""))
  }
  
  #NUMERIC
  if (is.numeric(data[,var])) {
    #If not normal then median + Q1 and Q3
    if (normal==0) { 
      #Only NA then return nothing:
      if (all(is.na(data[,var]))) {return("-")}

      x=format(round(median(data[,var],na.rm=TRUE),virg), nsmall = virg)
      y=format(round(quantile(data[,var],na.rm=TRUE,type=2)[2],virg), nsmall = virg)
      z=format(round(quantile(data[,var],na.rm=TRUE,type=2)[4],virg), nsmall = virg)
      return(paste(x," (",y,"-",z,")",sep=""))
      
    } else {
      #Only NA then return nothing:
      if (all(is.na(data[,var]))) {return("-")}
      
      x=format(round(mean(data[,var],na.rm=TRUE),virg), nsmall = virg)
      y=format(round(sd(data[,var],na.rm=TRUE),virg), nsmall = virg)
      return(paste(x," \xb1 ",y,sep="")) #Only unicode string http://unicode.scarfboy.com/?s=U%2b00b1
    }
  }
}



##########################


#Check normality: returns 0 or 1 (1 implies normality)
#Uses the Shapiro test to evaluate normality.
# If no TRT is specified it only checks normality in the total population but if
#   a TRT is specified it will check normality in the global pop and for each TRT pop.
# Normality is returned if in the global pop and all TRT pops the variable is normally distributed.
##EX: #CTRL+SHIFT+C to comment/uncomment below
# data("mtcars")
# checknormality(data=mtcars,var="qsec")
# qqnorm(mtcars$qsec,main="QQ plot of data",pch=19)
# qqline(mtcars$qsec, col = 2)
# checknormality(data=mtcars,var="disp")
# qqnorm(mtcars$disp,main="QQ plot of data",pch=19)
# qqline(mtcars$disp, col = 2)
# checknormality(data=mtcars,var="am") #All factors will automatically return 0
checknormality <- function(data,var,TRT=NULL) {
  data=droplevels(data)
  if (is.null(TRT)) {
    if (is.numeric(data[,var])) {
      #Cannot use Shapiro with less than 3 values, returns not normal:
      if (length(table(data[,var])) < 3) { return(0)}
      if (class(try(shapiro.test(data[,var]),silent=TRUE))=="try-error") {stop(paste0("Shapiro test did not work, check if ",var," is supposed to be numeric or check the data."))}
      if (shapiro.test(data[,var])$p < 0.05) {return(0)} else {return(1)}
    } else {return(0)} #if not numeric thus factor
  } else {
    if (!is.numeric(data[,var])) {return(0)} #if not numeric then quali thus not normal
    #Normality in each TRT group:
    data[,TRT]=as.factor(data[,TRT])
    for (i in 1:length(levels(data[,TRT]))) {
      if (sum(table(data[,var],data[,TRT])[,i])<3) {return(0)} #sum and not length because it is per group and thus has 0 for some values
      if (class(try(shapiro.test(data[,var][which(data[,TRT]==levels(data[,TRT])[i])]),silent=TRUE))=="try-error") {stop(paste0("Shapiro test did not work, check if ",var," is supposed to be numeric or check the data."))}
      if (shapiro.test(data[,var][which(data[,TRT]==levels(data[,TRT])[i])])$p < 0.05) {return(0)}
    }
    #Normality for the global pop
    if (class(try(shapiro.test(data[,var]),silent=TRUE))=="try-error") {stop(paste0("Shapiro test did not work, check if ",var," is supposed to be numeric or check the data."))}
    if (shapiro.test(data[,var])$p < 0.05) {return(0)}
    
    #If all TRT groups and the global pop. are all normal then no return was triggered and the var is normal
    #normal in each TRT group and the global pop thus 1 is returned
    return(1)
  }
}


#############






#Function to get table with getresults and names
#Can be used on its own to get means and freq for the global pop but no FlexTable is returned
#FlexTable = Table that can be exported to Word
##EX: #CTRL+SHIFT+C to comment/uncomment below
# data("mtcars")
# mtcars$am=as.factor(mtcars$am)
# mtcars$cyl=as.factor(mtcars$cyl)
# mtcars$vs=as.factor(mtcars$vs)
# vars=c("mpg","cyl","hp","vs")
# names=c("Miles/(US) gallon" , "Number of cylinders" ,"Gross horsepower" , "V/S")
# nolevels=vars
# mtcars$vs=relevel(mtcars$vs,ref="1")
# #nolevels means that any variable with only 2 levels will be on a single line and the first level will be used, thus relevel needed for some
# FTdesc(data=mtcars,vars=vars,names=names,nolevels=vars,virg=c(1,2,3))
# #virg can be specified: it's the amount of decimals used for the numeric variable
# #The order corresponds to the order of each numeric variable in vars ; 3rd value of virg here is not used
# #The \t means a tabulation in a word document, used for all levels of a factor
# #reference is TRUE by default, any other value removes the ", n (%)" after the label for categorical variables
FTdesc <- function(data,vars,names,virg=NULL,nolevels=NULL,normal=NULL,reference=TRUE)
{
  data=droplevels(data)
  #Setting virg to the correct length:
  if (is.null(virg)) {
    virg=rep(0,length(vars))
  } else {
    #If virg length = vars length then AT LEAST all numeric variables have a value of virg specified
    if (length(virg) != length(vars)) {
      virg2=c()
      compt=0
      for (i in 1:length(vars)) {
        #If numeric :
        if (vars[i] %in% vars[rapply(data[,vars],function(x) is.numeric(x))]) { #List of numeric variables
          compt=compt+1
          #If fewer decimals are specified than the number of numeric variables by default 0 decimals will be used
          if (is.na(virg[compt])) {virg2=c(virg2,0)} else {virg2=c(virg2,virg[compt])}
          #If virg is longer the additional values are ignored
        }
        #Else if factor then 0 is specified which is never used
        else {virg2=c(virg2,0)}
         }
      virg=virg2
      }
    }
  
  
  #If normal is empty then use checknormality function
  if (is.null(normal)) {
    normal=c()
    for (i in 1:length(vars)) {
      normal=c(normal,checknormality(data=data,var=vars[i]))
    }
  }
  
  tab=data.frame(A= character(), B= character())
  tab2=data.frame(A= character(), B= character())
  if (length(vars) != length(names)) {stop("Varlist and names are not the same length.")}
  for (i in 1:length(vars)) {
    if (is.factor(data[,vars[i]])) {
      # If in nolevels and length = 2 then 1 line: (OR 1 level meaning 100% answered by this 1 level; EX: only YES answers)
      if (vars[i] %in% nolevels && (length(levels(data[,vars[i]]))==2 || length(levels(data[,vars[i]]))==1) ) { 
        if (reference == TRUE) {tab2=data.frame(A=c(paste0(names[i],", n (%)")),B=getresult(data=data,var=vars[i],levelj=levels(data[,vars[i]])[1],virg=virg[i]))} else {tab2=data.frame(A=c(paste0(names[i],"")),B=getresult(data=data,var=vars[i],levelj=levels(data[,vars[i]])[1],virg=virg[i]))}
        tab=rbind(tab,tab2)
      } else {
        #Add exception, if length levels = 0 then run code above and return "-"
        if (all(is.na(data[,vars[i]]))) {
          if (reference == TRUE ) {tab2=data.frame(A=c(paste0(names[i],", n (%)")),B=getresult(data=data,var=vars[i],levelj=levels(data[,vars[i]])[1],virg=virg[i]))} else {tab2=data.frame(A=c(paste0(names[i],"")),B=getresult(data=data,var=vars[i],levelj=levels(data[,vars[i]])[1],virg=virg[i]))}
          tab=rbind(tab,tab2)
        } else {
          #Else multiple rows, one for each level
          if (reference == TRUE) {tab2=data.frame(A=c(paste0(names[i],", n (%)")),B="")} else {tab2=data.frame(A=c(paste0(names[i],"")),B="")}
          tab=rbind(tab,tab2)
          for (j in levels(data[,vars[i]])) {
            tab2=data.frame(A=c(paste0("\t",j)),B=getresult(data=data,var=vars[i],levelj=j,virg=virg[i])) 
            tab=rbind(tab,tab2)
          }
        }
        
      }
    } else {
      tab2=data.frame(A=c(names[i]),B=getresult(data=data,var=vars[i],virg=virg[i],normal=normal[i]))
      tab=rbind(tab,tab2)
    }
  }
  return(tab)
}








#Function to get the pvalue with the corresponding test
##EX: #CTRL+SHIFT+C to comment/uncomment below
# data("mtcars")
# mtcars$am=as.factor(mtcars$am)
# table(mtcars$am) #TRT group
# mtcars$cyl=as.factor(mtcars$cyl)
# comptest(data=mtcars,var="cyl",grp="am")
# table(mtcars$cyl,mtcars$am) #effectifs < 5 donc fisher
# chisq.test(mtcars$cyl,mtcars$am,correct=FALSE)
# fisher.test(mtcars$cyl,mtcars$am)
# comptest(data=mtcars,var="hp",grp="am",nrpval=5)
# qqnorm(mtcars$hp) #checknormality returns 0, not normal, it is possible to force it with normal option
# comptest(data=mtcars,var="hp",grp="am",normal=1,nrpval=3)
comptest <- function(data,var,grp,normal=NULL,nrpval) {
  data=droplevels(data)
  #if no nrpval is specified then 4 is used by default
  if (missing(nrpval)) {nrpval=4}
  
  #Sets the lower than value using nrpval
  lowervalch=paste0("0.",paste(rep(0,nrpval-1),collapse=""),"1")
  lowerval=as.numchar(lowervalch)
  lowervalch=paste0("<",lowervalch)
  
  
  #If normal is empty then use checknormality function
  if (is.null(normal)) {
    normal=checknormality(data=data,var=var,TRT=grp)
  }
  
  
  if (is.factor(data[,var])) {
    
    #If only values in 1 group then return without test:
    if (sum(table(data[,grp],data[,var])[1,])==0 || sum(table(data[,grp],data[,var])[2,])==0) {
      return("-")
    }
    #If only 1 level then return without test:
    if (length(levels(droplevels(data[,var])))==1) {return("-")}
    
    #Fisher is used if any of the expected values is less than 5
    #(To check for expected you must first do a chisq.test then X$expected, thus warnings)
    if (any(suppressWarnings(chisq.test(data[,var],data[,grp],correct=FALSE))$expected<5)) {
      #If more than 10 categories send a warning that the fisher test might take very long to run
      if (dim(suppressWarnings(chisq.test(data[,var],data[,grp],correct=FALSE)$resid))[1] > 20) {warning(paste0("Fisher test started with >20 categories, this test might take very long",". For variable: ",var))}
      #If an error with the fisher test then the pvalue will be simulated by repeating the simulation 100k times
      if (class(try(fisher.test(data[,var],data[,grp]),silent = TRUE))=="try-error") {
        if (fisher.test(data[,var],data[,grp],simulate.p.value=TRUE,B=1e5)$p.value < lowerval) {return(lowervalch)} else {return(formatC(fisher.test(data[,var],data[,grp],simulate.p.value=TRUE,B=1e5)$p.value, format='f', digits=nrpval ))}
      }
      if (fisher.test(data[,var],data[,grp])$p.value < lowerval) {return(lowervalch)} else {return(formatC(fisher.test(data[,var],data[,grp])$p.value, format='f', digits=nrpval ))}
    } else {
      if (chisq.test(data[,var],data[,grp],correct=FALSE)$p.value < lowerval) {return(lowervalch)} else {return(formatC(chisq.test(data[,var],data[,grp],correct=FALSE)$p.value, format='f', digits=nrpval ))}
    }
    
  } else {
    #If only values in 1 group then return without test:
    if (sum(table(data[,grp],data[,var])[1,])==0 || sum(table(data[,grp],data[,var])[2,])==0) {
      return("-")
    }
    
    #normal==0 means that it is not normal
    if (normal==0) {
      
      #If only 2 groups wilcoxon
      if (length(names(table(data[,grp]))) < 3) {
        #If wilcox pvalue is NA then return -
        #The exact=FALSE if for when there are ties (if none then pvalue is a little different but barely)
        if (is.na(suppressWarnings(wilcox.test(data[which(data[,grp]==names(table(data[,grp]))[1]),][,var],data[which(data[,grp]==names(table(data[,grp]))[2]),][,var]))$p.value)) {return("-")}
        
        if (suppressWarnings(wilcox.test(data[which(data[,grp]==names(table(data[,grp]))[1]),][,var],data[which(data[,grp]==names(table(data[,grp]))[2]),][,var]))$p.value < lowerval) {return(lowervalch)} else {return(formatC( suppressWarnings(wilcox.test(data[which(data[,grp]==names(table(data[,grp]))[1]),][,var],data[which(data[,grp]==names(table(data[,grp]))[2]),][,var]))$p.value, format='f', digits=nrpval ))}
        #Else Kruskall Wallis
      } else {
        if (kruskal.test(data[,var] ~ data[,grp])$p.value < lowerval) {return(lowervalch)} else {return(return(formatC(kruskal.test(data[,var] ~ data[,grp])$p.value, format='f', digits=nrpval )))}
      }
    } else {
      
      #If only 2 groups student
      if (length(names(table(data[,grp]))) < 3) {
        #Bartlett test to check for equal variances for the following t.test
        #Possible with leveneTest() from car library (Apparently: Bartlett better for small pop sizes but Levene better in general)
        if (bartlett.test(data[,var],data[,grp])$p.value < 0.05) {
          if (t.test(data[which(data[,grp]==names(table(data[,grp]))[1]),][,var],data[which(data[,grp]==names(table(data[,grp]))[2]),][,var])$p.value < lowerval) {return(lowervalch)} else {return(formatC( t.test(data[which(data[,grp]==names(table(data[,grp]))[1]),][,var],data[which(data[,grp]==names(table(data[,grp]))[2]),][,var])$p.value, format='f', digits=nrpval ))}
        } else {
          if (t.test(data[which(data[,grp]==names(table(data[,grp]))[1]),][,var],data[which(data[,grp]==names(table(data[,grp]))[2]),][,var],var.equal=TRUE)$p.value < lowerval) {return(lowervalch)} else {return(formatC( t.test(data[which(data[,grp]==names(table(data[,grp]))[1]),][,var],data[which(data[,grp]==names(table(data[,grp]))[2]),][,var],var.equal=TRUE)$p.value, format='f', digits=nrpval ))}
        }
        #Else ANOVA
      } else {
        #if not equal variances use a kruskal.test instead
        if (bartlett.test(data[,var],data[,grp])$p.value < 0.05) {if (kruskal.test(data[,var] ~ data[,grp])$p.value < lowerval) {return(lowervalch)} else {return(return(formatC(kruskal.test(data[,var] ~ data[,grp])$p.value, format='f', digits=nrpval )))}}
        if (summary(aov(data[,var] ~ data[,grp]))[[1]][["Pr(>F)"]][1] < lowerval) {return(lowervalch)} else {return(return(formatC(summary(aov(data[,var] ~ data[,grp]))[[1]][["Pr(>F)"]][1], format='f', digits=nrpval )))}
      }
    }
  }
}


















###### function to create the comparative table with a treatment and tests:

#vars=list of variable names
#names=list of labels for the variables
#TRT= comparative group that needs to be specified else only a desc table for the total pop
#virg = decimals, order is the same order as the numeric variables are in vars
#nolevels = vector of variables with 2 levels but transformed to 1 row (ex: Sexe, Homme Femme -> 1 Line Homme)
#normal=NULL by default, simple Shapiro else vector: 0 = not normal (or for factors) and 1 = normality (same length as vars)
#nrpval=number of decimals to be used for the pvalue
#boldpval=0.05 by default, puts all pvalues lower than set number in bold
#reference=TRUE by default, any other value will remove the ", n (%)" after the label for categorical variables
#NOTE: All variables need to be numeric or factor (TRT variable too!) and no dates

##EX similar to FTdesc function
##EX: #CTRL+SHIFT+C to comment/uncomment below
# data("mtcars")
# mtcars$am=as.factor(mtcars$am)
# mtcars$cyl=as.factor(mtcars$cyl)
# mtcars$vs=as.factor(mtcars$vs)
# vars=c("mpg","cyl","hp","vs")
# names=c("Miles/(US) gallon" , "Number of cylinders" ,"Gross horsepower" , "V/S")
# nolevels=vars
# mtcars$vs=relevel(mtcars$vs,ref="1")
# levels(mtcars$am)=c("Not Automatic","Automatic")
# #nolevels means that any variable with only 2 levels will be on a single line and the first level will be used, thus relevel needed
# TableDesc(data=mtcars,vars=vars,names=names,TRT="am",nolevels=vars,virg=c(2,1,33),nrpval=5)
# #virg can be specified: it's the amount of decimals used for the numeric variable
# #The order corresponds to the order of each numeric variable in vars
# #Simple version:
# TableDesc(data=mtcars,vars=vars,names=vars,TRT="am") #note that without nolevels the variable vs is on 3 lines
TableDesc = function(data,vars,names,TRT=NULL,virg=NULL,nolevels=NULL,normal=NULL,nrpval,boldpval=0.05,reference=TRUE) {
  data=droplevels(data)
  #If not data.frame but tibble:
  if (class(data)[1]!="data.frame") {
    tempcols=colnames(data)
    data=data.frame(data)
    colnames(data)=tempcols
  }
  #If TRT is not factor it is set to a factor:
  if (!is.factor(data[,TRT]) && !is.null(TRT)) {
    data[,TRT]=as.factor(data[,TRT])
    warning(paste0(TRT," was not set a factor, verify in the header that the reference levels are correct."))
  }
  
  #Stop if any vars are not in the data
  if (any(vars %ni% names(data))) {stop(paste0("Some variables are not in the dataset: ",paste(vars[which(vars %ni% names(data))],collapse=" ; ")))}
  
  #Stop if any vars is not numeric or factor
  for (i in vars) {if (class(data[,i]) %ni% c("factor","numeric")) {stop(paste0("Variable is not numeric or factor : ",i))}}
  
  #4 chiffres apres la virgule par defaut
  if (missing(nrpval)) {nrpval=4}
  
  #Setting virg to the correct length:
  if (is.null(virg)) {
    virg=rep(0,length(vars))
  } else {
    #If virg length = vars length then AT LEAST all numeric variables have a value of virg specified
    if (length(virg) != length(vars)) {
      virg2=c()
      compt=0
      for (i in 1:length(vars)) {
        #If numeric :
        if (vars[i] %in% vars[rapply(data[,vars],function(x) is.numeric(x))]) { #List of numeric variables
          compt=compt+1
          #If fewer decimals are specified than the number of numeric variables by default 0 decimals will be used
          if (is.na(virg[compt])) {virg2=c(virg2,0)} else {virg2=c(virg2,virg[compt])}
          #If virg is longer the additional values are ignored
        }
        #Else if factor then 0 is specified which is never used
        else {virg2=c(virg2,0)}
      }
      virg=virg2
    }
  }
  
  #If normal is empty then use checknormality function for each variable
  #OR if normal is not the same length as vars
  if (is.null(normal) || length(normal)!=length(vars)) {
    if (length(normal)!=length(vars) && !is.null(normal)) {warning("Length of vector 'normal' was not equal to length of 'vars'; checknormality function was used.")}
    normal=c()
    for (i in 1:length(vars)) {
      normal=c(normal,checknormality(data=data,var=vars[i],TRT=TRT))
    }
  }
  
  #If names is not given the variable names are used as labels
  if (is.null(names)) {names=vars}
  
  #Use the FTdesc function to get the description of the general population
  tab=FTdesc(data=data,vars=vars,names=names,virg=virg,nolevels=nolevels,normal=normal,reference=reference)
  tab  
  colnames(tab)=c("",paste("All","\n(n=",dim(data)[1],")",sep=""))
  
  #If no TRT is specified:
  if ((TRT %ni% colnames(data)) && !is.null(TRT) ) {TRT=NULL} 
  if (is.null(TRT)) {
    #Reporting with officer package
    Table=tab
    colnames(Table)[1]=" "
    Ftab=flextable(Table)
    Ftab=border_remove(Ftab)
    Ftab <- bold(Ftab, part = "header")
    font.name="Times New Roman"
    font.size=10
    Ftab <- flextable::font(Ftab,fontname=font.name,part ="all")
    Ftab <- fontsize(Ftab,size=font.size,part ="all")
    Ftab <- fontsize(Ftab,size=11,part ="header")
    Ftab <- align_text_col(Ftab, align = "left", header = TRUE)
    Ftab <- flextable::align(Ftab,i=1,j=2:ncol(Table),align = "center", part = "header")
    Ftab <- flextable::align(Ftab,i=1:nrow(Table),j=2:ncol(Table),align = "center", part = "body")
    Ftab <- hline(Ftab,border = fp_border(width = 1), part = "header" )
    Ftab <- hline_top(Ftab, border = fp_border(width = 1), part = "header" )
    Ftab <- hline_bottom(Ftab, border = fp_border(width = 1), part = "body")
    Ftab <- autofit(Ftab)
    Ftab <- height_all(Ftab, height=0.1, part = "body")
    Ftab <- height_all(Ftab, height=0.3, part = "header")
    Ftab <- height_all(Ftab, height=0.1, part = "footer")
    Ftab <- width(Ftab,j=1, width = 2.5)
    # Ftab <- width(Ftab,j=2:(ncol(Table)-1), width = 1.1) #no TRT thus only All column
    Ftab <- width(Ftab,j=ncol(Table), width = 1.1) #was 0.6
    
    boldrows=suppressWarnings(which(as.numchar(gsub("<","",Table[,ncol(Table)]))<boldpval))  #number of row to bold pvalue <0.05 by default
    for (rownr in boldrows)  {
      Ftab <- bold(Ftab,i=rownr,j=ncol(Table))
    }
    
    return(Ftab)
  }
  
  #Adding the values for the different groups
  for (trtgrp in names(table(data[,TRT]))) {
    temp=c()
    
    for (i in 1:length(vars)) {
      if (is.factor(data[,vars[i]])) {
        
        #If in nolevels and 2 levels (OR 1 level meaning 100% answered)
        if (vars[i] %in% nolevels && (length(levels(data[,vars[i]]))==2 || length(levels(data[,vars[i]]))==1)) {
          temp=c(temp,B=getresult(data=data[which(data[,TRT]==trtgrp),],var=vars[i],levelj=levels(data[,vars[i]])[1],virg=virg[i]))
        } else {
          #Add condition if only NA then return "-" 
          if (length(levels(data[which(data[,TRT]==trtgrp),][,vars[i]]))==0) {
          #getresult will return "-" if only NA
          temp=c(temp,B=getresult(data=data[which(data[,TRT]==trtgrp),],var=vars[i],levelj=levels(data[,vars[i]])[1],virg=virg[i]))
          } else {
            #Add empty line and add result for each level
            temp=c(temp,B="")
            for (j in levels(data[,vars[i]])) {
              temp=c(temp,B=getresult(data=data[which(data[,TRT]==trtgrp),],var=vars[i],levelj=j,virg=virg[i]))
            }
          }
          
        }
        
      } else {
        #If numeric then simple getresult function
        temp=c(temp,B=getresult(data=data[which(data[,TRT]==trtgrp),],var=vars[i],virg=virg[i],normal=normal[i]))
      }
    }
    
    tab=cbind(tab,temp)
  }
  
  
  
  temp=c()
  compt=0
  #Create progressbar
  pb <- txtProgressBar(min = 0, max = length(vars), style = 3)
  for (i in vars) {
    compt=compt+1
    if (is.factor(data[,i])) {
      
      #Check if any of the groups has only NA values, if so
      tempnotest=0
      for (j in 1:length(levels(data[,TRT]))) {
        if (all(is.na(data[,i][which(data[,TRT]==names(table(data[,TRT]))[j])]))) {tempnotest=1}
      }
      
      #If one group only has NA values return "-":
        if (tempnotest==0) {
          temp=c(temp,comptest(data=data,var=i,grp=TRT,nrpval=nrpval))
        } else { temp=c(temp,"-")}
          
      
      # If nolevels and 1 or 2 levels in total then 1 line:
      if (!(i %in% nolevels && (length(levels(data[,i]))==2 || length(levels(data[,i]))==1))) {
        for (j in levels(data[,i])) {
          temp=c(temp,"")
        }
      }
      
    } else {
      temp=c(temp,comptest(data=data,var=i,grp=TRT,normal=normal[compt],nrpval=nrpval))
    }
    
    # update progress bar
    setTxtProgressBar(pb, compt)
  }
  
  #Close progress bar
  close(pb)
  
  temp
  
  
  #Headers for the columns, to change the names you must directly change the level names of the TRT variable
  tab=cbind(tab,temp)
  headers=c()
  for (trtgrp in names(table(data[,TRT]))) {
    headers=c(headers,paste0(trtgrp,"\n(n=",dim(data[which(data[,TRT]==trtgrp),])[1],")"))
  }
  colnames(tab)=c("",paste("All","\n(n=",dim(data)[1],")",sep="")
                  ,headers,"p")
  
  
  
  #Reporting with officer package
  Table=tab
  colnames(Table)[1]=" "
  Ftab=flextable(Table)
  Ftab=border_remove(Ftab)
  Ftab <- bold(Ftab, part = "header")
  font.name="Times New Roman"
  font.size=10
  Ftab <- flextable::font(Ftab,fontname=font.name,part ="all")
  Ftab <- fontsize(Ftab,size=font.size,part ="all")
  Ftab <- fontsize(Ftab,size=11,part ="header")
  Ftab <- align_text_col(Ftab, align = "left", header = TRUE)
  Ftab <- flextable::align(Ftab,i=1,j=2:ncol(Table),align = "center", part = "header")
  Ftab <- flextable::align(Ftab,i=1:nrow(Table),j=2:ncol(Table),align = "center", part = "body")
  Ftab <- hline(Ftab,border = fp_border(width = 1), part = "header" )
  Ftab <- hline_top(Ftab, border = fp_border(width = 1), part = "header" )
  Ftab <- hline_bottom(Ftab, border = fp_border(width = 1), part = "body")
  Ftab <- autofit(Ftab)
  Ftab <- height_all(Ftab, height=0.1, part = "body")
  Ftab <- height_all(Ftab, height=0.3, part = "header")
  Ftab <- height_all(Ftab, height=0.1, part = "footer")
  Ftab <- width(Ftab,j=1, width = 2.5)
  Ftab <- width(Ftab,j=2:(ncol(Table)-1), width = 1.1)
  Ftab <- width(Ftab,j=ncol(Table), width = 0.6)
  
  boldrows=suppressWarnings(which(as.numchar(gsub("<","",Table[,ncol(Table)]))<boldpval))  #number of row to bold pvalue <0.05 by default
  for (rownr in boldrows)  {
    Ftab <- bold(Ftab,i=rownr,j=ncol(Table))
  }
  
  return(Ftab)
  
}



#Function to make a simple table with the amount of missing values per variable:
##EX: #CTRL+SHIFT+C to comment/uncomment below
# mtcars$hp[1:5]=NA
# tablemissing(data=mtcars,vars=c("hp","cyl"),names=vars)
# # #cyl has no missing values thus not added to the table
# # #Adding a TRT grouping variable will specificy the missing values per level of set variable
# mtcars$am=as.factor(mtcars$am)
# levels(mtcars$am)=c("Automatic","Non Automatic")
# tablemissing(data=mtcars,vars=c("hp","cyl"),names=vars,TRT="am")
tablemissing = function (data,vars,names=vars,TRT=NULL) {
  data=droplevels(data)
  #If TRT is not factor it is set to a factor:
  if (!is.factor(data[,TRT]) && !is.null(TRT)) {
    data[,TRT]=as.factor(data[,TRT])
    warning(paste0(TRT," was not set a factor, verify that the reference levels are correct."))
  }
  
  #Stop if any vars are not in the data
  if (any(vars %ni% names(data))) {stop(paste0("Some variables are not in the dataset: ",paste(vars[which(vars %ni% names(data))],collapse=" ; ")))}
  #Stop if any vars is not numeric or factor
  for (i in vars) {if (class(data[,i]) %ni% c("factor","numeric")) {stop(paste0("Variable is not numeric or factor : ",i))}}
  
  #If names is not given the variable names are used as labels
  if (is.null(names)) {names=vars}
  
  if (is.null(TRT)) {
    nnames=c()
    miss=c()
    for (i in 1:length(vars)) {
      if (any(is.na(data[,vars[i]]))) {
        nnames=c(nnames,names[i])
        miss=c(miss,length(which(is.na(data[,vars[i]])==TRUE)))
      }
    }
    if (is.null(nnames)) {stop("No missing data was found.")}
    df=cbind(nnames,miss)
    df=data.frame(df)
    colnames(df)=c("Variable","Missing")
    return(bold(hline(hline_top(hline_bottom(flextable::width(fontsize(flextable::font(theme_booktabs(flextable(df)), fontname = "Times New Roman", part = "all"), size = 10, part = "body"),width=dim_pretty(fontsize(flextable::font(theme_booktabs(flextable(df)), fontname = "Times New Roman", part = "all"), size = 10, part = "body"))$widths), border = fp_border(width = 1), part = "body"), border = fp_border(width = 1), part = "header" ),border = fp_border(width = 1), part = "header" ), part = "header"))
  } else {
    matrixvector=c()
    for (i in 1:length(vars)) {
      if (any(is.na(data[,vars[i]]))) {
        matrixvector=c(matrixvector,names[i],length(which(is.na(data[,vars[i]])==TRUE)))
        for (j in 1:length(levels(data[,TRT]))) {
          matrixvector=c(matrixvector,length(which(is.na(data[which(data[,TRT]==levels(data[,TRT])[j]),][,vars[i]])==TRUE)))
        }
      }
    }
    if (is.null(matrixvector)) {stop("No missing data was found.")}
    df=matrix(matrixvector,ncol=2+length(levels(data[,TRT])),byrow=TRUE)
    levelsn=c()
    for (lvl in levels(data[,TRT])) {levelsn=c(levelsn,paste0("\n(n=",dim(data[which(data[,TRT]==lvl),])[1],")"))}
    df=data.frame(df)
    colnames(df)=c("Variable",paste0("Missing \n Total","\n(n=",dim(data)[1],")"),paste0("Missing \n ",paste0(levels(data[,TRT]),levelsn)))
    return(width(width(bold(hline(hline_top(hline_bottom(flextable::width(fontsize(flextable::font(theme_booktabs(flextable(df)), fontname = "Times New Roman", part = "all"), size = 10, part = "body"),width=dim_pretty(fontsize(flextable::font(theme_booktabs(flextable(df)), fontname = "Times New Roman", part = "all"), size = 10, part = "body"))$widths), border = fp_border(width = 1), part = "body"), border = fp_border(width = 1), part = "header" ),border = fp_border(width = 1), part = "header" ), part = "header"),j=2:(2+length(levels(data[,TRT]))),width=1.1),j=1,width=3.1))
  }
}



###############################################################################
###############################################################################
###############################################################################
##EX: #CTRL+SHIFT+C to comment/uncomment below
# data("mtcars")
# mtcars$am=as.factor(mtcars$am)
# mtcars$cyl=as.factor(mtcars$cyl)
# mtcars$vs=as.factor(mtcars$vs)
# vars=c("cyl","hp","vs","am")
# names=c("Number of cylinders", "Gross horsepower", "V/S" ,"Automatic"  )
# mod=lm(mpg~cyl+hp+vs+am,data=mtcars)
# report.lm(mod,data=mtcars,labels=names,nrvirg=c(3,4,1,2),nrpval=3,reference=TRUE)
# report.lm(mod,data=dat,nrvirg=c(3,4),nrpval=5,SD=TRUE,reference=FALSE)
# #nrvirg= specify for each variable the number of decimals; if not specified 2 is used by default
# #SD=FALSE by default, it TRUE it will not use 95% CI but Beta +/- SD

report.lm=function (mod, data, labels,nrvirg=NA,nrpval,SD=FALSE,reference=TRUE){
  

  
  #Warning if there is a : in coeff names (interaction term) then might not compute
  if (any(sapply(names(mod$coefficients),function(x) grepl(":",x)))) {warning("Variable names contained ':', if the model contains an interaction term the function might not compute and if it does make sure to verify with summary(mod) that it is correct!")}
  
  #Stopping conditions
  if (missing(mod)) {
    stop("Function needs to have a model")
  }
  # Coef
  co <- mod$coef[-1]
  nvar <- length(co)
  varlist=attr(mod$terms, "term.labels")
  
  
  
  
  if (class(try(car::Anova(mod, test.statistic="Wald"),silent=TRUE))[1]=="try-error") {stop(paste0("Anova Wald statistics could not be used"))}
  ano=car::Anova(mod, test.statistic="Wald")
  ano=ano[-dim(ano)[1],]
  
  modxlevels=mod$xlevels
  
  
  #Par defaut 4 chiffres apres la virgule pour les pvalues
  if (missing(nrpval)) {
    nrpval=4
  }
  
  if (nrpval<2) {nrpval=2}
  
  if (missing(labels)) {labels=varlist}
  
  #labels not same as length(mod$coefficients) then use mod$coefficients instead
  if (length(varlist) != length(labels)) {
    warning("labels length was not the same as the number of variables in the model")
    labels=varlist
  }
  
  
  
  #Par defaut 2 chiffres apres la virgule pour les HR/Beta et CI
  Pnames=c()
  tempnrvirg=c()
  for (i in 1:length(varlist)) {
    
    if (varlist[i] %in% names(modxlevels)) {
      if (varlist[i] %in% rownames(ano)[which(ano$Df>1)]) {
        Pnames=c(Pnames,varlist[i])
        if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
        for (j in 1:(ano$Df[i]-1)) {
          Pnames=c(Pnames,varlist[i])
          if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
        }
      } else {
        Pnames=c(Pnames,varlist[i])
        if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
      }
    } else {
      Pnames=c(Pnames,varlist[i])
      if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
    }
  }
  Pnames
  nrvirg=tempnrvirg
  
  
  
  
  
  #pvalues
  pval = summary(mod)$coefficients[,4][-1]
  pvalf = format(round(pval, nrpval), nsmall = nrpval, scientific = F,
                 digits = nrpval)
  
  lowervalch=paste0("0.",paste(rep(0,nrpval-1),collapse=""),"1")
  lowerval=as.numchar(lowervalch)
  pvalf[pval < lowerval] = paste0("<",lowervalch)
  
  
  
  li <- format(round(confint(mod)[,1][-1], 10), nsmall = 10, scientific = F,digits = 10)
  li=gsub(" ", "", li, fixed = TRUE)
  ls <- format(round(confint(mod)[,2][-1], 10), nsmall = 10, scientific = F,digits = 10)
  ls=gsub(" ", "", ls, fixed = TRUE)
  li2=as.numchar(li)
  ls2=as.numchar(ls)
  co2=as.numchar(co)
  std2=as.numchar(summary(mod)$coef[,2][-1])
  li=c()
  ls=c()
  co=c()
  std=c()
  for (i in 1:length(li2)) {
    li=c(li,as.character(format(round(li2[i], nrvirg[i]), nsmall = nrvirg[i], scientific = F,digits = nrvirg[i])))
    ls=c(ls,as.character(format(round(ls2[i], nrvirg[i]), nsmall = nrvirg[i], scientific = F,digits = nrvirg[i])))
    co=c(co,as.character(format(round(co2[i], nrvirg[i]), nsmall = nrvirg[i], scientific = F,digits = nrvirg[i])))
    std=c(std,as.character(format(round(std2[i], nrvirg[i]), nsmall = nrvirg[i], scientific = F,digits = nrvirg[i])))
  }
  li
  ls
  co
  ci = paste0("[", li, ";", ls, "]")
  if (SD==TRUE) {co=paste0(co," \xb1 ",std)}
  
  
  
  Pdf=data.frame(Pnames,Beta = co, '95% CI' = ci,
                 'Pval' = pvalf,check.names=FALSE)
  
  
  
  Cnames=c()
  for (i in 1:length(varlist)) {
    
    if (varlist[i] %in% names(modxlevels)) {
      if (varlist[i] %in% rownames(ano)[which(ano$Df>1)]) {
        if (reference==TRUE) {Cnames=c(Cnames,paste0(labels[i]," (ref. ",modxlevels[varlist[i]][[1]][1],")"))} else {Cnames=c(Cnames,paste0(labels[i]))}
        
        for (j in 1:(ano$Df[i])) {
          Cnames=c(Cnames,paste0("\t",modxlevels[varlist[i]][[1]][j+1]))
        }
      } else {
        if (reference==TRUE) {Cnames=c(Cnames,paste0(labels[i]," (ref. ",modxlevels[varlist[i]][[1]][1],")"))} else {Cnames=c(Cnames,paste0(labels[i]))}
        
      }
    } else {
      Cnames=c(Cnames,labels[i])
    }
  }
  Cnames
  
  Variables=c()
  Beta=c()
  CI=c()
  Pvalue=c()
  compt=0
  compthr=0
  
  for (i in unique(Pdf$Pnames)) {
    compt=compt+1
    compthr=compthr+1
    
    if (i %in% rownames(ano)[which(ano$Df>1)]) {
      Variables=c(Variables,Cnames[compt])
      Beta=c(Beta,"")
      CI=c(CI,"")
      
      temppval=format(round(ano[,4][which(rownames(ano)==i)], nrpval), nsmall = nrpval, scientific = F,
                      digits = nrpval)
      if (as.numchar(temppval)< lowerval) {temppval = paste0("<",lowervalch)}
      
      Pvalue=c(Pvalue,temppval)
      
      
      for (j in 1:(ano$Df[which(rownames(ano)==i)])) {
        compt=compt+1
        if (j!=1) {compthr=compthr+1}
        Variables=c(Variables,Cnames[compt])
        Beta=c(Beta,co[compthr])
        CI=c(CI,ci[compthr])
        Pvalue=c(Pvalue,"")
      }
      
    } else {
      Variables=c(Variables,Cnames[compt])
      Beta=c(Beta,co[compthr])
      CI=c(CI,ci[compthr])
      Pvalue=c(Pvalue,pvalf[compthr])
    }
  }
  
  
  # Final data frame with all information
  results = data.frame(Variables, Beta = Beta, '95% CI' = CI,
                       'p value' = Pvalue,check.names=FALSE)
  rownames(results) = NULL
  if (SD==TRUE) {
    results=results[,-3]
    colnames(results)[2]="Beta \xb1 SD"
  }
  
  Table=results
  colnames(Table)[1]=" "
  Ftab=flextable(Table)
  Ftab=border_remove(Ftab)
  Ftab <- bold(Ftab, part = "header")
  font.name="Times New Roman"
  font.size=10
  Ftab <- flextable::font(Ftab,fontname=font.name,part ="all")
  Ftab <- fontsize(Ftab,size=font.size,part ="all")
  Ftab <- fontsize(Ftab,size=11,part ="header")
  Ftab <- align_text_col(Ftab, align = "left", header = TRUE)
  Ftab <- flextable::align(Ftab,i=1,j=2:ncol(Table),align = "center", part = "header")
  Ftab <- flextable::align(Ftab,i=1:nrow(Table),j=2:ncol(Table),align = "center", part = "body")
  Ftab <- hline(Ftab,border = fp_border(width = 1), part = "header" )
  Ftab <- hline_top(Ftab, border = fp_border(width = 1), part = "header" )
  Ftab <- hline_bottom(Ftab, border = fp_border(width = 1), part = "body")
  Ftab <- autofit(Ftab)
  Ftab <- height_all(Ftab, height=0.1, part = "body")
  Ftab <- height_all(Ftab, height=0.3, part = "header")
  Ftab <- height_all(Ftab, height=0.1, part = "footer")
  Ftab <- width(Ftab,j=1, width = 2.5)
  Ftab <- width(Ftab,j=2:(ncol(Table)-1), width = 1.1)
  Ftab <- width(Ftab,j=ncol(Table), width = 0.6)
  
  boldrows=suppressWarnings(which(as.numchar(gsub("<","",Table[,ncol(Table)]))<0.05))  #number of row to bold pvalue <0.05 by default
  for (rownr in boldrows)  {
    Ftab <- bold(Ftab,i=rownr,j=ncol(Table))
      }
  
  return(Ftab)
}
###############################################################################






###############################################################################
###############################################################################
###############################################################################
##### Function to get a sas dataset into a useable R dataset for DESC functions 
#sastor does the following:
# -Return a list of 3 elements; the dataset, the list of variables and the labels of all the variables
# -removes certain character strings from the variable labels if specified
# -set most variables to the correct numeric/factor format (not 100% accuracy)
# NOT CORRECT THUS COMMENTED:-Those with only 1 level get their NA values set to 0 (meaning only YES responses for example)
# -All those with the levels NON/OUI will be releveled to OUI as reference level
# -If levels contain Oui and Non (+ others) then Oui will be put al reference
# -droplevels on all factors
##EX:
# #remlabel allows to subtract a character string from all the label names of the variables
# listd=sastor(bilanbio,remlabel=c("Visite d'inclusion M0 : ","Visite de Suivi M1 : "))
# data=listd[[1]] #sastor return list of 3 elements; first the dataset
# vars=listd[[2]] #then the list of variables
# names=listd[[3]] #then the labels for all the variables
# 
# #Sometimes not all are correctly set to numeric/factor:
# data$HEMHG=as.numchar(data$HEMHG)
# data$HEMLEU=as.numchar(data$HEMLEU)
# data$HEMEOS=as.numchar(data$HEMEOS)
# data$BIOKA=as.numchar(data$BIOKA)
# data$MYCIGET=as.numchar(data$MYCIGET)
# data$MYCIGES=as.numchar(data$MYCIGES)
# data$UAIGG=as.numchar(data$UAIGG)
# 
# #Never use the SUBJID and the TRT if they are used afterwards
# #Also make sure the TRT variable is a factor
# names=names[-which(colnames(data) %in% c("SUBJID","TRT"))]
# vars=colnames(data)[-which(colnames(data) %in% c("SUBJID","TRT"))]
# TableDesc(data=data,vars=vars,names=names,TRT="TRT")
# FTable3=TableDesc(data=data,vars=vars,names=names,TRT="TRT",nolevels=vars)
# FTable3
sastor=function (data,remlabel=NULL){
  
#Before setting to data.frame get list of all labels:
  names=lapply(data,function(x) attr(x, "label"))
  
#Pour les variables sans labels on remet les nom des variables
  names[which(lapply(names,function(x) length(x))>1)]=names(names[which(lapply(names,function(x) length(x))>1)])
  names[which(lapply(names,function(x) is.null(x))==TRUE)]=names(names[which(lapply(names,function(x) is.null(x))==TRUE)])
  names=unlist(names,use.names=FALSE)
  
#Enlever du label de SAS la visite
  if (!is.null(remlabel)) {
    for (i in remlabel) {
    names=gsub(i,"",names)
    }
  }
  
  chars=lapply(data,function(x) attr(x, "format"))
  
#If format is null then numeric else as factor:
  #Remove format==TIME
  chars[which(lapply(chars,function(x) (x=="BEST" || is.null(x)))!=TRUE)]="as_factor"
  chars[which(lapply(chars,function(x) (x=="BEST" || is.null(x)))==TRUE)]="as.numchar"
  chars=unlist(chars,use.names=FALSE)
  
#set all the variables to as.factor or as.numeric
  data=data.frame(data)
  vars=colnames(data)
  for (i in 1:length(chars)) {
    if (!inherits(data[,i],"Date")) {data[,i]=eval(parse(text=paste0(chars[i],"(","data[,'",vars[i],"']",")")))}
    
  }
# #those with only 1 level get all their NA values set to 0
#   tempset=data[,which(rapply(data,function(x) length(levels(x)))==1)]
#   for (i in colnames(tempset)) {
#     data[,i]=as.character(data[,i])
#     data[,i][which(is.na(data[,i]))]=0
#     data[,i]=as.factor(data[,i])
#   }
# #Les variables avec seulement 1 réponse qui restent: (que des 0)
#   colnames(data[,which(rapply(data,function(x) length(levels(x)))==1)])
  
#Relevel tous les facteurs a 2 niveaux Oui et Non dans ref=Oui and 1 and 0 to ref 1
  tempset=data[,which(rapply(data,function(x) length(levels(x)))==2)]
  for (i in colnames(tempset)) {
    if (all(levels(data[,i])==c("Non","Oui"))) {
      data[,i]=relevel(data[,i],ref="Oui")
    }
    if (all(levels(data[,i])==c("NON","OUI"))) {
      data[,i]=relevel(data[,i],ref="OUI")
    }
    if (all(levels(data[,i])==c("0","1"))) {
      data[,i]=relevel(data[,i],ref="1")
    }
  }
  #Relevel if levels contain Oui or Non (makes the above useless)
  for (i in colnames(data)) {
    if (all(c("Non","Oui") %in% levels(data[,i]))) {
      data[,i]=relevel(data[,i],ref="Oui")
    }
    if (all(c("NON","OUI") %in% levels(data[,i]))) {
      data[,i]=relevel(data[,i],ref="OUI")
    }
    if (all(c("NO","YES") %in% levels(data[,i]))) {
      data[,i]=relevel(data[,i],ref="YES")
    }
    if (all(c("No","Yes") %in% levels(data[,i]))) {
      data[,i]=relevel(data[,i],ref="Yes")
    }
  }
#Droplevels for all factors:
  tempset=data[,which(rapply(data,function(x) is.factor(x)))]
  for (i in colnames(tempset)) {
    data[,i]=droplevels(data[,i])
  }
  
  return(list(data,vars,names))
}
###############################################################################























##### Univariate results of a Cox model:
# llibrary(survival)
# data(cancer)
# cancer$sex=as.factor(cancer$sex)
# cancer$ph.ecog=as.factor(cancer$ph.ecog)
# report.coxunivariate(data=cancer,vars=c("age","ph.ecog","ph.karno","wt.loss"),
#     names=c("Age","Ph Ecog","Ph Karno","Weight loss"),time="time",status="status",TRT="sex")

report.coxunivariate <- function(data,vars,names,time,status,nrpval=4,TRT=NULL)
{
  if (missing(time) || missing(status)) {stop("Both time and status need to be specified")}
  
  
  surv=paste0("Surv(",time,",",status,")")
  
  lowervalch=paste0("0.",paste(rep(0,nrpval-1),collapse=""),"1")
  lowerval=as.numchar(lowervalch)
  lowervalch=paste0("<",lowervalch)
  
  if (!is.null(TRT)) {
    treatmentp=paste0(TRT,"+")
  } else {treatmentp=NULL}
  
  cp=c()
  czph=c()
  for (loop.var in vars) {
    
    
    #invisible(capture.output) is used here to prevent any prints from the model in the console
    invisible(capture.output(mod <- coxph(formula(paste(surv,"~",treatmentp,paste(loop.var,collapse="+")))
                                          ,data=data)))
    
    #Either use the global pvalue (Anova or anova for stratified models) or the pvalue from the model
    if (class(try(car::Anova(mod, test.statistic="Wald"),silent=TRUE))[1]=="try-error") {stop(paste0("Anova Wald statistics could not be used for var: ",loop.var))}
    if (any(is.na(mod$coefficients))) {stop(paste0("Convergence error for variable: ",loop.var))}
    ano <- car::Anova(mod, test.statistic="Wald")
    anoDf <- ano$Df
    if (anoDf[which(rownames(ano)==loop.var)]>1) {
      pval <- ano$'Pr(>Chisq)'[which(rownames(ano)==loop.var)]
      if (pval < lowerval) {pval=lowervalch} else {pval=formatC(pval, format='f', digits=nrpval )}
      zph=cox.zph(mod)$table[which(rownames(cox.zph(mod)$table)=="GLOBAL"),3]
      if (zph < lowerval) {zph=lowervalch} else {zph=formatC(zph, format='f', digits=nrpval )}
    } else {
      level <- 0.95
      z <- abs(qnorm((1 - level)/2))
      se <-sqrt(diag(mod$var))
      pval <- 1 -pchisq((mod$coef/se)^2,1)
      pval=pval[which(grepl( loop.var,names(mod$coeff)))]
      if (pval < lowerval) {pval=lowervalch} else {pval=formatC(pval, format='f', digits=nrpval )}
      zph=cox.zph(mod)$table[which(grepl(loop.var,rownames(cox.zph(mod)$table))),3]
      if (zph < lowerval) {zph=lowervalch} else {zph=formatC(zph, format='f', digits=nrpval )}
    }
    cp=c(cp,pval)
    czph=c(czph,zph)
  }
  

  
  #Keep only var with pvalue below threshold, default = 0.2
  return(data.frame(vars,cp,czph))
}
#####################################################################








###############################################################################
###############################################################################
###############################################################################
##EX: #CTRL+SHIFT+C to comment/uncomment below
# data("Arrests")
# Arrests$employed=as.factor(Arrests$employed)
# Arrests$sex=as.factor(Arrests$sex)
# Arrests$checks=as.factor(Arrests$checks)
# levels(Arrests$checks)=c("0","1","2","3 ou +","3 ou +","3 ou +","3 ou +")
# Arrests$age=as.numeric(as.character(Arrests$age))
# names=c("Age", "Sexe","Nr. of checks")
# mod=glm(employed~age+sex+checks,data=Arrests,family="binomial")
# report.logi(mod,data=Arrests,labels=names,nrvirg=c(2,4,1,2),nrpval=3,reference=TRUE)

report.logi=function (mod, data, labels,nrvirg=NA,nrpval,reference=TRUE){
  #Warning if there is a : in coeff names (interaction term) then might not compute
  if (any(sapply(names(mod$coefficients),function(x) grepl(":",x)))) {warning("Variable names contained ':', if the model contains an interaction term the function might not compute and if it does make sure to verify with summary(mod) that it is correct!")}
  
  #Stopping conditions
  if (missing(mod)) {
    stop("Function needs to have a model")
  }
  # Coef
  co <- mod$coef[-1]
  nvar <- length(co)
  varlist=attr(mod$terms, "term.labels")
  
  #Exp estimate for OR
  co <- exp(co)
  
  if (class(try(car::Anova(mod, test.statistic="Wald"),silent=TRUE))[1]=="try-error") {stop(paste0("Anova Wald statistics could not be used"))}
  ano=car::Anova(mod, test.statistic="Wald")

  modxlevels=mod$xlevels
  
  
  #Par defaut 4 chiffres apres la virgule pour les pvalues
  if (missing(nrpval)) {
    nrpval=4
  }
  
  if (nrpval<2) {nrpval=2}
  
  if (missing(labels)) {labels=varlist}
  
  #labels not same as length(mod$coefficients) then use mod$coefficients instead
  if (length(varlist) != length(labels)) {
    warning("labels length was not the same as the number of variables in the model")
    labels=varlist
  }
  
  
  
  
  
  #Par defaut 2 chiffres apres la virgule pour les HR/Beta et CI
  Pnames=c()
  tempnrvirg=c()
  for (i in 1:length(varlist)) {
    
    if (varlist[i] %in% names(modxlevels)) {
      if (varlist[i] %in% rownames(ano)[which(ano$Df>1)]) {
        Pnames=c(Pnames,varlist[i])
        if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
        for (j in 1:(ano$Df[i]-1)) {
          Pnames=c(Pnames,varlist[i])
          if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
        }
      } else {
        Pnames=c(Pnames,varlist[i])
        if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
      }
    } else {
      Pnames=c(Pnames,varlist[i])
      if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
    }
  }
  Pnames
  nrvirg=tempnrvirg
  
  
  
  
  
  #pvalues
  pval = summary(mod)$coefficients[,4][-1]
  pvalf = format(round(pval, nrpval), nsmall = nrpval, scientific = F,
                 digits = nrpval)
  
  lowervalch=paste0("0.",paste(rep(0,nrpval-1),collapse=""),"1")
  lowerval=as.numchar(lowervalch)
  pvalf[pval < lowerval] = paste0("<",lowervalch)
  
  
  
  li <- format(round(exp(confint.default(mod))[,1][-1], 10), nsmall = 10, scientific = F,digits = 10)
  li=gsub(" ", "", li, fixed = TRUE)
  ls <- format(round(exp(confint.default(mod))[,2][-1], 10), nsmall = 10, scientific = F,digits = 10)
  ls=gsub(" ", "", ls, fixed = TRUE)
  li2=as.numchar(li)
  ls2=as.numchar(ls)
  co2=as.numchar(co)
  std2=as.numchar(summary(mod)$coef[,2][-1])
  li=c()
  ls=c()
  co=c()
  std=c()
  for (i in 1:length(li2)) {
    li=c(li,as.character(format(round(li2[i], nrvirg[i]), nsmall = nrvirg[i], scientific = F,digits = nrvirg[i])))
    ls=c(ls,as.character(format(round(ls2[i], nrvirg[i]), nsmall = nrvirg[i], scientific = F,digits = nrvirg[i])))
    co=c(co,as.character(format(round(co2[i], nrvirg[i]), nsmall = nrvirg[i], scientific = F,digits = nrvirg[i])))
  }
  li
  ls
  co
  ci = paste0("[", li, ";", ls, "]")
  
  
  
  Pdf=data.frame(Pnames,Beta = co, '95% CI' = ci,
                 'Pval' = pvalf,check.names=FALSE)
  
  
  
  Cnames=c()
  for (i in 1:length(varlist)) {
    
    if (varlist[i] %in% names(modxlevels)) {
      if (varlist[i] %in% rownames(ano)[which(ano$Df>1)]) {
        if (reference==TRUE) {Cnames=c(Cnames,paste0(labels[i]," (ref. ",modxlevels[varlist[i]][[1]][1],")"))} else {Cnames=c(Cnames,paste0(labels[i]))}
        for (j in 1:(ano$Df[i])) {
          Cnames=c(Cnames,paste0("\t",modxlevels[varlist[i]][[1]][j+1]))
        }
      } else {
        if (reference==TRUE) {Cnames=c(Cnames,paste0(labels[i]," (ref. ",modxlevels[varlist[i]][[1]][1],")"))} else {Cnames=c(Cnames,paste0(labels[i]))}
      }
    } else {
      Cnames=c(Cnames,labels[i])
    }
  }
  Cnames
  
  Variables=c()
  Beta=c()
  CI=c()
  Pvalue=c()
  compt=0
  compthr=0
  
  for (i in unique(Pdf$Pnames)) {
    compt=compt+1
    compthr=compthr+1
    
    if (i %in% rownames(ano)[which(ano$Df>1)]) {
      Variables=c(Variables,Cnames[compt])
      Beta=c(Beta,"")
      CI=c(CI,"")
      
      temppval=format(round(ano[,3][which(rownames(ano)==i)], nrpval), nsmall = nrpval, scientific = F,
                      digits = nrpval)
      if (as.numchar(temppval)< lowerval) {temppval = paste0("<",lowervalch)}
      
      Pvalue=c(Pvalue,temppval)
      
      
      for (j in 1:(ano$Df[which(rownames(ano)==i)])) {
        compt=compt+1
        if (j!=1) {compthr=compthr+1}
        Variables=c(Variables,Cnames[compt])
        Beta=c(Beta,co[compthr])
        CI=c(CI,ci[compthr])
        Pvalue=c(Pvalue,"")
      }
      
    } else {
      Variables=c(Variables,Cnames[compt])
      Beta=c(Beta,co[compthr])
      CI=c(CI,ci[compthr])
      Pvalue=c(Pvalue,pvalf[compthr])
    }
  }
  
  
  # Final data frame with all information
  results = data.frame(Variables, OR = Beta, '95% CI' = CI,
                       'p value' = Pvalue,check.names=FALSE)
  rownames(results) = NULL
  
  Table=results
  colnames(Table)[1]=" "
  Ftab=flextable(Table)
  Ftab=border_remove(Ftab)
  Ftab <- bold(Ftab, part = "header")
  font.name="Times New Roman"
  font.size=10
  Ftab <- flextable::font(Ftab,fontname=font.name,part ="all")
  Ftab <- fontsize(Ftab,size=font.size,part ="all")
  Ftab <- fontsize(Ftab,size=11,part ="header")
  Ftab <- align_text_col(Ftab, align = "left", header = TRUE)
  Ftab <- flextable::align(Ftab,i=1,j=2:ncol(Table),align = "center", part = "header")
  Ftab <- flextable::align(Ftab,i=1:nrow(Table),j=2:ncol(Table),align = "center", part = "body")
  Ftab <- hline(Ftab,border = fp_border(width = 1), part = "header" )
  Ftab <- hline_top(Ftab, border = fp_border(width = 1), part = "header" )
  Ftab <- hline_bottom(Ftab, border = fp_border(width = 1), part = "body")
  Ftab <- autofit(Ftab)
  Ftab <- height_all(Ftab, height=0.1, part = "body")
  Ftab <- height_all(Ftab, height=0.3, part = "header")
  Ftab <- height_all(Ftab, height=0.1, part = "footer")
  Ftab <- width(Ftab,j=1, width = 2.5)
  Ftab <- width(Ftab,j=2:(ncol(Table)-1), width = 1.1)
  Ftab <- width(Ftab,j=ncol(Table), width = 0.6)
  
  boldrows=suppressWarnings(which(as.numchar(gsub("<","",Table[,ncol(Table)]))<0.05))  #number of row to bold pvalue <0.05 by default
  for (rownr in boldrows)  {
    Ftab <- bold(Ftab,i=rownr,j=ncol(Table))
  }
  
  return(Ftab)
}
###############################################################################







###############################################################################
###############################################################################
###############################################################################
##EX: #CTRL+SHIFT+C to comment/uncomment below
# library(survival)
# data("cancer")
# cancer$ph.ecog=as.factor(cancer$ph.ecog)
# cancer$sex=as.factor(cancer$sex)
# cancer$age=as.numeric(as.character(cancer$age))
# names=c("Age (année)","Sexe", "ph.ecog")
# mod=coxph(Surv(time,status)~age+sex+ph.ecog,data=cancer)
# report.cox(mod,data=cancer,labels=names,nrvirg=c(3,4,1,2),nrpval=3,reference=TRUE)

report.cox=function (mod, data, labels,nrvirg=NA,nrpval,reference=TRUE){
  
  #Warning if there is a : in coeff names (interaction term) then might not compute
  if (any(sapply(names(mod$coefficients),function(x) grepl(":",x)))) {warning("Variable names contained ':', if the model contains an interaction term the function might not compute and if it does make sure to verify with summary(mod) that it is correct!")}
  
  #Stopping conditions
  if (missing(mod)) {
    stop("Function needs to have a model")
  }
  # Coef
  co <- mod$coef
  nvar <- length(co)
  varlist=attr(mod$terms, "term.labels")
  
  #Exp estimate for OR
  co <- exp(co)
  
  if (class(try(car::Anova(mod, test.statistic="Wald"),silent=TRUE))[1]=="try-error") {stop(paste0("Anova Wald statistics could not be used"))}
  ano=car::Anova(mod, test.statistic="Wald")
  
  modxlevels=mod$xlevels
  
  
  #Par defaut 4 chiffres apres la virgule pour les pvalues
  if (missing(nrpval)) {
    nrpval=4
  }
  
  if (nrpval<2) {nrpval=2}
  
  if (missing(labels)) {labels=varlist}
  
  #labels not same as length(mod$coefficients) then use mod$coefficients instead
  if (length(varlist) != length(labels)) {
    warning("labels length was not the same as the number of variables in the model")
    labels=varlist
  }
  
  
  
  
  
  #Par defaut 2 chiffres apres la virgule pour les HR/Beta et CI
  Pnames=c()
  tempnrvirg=c()
  for (i in 1:length(varlist)) {
    
    if (varlist[i] %in% names(modxlevels)) {
      if (varlist[i] %in% rownames(ano)[which(ano$Df>1)]) {
        Pnames=c(Pnames,varlist[i])
        if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
        for (j in 1:(ano$Df[i]-1)) {
          Pnames=c(Pnames,varlist[i])
          if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
        }
      } else {
        Pnames=c(Pnames,varlist[i])
        if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
      }
    } else {
      Pnames=c(Pnames,varlist[i])
      if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
    }
  }
  Pnames
  nrvirg=tempnrvirg
  
  
  
  
  
  #pvalues
  pval = summary(mod)$coefficients[,5]
  pvalf = format(round(pval, nrpval), nsmall = nrpval, scientific = F,
                 digits = nrpval)
  
  lowervalch=paste0("0.",paste(rep(0,nrpval-1),collapse=""),"1")
  lowerval=as.numchar(lowervalch)
  pvalf[pval < lowerval] = paste0("<",lowervalch)
  
  
  
  li <- format(round(exp(confint(mod))[,1], 10), nsmall = 10, scientific = F,digits = 10)
  li=gsub(" ", "", li, fixed = TRUE)
  ls <- format(round(exp(confint(mod))[,2], 10), nsmall = 10, scientific = F,digits = 10)
  ls=gsub(" ", "", ls, fixed = TRUE)
  li2=as.numchar(li)
  ls2=as.numchar(ls)
  co2=as.numchar(co)
  std2=as.numchar(summary(mod)$coef[,2])
  li=c()
  ls=c()
  co=c()
  for (i in 1:length(li2)) {
    li=c(li,as.character(format(round(li2[i], nrvirg[i]), nsmall = nrvirg[i], scientific = F,digits = nrvirg[i])))
    ls=c(ls,as.character(format(round(ls2[i], nrvirg[i]), nsmall = nrvirg[i], scientific = F,digits = nrvirg[i])))
    co=c(co,as.character(format(round(co2[i], nrvirg[i]), nsmall = nrvirg[i], scientific = F,digits = nrvirg[i])))
  }
  li
  ls
  co
  ci = paste0("[", li, ";", ls, "]")
  
  
  
  Pdf=data.frame(Pnames,Beta = co, '95% CI' = ci,
                 'Pval' = pvalf,check.names=FALSE)
  
  
  
  Cnames=c()
  for (i in 1:length(varlist)) {
    
    if (varlist[i] %in% names(modxlevels)) {
      if (varlist[i] %in% rownames(ano)[which(ano$Df>1)]) {
        if (reference==TRUE) {Cnames=c(Cnames,paste0(labels[i]," (ref. ",modxlevels[varlist[i]][[1]][1],")"))} else {Cnames=c(Cnames,paste0(labels[i]))}
        for (j in 1:(ano$Df[i])) {
          Cnames=c(Cnames,paste0("\t",modxlevels[varlist[i]][[1]][j+1]))
        }
      } else {
        if (reference==TRUE) {Cnames=c(Cnames,paste0(labels[i]," (ref. ",modxlevels[varlist[i]][[1]][1],")"))} else {Cnames=c(Cnames,paste0(labels[i]))}
      }
    } else {
      Cnames=c(Cnames,labels[i])
    }
  }
  Cnames
  
  Variables=c()
  Beta=c()
  CI=c()
  Pvalue=c()
  compt=0
  compthr=0
  
  for (i in unique(Pdf$Pnames)) {
    compt=compt+1
    compthr=compthr+1
    
    if (i %in% rownames(ano)[which(ano$Df>1)]) {
      Variables=c(Variables,Cnames[compt])
      Beta=c(Beta,"")
      CI=c(CI,"")
      
      temppval=format(round(ano[,3][which(rownames(ano)==i)], nrpval), nsmall = nrpval, scientific = F,
                      digits = nrpval)
      if (as.numchar(temppval)< lowerval) {temppval = paste0("<",lowervalch)}
      
      Pvalue=c(Pvalue,temppval)
      
      
      for (j in 1:(ano$Df[which(rownames(ano)==i)])) {
        compt=compt+1
        if (j!=1) {compthr=compthr+1}
        Variables=c(Variables,Cnames[compt])
        Beta=c(Beta,co[compthr])
        CI=c(CI,ci[compthr])
        Pvalue=c(Pvalue,"")
      }
      
    } else {
      Variables=c(Variables,Cnames[compt])
      Beta=c(Beta,co[compthr])
      CI=c(CI,ci[compthr])
      Pvalue=c(Pvalue,pvalf[compthr])
    }
  }
  
  
  # Final data frame with all information
  results = data.frame(Variables, HR = Beta, '95% CI' = CI,
                       'p value' = Pvalue,check.names=FALSE)
  rownames(results) = NULL
  
  Table=results
  colnames(Table)[1]=" "
  Ftab=flextable(Table)
  Ftab=border_remove(Ftab)
  Ftab <- bold(Ftab, part = "header")
  font.name="Times New Roman"
  font.size=10
  Ftab <- flextable::font(Ftab,fontname=font.name,part ="all")
  Ftab <- fontsize(Ftab,size=font.size,part ="all")
  Ftab <- fontsize(Ftab,size=11,part ="header")
  Ftab <- align_text_col(Ftab, align = "left", header = TRUE)
  Ftab <- flextable::align(Ftab,i=1,j=2:ncol(Table),align = "center", part = "header")
  Ftab <- flextable::align(Ftab,i=1:nrow(Table),j=2:ncol(Table),align = "center", part = "body")
  Ftab <- hline(Ftab,border = fp_border(width = 1), part = "header" )
  Ftab <- hline_top(Ftab, border = fp_border(width = 1), part = "header" )
  Ftab <- hline_bottom(Ftab, border = fp_border(width = 1), part = "body")
  Ftab <- autofit(Ftab)
  Ftab <- height_all(Ftab, height=0.1, part = "body")
  Ftab <- height_all(Ftab, height=0.3, part = "header")
  Ftab <- height_all(Ftab, height=0.1, part = "footer")
  Ftab <- width(Ftab,j=1, width = 2.5)
  Ftab <- width(Ftab,j=2:(ncol(Table)-1), width = 1.1)
  Ftab <- width(Ftab,j=ncol(Table), width = 0.6)
  
  boldrows=suppressWarnings(which(as.numchar(gsub("<","",Table[,ncol(Table)]))<0.05))  #number of row to bold pvalue <0.05 by default
  for (rownr in boldrows)  {
    Ftab <- bold(Ftab,i=rownr,j=ncol(Table))
  }
  
  return(Ftab)
}
###############################################################################












###############################################################################
###############################################################################
###############################################################################
# EX: #CTRL+SHIFT+C to comment/uncomment below
# data("mtcars")
# mtcars$vs=as.factor(mtcars$vs)
# mtcars$cyl=as.factor(mtcars$cyl)
# mtcars$hp=as.numeric(as.character(mtcars$hp))
# mtcars$cyl=relevel(mtcars$cyl,ref="6")
# names=c("Horsepower","VS", "Cylinder")
# mod=glm(carb~hp+vs+cyl,data=mtcars,family="poisson")
# report.poisson(mod,data=mtcars,labels=names,nrvirg=c(3,4,1,2),nrpval=3,reference=TRUE)

report.poisson=function (mod, data, labels,nrvirg=NA,nrpval,reference=TRUE){
  
  #Warning if there is a : in coeff names (interaction term) then might not compute
  if (any(sapply(names(mod$coefficients),function(x) grepl(":",x)))) {warning("Variable names contained ':', if the model contains an interaction term the function might not compute and if it does make sure to verify with summary(mod) that it is correct!")}
  
  #Stopping conditions
  if (missing(mod)) {
    stop("Function needs to have a model")
  }
  # Coef
  co <- summary(mod)$coef[-1]
  nvar <- length(co)
  varlist=attr(mod$terms, "term.labels")
  
  
  #Exp estimate for RR
  co <- exp(co)
  
  
  if (class(try(car::Anova(mod, test.statistic="Wald"),silent=TRUE))[1]=="try-error") {stop(paste0("Anova Wald statistics could not be used"))}
  ano=car::Anova(mod, test.statistic="Wald")
  
  modxlevels=mod$xlevels
  
  
  #Par defaut 4 chiffres apres la virgule pour les pvalues
  if (missing(nrpval)) {
    nrpval=4
  }
  
  if (nrpval<2) {nrpval=2}
  
  if (missing(labels)) {labels=varlist}
  
  #labels not same as length(mod$coefficients) then use mod$coefficients instead
  if (length(varlist) != length(labels)) {
    warning("labels length was not the same as the number of variables in the model")
    labels=varlist
  }
  
  
  
  
  
  #Par defaut 2 chiffres apres la virgule pour les HR/Beta et CI
  Pnames=c()
  tempnrvirg=c()
  for (i in 1:length(varlist)) {
    
    if (varlist[i] %in% names(modxlevels)) {
      if (varlist[i] %in% rownames(ano)[which(ano$Df>1)]) {
        Pnames=c(Pnames,varlist[i])
        if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
        for (j in 1:(ano$Df[i]-1)) {
          Pnames=c(Pnames,varlist[i])
          if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
        }
      } else {
        Pnames=c(Pnames,varlist[i])
        if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
      }
    } else {
      Pnames=c(Pnames,varlist[i])
      if (is.na(nrvirg[i])) {tempnrvirg=c(tempnrvirg,2)} else {tempnrvirg=c(tempnrvirg,nrvirg[i])}
    }
  }
  Pnames
  nrvirg=tempnrvirg
  
  
  
  
  
  #pvalues
  pval = summary(mod)$coefficients[,4][-1]
  pvalf = format(round(pval, nrpval), nsmall = nrpval, scientific = F,
                 digits = nrpval)
  
  lowervalch=paste0("0.",paste(rep(0,nrpval-1),collapse=""),"1")
  lowerval=as.numchar(lowervalch)
  pvalf[pval < lowerval] = paste0("<",lowervalch)
  
  
  
  li <- format(round(exp(confint.default(mod))[,1][-1], 10), nsmall = 10, scientific = F,digits = 10)
  li=gsub(" ", "", li, fixed = TRUE)
  ls <- format(round(exp(confint.default(mod))[,2][-1], 10), nsmall = 10, scientific = F,digits = 10)
  ls=gsub(" ", "", ls, fixed = TRUE)
  li2=as.numchar(li)
  ls2=as.numchar(ls)
  co2=as.numchar(co)
  std2=as.numchar(summary(mod)$coef[,2][-1])
  li=c()
  ls=c()
  co=c()
  for (i in 1:length(li2)) {
    li=c(li,as.character(format(round(li2[i], nrvirg[i]), nsmall = nrvirg[i], scientific = F,digits = nrvirg[i])))
    ls=c(ls,as.character(format(round(ls2[i], nrvirg[i]), nsmall = nrvirg[i], scientific = F,digits = nrvirg[i])))
    co=c(co,as.character(format(round(co2[i], nrvirg[i]), nsmall = nrvirg[i], scientific = F,digits = nrvirg[i])))
  }
  li
  ls
  co
  ci = paste0("[", li, ";", ls, "]")
  
  
  
  Pdf=data.frame(Pnames,Beta = co, '95% CI' = ci,
                 'Pval' = pvalf,check.names=FALSE)
  
  
  
  Cnames=c()
  for (i in 1:length(varlist)) {
    
    if (varlist[i] %in% names(modxlevels)) {
      if (varlist[i] %in% rownames(ano)[which(ano$Df>1)]) {
        if (reference==TRUE) {Cnames=c(Cnames,paste0(labels[i]," (ref. ",modxlevels[varlist[i]][[1]][1],")"))} else {Cnames=c(Cnames,paste0(labels[i]))}
        for (j in 1:(ano$Df[i])) {
          Cnames=c(Cnames,paste0("\t",modxlevels[varlist[i]][[1]][j+1]))
        }
      } else {
        if (reference==TRUE) {Cnames=c(Cnames,paste0(labels[i]," (ref. ",modxlevels[varlist[i]][[1]][1],")"))} else {Cnames=c(Cnames,paste0(labels[i]))}
      }
    } else {
      Cnames=c(Cnames,labels[i])
    }
  }
  Cnames
  
  Variables=c()
  Beta=c()
  CI=c()
  Pvalue=c()
  compt=0
  compthr=0
  
  for (i in unique(Pdf$Pnames)) {
    compt=compt+1
    compthr=compthr+1
    
    if (i %in% rownames(ano)[which(ano$Df>1)]) {
      Variables=c(Variables,Cnames[compt])
      Beta=c(Beta,"")
      CI=c(CI,"")
      
      temppval=format(round(ano[,3][which(rownames(ano)==i)], nrpval), nsmall = nrpval, scientific = F,
                      digits = nrpval)
      if (as.numchar(temppval)< lowerval) {temppval = paste0("<",lowervalch)}
      
      Pvalue=c(Pvalue,temppval)
      
      
      for (j in 1:(ano$Df[which(rownames(ano)==i)])) {
        compt=compt+1
        if (j!=1) {compthr=compthr+1}
        Variables=c(Variables,Cnames[compt])
        Beta=c(Beta,co[compthr])
        CI=c(CI,ci[compthr])
        Pvalue=c(Pvalue,"")
      }
      
    } else {
      Variables=c(Variables,Cnames[compt])
      Beta=c(Beta,co[compthr])
      CI=c(CI,ci[compthr])
      Pvalue=c(Pvalue,pvalf[compthr])
    }
  }
  
  
  # Final data frame with all information
  results = data.frame(Variables, "Rate Ratio" = Beta, '95% CI' = CI,
                       'p value' = Pvalue,check.names=FALSE)
  rownames(results) = NULL
  
  Table=results
  colnames(Table)[1]=" "
  Ftab=flextable(Table)
  Ftab=border_remove(Ftab)
  Ftab <- bold(Ftab, part = "header")
  font.name="Times New Roman"
  font.size=10
  Ftab <- flextable::font(Ftab,fontname=font.name,part ="all")
  Ftab <- fontsize(Ftab,size=font.size,part ="all")
  Ftab <- fontsize(Ftab,size=11,part ="header")
  Ftab <- align_text_col(Ftab, align = "left", header = TRUE)
  Ftab <- flextable::align(Ftab,i=1,j=2:ncol(Table),align = "center", part = "header")
  Ftab <- flextable::align(Ftab,i=1:nrow(Table),j=2:ncol(Table),align = "center", part = "body")
  Ftab <- hline(Ftab,border = fp_border(width = 1), part = "header" )
  Ftab <- hline_top(Ftab, border = fp_border(width = 1), part = "header" )
  Ftab <- hline_bottom(Ftab, border = fp_border(width = 1), part = "body")
  Ftab <- autofit(Ftab)
  Ftab <- height_all(Ftab, height=0.1, part = "body")
  Ftab <- height_all(Ftab, height=0.3, part = "header")
  Ftab <- height_all(Ftab, height=0.1, part = "footer")
  Ftab <- width(Ftab,j=1, width = 2.5)
  Ftab <- width(Ftab,j=2:(ncol(Table)-1), width = 1.1)
  Ftab <- width(Ftab,j=ncol(Table), width = 0.6)
  
  boldrows=suppressWarnings(which(as.numchar(gsub("<","",Table[,ncol(Table)]))<0.05))  #number of row to bold pvalue <0.05 by default
  for (rownr in boldrows)  {
    Ftab <- bold(Ftab,i=rownr,j=ncol(Table))
  }
  
  return(Ftab)
}
###############################################################################











###############################################################################
# Varlist function to get the variables used in a model:
###############################################################################
fvarlist=function (mod){
return(attr(mod$terms, "term.labels"))
}
###############################################################################

















#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

############## jdk.desc (cf. ClinReport package on CRAN)

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

jdk.desc <- function(varlist,group=NULL,data,names,spaces=TRUE) {
  data=droplevels(data)
  temp <- data
  temp$int <- 1
  totonly=FALSE
  if (is.null(group)) {
    temp$group=c(rep(1,floor(nrow(temp)/2)),rep(2,ceiling(nrow(temp)/2)))
    temp$group=as.factor(temp$group)
    totonly=TRUE
    group="group"
  }
  desc <- NULL
  
  #If names is not given the variable names are used as labels
  if (is.null(names)) {names=varlist}
  
  if (length(names) != length(varlist)) {stop("Varlist and names are not the same length.")}
  for (i in 1:length(varlist)) {
    if (group == varlist[i]) {next}
    if (is.factor(data[,varlist[i]])) {
      a <- jdk.report.quali(data=temp,y=varlist[i],x1=group,x2="int",total=TRUE)
      a$Factor <- NULL
      a <- cbind(Variable=names[i],a)
      colnames(a)[2]="Statistics/levels"
      #t <- data.frame(Variable=NA,"Statistics/levels"=NA,setNames(as.list(c(NA,NA,NA)),colnames(a)[3:ncol(a)]),check.names=FALSE)
    } else {
      a2 <- jdk.report.quali(data=temp,y=varlist[i],x1=group,x2="int",total=T)
      a2$Factor <- NULL
      a <- report.quanti(data=temp,y=varlist[i],x1=group,total=T)
      colnames(a)=colnames(a2)
      colnames(a)[1]="Statistics/levels"
      a <- cbind(Variable=names[i],a)
    }
    empty <- apply(a[1,],2,function(x) x=" ")
    #colnames(a) <- colnames(t)
    if (spaces==TRUE) {desc <- suppressWarnings(rbind(desc,empty,a))} else {desc <- suppressWarnings(rbind(desc,a))}
  }
  if (spaces==TRUE) {
    desc=desc[-1,]
    desc[,1]=as.character(desc[,1])
    desc[,1][which(is.na(desc[,1]))]=" "
  }
  
  #if (length(unique(colnames(desc))) != length(colnames(desc))) {colnames(desc[3:ncol(desc)])=levels(data[,group])}
  
  if (totonly) {
    desc=desc[,-c(3,4)]
    Table=desc
    desc=merge_v(flextable(desc), j = 1, part = "body")
    
    Ftab=desc
    Ftab=border_remove(Ftab)
    Ftab <- bold(Ftab, part = "header")
    font.name="Times New Roman"
    font.size=10
    Ftab <- flextable::font(Ftab,fontname=font.name,part ="all")
    Ftab <- fontsize(Ftab,size=font.size,part ="all")
    Ftab <- fontsize(Ftab,size=11,part ="header")
    Ftab <- align_text_col(Ftab, align = "left", header = TRUE)
    Ftab <- flextable::align(Ftab,i=1,j=2:ncol(Table),align = "center", part = "header")
    Ftab <- flextable::align(Ftab,i=1:nrow(Table),j=2:ncol(Table),align = "center", part = "body")
    Ftab <- hline(Ftab,border = fp_border(width = 1), part = "header" )
    Ftab <- hline_top(Ftab, border = fp_border(width = 1), part = "header" )
    Ftab <- hline_bottom(Ftab, border = fp_border(width = 1), part = "body")
    Ftab <- autofit(Ftab)
    Ftab <- height_all(Ftab, height=0.1, part = "body")
    Ftab <- height_all(Ftab, height=0.3, part = "header")
    Ftab <- height_all(Ftab, height=0.1, part = "footer")
    # Ftab <- width(Ftab,j=1, width = 2.5)
    # Ftab <- width(Ftab,j=2:(ncol(Table)-1), width = 1.1)
    # Ftab <- width(Ftab,j=ncol(Table), width = 0.6)
    desc=Ftab
    
    return(desc)
  }
  
  Table=desc
  desc=merge_v(flextable(desc), j = 1, part = "body")
  
  Ftab=desc
  Ftab=border_remove(Ftab)
  Ftab <- bold(Ftab, part = "header")
  font.name="Times New Roman"
  font.size=10
  Ftab <- flextable::font(Ftab,fontname=font.name,part ="all")
  Ftab <- fontsize(Ftab,size=font.size,part ="all")
  Ftab <- fontsize(Ftab,size=11,part ="header")
  Ftab <- align_text_col(Ftab, align = "left", header = TRUE)
  Ftab <- flextable::align(Ftab,i=1,j=2:ncol(Table),align = "center", part = "header")
  Ftab <- flextable::align(Ftab,i=1:nrow(Table),j=2:ncol(Table),align = "center", part = "body")
  Ftab <- hline(Ftab,border = fp_border(width = 1), part = "header" )
  Ftab <- hline_top(Ftab, border = fp_border(width = 1), part = "header" )
  Ftab <- hline_bottom(Ftab, border = fp_border(width = 1), part = "body")
  Ftab <- autofit(Ftab)
  Ftab <- height_all(Ftab, height=0.1, part = "body")
  Ftab <- height_all(Ftab, height=0.3, part = "header")
  Ftab <- height_all(Ftab, height=0.1, part = "footer")
  # Ftab <- width(Ftab,j=1, width = 2.5)
  # Ftab <- width(Ftab,j=2:(ncol(Table)-1), width = 1.1)
  # Ftab <- width(Ftab,j=ncol(Table), width = 0.6)
  desc=Ftab
  
  return(desc)
}








jdk.report.quali <- function (data, y = NULL, x1 = NULL, x2 = NULL, x1.name = "Treatment", 
                              x2.name = "Factor", variable.name = "Levels", total = F, 
                              round = 0, order = NULL, add.N = F, id = NULL, percent = T, 
                              frequencie = T, at.row = NULL) 
{
  data[, y] = as.factor(data[, y])
  if (!is.null(x1)) {
    data[, x1] = as.factor(data[, x1])
  }
  if (!is.null(x2)) {
    data[, x2] = as.factor(data[, x2])
  }
  data[, y] = addNA(data[, y])
  freq = data.frame(table(data[, y], data[, x1], data[, x2]))
  n = data.frame(table(data[, x1], data[, x2]))
  freq$Var1 = addNA(freq$Var1)
  freq_list = split(freq, freq$Var1)
  for (i in 1:length(levels(data[, y]))) {
    if (percent & frequencie) {
      stat = round(100 * (freq_list[[i]]$Freq/n$Freq), 
                   round)
      freq_list[[i]]$stat = paste0(freq_list[[i]]$Freq, 
                                   " (", format(stat, nsmall = round), "%)")
    }
    if (!percent & frequencie) {
      freq_list[[i]]$stat = freq_list[[i]]$Freq
    }
    if (percent & !frequencie) {
      stat = round(100 * (freq_list[[i]]$Freq/n$Freq), 
                   round)
      freq_list[[i]]$stat = stat
    }
    if (!percent & !frequencie) {
      stop("percent and freq arguments cannot be both set to FALSE")
    }
    freq_list[[i]]$Freq = NULL
  }
  m = reshape2::melt(freq_list, id.vars = c("Var1", "Var2", "Var3"), 
                     measure.vars = c("stat"))
  # m$value = gsub(" ", "", m$value)
  freq = reshape2::dcast(m, Var3 + Var1 ~ Var2)
  colnames(freq)[1] = x2.name
  colnames(freq)[2] = variable.name
  freq[, 2] = as.character(freq[, 2])
  freq[is.na(freq[, 2]), 2] = paste0("Missing")
  if (total) {
    freq.tot = report.quali.1(data = data, y = y, x1 = x2)
    freq.tot = reshape2::melt(freq.tot, measure.vars = colnames(freq.tot)[-1], 
                              variable.name = x2.name, value.name = "Total")
    freq = data.frame(freq, Total = freq.tot[freq.tot[, 2] != 
                                               "Statistics", "Total"])
    if (length(which(colnames(freq) == "NA")) > 0) {
      freq = freq[, -which(colnames(freq) == "NA")]
    }
    # freq[, -c(1, 2)] = apply(freq[, -c(1, 2)], 2, function(x) gsub(" ", 
    #                                                                "", x))
    if (!is.null(order)) {
      freq[, variable.name] = factor(freq[, variable.name], 
                                     levels = order)
      freq = freq[order(freq[, x2.name], freq[, variable.name]), 
                  ]
    }
  }
  if (add.N & !is.null(id)) {
    N = tapply(data[, id], data[, x1], function(x) length(unique(x)))
    colnames(freq)[3:(3 + length(levels(data[, x1])) - 1)] = paste(levels(data[, 
                                                                               x1]), "(N=", N, ")", sep = "")
    if (total) {
      colnames(freq)[colnames(freq) == "Total"] = paste("Total", 
                                                        "(N=", sum(N), ")", sep = "")
    }
  }
  N = apply(ftable(data[, x1], data[, x2]), 1, max)
  colnames(freq)[-c(1, 2)] = paste0(colnames(freq)[-c(1, 2)], 
                                    " (N=", N, ")")
  if (!is.null(at.row)) {
    freq = spacetable(freq, at.row = at.row)
  }
  if (total) {
    colnames(freq)[length(colnames(freq))] = paste("Total"," (N=", sum(N), ")", sep = "")
  }
  freq
}




# TODO: Add comment
# 
###############################################################################

#' @export

report.quali.1=function(data,y,x1,x1.name="Treatment",variable.name="Levels",
                        total=F,add.N=F,id=NULL)
{
  
  temp=data
  temp$int=1
  freq=report.quali(temp,y,x1,x1.name="Treatment",
                    x2="int",variable.name=variable.name,total=total,add.N=add.N,id=id)[,-1]
  
  
  freq$Statistics=rep("n(%)",nrow(freq))
  
  freq=freq[,c(1,ncol(freq),(2:(ncol(freq)-1)))]
  
  colnames(freq)[1]=variable.name
  
  #	N=tapply(data[,x1],data[,x1],function(x)length(x))
  #	
  #	colnames(freq)[-c(1,2)]=paste0(colnames(freq)[-c(1,2)]," (N=",N,")")
  #	
  freq
  
  
}

# test
# library(ggplot2)
# library(reshape2)
# data(mpg)
#report.quali.1(data=as.data.frame(mpg),y="cyl",x1="drv")








# TODO: Add comment
# 
###############################################################################

#' Returns "qualitative" statistics (frequencies and percentages) of one factor, splited 
#' according to the levels of two factors.
#' 
#'
#' @param data a data.frame object
#' @param y Character indicating a factor in the data
#' @param x1 Character indicating a factor in the data (levels will be displayed in columns)
#' @param x2 Character indicating a factor in the data (levels will be displayed in lines)
#' @param total Boolean which indicates if a column Total should added or not
#' @param round  Number to indicate how to round statistic values
#' @param order Character vector to reorder the levels of \code{x1}
#' @param add.N Boolean to add sample size in column headers
#' @param id Vector to identify the sample size in each group
#' @param percent Boolean to indicate if percentages should be reported
#' @param frequencie Boolean to indicate if frequencies should be reported
#' 
#' @description
#' Compute and report frequencies and percentages by levels of \code{y} and by levels of \code{x1}
#' and \code{x2}.
#' 
#' 
#' @details
#' This function computes and reports qualitative statistics by levels of \code{y} and by levels of \code{x1}
#' and \code{x2}. Run the example to show the results. If \code{total=T}, the last column is the statistics
#' performed overall levels of \code{x1} for each levels of \code{x2}. 

#' @return  
#' Data frame object. The first column is the levels
#' of the variable \code{x2}. The second column list the levels of \code{y}. The other columns
#' are the levels of the variable x1 (and column Total, if \code{total=T}).   
#' 
#' @seealso \code{\link{report.quali}} 

#' @examples
#' \dontrun{
#' data(data)
#' report.quali(data=data,y="clinical_cure",x1="GROUP",x2="timepoint_num")
#' }
#' 

#' @export

report.quali=function(data,y=NULL,x1=NULL,x2=NULL,x1.name="Treatment",
                      x2.name="Factor",
                      variable.name="Levels",total=F,
                      round=0,order=NULL,add.N=F,id=NULL,percent=T,frequencie=T,at.row=NULL)
{
  
  # transform into factor (just in case)
  
  data[,y]=as.factor(data[,y])
  if(!is.null(x1))
  {
    data[,x1]=as.factor(data[,x1])
  }
  
  if(!is.null(x2))
  {
    data[,x2]=as.factor(data[,x2])
  }
  
  # add NA as category
  # (to be more easily handled)
  
  data[,y]=addNA(data[,y])
  
  
  # Compute frequency and total sample size for percentage
  
  freq=data.frame(table(data[,y],
                        data[,x1],data[,x2]))
  
  n=data.frame(table(data[,x1],data[,x2]))
  
  # add NA as category on frequency too
  
  freq$Var1=addNA(freq$Var1)
  
  # add percent to frequencies
  
  freq_list=split(freq,freq$Var1)
  
  for(i in 1:length(levels(data[,y])))
  {
    if(percent & frequencie)
    {
      stat=round(100*(freq_list[[i]]$Freq/n$Freq),round)
      
      freq_list[[i]]$stat=paste0(freq_list[[i]]$Freq," (",
                                 format(stat,nsmall = round),"%)")
    }
    
    if(!percent & frequencie) #only frequencies are reported
    {
      freq_list[[i]]$stat=freq_list[[i]]$Freq
    }
    
    if(percent & !frequencie) #only percentages are reported
    {
      stat=round(100*(freq_list[[i]]$Freq/n$Freq),round)
      freq_list[[i]]$stat=stat
    }
    
    if(!percent & !frequencie) #only percentages are reported
    {
      stop("percent and freq arguments cannot be both set to FALSE")
    }
    
    freq_list[[i]]$Freq=NULL
  }
  
  # reshape list to data using melt and dcast
  
  m=reshape2::melt(freq_list,id.vars=c("Var1","Var2","Var3"),
                   measure.vars=c( "stat"))
  
  # m$value=gsub(" ","",m$value)
  
  freq=reshape2::dcast(m,Var3+Var1~Var2)
  
  colnames(freq)[1]=x2.name
  colnames(freq)[2]=variable.name
  freq[,2]=as.character(freq[,2])
  #	freq[is.na(freq[,2]),2]=paste0("Missing"," (",y,")")
  freq[is.na(freq[,2]),2]=paste0("Missing")
  
  # Add column Total if requested
  
  if(total)
  {
    freq.tot=report.quali.1(data=data,y=y,
                            x1=x2)
    
    
    freq.tot=reshape2::melt(freq.tot,measure.vars=colnames(freq.tot)[-1],
                            variable.name=x2.name,value.name="Total")
    
    freq=data.frame(freq,Total=freq.tot[freq.tot[,2]!="Statistics",
                                        "Total"])
    
    if(length(which(colnames(freq)=="NA"))>0)
    {
      freq=freq[,-which(colnames(freq)=="NA")]
    }
    
    # freq[,-c(1,2)]=apply(freq[,-c(1,2)],2,function(x)gsub(" ","",x))
    
    #		colnames(freq)[2]=""
    
    if(!is.null(order))
    {
      freq[,variable.name]=factor(freq[,variable.name],levels=order)
      freq=freq[order(freq[,x2.name],freq[,variable.name]),]
    }
    
  }
  
  # Add sample size in column header if requested
  # using subject identification number
  # works only if id is filled
  
  if(add.N & !is.null(id))
  {
    N=tapply(data[,id],data[,x1],function(x) length(unique(x)))
    colnames(freq)[3:(3+length(levels(data[,x1]))-1)]=paste(levels(data[,x1]),
                                                            "(N=",N,")",sep="")
    
    if(total)
    {
      colnames(freq)[colnames(freq)=="Total"]=paste("Total",
                                                    "(N=",sum(N),")",sep="")
    }
  }
  
  
  
  N=apply(ftable(data[,x1],data[,x2]),1,max)
  colnames(freq)[-c(1,2)]=paste0(colnames(freq)[-c(1,2)]," (N=",N,")")
  
  
  if(!is.null(at.row))
  {
    freq=spacetable(freq,at.row=at.row)
  }
  
  freq
}











###############################################################################

#' Returns "quantitative" statistics of one numerical variable, splited 
#' according to the levels of two factors.
#' 
#'
#' @param data Data.frame object
#' @param y Character indicating a numerical vector in the data frame passed to \code{data} argument
#' @param x1 Character indicating a factor in the data (levels will be displayed in columns)
#' @param x2 Character indicating a factor in the data (levels will be displayed in lines)
#' @param total Logical which indicates if a column Total should added or not
#' @param list_stat Character vector with value(s) in:
#'  "N","mean","sd","median","mad","q1","q3","min","max","missing","geomean","geomean1" 
#' @param regroup Logical to indicate if statistics should be regrouped with the regroup_list argument
#' @param regroup_list Character vector in: 
#' "mean_sd","median_mad","q1_q3","min_max" to indicate which statistics are to be regrouped (see Details below)
#' @param order_stat Numerical vector to indicate the order in which to report statistics
#' @param label_stat  Character vector to indicate the label of the statistics
#' @param round Numeric to indicate how to round statistics
#' @param total Logical to indicate if a "Total" column should be added
#' @param scientific Logical to indicate if statistics should be displayed in scientific notations
#' @param digits Numeric (used if scientifc=T) to indicate how many digits to use in scientific notation


#' @description
#' \code{report.quanti} 
#' Returns quantitative descriptive statistics such as mean, median, standard deviation etc...
#' 
#' 
#' @details
#' This function computes and reports quantitative statistics on \code{y} by level of \code{x1}
#' and \code{x2}. 
#' You can run the example to show the results. If \code{total=T}, the last column is the statistics
#' performed overall levels of \code{x1} for each levels of \code{x2}. 
#' Quantiles are calculated using type 3 (SAS) algorithms.
#' "geomean" compute the geometric mean defined as exp(mean(log(y)))
#' "geomean1" compute the geometric mean defined as exp(mean(log(y+1))-1)
#' for the case where y is equal 0 for some values.
#' 
#' \code{N} returns the number of observations (including NA values)
#' \code{n} returns the number of observations (excluding NA values)

#' @return  
#' A data frame object. The first column is the levels
#' of the variable \code{x2}. The second column list the statistics. The other columns
#' are the levels of the variable \code{x1} (and column Total, if \code{total=T}).   
#' 
#' @seealso \code{\link{report.quali}} 

#' @examples
#' \dontrun{
#' data(data)
#' 
#' # Quantitative statistics with no factor
#' 
#' report.quanti(data=data,y="DEMEANOUR_num")
#' 
#' # Quantitative statistics with one factor
#' 
#' report.quanti(data=data,y="DEMEANOUR_num",x1="GROUP")
#' 
#' # One factor with total column
#' 
#' report.quanti(data=data,y="DEMEANOUR_num",x1="GROUP",total=T)
#' 
#' # Quantitative statistics with two factors
#' 
#' report.quanti(data=data,y="DEMEANOUR_num",x1="GROUP",x2="TIMEPOINT")
#' 
#' # Quantitative statistics with two factors and a total column
#' 
#' report.quanti(data=data,y="DEMEANOUR_num",x1="GROUP",x2="TIMEPOINT",total=T)
#' 
#' 
#' # Quantitative statistics with two factors and only 3 statistics
#' # without regrouping statistics and with new names for statistics
#' 
#' report.quanti(data=data,y="DEMEANOUR_num",x1="GROUP",x2="TIMEPOINT",total=F,
#' list_stat=c("N","mean","sd","geomean1"),
#' label_stat=c("N","M","Standard deviation","Geo M"),
#' regroup=F,order_stat=c("N"=1,"mean"=2,"sd"=3,"geomean1"=4))
#' 
#' #' # Quantitative statistics with spacing rows
#' 
#' report.quanti(data=data,y="DEMEANOUR_num",x1="GROUP",x2="TIMEPOINT",total=T,at.row=6)

#' }
#' 
#' @export




report.quanti=function(data,y,x1=NULL,x2=NULL,list_stat=c("N",
                                                          "mean","sd","median","mad","q1","q3","min",
                                                          "max","missing"),regroup=T,regroup_list=c("mean_sd",
                                                                                                    "median_mad","q1_q3","min_max"),
                       order_stat=c("N"=1,"mean_sd"=2,
                                    "median_mad"=3,"q1_q3"=4,
                                    "min_max"=5,"missing"=6),
                       label_stat=c("N","Mean (SD)","Median (MAD)",
                                    "[Q1;Q3]","[Min;Max]","Missing"),
                       variable.name=y,
                       round=2,
                       total=F,scientific=F,digits=NULL,at.row=NULL)
{
  
  # TODO: add.N=T and N parameters to add N in columns
  
  # TODO add case if x1 and x2 are NULL
  
  # TODO add test:
  
  
  # if regroup=T and length(regroup_list)>1 test if
  
  ## median_mad then median and mad in list_stat
  ## q1_q3 then q1 and q3 in list_stat
  ## min_max then min and max in list_stat
  
  # TODO add test_that tests
  
  
  ################################
  # Check requirement
  ################################
  
  if(regroup==T & length(regroup_list)>1)
  {
    if("%in%"("mean_sd",regroup_list))
    {
      if(!"%in%"("mean",list_stat) | !"%in%"("sd",list_stat))
      {
        stop("mean and sd should be in list_stat")
      }
    }	
  }
  
  
  if(length(unique(order_stat))!=length(order_stat))
  {
    stop("length(unique(order_stat)) should be equal to length(order_stat)")
  }
  
  if(!is.logical(total))
  {
    stop("Argument total argument must be logical")
  }
  
  
  if(!is.numeric(digits) & !is.null(digits))
  {
    stop("Argument digits must be numeric")
  }
  
  
  if(length(order_stat)!=length(label_stat))
  {
    stop("length(order_stat) should be equal to length(label_stat)")
  }
  
  
  if(any(!"%in%"(names(order_stat),c(list_stat,regroup_list))))
  {
    stop("All names in order_stat should be in c(list_stat,regroup_list)")
  }
  
  
  if(regroup==T & length(regroup_list)==0)
  { 
    stop("If regroup=T, regroup_list should not be empty")
  }
  
  if(!is.numeric(data[,y]))
  {
    
    if(length(which(is.na(data[,y])))==length(data[,y]))
    {
      data[,y]=as.numeric(data[,y])
    }else
    {
      stop(paste(y,"=y should be a numeric variable"))
    }
  }
  
  if(!is.numeric(round))
  {
    stop(paste("round should be numeric"))
  }
  
  
  ################################
  # start function
  ################################
  
  # group data by factors
  
  if(!is.null(x2) & !is.null(x1))
  {
    by_GROUP=data %>% group_by_(x1)%>% group_by_(x2,add=T)	
  }
  
  
  if(is.null(x2) & !is.null(x1))
  {
    by_GROUP=data %>% group_by_(x1)
  }
  
  
  if(is.null(x2) & is.null(x1))
  {
    by_GROUP=data
  }
  
  # define statistics
  
  N=as.formula(paste0("~","length(.)"))
  n=as.formula(paste0("~","length(",y,",na.rm=T)"))
  mean=as.formula(paste0("~","mean(",y,",na.rm=T)"))
  sd=as.formula(paste0("~","sd(",y,",na.rm=T)"))
  median=as.formula(paste0("~","median(",y,",na.rm=T)"))
  mad=as.formula(paste0("~","mad(",y,",na.rm=T)"))
  q1=as.formula(paste0("~","quantile(",y,",na.rm=T,0.25,type =3)"))
  q3=as.formula(paste0("~","quantile(",y,",na.rm=T,0.75,type =3)"))
  min=as.formula(paste0("~","suppressWarnings(min(",y,",na.rm=T))"))
  max=as.formula(paste0("~","suppressWarnings(max(",y,",na.rm=T))"))
  missing=as.formula(paste0("~","length(",y,"[is.na(",y,")])"))
  
  geomean_func=function(x) exp(mean(log(x),na.rm=T))
  geomean_func1=function(x) exp(mean(log(x+1),na.rm=T)-1)
  
  
  geomean=as.formula(paste0("~","geomean_func(",y,")"))
  geomean1=as.formula(paste0("~","geomean_func1(",y,")"))
  
  #	select statistics
  
  stat_list=c("N"=N,"n"=n,"mean"=mean,
              "sd"=sd,"median"=median,"mad"=mad,
              "q1"=q1,"q3"=q3,"min"=min,"max"=max,
              "missing"=missing,
              "geomean"=geomean,
              "geomean1"=geomean1)
  
  stat_list=stat_list["%in%"(names(stat_list),list_stat)]
  
  # compute statistics
  
  stat=data.frame(by_GROUP %>% summarise_at(.funs=funs_(dots=stat_list),.vars=y))
  
  # format outputs
  
  list_stat2=list_stat[!"%in%"(list_stat,c("N","missing"))]
  stat[,list_stat2]=format(round(stat[,list_stat2],round),nsmall=round,scientific=scientific,digits=digits)
  
  # Regroup stat
  
  if(regroup)
  {
    if("%in%"("mean_sd",regroup_list))
    {
      stat$mean_sd=paste0(stat$mean,"(",stat$sd,")")
      stat$mean=NULL
      stat$sd=NULL
    }
    
    if("%in%"("median_mad",regroup_list))
    {
      stat$median_mad=paste0(stat$median,"(",stat$mad,")")
      stat$median=NULL
      stat$mad=NULL
    }
    
    if("%in%"("min_max",regroup_list))
    {
      stat$min_max=paste0("[",stat$min,";",stat$max,"]")
      stat$min=NULL
      stat$max=NULL
    }
    
    if("%in%"("q1_q3",regroup_list))
    {
      stat$q1_q3=paste0("[",stat$q1,";",stat$q3,"]")
      stat$q1=NULL
      stat$q3=NULL
    }
  }
  
  # set as factors (for melt droppings message)
  
  if(!is.null(x2) & !is.null(x1))
  {
    stat[,-c(1,2)]=apply(stat[,-c(1,2)],2,as.factor)
  }
  
  if(!is.null(x2) & is.null(x1))
  {
    stat[,-c(1)]=apply(stat[,-c(1)],2,as.factor)
  }
  
  if(is.null(x2) & is.null(x1))
  {
    stat=apply(stat,2,as.factor)
  }
  
  
  # reshape the damn thing
  
  measure.var=colnames(stat)["%in%"(colnames(stat),c(list_stat,
                                                     regroup_list))]
  
  m=reshape2::melt(stat,id.vars=c(x1,x2),measure.var=measure.var,
                   variable.name="Statistics")
  
  if(!is.null(x2) & !is.null(x1))
  {
    form=as.formula(paste0(x2,"+","Statistics~",x1))
  }
  
  if(is.null(x2) & !is.null(x1))
  {
    form=as.formula(paste0("Statistics~",x1))		
  }
  
  if(is.null(x2) & is.null(x1))
  {
    form=as.formula(paste0("Statistics~."))		
  }
  
  if(!is.null(x2) | !is.null(x1))
  {
    stat2=reshape2::dcast(m,form,value.var="value")
    stat2$Statistics=factor(stat2$Statistics,levels=names(order_stat)[order(order_stat)])
  }else
  {
    stat2=m
    stat2$Statistics=rownames(m)
    rownames(stat2)=NULL
    stat2$Statistics=factor(stat2$Statistics,levels=names(order_stat)[order(order_stat)])
    stat2=stat2[,c(2,1)]
    colnames(stat2)[2]=y
  }
  
  
  if(!is.null(x2))
  {
    stat2=stat2[order(stat2[,x2],stat2$Statistics),]
  }else
  {
    stat2=stat2[order(stat2$Statistics),]
  }
  
  if(is.null(x1) & is.null(x2))
  {
    colnames(stat2)[2]=variable.name
  }
  
  levels(stat2$Statistics)[levels(stat2$Statistics)=="N"]="N"
  levels(stat2$Statistics)[levels(stat2$Statistics)=="mean_sd"]="Mean (SD)"
  levels(stat2$Statistics)[levels(stat2$Statistics)=="median_mad"]="Median (MAD)"
  levels(stat2$Statistics)[levels(stat2$Statistics)=="min_max"]="[Min;Max]"
  levels(stat2$Statistics)[levels(stat2$Statistics)=="missing"]="Missing"
  levels(stat2$Statistics)[levels(stat2$Statistics)=="q1_q3"]="[Q1;Q3]"
  
  #	stat2
  
  # add Total, if requested
  
  if(total)
  {
    if(is.null(x2))
    {
      data$intercept=1
      temp=report.quanti(data=data,y=y,x1="intercept")
      stat2=cbind(stat2,Total=temp[,2])
    }
    
    if(!is.null(x2))
    {
      data$intercept=1
      temp=report.quanti(data=data,y=y,x1=x2,x2="intercept")
      temp$intercept=NULL
      temp=reshape2::melt(temp,id.vars="Statistics")
      stat2=cbind(stat2,Total=temp[,3])
    }
    
  }
  
  if(!is.null(at.row))
  {
    stat2=spacetable(stat2,at.row=at.row)
  }
  
  return(stat2)
}


#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################








# 
# 
# 
# #################################
# #Reporting example:
# 
# doc <- read_docx(path = "I:/10-methodologiste/JOE/ELITE/PRGM/template.docx")
# 
# doc <- body_add_break(doc)
# doc <- body_add_par(doc,value = "Sommaire", style = "header") 
# doc <- body_add_toc(doc,level = 2)
# doc <- body_add_break(doc)
# 
# doc <- body_add_par(doc,value = "Tableau 1. Caract?ristiques d?mographiques et ant?c?dents de la population ", style = "heading 1")
# doc <- body_add_flextable(doc,value = FTable10) 
# doc <- body_add_break(doc)
# 
# doc <- body_add_par(doc,value = "Tableau 2. Caract?ristiques d?mographiques et ant?c?dents de la population ", style = "heading 1")
# doc <- body_add_par(doc,"", style = "Normal")
# doc <- body_add_flextable(doc,value = FTable2) 
# doc <- body_add_break(doc)
# 
# doc <- body_add_par(doc,value = "Graphique22222. Caract?ristiques d?mographiques et ant?c?dents de la population ", style = "heading 1")
# doc <- body_add_par(doc,"", style = "Normal")
# filename <- tempfile(fileext = ".emf")
# emf(file = filename, width = 6, height = 6)
# print(plot1$plot) #$plot for ggsurvplot type
# dev.off()
# doc <- body_add_img(doc,src = filename, width = 6, height = 6)
# doc <- body_add_break(doc)
# 
# doc <- body_add_par(doc,value = "333333 Caract?ristiques d?mographiques et ant?c?dents de la population ", style = "heading 1")
# doc <- body_add_par(doc,"", style = "Normal")
# filename <- tempfile(fileext = ".emf")
# emf(file = filename, width = 6, height = 6)
# print(plot(1,1)) #simple print for any ggplot or normal plot
# dev.off()
# doc <- body_add_img(doc,src = filename, width = 6, height = 6)
# doc <- body_add_break(doc)
# 
# print(doc, target = "I:/10-methodologiste/JOE/NEBULAMB/PRGM/test.docx")