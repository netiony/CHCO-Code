# ---------------------------------------------------------------------------- #
# PURPOSE
# Function to create a table of the number of records with non-missing
# values for a list of analysis variables at each study visit.  This
# is useful for Longitudinal Datasets with lots of variables to get
# a picture of the amount of information available at different visits.
#
# At each visit the difference between the number of non-missing values
# for a variable and the total number of participants in the study is
# the sum of two numbers:
#     number of records with missing values for the variable at the visit
#     +
#     number of participants who did not attend the visit
#
#
# INPUTS
# dframe:        dataframe containing the variables of interest
#
# visit.var:     a factor variable indicating the visit
#
# analysis.vars: a character vector containing the names of the variables in
#                dframe which are to be tabulated for counts of non-missing
#                values at each visit time (visit.var)
#
# OUTPUTS
#
# A nvar x nvisit table of the number of non-missing values for each of the
# variables in analysis.vars by visit
#
#
# COMMENTS:
# This function requires that the variable visit.var be passed as a factor
# variable.
#
# ---------------------------------------------------------------------------- #
MissTable <- function(dframe, visit.var, analysis.vars){
  list.nonmiss <- by(dframe, visit.var, notNA.byvis, analysis.vars)
  matx.nonmiss <- matrix(unlist(list.nonmiss,recursive=F), nrow=length(list.nonmiss[[1]]), ncol=length(list.nonmiss))
  rownames(matx.nonmiss) <- names(list.nonmiss[[1]])
  colnames(matx.nonmiss) <- names(list.nonmiss)
  return(matx.nonmiss)
}

notNA.byvis <- function(dframe, vars){
   if(dim(dframe)[1] > 1) apply(!apply(dframe[,vars],2,is.na),2,sum)
   else if(dim(dframe)[1]==1) 1*(!apply(dframe[,vars],2,is.na))
}
