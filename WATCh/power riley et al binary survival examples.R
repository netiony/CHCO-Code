# Estimate R2_csadj

###############################################################
# Criterion 1 - sample size to ensure shrinkage factor >= 0.9 #
###############################################################

# p is number of candidate predictors
p_1 <- 20
# Svh is shrinkage factor, should be 0.9
Svh_1 <- 0.9
# R2_csadj is anticipated value of adjusted Cox-Snell R2 (which needs to be estimated per section 3 of Riley et al)
R2_csadj_1 <- 0.1
# calculate sample size needed for targeted expected shrinkage
denom_1 <- (Svh_1 - 1)*(log(1-(R2_csadj_1/Svh_1)))
n_crit_1 <- p_1/denom_1

#########################################################################################
# Criterion 2 - ensure small absolute difference in apparent and adjusted R2_Nagelkerke #
#########################################################################################

# we calculate required shrinkage factor Svh needed
# if that shrinkage factor is larger than what is used in criterion 1, need to recalculate sample size for
# Criterion 1 using the new shrinkage factor

# R2_csadj is anticipated value of adjusted Cox-Snell R2 (which needs to be estimated per section 3 of Riley et al)
R2_csadj_2 <- 0.1
# delta is difference between apparent and adjusted R2_Nagelkerke, usually 0.05
delta_2 <- 0.05
# R2_csadj_max is max value of adjusted Cox-Snell R2 (see eq 23 of Riley et al)
# N is expected sample size; E is expected number of events
N_2 <- 100
E_2 <- 50
out_prop_2 <- E_2/N_2
ln_L_null_2 <- E_2*log(out_prop_2) + (N_2 - E_2)*log(1-(E_2/N_2))
R2_csadj_max_2 <- 1 - (exp((2*ln_L_null_2)/N_2))
# calculate required shrinkage factor
Svh_req_2 <- R2_csadj_2 / (R2_csadj_2 + (delta_2*R2_csadj_max_2))

############################################################################################
# Criterion 3 - BINARY OUTCOME - ensure precise estimate of overall risk (model intercept) #
############################################################################################

# delta_3 is the margin of error, recommend 0.05
delta_3 <- 0.05
# phi_3 is true outcome proportion - use 0.5 as this leads to the largest required sample size
phi_3 <- 0.5
# calculate required N
n_crit_3 <- ((1.96/delta_3)^2)*phi_3*(1-phi_3)

#############################################################################
# Criterion 3 - SURVIVAL OUTCOME - ensure precise estimate of overall risk  #
#############################################################################

# consider precision of estimated cumulative incidence (outcome risk) at a time point of interest (T)
# here we need to ensure that the lower and upper bounds of the 95% CI for the event rate
# are <= delta (i.e., 0.05) from the true value
# so the width of the CI should be <= 0.1

# lambda_hat is the estimate event rate (number of events per person year)
lambda_hat_3 <- 0.1
# T is total person years of follow-up at time point of interest
T_3 <- 10
# N_3 is total person years of follow-up in the study
N_3 <- 2366

# calculate lower and upper bound of 95% CI
a_3 <- 1.96*sqrt(lambda_hat_3/N_3)
b_upper_3 <- lambda_hat_3 + a_3
b_lower_3 <- lambda_hat_3 - a_3
lower_3 <- 1 - exp(-(b_lower_3)*T_3)
upper_3 <- 1 - exp(-(b_upper_3)*T_3)
ci_width_3 <- (upper_3 - lower_3)
  

