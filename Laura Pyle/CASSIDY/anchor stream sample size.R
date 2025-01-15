# Example: T1D at SCH (overall)
N_tot <- 1750073
sigma_2_p <- 0.001
p <- .00228
phi_p <- p*2

# try worked example
N_tot <- 250
sigma_2_p <- 0.000625
p <- .15
phi_p <- p*2


num <- N_tot*sigma_2_p
denom <- p*(1-phi_p)
Psi <- ((num/denom) + 1) ^ -1
N_sample <- Psi*N_tot