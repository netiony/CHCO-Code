p <- 6
n <- 60
R2 <- 0.42

alpha_hat <- 113
sigma <- 40

a <- R2*(n-p-1)+p
b <- log(1-(a/(n-1)))
c <- (p-2)/(n*b)
Sc <- 1+c

delta <- (p*(1-R2))/(n-1)

rad <- (sigma*sigma*(1-R2))/n
CI_lo_int <- alpha_hat - (1.96*sqrt(rad))
CI_hi_int <- alpha_hat + (1.96*sqrt(rad))
