function (run) {
  #REQUIRED USER SPECIFICATIONS PORTION
  alpha = 0.05
  power = 0.80
  betaxz= 1
  sigmasq = 1
  xzvec = c(80.34, 47, 85.95, 80.22, 58, 90.14, 72, 90.75, 73, 91.2, 
            55, 89.97, 88.66, 58, 89.81, 52, 77.25, 54, 93, 55, 45, 90, 85, 
            31, 79.38, 87.7, 65, 86.82, 53, 85.82, 61, 88.75, 59, 91.17, 79, 
            54, 84.12, 68, 88.45, 46, 82.59, 58, 90.02, 48, 64, 85.9, 83, 91.85, 
            51, 95.23, 84.15, 88.33, 52, 89.45, 51, 79.22, 57, 90.64, 54, 84.56, 
            66, 87.98, 68, 81.87, 56, 95.01, 89.5, 86.91, 54, 85.52, 54, 50, 
            87.59, 54, 67, 43, 82.15, 88.97, 62, 87.93, 77, 87.5, 53, 42, 85.68, 
            76, 86.9, 51, 90.96, 88.99, 89.65, 35, 87.93, 36, 80.43, 47, 94.34, 
            32, 88.3, 85, 84, 51, 86.82,  55, 83.64, 50, 90.16, 59, 87.73, 
            90.95, 45, 62,86.44, 72, 80.45, 75, 88.14, 64, 82.02, 55, 80.92, 
            61, 88.66, 61, 73)
  #END OF REQUIRED USER SPECIFICATION
  
  xz = matrix(xzvec, length(xzvec)/2,2,byrow = TRUE) 
  xe = cbind(xz, xz[,1]* xz[,2])
  n = nrow(xe)
  cmean = apply(xe, 2,mean)
  xc = xe-matrix(rep(cmean, 40),40,3,byrow = TRUE)
  h = matrix(rep(0,9),3,3)
  hh = h %x% h
  for (i in 1:n)
  {
    xci = xc[i,,drop = FALSE]
    h = h + (t(xci)%*%(xci))
    hh = hh + (t(xci)%*%(xci))%x%(t(xci)%*%(xci))
  }
  sigm = h/n
  psi = hh/n
  isigm = solve(sigm)
  muw = 1/isigm[3,3]
  vw = isigm[3,,drop = FALSE]
  varw = (muw^4)*(((vw)%x%(vw))%*%psi%*%(t(vw)%x%t(vw))-muw^(-2)) 
  lfxz = betaxz* sqrt(muw/sigmasq)
  print(c(alpha, power, n,betaxz, sigmasq),digits = 4)
  print(c(muw, varw, lfxz),digits = 4)
  
  #for numerical integration
  numint = 1000
  coevec = c(1,rep(c(4,2),numint/2-1),4,1)
  int = qnorm(0.999995,0,1)
  interval = 2* int/numint
  zvec = interval*seq(0,numint) + (-int)
  wzpdf = (interval/3)* coevec* dnorm(zvec, 0,1)
  
  #st approach
  stpower = 0
  m = 5
  while (stpower < power){
    m = m+1
    tcrit = qt(1-alpha/2,m-4,0)
    stpower = 1-pt(tcrit, m-4,sqrt(m)*lfxz)
  }
  nst = m
  
  #nt approach
  ntpower = 0
  m = max(nst-10,5)
  while (ntpower < power){
    m = m+1
    tcrit = qt(1-alpha/2,m-4,0)
    wvec = sqrt(varw/(m-1))* zvec + muw
    wvec = wvec*(wvec>0)
    ntpower = 1-sum(wzpdf*pt(tcrit, m-4,betaxz*sqrt((m-1)*wvec/sigmasq)))
  }
  nnt = m
  
  #recalculate ntpower for sample size nst
  m = nst
  tcrit = qt(1-alpha/2,m-4,0)
  wvec = sqrt(varw/(m-1))* zvec + muw
  wvec = wvec*(wvec>0)
  
  ntpower_nst =  1-sum(wzpdf*pt(tcrit, m-4,betaxz*sqrt((m-1)*wvec/sigmasq))) 
  dn = nnt-nst
  dntpower = ntpower-ntpower_nst
  print(c(nst, stpower, nnt, ntpower),digits = 4)
  print(c(dn, ntpower_nst, dntpower),digits = 4)
}

# original code below
function () {
  #REQUIRED USER SPECIFICATIONS PORTION
  alpha <- 0.05
  power <- 0.90
  betaxz <- 1
  sigmasq <- 16
  xzvec <- c( 0.11, -1.02, 0.58, -0.46, 0.27, 0.51, 0.64, 0.35, 0.76, 0.48, 0.98, 0.26, -0.76,
            0.06, -0.18, 0.15, 0.78, 0.70, 0.18,
            0.47, -0.58, 0.91, 0.28, 1.18, 1.14, 1.43, 0.83, -0.86, -0.78, 0.17, 0.61, -0.17, 0.08, 0.74,
            -0.67, -1.70, 1.52, 0.32, 0.18,
            0.85, 0.04, 2.06, 1.08, -0.31, -0.15, -0.62, -0.50, 0.79, -0.30, -0.02, 0.60, 0.56, -0.49,
            0.60, 0.87, 0.34, -0.29, -0.66, -1.04, 1.30, 0.14, -1.35,
            -1.12, -0.79, 0.74, 1.68, -0.69, -1.44, -0.80, -1.01, -3.21, -1.91, -0.42, -0.49, 2.79,
            2.35, -0.47, -0.96, -0.77, -1.58)
  #END OF REQUIRED USER SPECIFICATION
  xz = matrix(xzvec, length(xzvec)/2,2,byrow = TRUE) 
  xe = cbind(xz, xz[,1]* xz[,2])
  n = nrow(xe)
  cmean = apply(xe, 2,mean)
  xc = xe-matrix(rep(cmean, 40),40,3,byrow = TRUE)
  h = matrix(rep(0,9),3,3)
  hh = h %x% h
  for (i in 1:n)
  {
    xci = xc[i,,drop = FALSE]
    h = h + (t(xci)%*%(xci))
    hh = hh + (t(xci)%*%(xci))%x%(t(xci)%*%(xci))
  }
  sigm = h/n
  psi = hh/n
  isigm = solve(sigm)
  muw = 1/isigm[3,3]
  vw = isigm[3,,drop = FALSE]
  varw = (muw^4)*(((vw)%x%(vw))%*%psi%*%(t(vw)%x%t(vw))-muw^(-2)) 
  lfxz = betaxz* sqrt(muw/sigmasq)
  print(c(alpha, power, n,betaxz, sigmasq),digits = 4)
  print(c(muw, varw, lfxz),digits = 4)
  
  #for numerical integration
  numint = 1000
  coevec = c(1,rep(c(4,2),numint/2-1),4,1)
  int = qnorm(0.999995,0,1)
  interval = 2* int/numint
  zvec = interval*seq(0,numint) + (-int)
  wzpdf = (interval/3)* coevec* dnorm(zvec, 0,1)
  
  #st approach
  stpower = 0
  m = 5
  while (stpower < power){
    m = m+1
    tcrit = qt(1-alpha/2,m-4,0)
    stpower = 1-pt(tcrit, m-4,sqrt(m)*lfxz)
  }
  nst = m
  
  #nt approach
  ntpower = 0
  m = max(nst-10,5)
  while (ntpower < power){
    m = m+1
    tcrit = qt(1-alpha/2,m-4,0)
    wvec = sqrt(varw/(m-1))* zvec + muw
    wvec = wvec*(wvec>0)
    ntpower = 1-sum(wzpdf*pt(tcrit, m-4,betaxz*sqrt((m-1)*wvec/sigmasq)))
  }
  nnt = m
  
  #recalculate ntpower for sample size nst
  m = nst
  tcrit = qt(1-alpha/2,m-4,0)
  wvec = sqrt(varw/(m-1))* zvec + muw
  wvec = wvec*(wvec>0)
  
  ntpower_nst =  1-sum(wzpdf*pt(tcrit, m-4,betaxz*sqrt((m-1)*wvec/sigmasq))) 
  dn = nnt-nst
  dntpower = ntpower-ntpower_nst
  print(c(nst, stpower, nnt, ntpower),digits = 4)
  print(c(dn, ntpower_nst, dntpower),digits = 4)
}
