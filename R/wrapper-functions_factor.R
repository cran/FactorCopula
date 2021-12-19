# Wrapper functions for mle for 
# (1) one-factor copula model.
# (2) two-factor copula model.
# (3) M2 statistic for one-factor and two-factor copula models.
#------------------------------------------------------------
# MLE Estimation of one-factor model
# Input:
# 1) continuous: continuous data.
# 2) ordinal: ordinal data.
# 3) count: count data.
# 4) copF1: copula names for the first factor.
# 5) gl: Gauss legendre quadrature points and weights.

# Outputs: 
# 1) "cutpoints": cutpoints for ordinal data.
# 2) "uni.count.est": univariate estimates for count variables.
# 3) "loglik": log-likelihood value.
# 4) "cpar": estimates copula parameters.
# 5) "taus": estimated Kendall taus.
# 6) "SEs": Standard errors of the estimated taus.

mle1factor=function(continuous=NULL, ordinal=NULL, count=NULL, copF1, gl, hessian=F, print.level=0){
  #Gauss legendre points and weights
  gln<-gl$nodes
  glw<-gl$weights
  d = length(copF1)
  #--------------------------------------------------------
  #                 Continuous variables
  #-----------------                    -------------------
  u <- NULL
  if(!is.null(continuous)){
    n<-nrow(continuous)
    u<-apply(continuous,2,rank)/(n+1)
  }
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  cutpoints<-NULL
  if(!is.null(ordinal)){
    n<-nrow(ordinal)
    cutpoints<-cuts(ordinal)
  }
  #--------------------------------------------------------
  #                   Count variables
  #-------------------               ----------------------
  est_count<-NULL
  if(!is.null(count)){
    n<-nrow(count)
    d3<-ncol(count)
    est_count=matrix(NA,ncol = 2,nrow = d3 )
    for(j in 1:d3){
      a<-(var(count[,j])-mean(count[,j]))/mean(count[,j])^2
      b<-mean(count[,j])
      initial_values<-c(a ,b)
      est_count[j,]<-nlm(nbinom.lik, initial_values, x=count[,j])$e
    }
  }
  #--------------------------------------------------------
  #--------------- structuring the parameters -------------
  theta_struc <- initPar(copF1)
  init_val <- unlist(as.relistable(theta_struc))
  #---------------------- Estimation ----------------------
  mod <- nlm(loglk_f1, init_val, udat=u, ydat=ordinal, countdat=count,
             copF1, cutp=cutpoints, estcount=est_count,
             gln, glw, theta_struc, param=T, hessian=hessian, print.level = print.level)
  #---------------Restructure the parameters---------------
  deltas=list(relist(mod$e, skeleton = theta_struc))
  #Converting back to copula parameters
  cpar=delta2cpar(deltas, copF1)
  # convert the copula paramters to taus
  taus = rep(NA, d)
  try(taus <- mapply(function(x,y) par2tau(x,y), 
              x=copF1, y=cpar$f1), T)
  #--------------------------------------------------------
  # calculate SE of taus using the Delta method
  if(hessian==TRUE){
    covstruc = (relist(1:length(mod$e), skeleton = theta_struc))
    covariance = rep(NA, d)
    for(i in 1:d){
      try(covariance[i]<-solve(mod$h)[covstruc[[i]]],T)
      covariance[which(lengths(theta_struc)!=2)] = NA
    }
    variance = (relist(NA, skeleton = theta_struc))
    try(variance <- relist(diag(solve(mod$h)), skeleton = theta_struc),T)
    SEtaus = rep(NA, d)
    try(SEtaus <- mapply(function(x,y,z, j) deltamethodSE(x, y, z, j), 
           x=copF1, y=deltas[[1]], z=variance, j = covariance), TRUE)
    
  }else{
    SEtaus=NULL
  }
  #--------------------------------------------------------
  list("cutpoints"=cutpoints,"uni.count.est"=est_count,
       "cpar"=cpar,"taus"=taus,"SEs"=unlist(SEtaus), "loglik"=-mod$m)
} 

mle2factor=function(continuous=NULL, ordinal=NULL, count=NULL, copF1, copF2, gl, hessian=F, print.level=0){
  #Gauss legendre points and weights
  gln<-gl$nodes
  glw<-gl$weights
  d = length(copF1)
  #--------------------------------------------------------
  #                 Continuous variables
  #-----------------                    -------------------
  u <- NULL
  if(!is.null(continuous)){
    n<-nrow(continuous)
    u<-apply(continuous,2,rank)/(n+1)
  }
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  cutpoints<-NULL
  if(!is.null(ordinal)){
    n<-nrow(ordinal)
    cutpoints<-cuts(ordinal)
  }
  #--------------------------------------------------------
  #                   Count variables
  #-------------------               ----------------------
  est_count<-NULL
  if(!is.null(count)){
    n<-nrow(count)
    d3<-ncol(count)
    est_count=matrix(NA,ncol = 2,nrow = d3 )
    for(j in 1:d3){
      a<-(var(count[,j])-mean(count[,j]))/mean(count[,j])^2
      b<-mean(count[,j])
      initial_values<-c(a ,b)
      est_count[j,]<-nlm(nbinom.lik, initial_values, x=count[,j])$e
    }
  }
  #--------------------------------------------------------
  #--------------- structuring the parameters -------------
  theta_struc <- initPar(copF1, copF2)
  init_val <- unlist(as.relistable(theta_struc))
  #---------------------- Estimation ----------------------
  mod <- nlm(loglk_f2, init_val, udat=u, ydat=ordinal, countdat=count,
             copF1=copF1, copF2=copF2, cutp=cutpoints, estcount=est_count,
             gln, glw, theta_struc, param=T, hessian=hessian, print.level = print.level)
  #---------------Restructure the parameters---------------
  deltas=relist(c(mod$e), skeleton = theta_struc)
  #Converting back to copula parameters
  cpar=delta2cpar(deltas, copF1, copF2)
  # convert the copula paramters to taus
  tausf1 = tausf2 = rep(NA, d)
  try(tausf1 <- mapply(function(x,y) par2tau(x,y), 
                x = copF1, y = cpar$f1),T)
  try(tausf2 <- mapply(function(x,y) par2tau(x,y), 
                x = copF2, y = cpar$f2),T)
  #--------------------------------------------------------
  
  if(hessian==TRUE){
    covstruc = (relist(1:length(mod$e), skeleton = theta_struc))
    covariance.f1 = covariance.f2 = rep(NA, d)
    for(i in 1:d){
      try(covariance.f1[i] <- solve(mod$h)[covstruc$f1[[i]]],T)
      covariance.f1[which(lengths(theta_struc$f1)!=2)] = NA
      try(covariance.f2[i] <- solve(mod$h)[covstruc$f2[[i]]],T)
      covariance.f2[which(lengths(theta_struc$f2)!=2)] = NA
    }
    
    variance = (relist(NA, skeleton = theta_struc))
    try(variance <- relist(diag(solve(mod$h)), skeleton = theta_struc),T)
    SEtausf1 = SEtausf2 = rep(NA, d)
    # calculate SE of taus using the Delta method
    try(SEtausf1 <- mapply(function(x,y,z, j) deltamethodSE(x, y, z, j), 
                         x=copF1, y=deltas$f1, z=variance$f1, j = covariance.f1), TRUE)
    try(SEtausf2 <- mapply(function(x,y,z, j) deltamethodSE(x, y, z, j), 
                           x=copF2, y=deltas$f2, z=variance$f2, j = covariance.f2), TRUE)
  }else{
    SEtausf1=NULL
    SEtausf2=NULL
  }
  #--------------------------------------------------------
  list("cutpoints"=cutpoints,"uni.count.est"=est_count,
       "cpar"=cpar,"taus"=c(tausf1, tausf2),
       "SEs"=unlist(c(SEtausf1,SEtausf2)) , "loglik"=-mod$m)
} 


mle2factor.bvn=function(continuous=NULL, ordinal=NULL, count=NULL, copF1, copF2, gl, SpC=NULL, print.level=0){
  #Gauss legendre points and weights
  gln<-gl$nodes
  glw<-gl$weights
  #--------------------------------------------------------
  #                 Continuous variables
  #-----------------                    -------------------
  u <- NULL
  if(!is.null(continuous)){
    n<-nrow(continuous)
    u<-apply(continuous,2,rank)/(n+1)
  }
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  cutpoints<-NULL
  if(!is.null(ordinal)){
    n<-nrow(ordinal)
    cutpoints<-cuts(ordinal)
  }
  #--------------------------------------------------------
  #                   Count variables
  #-------------------               ----------------------
  est_count<-NULL
  if(!is.null(count)){
    n<-nrow(count)
    d3<-ncol(count)
    est_count=matrix(NA,ncol = 2,nrow = d3 )
    for(j in 1:d3){
      a<-(var(count[,j])-mean(count[,j]))/mean(count[,j])^2
      b<-mean(count[,j])
      initial_values<-c(a ,b)
      est_count[j,]<-nlm(nbinom.lik, initial_values, x=count[,j])$e
    }
  }
  
  d1 <- ncol(continuous)
  d2 <- ncol(ordinal)
  d3 <- ncol(count)
  d <- d1 + d2 + d3
  
  #--------------------------------------------------------
  #--------------- structuring the parameters -------------
  theta_struc <- initPar(copF1, copF2)
  rr <- cor(cbind(continuous, ordinal, count))
  f2dat <- factanal(covmat=rr,factors=2)
  tem1 <- c(f2dat$loadings[,1])
  tem2 <- c(f2dat$loadings[,2])
  init_val <-c(tem1, tem2)
  init_val[ (d + SpC) ]=0
  #---------------------- Estimation ----------------------
  mod.bvn <- nlm(loglk_f2_bvn, init_val, udat=u, ydat=ordinal, countdat=count,
             copF1=copF1, copF2=copF2, cutp=cutpoints, estcount=est_count,
             gln, glw, SpC, theta_struc, param=T, print.level = print.level)
  #===== Re-parameterise
  deltas.bvn=relist(mod.bvn$e, skeleton = theta_struc)
  #Converting back to copula parameters
  cpar.bvn=delta2cpar(deltas.bvn, copF1,copF2)
  #cpar.bvn[[2]][[7]]=0
  # Thetas and deltas 
  theta=unlist(cpar.bvn[[1]])
  delta=unlist(cpar.bvn[[2]])
  
  #Rotation Rotation Rotation Rotation Rotation Rotation Rotation Rotation
  ga1=theta
  ga2=delta*sqrt(1-theta^2)
  loadmat=cbind(ga1,ga2)
  rmat=loadmat%*%t(loadmat)
  diag(rmat)=1
  #print(rmat)
  rot1=varimax(loadmat)
  a1=as.matrix(rot1$loadings)
  rmat1=a1%*%t(a1)
  diag(rmat1)=1
  #print(rmat1)
  #permute two columns of loadings after varimax and
  #then convert to partial correlations
  rtheta<-rot1$l[,1]
  rdelta<-rot1$l[,2]/sqrt(1-rtheta^2)
  tau1.bvn=2/pi*asin(rtheta)
  tau2.bvn=2/pi*asin(rdelta)
  #--------------------------------------------------------
  list("cutpoints"=cutpoints,"uni.count.est"=est_count,
       "cpar" = cpar.bvn, "taus" = c(tau1.bvn, tau2.bvn), "loglik"=-mod.bvn$m)
} 


# purpose:
# calulates the matrix with the cutpoints in the uniform scale
# input:
# dat: A data matrix where the number of rows corresponds to an
# individual's response and each column represents an item
# Number of ordinal categories for each item, coded as 0,...,(ncat-1).
# Currently supported are items that have the same number of categories.
# output:
# A matrix with the cutpoints (without the boundary cutpoints) in the
# uniform scale where the number of rows corresponds to an ordinal category
# and each column represents an item
cuts<-function(odat)
{
  n = nrow(odat)
  d = ncol(odat)
  ncat = length(table(odat))
  if (min(odat) == 1) 
    odat = odat - 1
  a = matrix(NA, ncat - 1, d)
  for (j in 1:d) {
    fr = table(c(odat[, j], 0:(ncat - 1)))
    pmf = (fr - 1)/n
    cdf = cumsum(pmf)
    a[, j] = cdf[-ncat]
  }
  a
}
#----------------------------------------------------------
# The negative binomial likelihood
nbinom.lik = function(param, x){ 
  gam <- param[1]
  mu <- param[2]
  if(gam<0 || mu<0 ) {return(1.e10)}
  logl <- sum(dnbinom(x, size=1/gam, mu=mu, log=TRUE))
  -logl
}

#-----------------------------------------------------------
# M2 statistic for one-factor and two-factor
# Input:
# 1) ncontinuous: discretized continuous data.
# 2) ordinal: ordinal data.
# 3) ncount: discretized count data to fewer categories.
# 4) cpar: estimated copula parameters from mle.
# 5) copF1: copula names for first factor.


# Output: 
# 1) M2 statistic.
# 2) Degrees of freedom (df).
# 3) p-value.

M2.1F=function(tcontinuous=NULL, ordinal=NULL, tcount=NULL, cpar, copF1, gl){
  #Gauss legendre nodes and weights -------------------------
  gln<-gl$nodes
  glw<-gl$weights
  #cbind all ordinal data -----------------------------------
  allordinal<-cbind(tcontinuous, ordinal, tcount)
  d<-ncol(allordinal)
  if(d!=length(copF1)){stop("variables and copF1 are not equal in length")}
  
  #cutpoints for allordinal data ----------------------------
  pnorma<-cuts(allordinal)
  
  #categories for allordinal data ---------------------------
  #lengK<-sapply(apply(allordinal, 2, table), length)
  
  # Unlisting the lengths -----------------------------------
  #Kj<-as.numeric(lengK)
  Kj=rep(NA, d)
  for(i in 1:d){
    Kj[i] = length(table(allordinal[,i]))
  }
  #copula functions for M2 ----------------------------------
  m2cops<-M2copulas(d, cop_name_f1 = copF1, cop_name_f2 = NULL)
  
  #----------------------------------------------------------
  #----------------------------------------------------------
  #change negative to positive copula parameters
  thetaF1=as.list(mapply( function(x, y) if(x=="1rjoe" || x=="2rjoe" 
                                            || x=="1rgum"|| x=="2rgum"){abs(y)}else{y},
                          x=copF1, y=cpar[[1]]))
  #----------------------- M2 statistic ---------------------
  sdat<-store.1fact(dat=allordinal, pnorma=pnorma, theta=thetaF1,
                    pcondcop=m2cops$pcondcop$f1, pconddotcop=unlist(m2cops$pconddotcop$f1),
                    dcop=m2cops$dcop$f1, gln=gln, param=F)
  gof_m2<-M2.1fact(dat=allordinal, dnorma=sdat$dnorma, prob1=sdat$prob1, 
                   cden=sdat$cden, fden=sdat$fden, fdotden=sdat$fdotden, 
                   theta=thetaF1, Kj=Kj, glw=glw)
  #----------------------------------------------------------
  #----------------------------------------------------------
  return(gof_m2)
}

# 6) copF2: provide the  copula names 
#           to estimate M2 for the two-factor copula model.

M2.2F=function(tcontinuous=NULL, ordinal=NULL, tcount=NULL, cpar, copF1, copF2, gl, SpC=NULL){
  #Gauss legendre nodes and weights -------------------------
  gln<-gl$nodes
  glw<-gl$weights
  #cbind all ordinal data -----------------------------------
  allordinal<-cbind(tcontinuous, ordinal, tcount)
  d<-ncol(allordinal)
  if(d!=length(copF1)){stop("variables and copF1 are not equal in length")}
  
  #cutpoints for allordinal data ----------------------------
  pnorma<-cuts(allordinal)
  
  #categories for allordinal data ---------------------------
  #lengK<-sapply(apply(allordinal, 2, table), length)
  
  # Unlisting the lengths -----------------------------------
  #Kj<-as.numeric(lengK)
  Kj=rep(NA, d)
  for(i in 1:d){
    Kj[i] = length(table(allordinal[,i]))
  }
  #copula functions for M2 ----------------------------------
  m2cops<-M2copulas(d, cop_name_f1 = copF1, cop_name_f2 = copF2)
  
  #----------------------------------------------------------
  #----------------------------------------------------------
  #change negative to positive copula parameters
  thetaF1=as.list(mapply( function(x, y) 
    if(x=="1rjoe" || x=="2rjoe" || x=="1rgum"|| x=="2rgum"){
      abs(y)}else{y}, x=copF1, y=cpar[[1]]))
  thetaF2=as.list(mapply( function(x, y) 
    if(x=="1rjoe" || x=="2rjoe" || x=="1rgum"|| x=="2rgum"){
      abs(y)}else{y}, x=copF2, y=cpar[[2]]))
  #----------------------- M2 statistic ---------------------
  sdat<-store.2fact(dat=allordinal, pnorma=pnorma,
                    theta=thetaF1, delta=thetaF2, 
                    pcondcop1=m2cops$pcondcop$f1,
                    pcondcop2=m2cops$pcondcop$f2,
                    pconddotcop1=unlist(m2cops$pconddotcop$f1),
                    pconddotcop2=unlist(m2cops$pconddotcop$f2),
                    dcop1=m2cops$dcop$f1, dcop2=m2cops$dcop$f2,
                    gln=gln, param=F)
  
  gof_m2<-M2.2fact(dat=allordinal, dnorma=sdat$dnorma, prob1=sdat$prob1,
                   cden2=sdat$cden2, fden2=sdat$fden2,
                   fdotden2=sdat$fdotden2, fbarden2=sdat$fbarden2,
                   theta=thetaF1, delta=thetaF2, Kj=Kj, glw=glw, SpC, iprint=T)
  #----------------------------------------------------------
  #----------------------------------------------------------
  return(gof_m2)
}


#-------------------------------------------
# simulations for one-factor copula
#-------------------------------------------
# Input:
# (1) n: sample size.
# (2) d1: number of continuous variables.
# (3) d2: number of ordinal variables.
# (4) categ: number categories for ordinal variables.
# (5) theta: copula parameters in a list form.
# (6) copF1: names of bivaraite copulas.
r1factor=function(n, d1=0, d2=0, categ, theta, copF1){
  #Quantile of conditional copulas
  cops=copulas(d1, d2, d3=0, cop_name_f1=copF1, cop_name_f2=NULL)
  #Simulate
  d=d1+d2
  #change negative to positive copula parameters
  thetaF1=as.list(mapply( function(x, y) 
    if(x=="1rjoe" || x=="2rjoe" || x=="1rgum"|| x=="2rgum"){
      list(abs(y))}else{list(y)}, x=copF1, y=theta))
  u=r1Fcopula(n, d, thetaF1, cops$qcondcop$f1)
  if(d1!=0){
    udat=u[,1:d1]
  }else{
    udat=NULL
  }
  #Converting to ordinal data
  if(d2!=0){
    udat.d2 = cbind(u[ ,(d1+1):(d1+d2)])
    ydat <- matrix(NA, ncol = d2, nrow = n )
    for(j in 1:d2){
      ydat[,j] = continuous2ordinal.sim(udat.d2[, j], categ[j])
    }
  }else{
    ydat=NULL
  }
  simulated=cbind(udat, ydat)
  return(simulated)
}

#-------------------------------------------
# simulations for two-factor copula
#-------------------------------------------
# Input:
# (1) n: sample size.
# (2) d1: number of continuous variables.
# (3) d2: number of ordinal variables.
# (4) categ: number categories for ordinal variables.
# (5) theta: copula parameters in a list form for 1st factor.
# (6) delta: copula parameters in a list form for 2nd factor.
# (7) copF1: names of bivaraite copulas for 1st factor.
# (8) copF2: names of bivaraite copulas for 2nd factor.

r2factor=function(n, d1=0, d2=0, categ, theta, delta, 
                  copF1, copF2){
  #Quantile of conditional copulas
  cops=copulas(d1, d2, d3=0, cop_name_f1=copF1, cop_name_f2=copF2)
  #Simulate
  d=d1+d2
  #change negative to positive copula parameters
  thetaF1=as.list(mapply( function(x, y) 
    if(x=="1rjoe" || x=="2rjoe" || x=="1rgum"|| x=="2rgum"){
      list(abs(y))}else{list(y)}, x=copF1, y=theta))
  thetaF2=as.list(mapply( function(x, y) 
    if(x=="1rjoe" || x=="2rjoe" || x=="1rgum"|| x=="2rgum"){
      list(abs(y))}else{list(y)}, x=copF2, y=delta))
  
  u=r2Fcopula(n, d, thetaF1, thetaF2, cops$qcondcop$f1, cops$qcondcop$f2)
  if(d1!=0){
    udat=u[,1:d1]
  }else{
    udat=NULL
  }
  #Converting to ordinal data
  if(d2!=0){
    udat.d2 = cbind(u[ ,(d1+1):(d1+d2)])
    ydat <- matrix(NA, ncol = d2, nrow = n )
    for(j in 1:d2){
      ydat[,j] = continuous2ordinal.sim(udat.d2[, j], categ[j])
    }
  }else{
    ydat=NULL
  }
  simulated=cbind(udat, ydat)
  return(simulated)
}

