# Wrapper functions for
# (1) MLE bi-factor copula model.
# (2) MLE second-order copula model.
# (3) M2 statistic for the bi-factor and second-order copula models.
# (4) Simulation code for Bi-factor and Second-order.
# (5) Vuong's test for bi-factor and second-order copula models.
#------------------------------------------------------------
# MLE Estimation of one-factor model
# Input:
# 1) y: ordinal data.
# 2) copnames1: copula names for the items with common factor.
# 3) copnames2: copula names for the items with group-specifc factors.
# 4) gl: Gauss legendre quadrature points and weights.
# 5) ngrp: number of groups.
# 6) grpsize: vector of size for each group.

# Outputs:
# 1) "cutpoints": cutpoints for ordinal data.
# 2) "loglik": log-likelihood value.
# 3) "taus": estimated Kendall taus.
# 4) "SEs": Standard errors of the estimated taus.
 
mleBifactor=function(y, copnames1, copnames2, gl, ngrp, grpsize, hessian=F, print.level=0){
  #Gauss legendre points and weights
  gln<-gl$nodes
  glw<-gl$weights
  nq<-length(glw)
  d<-length(copnames1)
  n<-nrow(y)
  cutpoints<-cuts(y)
  #--------------------------------------------------------
  #--------------- structuring the parameters -------------
  initpar=initval(copnames1,copnames2)
  #EST
  mle<-nlm(bifactorllk,initpar,ydat=y,cutp=cutpoints,
           copX0=copnames1,copXgY=copnames2,
           nq=nq,ngrp=ngrp, grpsize=grpsize, glw=glw, gln=gln, SpC=NULL, param=F,
           print.level=print.level, hessian=hessian)
  #--------------------------------------------------------
  # calculate SE of taus using the Delta method
  if(hessian==TRUE){
    SE = rep(NA, d*2)
      try(SEtaus <-  sqrt(diag(solve(mle$h))) , T)
  }else{
    SEtaus=NULL
  }
  #--------------------------------------------------------
  list("cutpoints"=cutpoints,"taus"=mle$e,
       "SEs"=SEtaus, "loglik"=-mle$m)
}


mleSecond_order=function(y, copnames1, copnames2, gl, ngrp, grpsize, hessian=F, print.level=0){
  #Gauss legendre points and weights
  gln<-gl$nodes
  glw<-gl$weights
  nq<-length(glw)
  d<-length(copnames2)
  n<-nrow(y)
  cutpoints<-cuts(y)
  #--------------------------------------------------------
  #--------------- structuring the parameters -------------
  initpar=initval(copnames1,copnames2)
  #EST
  mle<-nlm(nestedllk,p=initpar,ydat=y,cutp=cutpoints,
           copX0=copnames1,copXgY=copnames2,
           nq=nq,ngrp=ngrp,grpsize=grpsize,glw=glw, gln=gln,SpC=NULL,
           hessian=hessian, print.level=print.level)
  #--------------------------------------------------------
  # calculate SE of taus using the Delta method
  if(hessian==TRUE){
    SE = rep(NA, d+ngrp)
    try(SEtaus <-  sqrt(diag(solve(mle$h))) , T)
  }else{
    SEtaus=NULL
  }
  #--------------------------------------------------------
  list("cutpoints"=cutpoints,"taus"=mle$e,
       "SEs"=SEtaus, "loglik"=-mle$m)
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



#-----------------------------------------------------------
# M2 statistic for bi-factor and second-order copula models
# Input:
# 2) y: ordinal data.
# 4) cpar: estimated copula parameters from mle.
# 5) copnames1: copula names for common factor.
# 5) copnames2: copula names for group factor.


# Output:
# 1) M2 statistic.
# 2) Degrees of freedom (df).
# 3) p-value.

M2Bifactor=function(y,cpar, copnames1, copnames2, gl, ngrp, grpsize){
  #Gauss legendre nodes and weights -------------------------
  gln<-gl$nodes
  glw<-gl$weights
  nq<-length(gln)
  d<-ncol(y)

  #cutpoints for ordinal data ----------------------------
  cutpoints<-cuts(y)

  #copula functions for M2 ----------------------------------
  m2copnm<-M2copulas_structured(copnames1,copnames2)
  pcondX0<-m2copnm$pcondcopX0
  pcondXg<-m2copnm$pcondcopXg
  pconddotX0<-m2copnm$pcond2copX0
  pconddotXg<-m2copnm$pcond2copXg
  dcopXg<-m2copnm$dcopXg
  dcopX0<-m2copnm$dcopX0

  cpar=mapply( function(x, y)
    if(x=="1rjoe" || x=="2rjoe" || x=="1rgum"|| x=="2rgum"){
      abs(y)}else{y}, x=c(copnames1,copnames2), y=cpar)

  #----------------------- M2 statistic ---------------------
  sdat<-store.bifact(cutp=cutpoints,cpar[1:d],cpar[(d+1):(2*d)],
                     pcondX0,pcondXg,pconddotX0,pconddotXg,
                     dcopX0,dcopXg,nq,glw,gln)
  gof_m2<-M2.bifactor(dat=y,sdat$dnorma,sdat$prob1,sdat$cden2,sdat$fden2,
                            sdat$fdotden2,sdat$fbarden2,ngrp,grpsize,glw,SpC=NULL,iprint=F)
  #----------------------------------------------------------
  #----------------------------------------------------------
  return(gof_m2)
}


M2Second_order=function(y,cpar, copnames1, copnames2, gl, ngrp, grpsize){
  #Gauss legendre nodes and weights -------------------------
  gln<-gl$nodes
  glw<-gl$weights
  nq<-length(gln)
  d<-ncol(y)

  #cutpoints for ordinal data ----------------------------
  cutpoints<-cuts(y)

  #copula functions for M2 ----------------------------------
  m2copnm<-M2copulas_structured(copnames1,copnames2)
  pcondY<-m2copnm$pcondcopXg
  pconddotY<-m2copnm$pcond2copXg
  dcopY<-m2copnm$dcopXg
  dcopX<-m2copnm$dcopX0
  ddotcopX<-m2copnm$d2copX0

  cpar=mapply( function(x, y)
    if(x=="1rjoe" || x=="2rjoe" || x=="1rgum"|| x=="2rgum"){
      abs(y)}else{y}, x=c(copnames1,copnames2), y=cpar)

  #----------------------- M2 statistic ---------------------
  sdat<-store.nested(cutpoints,thX0Xg=cpar[1:ngrp],thXgYj=cpar[(1+ngrp):(d+ngrp)],
                     pcondY,pconddotY,dcopY,
                     dcopX,ddotcopX,nq,ngrp,glw,gln)
  gof_m2<-M2.nested(dat=y,sdat$dnorma,sdat$prob1,sdat$cden,sdat$pmf,
                    sdat$pmfdot,sdat$fden,sdat$fdotden,ngrp,
                    grpsize,glw,gln,SpC=NULL,iprint=F)
  #----------------------------------------------------------
  #----------------------------------------------------------
  return(gof_m2)
}


#-------------------------------------------
# simulations for bi-factor copula
#-------------------------------------------
# Input:
# (1) n: sample size.
# (3) d: number of ordinal variables.
# (4) categ: number categories for ordinal variables.
# (5) theta1: copula parameters for common factors.
# (5) theta2: copula parameters for group-specific factors.
# (6) copnames1: names of bivaraite copulas for common factors.
# (7) copnames2: names of bivaraite copulas for group-specific factors.
rBifactor=function(n, d, grpsize, categ, copnames1,copnames2, theta1, theta2){
  #Quantile of conditional copulas
  cops=copulas_bifactor(copnames1,copnames2)

  theta=mapply( function(x, y)
    if(x=="1rjoe" || x=="2rjoe" || x=="1rgum"|| x=="2rgum"){
      abs(y)}else{y}, x=c(copnames1,copnames2), y=c(theta1,theta2))

  u=simulation.bifactor(n , grpsize , cops$qcondcopX0, cops$pcondcopXgY,
                        theta)
  #Converting to ordinal data
    ydat <- matrix(NA, ncol = d, nrow = n )
    for(j in 1:d){
      ydat[,j] = continuous2ordinal.sim(u[, j], categ[j])
    }
  return(ydat)
}

#-------------------------------------------
# simulations for second-order factor copula
#-------------------------------------------
# Input:
# (1) n: sample size.
# (3) d: number of ordinal variables.
# (4) categ: number categories for ordinal variables.
# (5) theta1: copula parameters for common and group-specifc factor.
# (5) theta2: copula parameters for items with group factor.
# (6) copnames1: names of bivaraite copulas for common and group-specifc factor.
# (6) copnames2: names of bivaraite copulas for items with group factor.

rSecond_order=function(n, d, grpsize, categ, copnames1, copnames2, theta1, theta2){
  #Quantile of conditional copulas
  cops=copulas_nested(copnames1,copnames2)

  theta=mapply( function(x, y)
    if(x=="1rjoe" || x=="2rjoe" || x=="1rgum"|| x=="2rgum"){
      abs(y)}else{y}, x=c(copnames1,copnames2), y=c(theta1,theta2))

  u=simulation.nested(n, grpsize, cops$qcopX0Xg,
                      cops$qcondXgY, theta)

  #Converting to ordinal data
  ydat <- matrix(NA, ncol = d, nrow = n)
  for(j in 1:d){
    ydat[,j] = continuous2ordinal.sim(u[, j], categ[j])
  }
  return(ydat)
}



#Model selection algorithm for the bi-factor copula model
selectBifactor = function(y, grpsize, copnames, gl){
  gln <- gl$nodes
  glw <- gl$weights
  nq<-length(gln)
  n <- nrow(y)
  cutp <- cuts(y)
  d <- ncol(y)
  #number of groups
  ngrp=length(grpsize)

  selected.model=Select.bifactor(ydat = y, cutp=cutp , ngrp=ngrp,
                                 grpsize=grpsize,
                                 allcop=copnames, gln=gln, glw=glw,
                                 nq=nq, n=n, d=d,
                                 SpC = NULL)
  selected.model
}
#cop=c("bvn","bvt2","bvt3","bvt4","bvt5","bvt6",
#"bvt7","bvt8","bvt9","frk","gum","sgum")
#copulas.names=selectBifactor(ydat, grpsize,
#copnames=cop, gl, SpC = NULL)

#Model selection algorithm for the second-order copula model
selectSecond_order = function(y, grpsize, copnames, gl){

  gln <- gl$nodes
  glw <- gl$weights
  nq <- length(gln)
  n <- nrow(y)
  cutp <- cuts(y)
  d <- ncol(y)
  #number of groups
  ngrp=length(grpsize)

  selected.model=Select.nested(ydat=y,cutp=cutp, ngrp=ngrp,
                               grpsize=grpsize, allcop=copnames,
                               gl=gl,gln=gln, glw=glw,nq=nq,n=n,
                               d=d, SpC = NULL)
  selected.model
}

#cop=c("bvn","bvt2","bvt3","bvt4","bvt5","bvt6",
#"bvt7","bvt8","bvt9","frk","gum","sgum")
#copulas.names=select2ndOrder(ydat, grpsize,
#copnames=cop, gl, SpC = NULL)


#---------------------------------------- Vuong's test
#models: choose a number in [1,...,4]
#         1: M1:bifactor  , M2:bifactor
#         2: M1:nested    , M2:nested
#         3: M1:nested    , M2:bifactor
#         4: M1:bifactor  , M2:nested
#cpar.m1: vector of copula paramters for model 1, starting with copula parameters that link the 
# items with common factor (bifactor),
# or group factors with common facotr (second-order).
#cpar.m2: vector of copula paramters for model 2, starting with copula parameters that link the 
# items with common factor (bifactor),
# or group factors with common facotr (second-order).
# copnames.m1: vector of copula families names for model 1, starting with copulas that link the 
# items with common factor (bifactor),
# or group factors with common facotr (second-order).
# copnames.m2: vector of copula families names for model 2, starting with copulas that link the 
# items with common factor (bifactor),
# or group factors with common facotr (second-order).
#y: items response data
#grpsize: vector contains sizes of group items. 

vuong_structured = function(models, cpar.m1, copnames.m1, 
                 cpar.m2, copnames.m2, 
                 y, grpsize){
  d=ncol(y)
  ngrp=length(grpsize)
  
  if(models==1){
    copF0.M1 = copnames.m1[1:d]
    copFg.M1 = copnames.m1[(d+1):(2*d)]
    copF0.M2 = copnames.m2[1:d]
    copFg.M2 = copnames.m2[(d+1):(2*d)]
    out=vuong.bifactor(cpar.M1=cpar.m1,copF0.M1, copFg.M1, cpar.M2=cpar.m2, copF0.M2, copFg.M2, 
                       ordinal=y,ngrp=ngrp, grpsize=grpsize, nq=25, param=F)
  }else if(models==2){
    copF0.M1 = copnames.m1[1:ngrp]
    copFg.M1 = copnames.m1[(ngrp+1):(ngrp+d)]
    copF0.M2 = copnames.m2[1:ngrp]
    copFg.M2 = copnames.m2[(ngrp+1):(ngrp+d)]
    out=vuong.nested(cpar.M1=cpar.m1,copF0.M1, copFg.M1, cpar.M2=cpar.m2, copF0.M2, copFg.M2, 
                     ordinal=y,ngrp=ngrp, grpsize=grpsize, nq=25, param=F)
    
  }else if(models==3){
    copF0.M1 = copnames.m1[1:ngrp]
    copFg.M1 = copnames.m1[(ngrp+1):(ngrp+d)]
    copF0.M2 = copnames.m2[1:d]
    copFg.M2 = copnames.m2[(d+1):(2*d)]
    out=vuong.bifactor.nested(cpar.M1=cpar.m1,copF0.M1, copFg.M1, cpar.M2=cpar.m2, copF0.M2, copFg.M2, 
                              ordinal=y,ngrp=ngrp, grpsize=grpsize, nq=25, param=F)
    
  }else if(models==4){
    copF0.M2 = copnames.m2[1:ngrp]
    copFg.M2 = copnames.m2[(ngrp+1):(ngrp+d)]
    copF0.M1 = copnames.m1[1:d]
    copFg.M1 = copnames.m1[(d+1):(2*d)]
    out=vuong.nested.bifactor(cpar.M1=cpar.m1,copF0.M1, copFg.M1, cpar.M2=cpar.m2, copF0.M2, copFg.M2, 
                              ordinal=y,ngrp=ngrp, grpsize=grpsize, nq=25, param=F)
    
  }else{
    out="select a viable number of models in Vuong's test"
  }
  return(out)
}