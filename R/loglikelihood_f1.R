# one-factor copula for mixed data. 

# purpose: the density function for one-factor copula for mixed data.
# input: 
# udat: continuous data on uniform scale.
# ydat: discrete data on ordinal scale.
# countdat: count data.
# dcop: density function for continuous data
# pcondcop: conditional cdf copula conditioned on latent variables.
# cutp: the cutpoints, and is calculated using the above function.
# theta: copula parameters,the first parameters are for the continuous dat.
# and followed with the parameters for the ordinal and count data
# gln:  Gauss Legendre quadrature points.
# glw:  Gauss Legendre quadrature weights.
# param:  logical arg (T or F) for parameterization.
dfcopula_f1=function(udat,ydat,countdat,dcop,pcondcop,
                     cutp,estcount,th,gln,glw,param){
  nq<-length(gln)
#------------------------------------------------------
#Continuous::::::
  if(!is.null(udat)){
  d1<-ncol(udat)#number of continuous variables
  n<-nrow(udat)#number of observations
  th_u<-th[1:d1]#copula parameters for continuous
  dcop_u=dcop[1:d1]#copula functions for continuous 
  }else{
    d1<-0
  }
#------------------------------------------------------
#Ordinal::::::
  if(!is.null(ydat)){
  n<-nrow(ydat)#number of observations
  K<-nrow(cutp)+1#number of categories
  d2<-ncol(cutp)#number of items
  th_ord<-th[(d1+1):(d1+d2)]#copula parameters for ordinal
  #Copula functions for ordinal variables
  pcondcopOrd=pcondcop[(d1+1):(d1+d2)]
#------------------------------------------------------
#----- pmf for ordinal-----------------------------
  condcdf<-array(NA,dim=c(nq,K-1,d2))
  for(j in 1:d2){
    pcondcopOrdinal<-pcondcopOrd[[j]]
    for(k in 1:(K-1)){
    condcdf[,k,j]=pcondcopOrdinal(cutp[k,j],gln,th_ord[[j]],param) 
    }
  }
#empty arrays for the ordinal matrices-----------------
  arr0<-array(0,dim=c(nq,1,d2))
  arr1<-array(1,dim=c(nq,1,d2))
  condcdf<-abind(arr0,condcdf,arr1,along=2)
#----- pmf for ordinal---------------------------------
  pmf<-array(NA,dim=c(nq,K,d2))
  for(j in 1:d2){
    for(k in 1:K){
      pmf[,k,j]=condcdf[,k+1,j]-condcdf[,k,j] 
    }
  }
  }else{
    d2<-0
  }
#------------------------------------------------------
#Count::::::
  if(!is.null(countdat)){
  d3<-ncol(countdat)
  n<-nrow(countdat)
  Fy<-matrix(NA,ncol = d3,nrow = n)
  Fy_1<-matrix(NA,ncol = d3,nrow = n)
  for(j in 1:d3){
  c1_j<-estcount[j,1]
  c2_j<-estcount[j,2]
  Fy[,j]=pnbinom(countdat[,j],size=1/c1_j,mu=c2_j)
  Fy_1[,j]=pnbinom(countdat[,j]-1,size=1/c1_j,mu=c2_j)
  }
#Copula parameters for count
  th_count=th[(d1+d2+1):(d1+d2+d3)]
#Copula functions for count
  pcondcopCount=pcondcop[(d1+d2+1):(d1+d2+d3)]
  }
#------------------------------------------------------
# the product of the factor model----------------------
  fproduct = matrix(NA,n,nq)
  for(i in 1:n){
#------------------------------------------------------
#Continuous--------------------------------------------
    prod_u<-1
    if(!is.null(udat)){
    uvec=udat[i,]
    for(j in 1:d1){
      dcop_j=dcop_u[[j]]
      temp=dcop_j(uvec[j],gln,th_u[[j]],param)
      prod_u=prod_u*temp
    }
    }
#------------------------------------------------------
#Ordinal-----------------------------------------------
    prod_y<-1
    if(!is.null(ydat)){
    for(j in 1:d2)
    {
      res<-pmf[,ydat[i,j]+1,j]
      prod_y<-prod_y*res
    }
    }
#------------------------------------------------------
#Count-------------------------------------------------
    prod_c<-1
    if(!is.null(countdat)){
    for(j in 1:d3){
      th_count_j<-th_count[[j]]
      pcondcopCount_j<-pcondcopCount[[j]]
      condcdf_count<-(pcondcopCount_j(Fy[i,j],gln,th_count_j,param)-
                        pcondcopCount_j(Fy_1[i,j],gln,th_count_j,param))
      prod_c<-prod_c*condcdf_count
    }
    }
#------------------------------------------------------
    fproduct[i, ]=prod_u*prod_y*prod_c
  }
#------------------------------------------------------
  (fproduct%*%glw)
}

# purpose: calculates log-likelihood for the copula model for mixed data.
# input:
# th: copula parameters,the first parameters are for the continuous dat.
# and followed with the parameters for the ordinal and count data
# udat: continuous data on uniform scale.
# ydat: ordinal data.
# countdat: count data.
# copF1: the name of the copulas for the one-factor model.
# cutp: the cutpoints for ordinal items.
# estcount: est parameters for count data.
# gln:  Gauss Legendre quadrature points.
# glw:  Gauss Legendre quadrature weights.
# param:  logical arg (T or F) for parameterization.
loglk_f1= function(theta,udat=NULL,ydat=NULL,countdat=NULL,copF1,cutp=NULL,
                      estcount=NULL,gln,glw,theta_struc,param){
  d1=ifelse(is.null(udat),0,ncol(udat))
  d2=ifelse(is.null(ydat),0,ncol(ydat))
  d3=ifelse(is.null(countdat),0,ncol(countdat))
  cops=copulas(d1,d2,d3,cop_name_f1=copF1,cop_name_f2=NULL)
  cpar=relist(theta,skeleton=theta_struc)
  lk=dfcopula_f1(udat,ydat,countdat,dcop=cops$dcop$f1,
                   pcondcop=cops$pcondcop$f1,
                   cutp,estcount,th=cpar,gln,glw,param)
  if (any(lk <= 0) || any(is.nan(lk)) || any(is.infinite(lk)) )
  {return(1e+10)}
  loglik= log(lk)
  return(- sum(loglik))
} 

#nlm(loglk_f1,initial.param,udat=u,ydat=ordinal,
#    countdat=count,copF1=rep("bvn",d),
#    cutp=cutpoints,estcount=est.count,gln,glw,theta_struc=parameters,
#    param=T,print.level = 2)

# The likelihood for to be used in vuong statistics
lk_f1= function(theta,udat=NULL,ydat=NULL,countdat=NULL,copF1,cutp=NULL,
                   estcount=NULL,gln,glw,theta_struc,param=T)
{
  d1=ifelse(is.null(udat),0,ncol(udat))
  d2=ifelse(is.null(ydat),0,ncol(ydat))
  d3=ifelse(is.null(countdat),0,ncol(countdat))
  cops=copulas(d1,d2,d3,cop_name_f1=copF1,cop_name_f2=NULL)
  cpar=relist(theta,skeleton=theta_struc)
  lk=dfcopula_f1(udat,ydat,countdat,dcop=cops$dcop$f1,
                 pcondcop=cops$pcondcop$f1,
                 cutp,estcount,th=cpar,gln,glw,param)
  return(lk)
} 
