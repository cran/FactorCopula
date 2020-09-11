# This code is for two-factor copula for mixed data. 


# purpose: the density function for two-factor copula for mixed data.
# input: 
# udat: continuous data on uniform scale
# ydat: discrete data on ordinal scale
# dcop: density function for continuous data
# pcondcop: conditional cdf copula conditioned on latent variables
# pcondcop1: conditional cdf copula conditioned on latent variables for one-factor copula
# pcondcop2: conditional cdf copula conditioned on latent variables for two-factor copula
# cutp: the cutpoints, and is calculated using the above function
# theta: parameters for the mixed data,the first parameters are for 
# the first factor for continuous dat,
# and followed with the parameters for the first factor for ordinal data,
# and similarly for the second factor.
# nq: number of Gauss Legendre quadrature points.
dfcopula_f2=function(udat,ydat,countdat,d,dcopf1,dcopf2,
            pcondcopf1,pcondcopf2,cutp,estcount,th,gln,glw,param){
  nq<-length(gln)
  #------------------------------------------------------
  #Continuous::::::
  if(!is.null(udat)){
    d1<-ncol(udat)#number of continuous variables
    n<-nrow(udat)#number of observations
    #copula parameters for continuous
    th_u1=th[[1]][1:d1]#copula parameters  for 1st factor.
    th_u2=th[[2]][1:d1]#copula parameters  for 2nd factor.
    #copula functions for continuous 
    dcop_u_f1=dcopf1[1:d1];
    dcop_u_f2=dcopf2[1:d1]
    pcondcop_u_f1=pcondcopf1[1:d1];
    pcondcop_u_f2=pcondcopf2[1:d1]
  }else{
    d1<-0
  }
  #------------------------------------------------------
  #Ordinal::::::
  if(!is.null(ydat)){
    n<-nrow(ydat)#number of observations
    K<-nrow(cutp)+1#number of categories
    d2<-ncol(cutp)#number of items
    #copula parameters for ordinal
    th_y1=th[[1]][(d1+1):(d1+d2)]#copula parameters  for 1st factor.
    th_y2=th[[2]][(d1+1):(d1+d2)]#copula parameters  for 2nd factor.
    #Copula functions for ordinal variables
    pcondcop_y_f1=pcondcopf1[(d1+1):(d1+d2)];
    pcondcop_y_f2=pcondcopf2[(d1+1):(d1+d2)]
    #------------------------------------------------------
    #----- pmf for ordinal-----------------------------
    # empty arrays to save the conditional one-factor copula in condcdf  
    condcdf<-array(NA,dim=c(nq,K-1,d2))
    condcdf2<-array(NA,dim=c(nq,nq,K-1,d2))
    for(j in 1:d2){
      cop_yj1=pcondcop_y_f1[[j]]
      cop_yj2=pcondcop_y_f2[[j]]
      for(k in 1:(K-1)){
        condcdf[,k,j]<-cop_yj1(cutp[k,j],gln,th_y1[[j]],param) 
        for(q in 1:nq)
        { 
          condcdf2[,q,k,j]<-cop_yj2(condcdf[q,k,j],gln,th_y2[[j]],param)
        }
      }
    }
    # binding the array of condcdf2 between 0 & 1.
    arr0<-array(0,dim=c(nq,nq,1,d2))
    arr1<-array(1,dim=c(nq,nq,1,d2))
    condcdf2<-abind(arr0,condcdf2,arr1,along=3)
    # calculate the finite difference for the discrete data. 
    fden2<-array(NA,dim=c(nq,nq,K,d2))
    for(j in 1:d2){
      for(k in 1:K)
      { 
        fden2[,,k,j]<-condcdf2[,,k+1,j]-condcdf2[,,k,j] 
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
    th_c1= th[[1]][(d1+d2+1):(d)]#copula parameters  for 1st factor.
    th_c2= th[[2]][(d1+d2+1):(d)]#copula parameters  for 1st factor.
    #Copula functions for count
    pcondcop_c_f1=pcondcopf1[(d1+d2+1):(d1+d2+d3)]
    pcondcop_c_f2=pcondcopf2[(d1+d2+1):(d1+d2+d3)]
  }
  # empty arrays to save density for one-factor copula for 
  # cont. dat  in dfactor1 and  dfactor2
  # and the conditional one-factor copula in condcop.
  dfactor_uf1=array(NA,dim=c(nq));
  condcop_uf1=array(NA,dim=c(nq))
  dfactor_uf2=array(NA,dim=c(nq,nq))

  # empty arrays to save the conditional 1-factor and 2-factor copula for count:
  condcdf1_cij<- array(NA,dim=c(nq))
  condcdf1_cij_1<- array(NA,dim=c(nq))
  condcdf2_cij<- array(NA,dim=c(nq,nq))
  
  # an empty matrix to save the pre-final output of the density for two-factor 
  # copula for mixed data, and the final step is to multiply the matrix with
  # the second Gauss Legendre quadrature weight points.
  
  #------------------------------------------------------
  # the product of the 2-factor model--------------------
  fproduct = matrix(NA,n,nq)
  for(i in 1:n){ 
    #------------------------------------------------------
    #Continuous--------------------------------------------
    prod_u=1;
    if(!is.null(udat)){
    uvec=udat[i,]
    for(j in 1:d1){
      uvec_j=uvec[j]
      # the jth copula parameter
      th_u1j=th_u1[[j]];th_u2j=th_u2[[j]]
      # the jth copula function
      pcondcop_u1j=pcondcop_u_f1[[j]]
      dcopula_f1j=dcopf1[[j]];dcopula_f2j=dcopf2[[j]]
      #================================================
      #calcualting the conditional and density copula for 1st factor
      condcop_uf1<-pcondcop_u1j(uvec_j,gln,th_u1j,param)
      dfactor_uf1<-dcopula_f1j(uvec_j,gln,th_u1j,param)
      #================================================
      for(q in 1:nq){
        dfactor_uf2[,q] = dcopula_f2j(condcop_uf1[q],gln,th_u2j,param)*(dfactor_uf1[q])
      }
      temp =  dfactor_uf2
      prod_u = prod_u * temp
    }
    }
    #------------------------------------------------------
    #Ordinal-----------------------------------------------
    prod_y=1
    if(!is.null(ydat)){
      for(j in 1:d2){
        res=fden2[,,ydat[i,j]+1,j] 
        prod_y = prod_y * res 
      }
    }
    #------------------------------------------------------
    #Count-------------------------------------------------
    prod_c<-1
    if(!is.null(countdat)){
      prod_c=1
      for(j in 1:d3){
        cdat=countdat[i,j]
        th_c1j=th_c1[[j]];th_c2j=th_c2[[j]]
        pcondcop_c1j= pcondcop_c_f1[[j]]
        pcondcop_c2j= pcondcop_c_f2[[j]]
        condcdf1_cij<-pcondcop_c1j(Fy[i,j],gln,th_c1j,param)
        condcdf1_cij_1<-pcondcop_c1j(Fy_1[i,j],gln,th_c1j,param)
        for(q in 1:nq){  
          condcdf2_cij[,q]<-(pcondcop_c2j(condcdf1_cij[q],gln,th_c2j,param)-
                 pcondcop_c2j(condcdf1_cij_1[q],gln,th_c2j,param))
          }
        temp = condcdf2_cij
        prod_c = prod_c * temp
      }
    }
  #------------------------------------------------------
    fproduct[i, ]=as.vector((prod_u*prod_y*prod_c) %*% glw)
  }
  #------------------------------------------------------
  fproduct %*% glw
  #------------------------------------------------------
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## log likelihood function
loglk_f2 = function(theta,udat=NULL,ydat=NULL,countdat=NULL,copF1,copF2,cutp=NULL,
                       estcount=NULL,gln,glw,theta_struc,param=T){
  d1=ifelse(is.null(udat),0,ncol(udat))
  d2=ifelse(is.null(ydat),0,ncol(ydat))
  d3=ifelse(is.null(countdat),0,ncol(countdat))
  d=d1+d2+d3
  cops=copulas(d1,d2,d3,cop_name_f1=copF1,cop_name_f2=copF2)
  cpar=relist(theta,skeleton=theta_struc)
  lk=dfcopula_f2(udat,ydat,countdat,d,dcopf1=cops$dcop$f1,dcopf2=cops$dcop$f2,
                   pcondcopf1=cops$pcondcop$f1,pcondcopf2=cops$pcondcop$f2,
                   cutp,estcount,th=cpar,gln,glw,param)
  if(any(lk <= 0) || any(is.nan(lk)) || any(is.infinite(lk))){
    return(1e+10)
    }
  loglik= log(lk)
  return(- sum(loglik))
} 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## log likelihood function for bvn model
loglk_f2_bvn = function(theta,udat=NULL,ydat=NULL,countdat=NULL,copF1,copF2,cutp=NULL,
                    estcount=NULL,gln,glw,SpC,theta_struc,param=T){
  d1=ifelse(is.null(udat),0,ncol(udat))
  d2=ifelse(is.null(ydat),0,ncol(ydat))
  d3=ifelse(is.null(countdat),0,ncol(countdat))
  d=d1+d2+d3
  cops=copulas(d1,d2,d3,cop_name_f1=copF1,cop_name_f2=copF2)
   
  theta[ (d + SpC) ]=0
  cpar=relist(theta,skeleton=theta_struc)
  lk=dfcopula_f2(udat,ydat,countdat,d,dcopf1=cops$dcop$f1,dcopf2=cops$dcop$f2,
                 pcondcopf1=cops$pcondcop$f1,pcondcopf2=cops$pcondcop$f2,
                 cutp,estcount,th=cpar,gln,glw,param)
  if(any(lk <= 0) || any(is.nan(lk)) || any(is.infinite(lk))){
    return(1e+10)
  }
  loglik= log(lk)
  return(- sum(loglik))
} 
