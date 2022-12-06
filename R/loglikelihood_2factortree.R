factor2vine_lk=function(ydat,A,pcondcopf1,pcondcopf2,vinepcop,cutp,th,nq,gln,glw,param)
{
  tree_arcs = t(rbind(A[1,-1], diag(A[-1,-1])))
  n<-nrow(ydat)#number of observations
  K<-nrow(cutp)+1#number of categories
  d<-ncol(cutp)#number of items
  d1=d-1
  th_factor1<-th[(1):(d)]#copula pars for factor model; 1st tree
  th_factor2<-th[(1+d):(d*2)]#copula pars for factor model; 2nd tree
  th_vine<-th[(1+d*2):(2*d+d-1)]#copula pars for 2nd tree
  #------------------------------------------------------
  #----- pmf for ordinal-----------------------------
  # empty arrays to save the conditional one-factor copula in condcdf  
  condcdf<-array(NA,dim=c(nq,K-1,d))
  condcdf2<-array(NA,dim=c(nq,nq,K-1,d))
  for(j in 1:d){
    pcondcopf1j=pcondcopf1[[j]]
    pcondcopf2j=pcondcopf2[[j]]
    for(k in 1:(K-1)){
      condcdf[,k,j]<-pcondcopf1j(cutp[k,j],gln,th_factor1[j]) 
      for(q in 1:nq)
      { 
        condcdf2[,q,k,j]<-pcondcopf2j(condcdf[q,k,j],gln,th_factor2[j])
      }
    }
  }
  # binding the array of condcdf2 between 0 & 1.
  arr0<-array(0,dim=c(nq,nq,1,d))
  arr1<-array(1,dim=c(nq,nq,1,d))
  condcdf2<-abind(arr0,condcdf2,arr1,along=3)
  # calculate the finite difference for the discrete data. 
  pmf_factor<-array(NA,dim=c(nq,nq,K,d))
  for(j in 1:d){
    for(k in 1:K)
    { 
      pmf_factor[,,k,j]<-condcdf2[,,k+1,j]-condcdf2[,,k,j] 
    }
  }

  #---------------------------------------------------------        
  fproduct = matrix(NA,n,nq)
  for(i in 1:n)
  {
    prod_num=prod_denom=1
    for(j in 1:d1)
    {
      j_index=tree_arcs[j,]
      j1= j_index[1]
      j2= j_index[2]
      
      Fyj_1p=condcdf2[,,ydat[i,j1]+2,j1]
      Fyj_1m=condcdf2[,,ydat[i,j1]+1,j1]
      
      Fyj_2p=condcdf2[,,ydat[i,j2]+2,j2]
      Fyj_2m=condcdf2[,,ydat[i,j2]+1,j2]
      
      vinepcopj=vinepcop[[j]]
      ypp=vinepcopj( Fyj_1p, Fyj_2p , th_vine[j])
      ypm=vinepcopj( Fyj_1p, Fyj_2m , th_vine[j])
      ymp=vinepcopj( Fyj_1m, Fyj_2p , th_vine[j])
      ymm=vinepcopj( Fyj_1m, Fyj_2m , th_vine[j])
      
      fdiff = ypp - ypm - ymp + ymm
      prod_num<-prod_num*fdiff
    }
    
    #j_A=table(c(tree_arcs))-1
    #j_A=c(0,rep( 1:d  ,tt))#0 to start from 2nd val
    j_A=A[1,-1]
    for(j in 2:d1){
      j_index=j_A[j]
      res = pmf_factor[,,ydat[i,j_index]+1,j_index]
      prod_denom = prod_denom * res
    }
    
    fproduct[i, ] = as.vector((prod_num/prod_denom)%*% glw)
  }
  fproduct %*% glw
}


# purpose: log likelihood of factor tree for item response
loglk_f2vine= function(tau,ydat,A,cutp,nq,gln,glw,copnm,param=F)
{
  d=ncol(ydat)
  boundlimits=LUbound(copnm)
  if(sum(mapply(function(x,y) sum(x<y[1]|x>y[2]),x=as.list(tau),
                y=boundlimits))>0){return(1.e10)}
  theta=mapply(function(x,y) tau2par(x,y),x=copnm,y=tau)
  ntheta=theta
  #if(!is.null(SpC)){ #ARIS CHANGE HERE 
  #  ntheta[ (d + SpC) ]=0}
  
  cops=copulas_factortree(copnm[1:d],copnm[(d+1):(2*d)],copnm[(2*d+1):(2*d+d-1)])
  lk=factor2vine_lk(ydat,A,cops$pcondF1,cops$pcondF2,cops$copvine,cutp,
                    th=ntheta,nq=nq,gln,glw,param=F)
  if (any(lk <= 0) || any(is.nan(lk))) return(1.e10)
  loglik=sum(log(lk))
  return(-loglik)
} 


# purpose: likelihood of factor tree for item response
lk_f2vine= function(tau,ydat,A,cutp,nq,gln,glw,copnm,param=F)
{
  d=ncol(ydat)
  boundlimits=LUbound(copnm)
  if(sum(mapply(function(x,y) sum(x<y[1]|x>y[2]),x=as.list(tau),
                y=boundlimits))>0){return(1.e10)}
  theta=mapply(function(x,y) tau2par(x,y),x=copnm,y=tau)
  cops=copulas_factortree(copnm[1:d],copnm[(d+1):(2*d)],copnm[(2*d+1):(2*d+d-1)])
  lk=factor2vine_lk(ydat,A,cops$pcondF1,cops$pcondF2,cops$copvine,cutp,
                    th=theta,nq=nq,gln,glw,param=F)
  return(lk)
} 

