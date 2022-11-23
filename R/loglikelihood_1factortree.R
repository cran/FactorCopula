factor1vine_lk=function(ydat,A, pcondcop,vinepcop,cutp,th,nq,gln,glw,param)
{
  tree_arcs = t(rbind(A[1,-1], diag(A[-1,-1])))
  n<-nrow(ydat)#number of observations
  K<-nrow(cutp)+1#number of categories
  d<-ncol(cutp)#number of items
  d1=d-1
  th_factor<-th[(1):(d)]#copula pars for factor model; 1st tree
  th_vine<-th[(1+d):(d+d-1)]#copula pars for 2nd tree
  #--------------------------------------------------
  #----- pmf for ordinal-----------------------------
  condcdf<-array(NA,dim=c(nq,K-1,d))
  for(j in 1:d){
    pcondcopj<-pcondcop[[j]]
    for(k in 1:(K-1)){
      condcdf[,k,j]=pcondcopj(cutp[k,j],gln,th_factor[j]) 
    }
  }
  #empty arrays for the ordinal matrices-----------------
  arr0<-array(0,dim=c(nq,1,d))
  arr1<-array(1,dim=c(nq,1,d))
  condcdf<-abind(arr0,condcdf,arr1,along=2)
  #----- pmf        -------------------------------------
  pmf_factor<-array(NA,dim=c(nq,K,d))
  for(j in 1:d){
    for(k in 1:K){
      pmf_factor[,k,j]=condcdf[,k+1,j]-condcdf[,k,j] 
    }
  }
  
  
  
  #---------------------------------------------------------        
  #fproduct = matrix(NA,n,nq)
  #fdiff<-array(NA,dim=c(nq,lngth_combnK,d1))
  
  prod_num=prod_denom=1
  for(j in 1:d1)
  {
    j_index=tree_arcs[j,]
    j1= j_index[1]
    j2= j_index[2]
    
    ydatij1=ydat[,j1]
    ydatij2=ydat[,j2]
    
    Fyj_1p=condcdf[,ydatij1+2,j1]
    Fyj_1m=condcdf[,ydatij1+1,j1]
    Fyj_2p=condcdf[,ydatij2+2,j2]
    Fyj_2m=condcdf[,ydatij2+1,j2]
    
    vinepcopj=vinepcop[[j]]
    th_vinej=th_vine[j]
    ypp=vinepcopj( Fyj_1p, Fyj_2p , th_vinej)
    ypm=vinepcopj( Fyj_1p, Fyj_2m , th_vinej)
    ymp=vinepcopj( Fyj_1m, Fyj_2p , th_vinej)
    ymm=vinepcopj( Fyj_1m, Fyj_2m , th_vinej)
    
    fdiff = ypp - ypm - ymp + ymm
    prod_num<-prod_num*fdiff
  }
  
  
  #j_A=table(c(tree_arcs))-1
  #j_A=c(0,rep( 1:d  ,tt))#0 to start from 2nd val
  j_A=A[1,-1]
  for(j in 2:d1){
    j_ind=j_A[j]
    res=pmf_factor[,ydat[,j_ind]+1,j_ind]
    prod_denom = prod_denom * res
  }
  
  fproduct = prod_num/prod_denom
  t(fproduct) %*% glw
}

# purpose: calculates log likelihood for the copula for mixed data
# and therefore is used to estimated parametrs
# input:
# param: parameters(first params for cont. data, and then for the discrete dat)
# and same as above
loglk_f1vine= function(tau,ydat,A,cutp,nq,gln,glw,copnm,param=F)
{
  d=ncol(ydat)
  boundlimits=LUbound(copnm)
  if(sum(mapply(function(x,y) sum(x<y[1]|x>y[2]),x=as.list(tau),y=boundlimits))>0){return(1.e10)}
  theta=mapply(function(x,y) tau2par(x,y),x=copnm,y=tau)
  d=ncol(ydat)
  cops=copulas_factortree(copnm[1:d],0,copnm[(d+1):(d+d-1)])
  lk=factor1vine_lk(ydat,A,pcondcop=cops$pcondF1,vinepcop=cops$copvine,cutp,th=theta,nq=nq,gln,glw,param=F)
  if (any(lk <= 0) || any(is.nan(lk))) return(1.e10)
  loglik=sum(log(lk))
  return(-loglik)
}

#mle=mleF1V(y, A, cop=c(rep("gum",d+d-1)), gl=gauss.quad.prob(15), 
#       hessian = F, print.level = 2)

lk_f1vine=function(tau,ydat,A,cutp,nq,gln,glw,copnm,param=F)
{
  d=ncol(ydat)
  boundlimits=LUbound(copnm)
  if(sum(mapply(function(x,y) sum(x<y[1]|x>y[2]),x=as.list(tau),y=boundlimits))>0){return(1.e10)}
  theta=mapply(function(x,y) tau2par(x,y),x=copnm,y=tau)
  d=ncol(ydat)
  cops=copulas_factortree(copnm[1:d],0,copnm[(d+1):(d+d-1)])
  lk=factor1vine_lk(ydat,A,pcondcop=cops$pcondF1,vinepcop=cops$copvine,cutp,th=theta,nq=nq,gln,glw,param=F)
  lk
}
