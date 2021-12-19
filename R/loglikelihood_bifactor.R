cuts<-function(dat)
{ n<-nrow(dat)
d<-ncol(dat)
K<-length(table(dat[,1]))
prob<-matrix(NA,K,d)
for(j in 1:d)
{ prob[,j]<-table(dat[,j])/n }
a<-matrix(NA,K-1,d)
for(j in 1:d)
{ s<-0
for(i in 1:(K-1))
{ s<-s+prob[i,j]
a[i,j]<-s
}
}
a
}


bipmf<-function(ydat,thX0,thXg,cutp,pcondcopX0,pcondcopXg,nq,ngrp,grpsize,glw,gln,param){
  K=nrow(cutp)+1
  d=ncol(cutp)
  n=nrow(ydat)
  condcdfX0<-array(NA,dim=c(nq,K-1,d))
  for(j in 1:d){
    pcondX0j=pcondcopX0[[j]]
    for(k in 1:(K-1))
    {
      condcdfX0[,k,j]<-pcondX0j(cutp[k,j],gln,thX0[j],param)
    }
  }
  condcdfXg<-array(NA,dim=c(nq,nq,K-1,d))
  for(j in 1:d){
    pcondXgj=pcondcopXg[[j]]
    for(k in 1:(K-1))
    {
      for(u in 1:nq)
      {
        condcdfXg[u,,k,j]<-pcondXgj(condcdfX0[u,k,j],gln,thXg[j],param)
      }
    }
  }
  arr0<-array(0,dim=c(nq,nq,1,d))
  arr1<-array(1,dim=c(nq,nq,1,d))
  condcdfXg<-abind(arr0,condcdfXg,arr1,along=3)
  pmf<-array(NA,dim=c(nq,nq,K,d))
  for(j in 1:d){
    for(k in 1:K)
    {
      pmf[,,k,j]<-condcdfXg[,,k+1,j]-condcdfXg[,,k,j]
    }
  }
  pmfprod<-matrix(NA,n,nq)
  for(i in 1:n)
  {
    pmf.g=1;ind=0
    for(g in 1:ngrp){
      ind1=ind+1;ind2=ind+grpsize[g];ind=ind+grpsize[g];pmf.j<-1;
      for(j in ind1:ind2){
        temp<-pmf[,,ydat[i,j]+1,j]
        pmf.j<-pmf.j*temp
      }
      pmf.g<-pmf.g*as.vector(pmf.j%*%glw)
    }
    pmfprod[i, ]= pmf.g
  }
  return(pmfprod%*%glw)
}


bifactorllk<-function(tau,ydat,cutp,copX0,copXgY,nq,ngrp,grpsize,glw,gln,SpC=NULL,param=FALSE){
  d<-ncol(cutp)
  #adjusting the number of columns for SpC for the group-specific factor
  if (!is.null(SpC)){
    tau[d+SpC] = 0.95
  }
  copnm=c(copX0,copXgY)
  boundlimits=LUbound(copnm)
  if(sum(mapply(function(x,y) sum(x<y[1]|x>y[2]),x=as.list(tau),y=boundlimits))>0){return(1.e10)}
  
  cpar=mapply(function(x,y) tau2par(x,y),x=copnm,y=tau)
  thX0=cpar[1:d]
  thXg= cpar[(d+1):(2*d)]

  cops=copulas_bifactor(copX0,copXgY)
  lk<-bipmf(ydat,thX0,thXg ,cutp ,pcondcopX0=cops$pcondcopX0,
            pcondcopXg=cops$pcondcopXgY,nq,ngrp,grpsize,glw,gln,param)
  if(any(lk <= 0) || any(is.nan(lk)) || any(is.infinite(lk))){return(1e+10)}
  llk=sum(log(lk))
  -llk
}



