#log-likelihood function for the Nested factor copula models

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

nestedcop.ir= function(ydat,cutp,dcop,pcondcop,qcondcop,th.x0g,th.xgy,nq,ngrp,grpsize,glw,gln,param)
{
  K=nrow(cutp)+1
  d=ncol(cutp)
  n=nrow(ydat)

  mgrid=meshgrid(gln,gln,nargout=2)

  #creating G dependant quadrature points, where G is the number of groups.
  gln.grp.dep=list()
  lengthgln=length(mgrid$y)
  gln.dep.mgrid=rep(NA,lengthgln)
  for(g in 1:ngrp){
    qcond.jg=qcondcop[[g]]
    for(i in 1:lengthgln){
      gln.dep.mgrid[i]=qcond.jg(mgrid$y[i],mgrid$x[i],th.x0g[g])
    }
    gln.dep.mgrid=matrix(gln.dep.mgrid,ncol = nq)
    gln.grp.dep[[g]]=gln.dep.mgrid
  }


  condcdf<-array(NA,dim=c(nq,nq,K-1,ngrp,d))
  ind=0
  for(g in 1:ngrp){
  ind1=ind+1;ind2=ind+grpsize[g];ind=ind+grpsize[g];
  for(j in ind1:ind2){
  pcondxgy=pcondcop[[j]]
    for(k in 1:(K-1))
    {
      condcdf[,,k,g,j]<-pcondxgy(cutp[k,j],c(gln.grp.dep[[g]]),th.xgy[j],param)
    }
  }
  }

  arr0<-array(0,dim=c(nq,nq,1,ngrp,d))
  arr1<-array(1,dim=c(nq,nq,1,ngrp,d))
  condcdf<-abind(arr0,condcdf,arr1,along=3)
  pmf<-array(NA,dim=c(nq,nq,K,ngrp,d))
  for(j in 1:d){
    for(k in 1:K){
      pmf[,,k,,j]<-condcdf[,,k+1,,j]-condcdf[,,k,,j]
    }
  }

  #dfactor.x0g<-array(NA,dim=c(nq,nq,ngrp))
  #for(g in 1:ngrp){
  #  dcopx0xg=dcop[[g]]
  #  for(q in 1:nq){
  #    dfactor.x0g[q,,g]=dcopx0xg( c(mgrid$x),
  #                                c(gln.grp.dep[[g]]),th.x0g[g],param)
  #  }
  #}
  #dfactor.x0g<-array(1,dim=c(nq,nq,ngrp))

  fproduct=matrix(NA,n)
  for(i in 1:n)
  {
    ind=0;gprod = 1
    for(g in 1:ngrp){
      ind1=ind+1;ind2=ind+grpsize[g];ind=ind+grpsize[g];
      jprod=1;
      for(j in ind1:ind2){
        jprod=jprod*pmf[,,ydat[i,j]+1,g,j]
      }

      dfactor.x0gy=jprod#*dfactor.x0g[,,g]

      temp1=as.vector(glw %*% dfactor.x0gy )
      gprod = gprod*temp1
    }
    fproduct[i, ]= glw %*%gprod
  }
  return(fproduct)
}

nestedllk=function(tau,ydat,cutp,copX0,copXgY,nq,ngrp,grpsize,glw,gln,SpC=NULL,param=FALSE){
  copnm=c(copX0,copXgY)
  if (!is.null(SpC)){
    tau[ngrp+SpC] = 0.95
  }
  boundlimits=LUbound(copnm)#need source("Low-Upp-Cop-Boundaries.R")
  if(sum(mapply(function(x,y)sum(x<y[1]|x>y[2]),x=as.list(tau),y=boundlimits))>0){return(1.e10)}

  cpar=mapply(function(x,y) tau2par(x,y),x=copnm,y=tau)
  thlength=length(tau)
  th.x0g=cpar[1:ngrp]
  th.xgy=cpar[(ngrp+1):thlength]

  cops=copulas_nested(copX0,copXgY)
  lk = nestedcop.ir(ydat=ydat,cutp=cutp,dcop=cops$dcopX0Xg,
                    pcondcop=cops$pcondcopXgY,qcondcop=cops$qcopX0Xg,
                    th.x0g=th.x0g,th.xgy=th.xgy,nq=nq,ngrp=ngrp,grpsize=grpsize,
                    glw,gln,param)
  if (any(lk <= 0) || any(is.nan(lk)) || any(is.infinite(lk)))  {return(1e+10)}
  llk=sum(log(lk))
  -llk
}
