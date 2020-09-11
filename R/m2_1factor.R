# M2 statistics: provides values for (1) M2 statistic
# (2) degrees of freedom (3) p-value.
# This code can be applied on distinct number of categories.

#Input:::::::::
# dat: ordinal data (includes discretized continuous).
# pnorma: cutpoints for ordinal data.
# theta: copula parameters.
# pcondcop: condtional copula cdfs in list form.
# pconddotcop: derivatives of pcondcop w.r.t theta.
# dcop: copula densities.
# gln: Guass legendre points.
# param: T or F for (param)eterization. Default is FALSE.
# outputs:::::: list of values of the following
# fden,pmf,dnorma,fdotden,prob1 --- See below.

store.1fact<-function(dat,pnorma,theta,pcondcop,pconddotcop,dcop,gln,param=F){ 
  d=ncol(pnorma)
  K=nrow(pnorma)+1
  nq<-length(gln)
  dnorma=dnorm(qnorm(pnorma))
  fden<-array(NA,dim=c(nq,K,d))
  condcdf<-array(NA,dim=c(nq,K-1,d))
  # length of theta in a list form
  ddot=length(pconddotcop)
  fdotden<-array(NA,dim=c(nq,K,ddot))
  conddotcdf<-array(NA,dim=c(nq,K-1,ddot))
  cden<-array(NA,dim=c(nq,K-1,d))
  for(j in 1:d){ 
    pcondcopula=pcondcop[[j]]
    dcopula=dcop[[j]]
    for(k in 1:(K-1)){ 
      condcdf[,k,j]<-pcondcopula(pnorma[k,j],gln,theta[[j]],param)
      cden[,k,j]<-dcopula(pnorma[k,j],gln,theta[[j]],param)
    }
  }
  #modify for two-parameter families.
  dotj = lengths(theta)
  dotpnorma=NULL
  for(i in 1:length(dotj)){
    for(j in 1:dotj[i]){
      dotpnorma=cbind(dotpnorma,pnorma[,i])
    }
  }
  dottheta=rep(theta,times= dotj)
  for(j in 1:ddot){ 
    pconddotcopula=pconddotcop[[j]]
    for(k in 1:(K-1)){ 
      conddotcdf[,k,j]<-pconddotcopula(dotpnorma[k,j],gln,dottheta[[j]])
    }
  }
  arr0<-array(0,dim=c(nq,1,d))
  arr1<-array(1,dim=c(nq,1,d))
  condcdf<-abind(arr0,condcdf,arr1,along=2)
  for(j in 1:d){
    for(k in 1:K){ 
      fden[,k,j]<-condcdf[,k+1,j]-condcdf[,k,j]
    }
  }
  arr0dot<-array(0,dim=c(nq,1,ddot))
  conddotcdf<-abind(arr0dot,conddotcdf,arr0dot,along=2)
  for(j in 1:ddot){
    for(k in 1:K){ 
      fdotden[,k,j]<-conddotcdf[,k+1,j]-conddotcdf[,k,j]
    }
  }
  allpnorma=rbind(0,pnorma,1)
  prob1<-matrix(NA,K,d)
  for(j in 1:d){
    for(k in 1:K){
      prob1[k,j]<-allpnorma[k+1,j]-allpnorma[k,j]
    }
  }
  list(fden=fden,dnorma=dnorma,cden=cden,fdotden=fdotden,prob1=prob1)
}

# all the estimated probabilities univariate and bivariate
# to derive \dot\pi_2 apply this function to all bivariate margins and stuck
# the results
pihat.1fact<-function(prob1,fden,Kj,glw)  #theta a vector with 2 parameters
{ 
  dims<-dim(fden)
  d<-dims[3]
  res2<-prob1[-1,]
  if( sum(res2 == 0) > 0 ){
    res2[res2==0]=NA
    res2=as.numeric(na.omit(c(res2)))
  }else{
    res2=res2
  }
  for(i1 in 1:(d-1)){
    for(i2 in (i1+1):d){
      Ki1= Kj[i1]
      Ki2= Kj[i2]
      res1<-NULL
      for(j1 in 2:Ki1){
        for(j2 in 2:Ki2){
          res1<-c(res1,(fden[,j1,i1]*fden[,j2,i2])%*%glw)
        }
      }
      res2<-c(res2,res1)
    }
  }
  res2
}


# the bivariate pairs
# d the dimension (numeric)
bivpairs.1fact<-function(d){
  res<-NULL
  for(id1 in 1:(d-1)){
    for(id2 in (id1+1):d){
      res<-rbind(res,c(id1,id2))
    }
  }
  res
}
# all the observed probabilities univariate and bivariate
# dat the response data
# alpha a matrix with the cutpoints for all the responses
# bpairs the bivariate pairs
piobs.1fact<-function(dat,bpairs,Kj){ 
  n<-nrow(dat)
  d<-ncol(dat)
  Kj1=Kj-1
  res<-NULL
  for(j in 1:d ){
    nKj1= Kj1[j]  
      for(jj in 1:nKj1){
        res<-c(res,sum(dat[,j] == jj )/n) 
      }       
  }
  for(i in 1:nrow(bpairs)){
    y1=bpairs[i,1]
    y2=bpairs[i,2]
    Kjy1= Kj1[y1]
    Kjy2= Kj1[y2]
    for(j1 in 1:Kjy1){
      for(j2 in 1:Kjy2){ 
        res<-c(res,sum(dat[,y1]== j1 & dat[,y2]== j2 )/n) 
      }
    }
  }
  res
}
################################################################################
#                              functions for                                   #
#                                  \Xi_2                                       #
################################################################################
# i is the indexing and j is the value
# cov(Pr(y_{i1}=j1),Pr(y_{i2}=j2))
# for binary there is a problem with the confusion in the cutpoints vector
# or matrix to be resolved later
cov2.1fact<-function(i1,j1,i2,j2,prob1,fden,glw)
{ if(i1==i2 & j1!=j2) { -prob1[j1,i1]*prob1[j2,i2]
} else {
  if(i1==i2 & j1==j2) { temp<-prob1[j1,i1]
  temp*(1-temp)
  } else {
    prob2<-(fden[,j1,i1]*fden[,j2,i2])%*%glw
    prob2-prob1[j1,i1]*prob1[j2,i2]}}
}

# i is the indexing and j is the value
# cov(Pr(y_{i1}=j1),Pr(y_{i2}=j2,y_{i3}=j3)
cov3.1fact<-function(i1,j1,i2,j2,i3,j3,prob1,fden,glw)
{ if((i1==i2 & j1!=j2) | (i1==i3 & j1!=j3) )
{ prob2<-(fden[,j2,i2]*fden[,j3,i3])%*%glw
-prob1[j1,i1]*prob2
} else {
  if((i1==i2) | (i1==i3))
  { prob2<-(fden[,j2,i2]*fden[,j3,i3])%*%glw
  prob2*(1-prob1[j1,i1])
  } else {
    temp<-fden[,j2,i2]*fden[,j3,i3]
    prob3<-(fden[,j1,i1]*temp)%*%glw
    prob2<-temp%*%glw
    prob1<-prob1[j1,i1]
    prob3-prob1*prob2
  }}
}

ndistinct.1fact=function(i1,i2,i3,i4)
{ tem=unique(c(i1,i2,i3,i4))
length(tem)
}

#i is the indexing and j is the value
#cov(Pr(y_{i1}=j1,y_{i2}=j2),Pr(y_{i3}=j3,y_{i4}=j4)
cov4.1fact<-function(i1,j1,i2,j2,i3,j3,i4,j4,fden,glw)
{ nd=ndistinct.1fact(i1,i2,i3,i4)
if(nd==4)
{ temp1<-fden[,j1,i1]*fden[,j2,i2]
temp2<-fden[,j3,i3]*fden[,j4,i4]
prob4<-(temp1*temp2)%*%glw
prob21<-temp1%*%glw
prob22<-temp2%*%glw
prob4-prob21*prob22
} else { if(nd==3)
{ if((i1==i3 & j1==j3) | (i2==i3 & j2==j3))
{ prob3<-(fden[,j1,i1]*fden[,j2,i2]*fden[,j4,i4])%*%glw
prob21<-(fden[,j1,i1]*fden[,j2,i2])%*%glw
prob22<-(fden[,j3,i3]*fden[,j4,i4])%*%glw
prob3-prob21*prob22
} else {
  if((i1==i4 & j1==j4) | (i2==i4 & j2==j4))
  { temp<-fden[,j1,i1]*fden[,j2,i2]
  prob3<-(temp*fden[,j3,i3])%*%glw
  prob21<-temp%*%glw
  prob22<-(fden[,j3,i3]*fden[,j4,i4])%*%glw
  prob3-prob21*prob22
  } else {
    if((i1==i3 & j1!=j3) | (i1==i4 & j1!=j4) | (i2==i3 & j2!=j3) |(i2==i4 & j2!=j4))
    { prob21<-(fden[,j1,i1]*fden[,j2,i2])%*%glw
    prob22<-(fden[,j3,i3]*fden[,j4,i4])%*%glw
    -prob21*prob22
    }}}
} else { if(nd==2)
{ if (j1!=j3 | j2!=j4)
{ prob21<-(fden[,j1,i1]*fden[,j2,i2])%*%glw
prob22<-(fden[,j3,i3]*fden[,j4,i4])%*%glw
-prob21*prob22
} else {
  prob2<-(fden[,j1,i1]*fden[,j2,i2])%*%glw
  prob2*(1-prob2) }
}}}
}


# the \Xi_2 matrix
Xi2.1fact<-function(prob1,fden,Kj,glw)
{ 
  dims<-dim(fden)
  d=dims[3]
  d1=d-1
  rows<-NULL
  
  for(i1 in 1:d){
    Kji1= Kj[i1]
    for(j1 in 2:Kji1){
      res<-NULL
      for(i2 in 1:d){
        Kji2= Kj[i2]
        for(j2 in 2:Kji2){
          res<-c(res,cov2.1fact(i1,j1,i2,j2,prob1,fden,glw))
        }
      }
      for(i2 in 1:d1){
        Kji2= Kj[i2]
        for(i3 in (i2+1):d){
          Kji3= Kj[i3]
          for(j2 in 2:Kji2){
            for(j3 in 2:Kji3){
              res<-c(res,cov3.1fact(i1,j1,i2,j2,i3,j3,prob1,fden,glw))
            }
          }
        }
      }
      rows<-rbind(rows,res)
    }
  }
  for(i1 in 1:d1){
    for(i2 in (i1+1):d){
      Kji1= Kj[i1]
      for(j1 in 2:Kji1){
        Kji2= Kj[i2]
        for(j2 in 2:Kji2){
          res<-NULL
          for(i3 in 1:d){
            Kji3= Kj[i3]
            for(j3 in 2:Kji3){
              res<-c(res,cov3.1fact(i3,j3,i1,j1,i2,j2,prob1,fden,glw))
            }
          }
          for(i3 in 1:d1){
            for(i4 in (i3+1):d){
              Kji3= Kj[i3]
              for(j3 in 2:Kji3){
                Kji4= Kj[i4]
                for(j4 in 2:Kji4){
                  res<-c(res,cov4.1fact(i1,j1,i2,j2,i3,j3,i4,j4,fden,glw))
                }
              }
            }
          }
          rows<-rbind(rows,res)
        }
      }
    }
  }
  rows
}


################################################################################
#                              functions for \Delta_2                          #
################################################################################

# The derivative of \Pr_{i1 i2, y1 y2} wrt a_{cutp i}  (three functions)
der.bivprob.cutp.1fact<-function(i1,i2,y1,y2,i,cutp,dnorma,cden,fden,glw){ 
  if(i1==i & (y1-1)==cutp){
    (dnorma[cutp+1,i]*cden[,cutp+1,i]*fden[,y2,i2])%*%glw
    } 
  else {
    if(i2==i & (y2-1)==cutp){
      (dnorma[cutp+1,i]*cden[,cutp+1,i]*fden[,y1,i1])%*%glw
      }
    else {
      if(i1==i & cutp==(y1-2)){
        (-dnorma[cutp+1,i]*cden[,cutp+1,i]*fden[,y2,i2])%*%glw
        } 
      else {
        if(i2==i & cutp==(y2-2)){
          (-dnorma[cutp+1,i]*cden[,cutp+1,i]*fden[,y1,i1])%*%glw
        } 
        else {0}
      }
    }
  }
}

# finally i stuck all the possible cases together
all.der.bivprob.cutp.1fact<-function(dnorma,cden,fden,Kj,glw){ 
  dims<-dim(fden)
  d=dims[3]
  Kj2=Kj-2
  cols<-NULL
  for(i in 1:d){
    Kj2i= Kj2[i]
    for(cutp in 0:Kj2i){
      res<-NULL
      for(i1 in 1:(d-1)){
        for(i2 in (i1+1):d){
          Kj2i1= Kj[i1]
          for(y1 in 2:Kj2i1){
            Kj2i2= Kj[i2]
            for(y2 in 2:Kj2i2){
              res<-c(res,der.bivprob.cutp.1fact(i1,i2,y1,y2,i,cutp,dnorma,cden,fden,glw))
            }
          }
        }
      }
      cols<-cbind(cols,res)
    }
  }
  cols
}


################################################################################
# The derivative of \Pr_{i1 i2, y1 y2} wrt theta_i  (three functions)
# first I calculate the inner derivative
# i denotes the index for the copula vector
der.bivprob.theta.1fact<-function(i1,i2,y1,y2,i,cmsm.ddotj,fden,fdotden,glw){
  if(i1==i){ 
    (fden[,y2,i2]*fdotden[,y1,i+cmsm.ddotj])%*%glw
    }
  else
    {
      if(i2==i){
        (fden[,y1,i1]*fdotden[,y2,i+cmsm.ddotj])%*%glw
        }
      else
        {
         0
        }
    }
  }

# finally I stuck everything together
all.der.bivprob.theta.1fact<-function(fden,fdotden,theta,Kj,glw)
{
  dims<-dim(fden)
  d=dims[3]
  # length of parameters (number of der.cond.cop. wrt parameters)
  ddot=dim(fdotden)[3]
  ddotj=lengths(theta)-1
  cmsm.ddot=c(0,cumsum(ddotj)[-length(ddotj)])
  cols<-NULL
  for(i in 1:d){
    for(idot in 0:ddotj[i]){
      cmsm.ddotj = cmsm.ddot[i]+idot
      res<-NULL
      for(i1 in 1:(d-1)){
        for(i2 in (i1+1):d){
          Kji1= Kj[i1]
          for(y1 in 2:Kji1){
            Kji2= Kj[i2]
            for(y2 in 2:Kji2){
              res<-c(res,der.bivprob.theta.1fact(i1,i2,y1,y2,i,cmsm.ddotj,fden,fdotden,glw))
            }
          }
        }
      }
      cols<-cbind(cols,res)
    }
  }
  cols
}

################################################################################
# all the derivatives of the univariate probabilities wrt to the cutpoints
all.der.uprob.1fact<-function(dnorma,theta){
  d=ncol(dnorma)
  ddot=sum(lengths(theta))
  r=nrow(dnorma)
  neg.dnorma=-dnorma
  rows.result= d*r
  result = matrix( 0 , ncol =  d*r+ddot , nrow =  d*r)
  diag(result)=neg.dnorma
  ndnorma=dnorma[-1,]
  subdiagonal=as.numeric(rbind(ndnorma,0))[-rows.result]
  diag(result[-nrow(result),-1])=subdiagonal
  if(sum(dnorma==0)>0){
    loc.zero=as.numeric(which(as.numeric(dnorma)==0))
    result = result[ -loc.zero, -loc.zero ]
  }else{
    result=result
  }
  return(result)
}

################################################################################
# put everything together to  compose the \Delta_2 matrix
Delta2.1fact<-function(dnorma,cden,fden,fdotden,theta,Kj,glw){
  res1<-all.der.uprob.1fact(dnorma,theta)
  res2<-all.der.bivprob.cutp.1fact(dnorma,cden,fden,Kj,glw)
  res3<-all.der.bivprob.theta.1fact(fden,fdotden,theta,Kj,glw)
  rbind(res1,cbind(res2,res3))
}

################################################################################
# R code for orthogonal complement (used by one of my students)
orthcomp <- function(x,tol=1.e-12){ 
  s <- nrow(x)
  q <- ncol(x)
  if (s<=q) { return('error, matrix must have more columns than rows') }
  x.svd <- svd(x,s)
  if (sum(abs(x.svd$d)<tol)!=0)
  { return('error, matrix must full column rank') }
  return(x.svd$u[,(q+1):s])
}

# the M_2 statistic
#Input: 
# dat: ordinal data.
# dnorma,prob1,cden,fden,fdotden: obtained from sdat function above.
# theta: copula parameters.
# Kj: number of categories for each variables. Vector form.
# glw: Gauss legendre weights.

M2.1fact<-function(dat,dnorma,prob1,cden,fden,fdotden,theta,Kj,glw,iprint=F){ 
  dims<-dim(fden)
  d=dims[3]
  n=nrow(dat)
  bpairs<-bivpairs.1fact(d)
  bm<-nrow(bpairs)
  #-------------------------------------------------------
  V<-Xi2.1fact(prob1,fden,Kj,glw)
  #-------------------------------------------------------
  if(iprint){ 
    cat("\ndim(V): dimension of Xi_2 matrix\n")
    print(dim(V))
  }
  #-------------------------------------------------------
  delta<-Delta2.1fact(dnorma,cden,fden,fdotden,theta,Kj,glw)
  #-------------------------------------------------------
  if(iprint)
  { cat("\ndim(delta): dimension of Delta_2 matrix\n")
    print(dim(delta))
  }
  #-------------------------------------------------------
  oc.delta<-orthcomp(delta)
  inv<-solve(t(oc.delta)%*%V%*%oc.delta)
  C2<-oc.delta%*%inv%*%t(oc.delta)
  pi.r<-pihat.1fact(prob1,fden,Kj,glw)
  p.r<-piobs.1fact(dat,bpairs,Kj)
  #-------------------------------------------------------
  stat<-nrow(dat)*t(p.r-pi.r)%*%C2%*%(p.r-pi.r)
  dof<-nrow(delta)-ncol(delta)
  pvalue<-pchisq(stat,dof, lower.tail = FALSE)
  #-------------------------------------------------------
  #discrepancies
  Kj1=Kj-1
  K1d=1:(sum(Kj1))
  bm<-nrow(bpairs)
  K=Kj1[bpairs[,1]]*Kj1[bpairs[,2]]
  val=abs(p.r[-K1d]-pi.r[-K1d])
  cmsm =cumsum(K)
  lst=c()
  for(i in 1:bm){
    lst[[i]]= val[(cmsm[i]-K[i]+1):cmsm[i]]
  }
  dg=round(sapply(lst,max)*n)
  #-------------------------------------------------------
  #Return output:
  list("M_2"=stat,"df"=dof,"p-value"=pvalue,"dg"=dg)
  #-------------------------------------------------------
}

