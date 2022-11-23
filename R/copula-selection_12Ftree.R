# Models selction algorithm for common 1-2-factor copula models
#library(FactorCopula)
#library(foreach)
#library(doMC)
##registerDoMC(3)

#copula selection for 1 factor, one copula family for
#all items with latent factor
copselect1Factor= function(y, allcop, gl, d){
  #====== The possible copulas for all variables
  length.allcop = length(allcop)
  copMatF1=matrix(allcop,nrow = length.allcop,ncol = d)
  coplngthF1=length.allcop
  
  j<-NA
  inner1<-matrix(NA,coplngthF1,d+1)
  for(j in 1:coplngthF1) #, .combine = rbind )%dopar%
  {
    print(c("F1",j))
    copF1<-copMatF1[j,]
    modF1=list(estimate=rep(NA,d),minimum=NA)
    try(modF1<-mle1factor(continuous=NULL,ordinal=y,
                          count=NULL, copF1, gl, hessian = F, print.level = 0),TRUE)
    inner1[j,]=c(modF1$loglik,copF1)
  }
  index1<-which.max(cbind(as.numeric(noquote(rbind(inner1))[,1])))
  copF1<-inner1[index1,-1]
  
  return(list("1st factor" = copF1))
}

#copula selection for 2 factor, one copula family for
#all items with first latent factor, and one
#copula family for all items with second latent factor,
copselect2Factor= function(y, allcop, gl, d){
  #====== The possible copulas for all variables
  length.allcop = length(allcop)
  copMatF1=copMatF2=matrix(allcop,nrow = length.allcop,ncol = d)
  coplngthF1=coplngthF2=length.allcop
  j<-NA
  inner1<-matrix(NA,coplngthF1,d+1)
  for(j in 1:coplngthF1) #, .combine = rbind )%dopar%
  {
    print(c("F1",j))
    copF1<-copMatF1[j,]
    modF1=list(estimate=rep(NA,d),minimum=NA)
    try(modF1<-mle2factor(NULL,y,NULL, copF1 =copF1, copF2 = rep("bvn",d), 
                          gl, hessian = F, print.level = 0),TRUE)
    inner1[j,]=c(modF1$loglik,copF1)
  }
  index1<-which.max(cbind(as.numeric(noquote(rbind(inner1))[,1])))
  copF1<-inner1[index1,-1]
  j<-NA
  inner2<-matrix(NA,coplngthF1,d+1)
  for(j in 1:coplngthF2)#, .combine = rbind )%dopar%
  {
    print(c("F2",j))
    copF2<-copMatF2[j,]
    modF2=list(estimate=rep(NA,d),minimum=NA)
    try(modF2<-mle2factor(NULL,y,NULL, copF1 =copF1, copF2 = copF2, gl, hessian = F, print.level = 0),TRUE)
    inner2[j,]=c(modF2$loglik,copF2)
  }
  index2<-which.max(cbind(as.numeric(noquote(rbind(inner2))[,1])))
  copF2<-inner2[index2,-1]
  
  return(list("1st factor" = copF1, "2nd factor" = copF2))
}


copselect1Ftree= function(y, A, f1copnames, vinecopnames, gl, d){

  coplngth = length(vinecopnames)
  copMat=matrix(vinecopnames,nrow = coplngth,ncol = d-1)
  j<-NA
  inner1<-matrix(NA,coplngth,3*d-1)
  for(j in 1:coplngth) #, .combine = rbind )%dopar%
  {
    print(c("1Ftree",j))
    cop1Ftree<-copMat[j,]
    copulas.all=c(f1copnames,cop1Ftree)
    mod1Ftree=list(estimate=rep(NA,d),minimum=NA)
    try(mod1Ftree<-mle1FactorTree(y,A, cop=copulas.all, gl, 
                              hessian = F, print.level = 0),TRUE)
      
    inner1[j,]=c(mod1Ftree$loglik,mod1Ftree$taus,cop1Ftree)
  }
  index1<-which.max(cbind(as.numeric(noquote(rbind(inner1))[,1])))
  taus1Ftree <- as.numeric(inner1[index1,(2):(2*d)])
  cop1Ftree<-inner1[index1,-(1:(2*d))]
  model.lglk<- as.numeric(inner1[index1,1])
  
  return(list("Vine tree copulas" = cop1Ftree,
              "Log-likelihood" = model.lglk, "estimated taus" = taus1Ftree))
}

copselect2Ftree=function(y, A, f1copnames,f2copnames, vinecopnames, gl, d){
#  d <- ncol(y)
  coplngth = length(vinecopnames)
  copMat=matrix(vinecopnames,nrow = coplngth,ncol = d-1)
  j<-NA
  inner1<-matrix(NA,coplngth,4*d-1)
  for(j in 1:coplngth) #, .combine = rbind )%dopar%
  {
    print(c("2Ftree",j))
    cop2Ftree<-copMat[j,]
    copulas.all=c(f1copnames,f2copnames,cop2Ftree)
    mod2Ftree=list(estimate=rep(NA,d),minimum=NA)
    try(mod2Ftree<-mle2FactorTree(y,A, cop=copulas.all, gl,SpC=NULL,
                              hessian = F, print.level = 0),TRUE)
    
    inner1[j,]=c(mod2Ftree$loglik,mod2Ftree$taus,cop2Ftree)
  }
  index1<-which.max(cbind(as.numeric(noquote(rbind(inner1))[,1])))
  #taus2Ftree <- as.numeric(inner1[index1,(2):(2*d)])
  #cop2Ftree<-inner1[index1,-(1:(2*d))]
  taus2Ftree <- as.numeric(inner1[index1, (2):(3 * d)])
  cop2Ftree <- inner1[index1, -(1:(3 * d))]
  model.lglk<- as.numeric(inner1[index1,1])
  
  return(list("Vine tree copulas" = cop2Ftree,
              "Log-likelihood" = model.lglk, "estimated taus" = taus2Ftree))
}




#--------
#y = as.matrix(read.table("ptsd.txt"))
#dim(y)
#apply(y,2,table)
#d=ncol(y)
#n=nrow(y)
#nq=15
#gl=gauss.quad.prob(nq)
#gln=gl$n
#glw=gl$w

#allcop=c("bvn","bvt2","bvt3","bvt4","bvt5","bvt6",
#         "bvt7","bvt8","bvt9","gum","rgum")

#copulas.1f=Select.1factor(ydat=y, allcop, gln, glw,
#                          nq, n, d)

#copulas.2f=Select.2factor(ydat=y, allcop, gln, glw,
#                          nq, n, d)

