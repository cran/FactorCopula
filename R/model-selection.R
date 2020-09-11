
# Model selection for the one-factor and two-factor copula model
# Input:
# 1) continuous: continuous data.
# 2) ordinal: ordinal data.
# 3) count: count data.
# 4) copnames: possible copula candidates.
# 5) gl: Gauss legendre quadrature points and weights.

# Output
# a vector of best copulas via AIC.

select1F=function(continuous=NULL, ordinal=NULL, count=NULL, copnamesF1, gl){
  
  #Gauss legendre points and weights
  gln<-gl$nodes
  glw<-gl$weights
  #Transforming the variables to U(0,1)
  #--------------------------------------------------------
  #                 Continuous variables
  #-----------------                    -------------------
  if(!is.null(continuous)){
    n<-nrow(continuous)
    u<-apply(continuous,2,rank)/(n+1)
    d1<-ncol(u)
  }else{
    d1<-0
  }
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  cutpoints<-NULL
  if(!is.null(ordinal)){
    n<-nrow(ordinal)
    cutpoints<-cuts(ordinal)
    d2<-ncol(cutpoints)
  }else{
    d2<-0
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
  }else{
    d3<-0
  }
  #number of variables
  d<-d1+d2+d3
  # The possible copula candidates for all variables
  #allcop<-list(copnamesF1)
  # creating a list of all copula candidates for the
  # 1st  factor
  coplstF1 <- copnamesF1 #rep(allcop,d)
  #calculating the length for the possible copulas for all variables
  coplngthF1<-as.numeric(sapply(coplstF1,length))
  #starting copulas
  copF1vec<-rep("frk",d)
  #========================================================
  for(j in 1:d){
    print(c("jF1",j))
    copd=coplngthF1[j]
    ii<-NA
    inner1<-matrix(NA, nrow = copd, ncol = 2*d+1)
      for( ii in 1:copd){
      print(c("iF1",ii))
      copF1vec[j]<-coplstF1[[j]][[ii]]
      #--------------------------------------------------------
      #--------------- structuring the parameters -------------
      theta_struc <- initPar(copF1vec)
      init_val <- unlist(as.relistable(theta_struc))
      modF1=list(estimate=rep(NA,d*2),minimum=NA)
      try(modF1 <- nlm(loglk_f1, init_val, udat=u, ydat=ordinal, 
                       countdat=count, copF1=copF1vec,
                       cutp=cutpoints, estcount=est_count,
                       gln, glw, theta_struc, param=T, 
                       hessian=F,print.level = 0),TRUE)
      #------- parameters
      deltas=list(relist(c(modF1$e), skeleton = theta_struc))
      #Converting back to copula parameters
      cpar=delta2cpar(deltas, copF1vec)
      # convert the copula paramters to taus
      taus = rep(NA, d)
      try(taus <- mapply(function(x,y) par2tau(x,y), 
                  x = copF1vec, y = cpar$f1 ), TRUE)
      
      #AIC
      aic.F1 = 2*modF1$minimum+2*length(modF1$e)
      inner1[ii, ] = c(aic.F1, copF1vec, taus)
    }
    index1<- which.min(cbind(as.numeric(noquote(rbind(inner1))[,1])))
    model.aic<- inner1[index1,1]
    copF1vec[j]<- inner1[index1, j+1]
    d.inner1=ncol(inner1)
    estimated.taus<- as.numeric(inner1[index1, (d+2):d.inner1])
  }
  
  return(list("1st factor" = copF1vec,
              "AIC" = model.aic, "estimated taus" = estimated.taus))
}


select2F=function(continuous=NULL, ordinal=NULL, count=NULL, copnamesF1, copnamesF2, gl){
#Gauss legendre points and weights
  gln<-gl$nodes
  glw<-gl$weights
#Transforming the variables to U(0,1)
  #--------------------------------------------------------
  #                 Continuous variables
  #-----------------                    -------------------
  if(!is.null(continuous)){
    n<-nrow(continuous)
    u<-apply(continuous,2,rank)/(n+1)
    d1<-ncol(u)
  }else{
    d1<-0
  }
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  cutpoints<-NULL
  if(!is.null(ordinal)){
    n<-nrow(ordinal)
    cutpoints<-cuts(ordinal)
    d2<-ncol(cutpoints)
  }else{
    d2<-0
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
  }else{
    d3<-0
  }
  #number of variables
  d<-d1+d2+d3
# The possible copula candidates for all variables
#allcop<-list(copnames)

# creating a list of all copula candidates for the
# 1st and 2nd factors
coplstF1 <- copnamesF1 #rep((allcop),d)
coplstF2 <- copnamesF2 #rep((allcop),d)

#calculating the length for the possible copulas for all variables
coplngthF1<-as.numeric(sapply(coplstF1,length))
coplngthF2<-as.numeric(sapply(coplstF2,length))

#starting copulas
copF1vec<-rep("frk",d)
copF2vec<-rep("frk",d)
#========================================================
for(j in 1:d){
  print(c("jF1",j))
  copd=coplngthF1[j]
  ii<-NA
  
  inner1 = matrix( NA, nrow = copd, ncol =  d+1)
  for( ii in 1:copd){
    print(c("iF1",ii))
    copF1vec[j]<-coplstF1[[j]][[ii]]
    #--------------------------------------------------------
    #--------------- structuring the parameters -------------
    theta_struc <- initPar(copF1vec, copF2vec)
    init_val <- unlist(as.relistable(theta_struc))
    #---------------------- Estimation ----------------------
    modF1=list(estimate=rep(NA,d*2),minimum=NA)
    try(modF1 <- nlm(loglk_f2, init_val, udat=u, ydat=ordinal, 
                     countdat=count,copF1=copF1vec, copF2=copF2vec, 
                     cutp=cutpoints, estcount=est_count,
                     gln, glw, theta_struc, param=T, 
                     hessian=F,print.level = 0),TRUE)
    
    #AIC
    aic.F1= 2*modF1$minimum+2*length(modF1$e)
    inner1[ii, ] = c(aic.F1,copF1vec)
  }
  index1=which.min(cbind(as.numeric(noquote(rbind(inner1))[,1])))
  copF1vec[j]<- inner1[index1,j+1]
}

for(j in 1:d){
  print(c("jF2",j))
  copd=coplngthF2[j]
  #empty arrays
  ii<-NA
  inner2 = matrix(NA, nrow = copd, ncol = 3*d+1)
  for( ii in 1:copd){
    print(c("iF2",ii))
    copF2vec[j]<-coplstF2[[j]][[ii]]
    #--------------------------------------------------------
    #--------------- structuring the parameters -------------
    theta_struc <- initPar(copF1vec, copF2vec)
    init_val <- unlist(as.relistable(theta_struc))
    modF2=list(estimate=rep(NA,d*2),minimum=NA)
    try(modF2 <- nlm(loglk_f2, init_val, udat=u, ydat=ordinal, 
                     countdat=count, copF1=copF1vec, copF2=copF2vec, 
                     cutp=cutpoints, estcount=est_count,
                     gln, glw, theta_struc, param=T, 
                     hessian=F,print.level = 0),TRUE)
    
    #------- parameters
    deltas=relist(c(modF2$e), skeleton = theta_struc)
    #Converting back to copula parameters
    cpar=delta2cpar(deltas, copF1vec, copF2vec)
    # convert the copula paramters to taus
    taus1 = rep(NA, d)
    try(taus1 <- mapply(function(x,y) par2tau(x,y), 
                x = copF1vec, y = cpar$f1), TRUE)
    taus2 = rep(NA, d)
    try(taus2 <- mapply(function(x,y) par2tau(x,y), 
                x = copF2vec, y = cpar$f2), TRUE)
    taus = c(taus1, taus2)
    
    #AIC
    aic.F2= 2*modF2$minimum+2*length(modF2$e)
    inner2[ii, ] = c(aic.F2, copF2vec, taus)
  }
  index2<- which.min(cbind(as.numeric(noquote(rbind(inner2))[,1])))
  model.aic<- inner2[index2,1]
  copF2vec[j]<- inner2[index2,j+1]
  d.inner2=ncol(inner2)
  estimated.taus<- as.numeric(inner2[index2,(2+d):d.inner2])
}

return(list("1st factor" = copF1vec, "2nd factor" = copF2vec,
             "AIC" = model.aic, "estimated taus" = estimated.taus))
}
