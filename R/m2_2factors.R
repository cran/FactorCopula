# This code is to calculate the M2 statistics 
# There is a modification in this code over 
# the code that Aris gave me (Sayed). These changes 
# are followed with this sign " #***** ". So,
# This code can be applied on 
# distinct number of categories of ordinal data.
################################################################################
################################################################################
#################################2-factor model M2##############################
################################################################################
################################################################################
#
#
#
# purpose:
# calulates the arrays fden2, cden2, fdotden2, fbarden2, etc for the computation
# of the M2 for the 2-factor model (see appendix of paper)
# we allow two different copula families, one for factor 1 and one for factor 2
# input:
# dat: A data matrix where the number of rows corresponds to an
# individual's response and each column represents an item
# Number of ordinal categories for each item, coded as 0,...,(ncat-1).
# Currently supported are items that have the same number of categories.
# theta: a vector with the copula parameters for the 1st factor
# delta: a vector with the copula parameters at the 2nd factor
# cutp: A matrix with the cutpoints (without the boundary cutpoints) in the
# uniform scale where the number of rows corresponds to an ordinal category
# and each column represents an item.
# pcondcop1: A function with the conditional copula cdf C_{2|1} for the 1st factor
# pconddotcop1: A function with the derivative with resect to the copula
# parameter of the conditional copula cdf C_{2|1} for the 1st factor
# dcop1: A function with the copula density c for the 1st factor
# pcondcop2: A function with the conditional copula cdf C_{2|1} for the 2nd factor
# pconddotcop2: A function with the derivative with resect to the copula
# parameter of the conditional copula cdf C_{2|1} for the 2nd factor
# dcop2: A function with the copula density c for the 2nd factor
# output:
# fden2: fden2 (see Appendix)
# dnorma: \phi[\Phi^{-1}(cutp)]
# cden2: cdens12 (see Appendix)
# fdotden2: fdotden2 (see Appendix)
# fbarden2: fbarden (see Appendix)
# prob1: the (K\times d) matrix with the univariate probabilities
store.2fact<-function(dat,pnorma,theta,delta,
                pcondcop1,pcondcop2,pconddotcop1,pconddotcop2,dcop1,dcop2,gln,param=F)
{ 
  nq<-length(gln)
  d=ncol(pnorma)
  K=nrow(pnorma)+1
  dnorma=dnorm(qnorm(pnorma))
  condcdf<-array(NA,dim=c(nq,K-1,d))
  cden<-array(NA,dim=c(nq,K-1,d))
  
  # First Factor ::::::::::::::::::::::::::::::::::::::::: First Factor
  for(j in 1:d)
  { 
    pcondcopula1=pcondcop1[[j]]
    dcopula1=dcop1[[j]]  
    for(k in 1:(K-1))
    { 
      condcdf[,k,j]<-pcondcopula1(pnorma[k,j],gln,theta[[j]],param)
      cden[,k,j]<-dcopula1(pnorma[k,j],gln,theta[[j]],param)
    }
  }
  
  
  #modify for two-parameter families.
  ddotF1=length(pconddotcop1)#length of copula functions provided
  #Check length of each element, 1 will be for one-par copulas
  # and 2 for two-par copulas.
  ddotF1j= lengths(c(theta))
  
  #Repeat parameters for two-param copulas.
  dottheta=rep(theta,times= ddotF1j)
  
  # Empty arrays to store values
  fdotden<-array(NA,dim=c(nq,K,ddotF1))
  conddotcdf<-array(NA,dim=c(nq,K-1,ddotF1))
  
  #repeat pnorma columns for two-paramter families if exist. 
  # Otherwise, it will not change.
  # The repetition is for two-par copulas, because there are two
  # derivatives functions wrt pars
  
  dotpnorma.theta=NULL
  for(i in 1:length(ddotF1j))
  {
    for(j in 1:ddotF1j[i])
    {
      dotpnorma.theta=cbind(dotpnorma.theta,pnorma[,i])
    }
  }
  
  
  for(j in 1:ddotF1)
  { 
    pconddotcopula1=pconddotcop1[[j]]
    for(k in 1:(K-1))
    { 
      conddotcdf[,k,j]<-pconddotcopula1(dotpnorma.theta[k,j],gln,dottheta[[j]])
    }
  }
  
  
  # Second Factor :::::::::::::::::::::::::::::::::::::::: Second Factor
  condcdf2<-array(NA,dim=c(nq,nq,K-1,d))
  cden2<-array(NA,dim=c(nq,nq,K-1,d))
  
  for(j in 1:d)
  { 
    pcondcopula2=pcondcop2[[j]]
    dcopula2=dcop2[[j]] 
    for(k in 1:(K-1))
    {
      for(u in 1:nq)
      { 
        condcdf2[,u,k,j]<-pcondcopula2(condcdf[u,k,j],gln,delta[[j]])
        temp<-dcopula2(condcdf[u,k,j],gln,delta[[j]])
        cden2[,u,k,j]<-temp*cden[u,k,j]
      }
    }
  }
  
  # Modify for two-parameter families.
  ddotF2=length(pconddotcop2)#length of copula functions provided
  #Check length of each element, 1 will be for one-par copulas
  # and 2 for two-par copulas.
  ddotF2j = lengths(delta)
  
  #Repeat parameters
  dotdelta=rep(delta,times= ddotF2j)
  
  dotpnorma.delta=NULL
  for(i in 1:length(ddotF2j))
  {
    for(j in 1:ddotF2j[i])
    {
      dotpnorma.delta=cbind(dotpnorma.delta,pnorma[,i])
    }
  }
  
  # Now, we repeat the 3 dimentional matrix
  # so it aligns with the two-par copulas
  # The two-par copulas will have two derivatives wrt each par.
  ddotF2j.1=ddotF2j-1
  cmsm.ddotF2=c(0,cumsum(ddotF2j.1)[-length(ddotF2j.1)])
  
  #Empty array to store values
  condcdf.2par=array(NA,dim=c(nq,K-1,ddotF2)) 
  for(i in 1:length(ddotF2j))
  {
    for(idot in 0:ddotF2j.1[i])
    {
      cmsmi=cmsm.ddotF2[i]
      condcdf.2par[,,idot+cmsmi+i]=condcdf[, ,i]
    }
  }
  
  
  
  # Repeat copulas
  dcopF2jdot=rep(dcop2,times=ddotF2j)
  
  #Empty arrays to store
  conddotcdf2<-array(NA,dim=c(nq,nq,K-1,ddotF2))
  for(j in 1:ddotF2)
  { 
    pconddotcopula2=pconddotcop2[[j]]
    for(k in 1:(K-1))
    {
      for(u in 1:nq)
      {
        conddotcdf2[,u,k,j]<-pconddotcopula2(condcdf.2par[u,k,j],gln,dotdelta[[j]])
      }
    }
  }
  
  #================================================================================
  # I need to repeat this loop according to the first factor families
  # if the first factor has 2-parameter copulas, then cden3 will
  # change accordingly
  
  ddotF1j.1=ddotF1j-1
  cmsm.ddotF1=c(0,cumsum(ddotF1j.1)[-length(ddotF1j.1)])
  
  #Empty array to store values
  condcdf.2parF1=array(NA,dim=c(nq,K-1,ddotF1)) 
  for(i in 1:length(ddotF1j))
  {
    for(idot in 0:ddotF1j.1[i])
    {
      cmsmi=cmsm.ddotF1[i]
      condcdf.2parF1[,,idot+cmsmi+i]=condcdf[, ,i]
    }
  }
  
  # Do the same again for conddotcdf.
  conddotcdf.2parF1<-array(NA,dim=c(nq,K-1,ddotF1))
  for(i in 1:length(ddotF1j))
  {
    for(idot in 0:ddotF1j.1[i])
    {
      cmsmi=cmsm.ddotF1[i]
      conddotcdf.2parF1[,,idot+cmsmi+i]=conddotcdf[, ,i]
    }
  }
  
  
  dcopF2jdotF1=rep(dcop2,times=ddotF1j)
  deltaF1=rep(delta,times=ddotF1j)
  
  cden3<-array(NA,dim=c(nq,nq,K-1,ddotF1))
  for(j in 1:ddotF1)
  { 
    dcopulaF2dot=dcopF2jdotF1[[j]] 
    for(k in 1:(K-1))
    {
      for(u in 1:nq)
      {
        temp<-dcopulaF2dot(condcdf.2parF1[u,k,j],gln,deltaF1[[j]])
        cden3[,u,k,j]<-temp*conddotcdf[u,k,j]
      }
    }
  }
  # Let 0 and 1/0 be the boundries.
  arr0<-array(0,dim=c(nq,nq,1,d))
  arr1<-array(1,dim=c(nq,nq,1,d))
  condcdf2<-abind(arr0,condcdf2,arr1,along=3)
  arr0dotF2<-array(0,dim=c(nq,nq,1,ddotF2))
  conddotcdf2<-abind(arr0dotF2,conddotcdf2,arr0dotF2,along=3)
  arr0dotF1<-array(0,dim=c(nq,nq,1,ddotF1))
  arr0dotF1<-array(0,dim=c(nq,nq,1,ddotF1))
  cden3<-abind(arr0dotF1,cden3,arr0dotF1,along=3)
  
  # Empty arrays to store the densities
  fden2<-array(NA,dim=c(nq,nq,K,d))
  fdotden2<-array(NA,dim=c(nq,nq,K,ddotF2))
  fbarden2<-array(NA,dim=c(nq,nq,K,ddotF1))
  for(j in 1:d)
  {
    for(k in 1:K)
    {
      fden2[,,k,j]<-condcdf2[,,k+1,j]-condcdf2[,,k,j]
    }
  }
  
  for(j in 1:ddotF2)
  {
    for(k in 1:K)
    {
      fdotden2[,,k,j]<-conddotcdf2[,,k+1,j]-conddotcdf2[,,k,j]
    }
  }
  
  for(j in 1:ddotF1)
  {
    for(k in 1:K)
    {
      fbarden2[,,k,j]<-cden3[,,k+1,j]-cden3[,,k,j]
    }
  }
  allcutp=rbind(0,pnorma,1)
  prob1<-matrix(NA,K,d)
  for(j in 1:d)
  {
    for(k in 1:K)
    {
      prob1[k,j]<-allcutp[k+1,j]-allcutp[k,j]
    }
  }
  list(fden2=fden2,dnorma=dnorma,cden2=cden2,fdotden2=fdotden2,
       fbarden2=fbarden2,prob1=prob1)
}


# purpose:
# To calculate all the estimated 2-factor univariate and bivariate probabilities
# input:
# prob1: the (K\times d) matrix with the univariate probabilities
# fden2: fden2 (see Appendix)
# output:
# a vector with all the estimated 2-factor univariate and bivariate probabilities
pihat.2fact<-function(prob1,fden2,Kj,glw)  #theta a vector with 2 parameters
{ 
  dims<-dim(fden2)
  d=dims[4]
  res2<-prob1[-1,] # the est. univariate probabilities
  
  #***** For distinct categories We delete the zeros as extra elements in the matrix.
  if( sum(res2 == 0) > 0 )
  {
    res2[res2==0]=NA
    res2=as.numeric(na.omit(c(res2)))     
  }else{
    res2=res2
  }
  
  for(i1 in 1:(d-1))
  {
    for(i2 in (i1+1):d)
    {
      
      Kji1= Kj[i1]  #*****
      Kji2= Kj[i2]  #*****
      
      res1<-NULL
      for(j1 in 2:Kji1)  #*****
      {
        for(j2 in 2:Kji2)  #*****
        { 
          res1<-c(res1,glw%*%((fden2[,,j1,i1]*fden2[,,j2,i2])%*%glw))
        }
      }
      res2<-c(res2,res1)
    }
  }
  res2
}

# the bivariate pairs
# d the dimension (numeric)
bivpairs.2fact<-function(d)
{
  res<-NULL
  for(id1 in 1:(d-1))
  {
    for(id2 in (id1+1):d)
    {
      res<-rbind(res,c(id1,id2))
    }
  }
  res
}

# all the observed probabilities univariate and bivariate
# dat the response data
# alpha a matrix with the cutpoints for all the responses
# bpairs the bivariate pairs
piobs.2fact<-function(dat,bpairs,Kj)
{ 
  n<-nrow(dat)
  d<-ncol(dat)
  
  Kj1=Kj-1
  res<-NULL
  for(j in 1:d )
  {
    nKj1= Kj1[j]  #******
    #X=1:k then for(jj in seq_along(X)) then dat[,j] == X[jj]
    for(jj in 1:nKj1)
    {
      res<-c(res,sum(dat[,j] == jj )/n) 
    }       
  }
  
  
  for(i in 1:nrow(bpairs))
  {
    
    y1=bpairs[i,1]
    y2=bpairs[i,2]
    
    Kjy1= Kj1[y1]
    Kjy2= Kj1[y2]
    
    #if(y1==1 ){K1 = 1}else{K1 = 4}   #******
    #if(y2==1 ){K2 = 1}else{K2 = 4}    #******
    
    
    #we can add extra arguments for different pairs
    #that have different number of categories.
    for(j1 in 1:Kjy1)
    {
      
      for(j2 in 1:Kjy2)
      { 
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
# purpose:
# to calculate the cov(Pr(y_{i1}=j1),Pr(y_{i2}=j2)) for the 2-factor model
# input:
# i_1: the  index for the 1st variable
# j_1: the value for the 1st variable
# i_2: the  index for the 2nd variable
# j_2: the value for the 2nd variable
# prob1: the (K\times d) matrix with the univariate probabilities
# fden2: fden2 (see Appendix)
# output:
# the covariance term
cov2.2fact<-function(i1,j1,i2,j2,prob1,fden2,glw)
{
  if(i1==i2 & j1!=j2)
  {
    -prob1[j1,i1]*prob1[j2,i2]
  }
  else
  {
    if(i1==i2 & j1==j2)
    {
      temp<-prob1[j1,i1]
      temp*(1-temp)
    }
    else
    {
      prob2<-glw%*%((fden2[,,j1,i1]*fden2[,,j2,i2])%*%glw)
      prob2-prob1[j1,i1]*prob1[j2,i2]
    }
  }
}



# purpose:
# to calculate cov(Pr(y_{i1}=j1),Pr(y_{i2}=j2,y_{i3}=j3) for the 2-factor model
# input:
# i_1: the  index for the 1st variable
# j_1: the value for the 1st variable
# i_2: the  index for the 2nd variable
# j_2: the value for the 2nd variable
# i_3: the  index for the 3rd variable
# j_3: the value for the 3rd variable
# prob1: the (K\times d) matrix with the univariate probabilities
# fden2: fden2 (see Appendix)
# output:
# the covariance term
cov3.2fact<-function(i1,j1,i2,j2,i3,j3,prob1,fden2,glw)
{
  if((i1==i2 & j1!=j2) | (i1==i3 & j1!=j3) )
  {
    prob2<-glw%*%((fden2[,,j2,i2]*fden2[,,j3,i3])%*%glw)
    -prob1[j1,i1]*prob2
  }
  else
  {
    if((i1==i2) | (i1==i3))
    {
      prob2<-glw%*%((fden2[,,j2,i2]*fden2[,,j3,i3])%*%glw)
      prob2*(1-prob1[j1,i1])
    }
    else
    {
      temp<-fden2[,,j2,i2]*fden2[,,j3,i3]
      prob3<-glw%*%((fden2[,,j1,i1]*temp)%*%glw)
      prob2<-glw%*%(temp%*%glw)
      prob1<-prob1[j1,i1]
      prob3-prob1*prob2
    }
  }
}


ndistinct.2fact=function(i1,i2,i3,i4)
{
  tem=unique(c(i1,i2,i3,i4))
  length(tem)
}


# purpose:
# to calculate cov(Pr(y_{i1}=j1,y_{i2}=j2),Pr(y_{i3}=j3,y_{i4}=j4) for the
# 2-factor model
# input:
# i_1: the  index for the 1st variable
# j_1: the value for the 1st variable
# i_2: the  index for the 2nd variable
# j_2: the value for the 2nd variable
# i_3: the  index for the 3rd variable
# j_3: the value for the 3rd variable
# i_4: the  index for the 3rd variable
# j_4: the value for the 3rd variable
# fden2: fden2 (see Appendix)
# output:
# the covariance term
cov4.2fact<-function(i1,j1,i2,j2,i3,j3,i4,j4,fden2,glw)
{
  nd=ndistinct.2fact(i1,i2,i3,i4)
  if(nd==4)
  {
    temp1<-fden2[,,j1,i1]*fden2[,,j2,i2]
    temp2<-fden2[,,j3,i3]*fden2[,,j4,i4]
    prob4<-glw%*%((temp1*temp2)%*%glw)
    prob21<-glw%*%(temp1%*%glw)
    prob22<-glw%*%(temp2%*%glw)
    prob4-prob21*prob22
  }
  else
  { if(nd==3)
  {
    if((i1==i3 & j1==j3) | (i2==i3 & j2==j3))
    {
      prob3<-glw%*%((fden2[,,j1,i1]*fden2[,,j2,i2]*fden2[,,j4,i4])%*%glw)
      prob21<-glw%*%((fden2[,,j1,i1]*fden2[,,j2,i2])%*%glw)
      prob22<-glw%*%((fden2[,,j3,i3]*fden2[,,j4,i4])%*%glw)
      prob3-prob21*prob22
    }
    else
    {
      if((i1==i4 & j1==j4) | (i2==i4 & j2==j4))
      {
        temp<-fden2[,,j1,i1]*fden2[,,j2,i2]
        prob3<-glw%*%((temp*fden2[,,j3,i3])%*%glw)
        prob21<-glw%*%(temp%*%glw)
        prob22<-glw%*%((fden2[,,j3,i3]*fden2[,,j4,i4])%*%glw)
        prob3-prob21*prob22
      }
      else
      {
        if((i1==i3 & j1!=j3) | (i1==i4 & j1!=j4) | (i2==i3 & j2!=j3) |(i2==i4 & j2!=j4))
        { 
          prob21<-glw%*%((fden2[,,j1,i1]*fden2[,,j2,i2])%*%glw)
          prob22<-glw%*%((fden2[,,j3,i3]*fden2[,,j4,i4])%*%glw)
          -prob21*prob22
        }
      }
    }
  }
    else 
    {
      if(nd==2)
      {
        if(j1!=j3 | j2!=j4)
        {
          prob21<-glw%*%((fden2[,,j1,i1]*fden2[,,j2,i2])%*%glw)
          prob22<-glw%*%((fden2[,,j3,i3]*fden2[,,j4,i4])%*%glw)
          -prob21*prob22
        }
        else
        {
          prob2<-glw%*%((fden2[,,j1,i1]*fden2[,,j2,i2])%*%glw)
          prob2*(1-prob2)
        }
      }
    }
  }
}



# purpose:
# to calculate the \Xi_2 matrix for the 2-factor model
# input:
# prob1: the (K\times d) matrix with the univariate probabilities
# fden2: fden2 (see Appendix)
# output:
# the \Xi_2 matrix
Xi2.2fact<-function(prob1,fden2,Kj,glw){ 
  dims<-dim(fden2)
  d=dims[4]
  d1=d-1
  rows<-NULL
  for(i1 in 1:d){
    Kji1= Kj[i1]   
    for(j1 in 2:Kji1){
      res<-NULL
      for(i2 in 1:d){
        Kji2= Kj[i2]  
        for(j2 in 2:Kji2){
          res<-c(res,cov2.2fact(i1,j1,i2,j2,prob1,fden2,glw))
        }
      }
      for(i2 in 1:d1){
        Kji2= Kj[i2]   
        for(i3 in (i2+1):d){
          Kji3= Kj[i3]  
          for(j2 in 2:Kji2){
            for(j3 in 2:Kji3){
              res<-c(res,cov3.2fact(i1,j1,i2,j2,i3,j3,prob1,fden2,glw))
            }
          }
        }
      }
      rows<-rbind(rows,res)
    }
  }
  
  for(i1 in 1:d1){
    Kji1= Kj[i1]    
    for(i2 in (i1+1):d){
      Kji2= Kj[i2]   
      for(j1 in 2:Kji1){
        for(j2 in 2:Kji2){
          res<-NULL
          for(i3 in 1:d){
            Kji3= Kj[i3]   
            for(j3 in 2:Kji3){
              res<-c(res,cov3.2fact(i3,j3,i1,j1,i2,j2,prob1,fden2,glw))
            }
          }
          for(i3 in 1:d1){
            Kji3= Kj[i3]   
            for(i4 in (i3+1):d){
              Kji4= Kj[i4]   
              for(j3 in 2:Kji3){
                for(j4 in 2:Kji4){
                  res<-c(res,cov4.2fact(i1,j1,i2,j2,i3,j3,i4,j4,fden2,glw))
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
# purpose:
# to calculate the derivative of \Pr_{i1 i2, y1 y2} wrt a_{cutp i}
# for the 2-factor model
# input:
# i_1: the  index for the 1st variable
# i_2: the  index for the 2nd variable
# y_1: the value for the 1st variable
# y_2: the value for the 2nd variable
# i: the index for the cutpoint
# cutp: A matrix with the cutpoints (without the boundary cutpoints) in the
# uniform scale where the number of rows corresponds to an ordinal category
# and each column represents an item.
# fden2: fden2 (see Appendix)
# dnorma: \phi[\Phi^{-1}(cutp)]
# cden2: cdens12 (see Appendix)
# output:
# the derivative of \Pr_{i1 i2, y1 y2} wrt a_{cutp i}
der.bivprob.cutp.2fact<-function(i1,i2,y1,y2,i,cutp,dnorma,cden2,fden2,glw){
  if(i1==i & (y1-1)==cutp){
    inner<-dnorma[cutp+1,i]*cden2[,,cutp+1,i]*fden2[,,y2,i2]
    glw%*%(inner%*%glw)
  }else{
    if(i2==i & (y2-1)==cutp){
      inner<-dnorma[cutp+1,i]*cden2[,,cutp+1,i]*fden2[,,y1,i1]
      glw%*%(inner%*%glw)
    }else{
      if(i1==i & cutp==(y1-2)){
        inner<--dnorma[cutp+1,i]*cden2[,,cutp+1,i]*fden2[,,y2,i2]
        glw%*%(inner%*%glw)
      }else{
        if(i2==i & cutp==(y2-2)){
          inner<--dnorma[cutp+1,i]*cden2[,,cutp+1,i]*fden2[,,y1,i1]
          glw%*%(inner%*%glw)
        }else{
          0
        }
      }
    }
  }
}


# purpose:
# to calculate the derivatives of \Pr_{i1 i2, y1 y2} wrt a_{cutp i} for all the
# bivariate probabilities and cutpoints for the 2-factor model
# input:
# fden2: fden2 (see Appendix)
# dnorma: \phi[\Phi^{-1}(cutp)]
# cden2: cdens12 (see Appendix)
# output:
# a matrix with the derivatives of \Pr_{i1 i2, y1 y2} wrt a_{cutp i} for all the
# bivariate probabilities and cutpoints
all.der.bivprob.cutp.2fact<-function(dnorma,cden2,fden2,Kj,glw){ 
  dims<-dim(fden2)
  d=dims[4]
  Kj2=Kj-2
  cols<-NULL
  for(i in 1:d){
    Kj2i= Kj2[i]  
    for(cutp in 0:Kj2i){
      res<-NULL
      for(i1 in 1:(d-1)){
        for(i2 in (i1+1):d){
          Kji1= Kj[i1]  
          for(y1 in 2:Kji1){
            Kji2= Kj[i2]  
            for(y2 in 2:Kji2){
              res<-c(res,der.bivprob.cutp.2fact(i1,i2,y1,y2,i,cutp,dnorma,cden2,fden2,glw))
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


# purpose:
# to calculate the derivative of \Pr_{i1 i2, y1 y2} wrt delta_i
# for the 2-factor model
# input:
# i_1: the  index for the 1st variable
# i_2: the  index for the 2nd variable
# y_1: the value for the 1st variable
# y_2: the value for the 2nd variable
# i: the index for the copula vector
# fden2: fden2 (see Appendix)
# fdotden2: fdotden2 (see Appendix)
# output:
# the derivative of \Pr_{i1 i2, y1 y2} wrt delta_i
der.bivprob.delta.2fact<-function(i1,i2,y1,y2,i,cmsm.ddotj,fden2,fdotden2,glw){
  if(i1==i){
    inner<-fden2[,,y2,i2]*fdotden2[,,y1,i+cmsm.ddotj]
    glw%*%(inner%*%glw)
  } else {
    if(i2==i)
    {
      inner<-fden2[,,y1,i1]*fdotden2[,,y2,i+cmsm.ddotj]
      glw%*%(inner%*%glw)
    } else {
      0
    }
  }
}

# purpose:
# to calculate the derivatives of \Pr_{i1 i2, y1 y2} wrt delta_i for all the
# bivariate probabilities and copula parameters at the 2nd factor
# for the 2-factor model
# input:
# fden: fden (see Appendix)
# fdotden: fdotden (see Appendix)
# output:
# a matrix with the derivatives of \Pr_{i1 i2, y1 y2} wrt delta_i  for all the
# bivariate probabilities and copula parameters at the 2nd factor
all.der.bivprob.delta.2fact<-function(fden2,fdotden2,delta,Kj,glw, SpC){ 
  dims<-dim(fden2)
  d=dims[4]
  # length of parameters (number of der.cond.cop. wrt parameters)
  ddotj=lengths(delta)-1
  cmsm.ddot=c(0,cumsum(ddotj)[-length(ddotj)])
  cols<-NULL
  
  #adjusting the number of columns if SpC is TRUE or FALSE.
  if (is.null(SpC)){
    dF2delta = 1:d 
  } else {
    dF2delta = (1:d)[-SpC]
  }
  
  for(i in dF2delta){
    for(idot in 0:ddotj[i]){
      cmsm.ddotj = cmsm.ddot[i]+idot
      res<-NULL
      for(i1 in 1:(d-1)){
        for(i2 in (i1+1):d){
          Kji1= Kj[i1]
          for(y1 in 2:Kji1){
            Kji2= Kj[i2] 
            for(y2 in 2:Kji2){
              res<-c(res,der.bivprob.delta.2fact(i1,i2,y1,y2,i,cmsm.ddotj,fden2,fdotden2,glw))
            }
          }
        } 
      }
      cols<-cbind(cols,res)
    }
  }
  cols
}
###############################################################################

# purpose:
# to calculate the derivative of \Pr_{i1 i2, y1 y2} wrt theta_i
# for the 2-factor model
# input:
# i_1: the  index for the 1st variable
# i_2: the  index for the 2nd variable
# y_1: the value for the 1st variable
# y_2: the value for the 2nd variable
# i: the index for the copula vector
# fden2: fden2 (see Appendix)
# fbarden2: fbarden (see Appendix)
# output:
# the derivative of \Pr_{i1 i2, y1 y2} wrt theta_i
der.bivprob.theta.2fact<-function(i1,i2,y1,y2,i,cmsm.ddotj,fden2,fbarden2,glw){
  if(i1==i){
    inner<-fden2[,,y2,i2]*fbarden2[,,y1,i+cmsm.ddotj]
    glw%*%(inner%*%glw)
    } 
  else {
    if(i2==i){
      inner<-fden2[,,y1,i1]*fbarden2[,,y2,i+cmsm.ddotj]
      glw%*%(inner%*%glw)
    }else{0}
  }
}

# purpose:
# to calculate the derivatives of \Pr_{i1 i2, y1 y2} wrt theta_i for all the
# bivariate probabilities and copula parameters at the 1st factor
# for the 2-factor model
# input:
# fden: fden (see Appendix)
# fdotden: fdotden (see Appendix)
# output:
# a matrix with the derivatives of \Pr_{i1 i2, y1 y2} wrt theta_i  for all the
# bivariate probabilities and copula parameters at the 1st factor
all.der.bivprob.theta.2fact<-function(fden2,fbarden2,theta,Kj,glw){
  dims<-dim(fden2)
  d=dims[4]
  # length of parameters (number of der.cond.cop. wrt parameters)
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
              res<-c(res,der.bivprob.theta.2fact(i1,i2,y1,y2,i,cmsm.ddotj,fden2,fbarden2,glw))
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
# purpose:
# to calculate all the derivatives of the univariate probabilities wrt to
# (a) the cutpoints  and (b) the copula parameters for the 2-factor model
# input:
# dnorma: \phi[\Phi^{-1}(cutp)]
# output: the matrix with derivatives of all univariate
# probabilities wrt to the cutpoints and copula parameters
################################################################################
# all the derivatives of the univariate probabilities wrt to the cutpoints
all.der.uprob.2fact<-function(dnorma,theta,delta, SpC){
  d=ncol(dnorma)
  ddot1=sum(lengths(theta))
  ddot2=sum(lengths(delta))
  ddot=ddot1+ddot2
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
  
  # Return the matrix with derivatives of univaraite margins
  # For two-factor copula model with BVN copulas (if SpC is TRUE) 
  # we exclude the last columnn that correspond to the nth 
  # copula that is set to independence.
  if ( is.null(SpC) ) {
    output = result
  }else {
    output = result[ , -SpC] 
  }
  return(output)
}

################################################################################
# purpose:
# to calculate  the Delta_2 matrix for the 2-factor model
# input:
# dnorma: \phi[\Phi^{-1}(cutp)]
# cden2: cdens12 (see Appendix)
# fden2: fden2 (see Appendix)
# fdotden2: fdotden2 (see Appendix)
# fbarden2: fbarden2 (see Appendix)
# output:
# the Delta_2 matrix
Delta2.2fact<-function(dnorma,cden2,fden2,fdotden2,fbarden2,theta,delta,Kj,glw, SpC)
{ #the derivatives of the univariate probabilities
  res1<-all.der.uprob.2fact(dnorma,theta,delta, SpC)
  #the derivatives of the bivariate probabilities wrt cutpoints
  res2<-all.der.bivprob.cutp.2fact(dnorma,cden2,fden2,Kj,glw)
  #the derivatives of the bivariate probabilities wrt copula parmaeters
  res3<-all.der.bivprob.theta.2fact(fden2,fbarden2,theta,Kj,glw)
  res4<-all.der.bivprob.delta.2fact(fden2,fdotden2,delta,Kj,glw, SpC)
  rbind(res1,cbind(res2,res3,res4))
}
################################################################################
# R code for orthogonal complement (used by one of my students)
orthcomp <- function(x,tol=1.e-12)
{ 
  s <- nrow(x)
  q <- ncol(x)
  if (s<=q) { return('error, matrix must have more columns than rows') }
  x.svd <- svd(x,s)
  if (sum(abs(x.svd$d)<tol)!=0)
  { return('error, matrix must full column rank') }
  return(x.svd$u[,(q+1):s])
}

# purpose:
# to calculate  the M_2 for the 2-factor model
# input:
# dat: A data matrix where the number of rows corresponds to an
# individual's response and each column represents an item
# Number of ordinal categories for each item, coded as 0,...,(ncat-1).
# Currently supported are items that have the same number of categories.
# dnorma: \phi[\Phi^{-1}(cutp)]
# prob1: the (K\times d) matrix with the univariate probabilities
# fden2: fden2 (see Appendix)
# fdotden2: fdotden2 (see Appendix)
# cden2: cdens12 (see Appendix)
# fbarden2: fbarden2 (see Appendix)
# iprint: debug option
# output:
# a list with
# "M_2"=M2 stat
# "df"=degrees of freedom
# "p-value"=pvalue of M2
# "diagnostic"=biavriate diagnostics
M2.2fact<-function(dat,dnorma,prob1,cden2,fden2,fdotden2,
                   fbarden2, theta, delta, Kj, glw, SpC, iprint){ 
  dims<-dim(fden2)
  d=dims[4]
  n=nrow(dat)
  bpairs<-bivpairs.2fact(d)
  #-------------------------------------------------------
  V<-Xi2.2fact(prob1,fden2,Kj,glw)
  #-------------------------------------------------------
  if(iprint){
    cat("\\ndim(V): dimension of Xi_2 matrix\\n")
    print(dim(V))
  }
  #-------------------------------------------------------
  delta<-Delta2.2fact(dnorma,cden2,fden2,fdotden2,fbarden2,theta,delta,Kj,glw, SpC)
  #-------------------------------------------------------
  if(iprint){ 
    cat("\\ndim(delta): dimension of Delta_2 matrix\\n")
    print(dim(delta))
  }
  #-------------------------------------------------------
  oc.delta<-orthcomp(delta)
  inv<-solve(t(oc.delta)%*%V%*%oc.delta)
  C2<-oc.delta%*%inv%*%t(oc.delta)
  pi.r<-pihat.2fact(prob1,fden2,Kj,glw)
  p.r<-piobs.2fact(dat,bpairs,Kj)
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


