# purpose:
# calulates the arrays fden2, cden2, fdotden2, fbarden2, etc for the computation
# of the M2 for the 2-factor model (see appendix of paper)
# input:

store.bifact<-function(cutp,thX0,thXg,
                       pcondX0,pcondXg,pconddotX0,pconddotXg,
                       dcopX0,dcopXg,nq,glw,gln){
  d=ncol(cutp)
  K=nrow(cutp)+1
  #Univariate probabilities (first-order marginals)
  dnorma=dnorm(qnorm(cutp))
  allcutp=rbind(0,cutp,1)
  prob1<-matrix(NA,K,d)
  for(j in 1:d){
    for(k in 1:K){
      prob1[k,j]<-allcutp[k+1,j]-allcutp[k,j]
    }
  }
  #Conditional copula "condcdf", derivative of conditional
  #copula "conddotcdf",and the copula density "cden" for the observed
  #variables and the common latent variable.
  condcdf<-array(NA,dim=c(nq,K-1,d))
  conddotcdf<-array(NA,dim=c(nq,K-1,d))
  cden<-array(NA,dim=c(nq,K-1,d))
  for(j in 1:d){
    pcondX0.j=pcondX0[[j]]
    pconddotX0.j=pconddotX0[[j]]
    dcopX0.j=dcopX0[[j]]
    for(k in 1:(K-1)){
      condcdf[,k,j]<-pcondX0.j(cutp[k,j],gln,thX0[j])
      conddotcdf[,k,j]<-pconddotX0.j(cutp[k,j],gln,thX0[j])
      cden[,k,j]<-dcopX0.j(cutp[k,j],gln,thX0[j])
    }
  }
  #Conditional copula "condcdf2", derivative of conditional
  #copula "conddotcdf2",and the copula density "cden2" for the observed
  #variables and the group-specifc latent variable.
  #cden3: derivatives of conditional copula wrt to
  #     first factor paramaeters. Here we use the chain rule
  #     i.e. c_{jX2}(C_{J|x_0},x_g) * Cdot_{j|X_0}
  condcdf2<-array(NA,dim=c(nq,nq,K-1,d))
  conddotcdf2<-array(NA,dim=c(nq,nq,K-1,d))
  cden2<-array(NA,dim=c(nq,nq,K-1,d))
  cden3<-array(NA,dim=c(nq,nq,K-1,d))
  for(j in 1:d){
    pcondXg.j=pcondXg[[j]]
    dcopXg.j=dcopXg[[j]]
    pconddotXg.j=pconddotXg[[j]]
    for(k in 1:(K-1))
    {
      for(u in 1:nq){
        condcdf2[,u,k,j]<-pcondXg.j(condcdf[u,k,j],gln,thXg[j])
        temp<-dcopXg.j(condcdf[u,k,j],gln,thXg[j])
        cden2[,u,k,j]<-temp*cden[u,k,j]
        cden3[,u,k,j]<-temp*conddotcdf[u,k,j]
        conddotcdf2[,u,k,j]<-pconddotXg.j(condcdf[u,k,j],gln,thXg[j])
      }
    }
  }
  # Here we calculate the finite differences
  # for all copula functions calculated above.
  arr0<-array(0,dim=c(nq,nq,1,d))
  arr1<-array(1,dim=c(nq,nq,1,d))
  condcdf2<-abind(arr0,condcdf2,arr1,along=3)
  conddotcdf2<-abind(arr0,conddotcdf2,arr0,along=3)
  cden3<-abind(arr0,cden3,arr0,along=3)
  fden2<-array(NA,dim=c(nq,nq,K,d))
  fbarden2<-array(NA,dim=c(nq,nq,K,d))
  fdotden2<-array(NA,dim=c(nq,nq,K,d))
  for(j in 1:d){
    for(k in 1:K){
      fden2[,,k,j]<-condcdf2[,,k+1,j]-condcdf2[,,k,j]
      fdotden2[,,k,j]<-conddotcdf2[,,k+1,j]-conddotcdf2[,,k,j]
      fbarden2[,,k,j]<-cden3[,,k+1,j]-cden3[,,k,j]
    }
  }
  list(fden2=fden2,dnorma=dnorma,cden2=cden2,fdotden2=fdotden2,
       fbarden2=fbarden2,prob1=prob1)
}

# purpose:
# To calculate all the observed univariate and bivariate probabilities
# input:
# dat: A data matrix where the number of rows corresponds to an
# individual's response and each column represents an item
# Number of ordinal categories for each item, coded as 0,...,(ncat-1).
# Currently supported are items that have the same number of categories.
# bpairs: the possible bivariate (pairwise) pairs
# output:
# a vector with all the observed
# univariate (pobs1) and bivariate probabilities (pobs2)
pobs1<-function(dat)
{
  n=nrow(dat)
  #excluding zero categories
  p1=c(apply(dat,2,table)[-1,]/n)
  p1
}

pobs2<-function(dat){
  d<-ncol(dat)
  n<-nrow(dat)
  K<-length(unique(dat[,1]))-1
  p2=NULL
  for(j1 in 1:(d-1)){
    for(j2 in (j1+1):d){
      p2=c(p2,c(table(as.factor(dat[,j2]),
                      as.factor(dat[,j1]))[-1,-1]/n))
    }
  }
  p2
}

#combining both observed probabilities
pobs<-function(dat){
  p1=pobs1(dat)
  p2=pobs2(dat)
  p=c(p1,p2)
  return(p)
}

bivpairs<-function(d)
{ res<-NULL
for(id1 in 1:(d-1))
{ for(id2 in (id1+1):d)
{ res<-rbind(res,c(id1,id2)) } }
res
}

#Model based univaraite ans bivaraite probabilities
pihat<-function(dat,ngrp,grpsize,prob1,fden2,glw)
{
  d=ncol(dat)
  n=nrow(dat)
  K<-length(unique(dat[,1]))

  #uninvariate probabilities
  pihat1=prob1[-1,]
  bpairs=bivpairs(d)

  res2=NULL
  for(i in 1:nrow(bpairs)){
    pairsit=bpairs[i,]
    grpsp=list()
    ind=1
    group1=rep(FALSE,ngrp)
    for(g in 1:ngrp){
      grpsp[[g]]=ind:cumsum(grpsize)[g]
      ind=ind+grpsize[g]
      if(pairsit[1]%in%grpsp[[g]]&&pairsit[2]%in%grpsp[[g]]){
        group1[g]=TRUE}
    }

    if(sum(group1)==1){
      res1<-NULL
      for(j1 in 2:K){
        for(j2 in 2:K){
          res1<-c(res1,glw%*%as.vector( glw%*%(fden2[,,j1,pairsit[1]]*
                                                   fden2[,,j2,pairsit[2]]) ) )
        }
      }
      res1
    }
    else{
      res1<-NULL
      for(j1 in 2:K){
        for(j2 in 2:K){
          res1<-c(res1, glw%*%as.vector((glw%*%fden2[,,j1,pairsit[1]])*
                                           (glw%*%fden2[,,j2,pairsit[2]])))
          #sumtemp1=rep(NA,25)
          #for(i in 1:25){sumtemp1[i]=sum(fden2[,i,j1,pairsit[1]]*glw)}
          #sumtemp2=rep(NA,25)
          #for(i in 1:25){sumtemp2[i]=sum(fden2[,i,j2,pairsit[2]]*glw)}
          #glw%*% (sumtemp1*sumtemp2)

        }
      }
      res1
    }
    res2=c(res2,res1)
  }
  return(c(pihat1,res2))
}

#============================= Covariance


cov2.bifactor<-function(i1,j1,i2,j2,prob1,fden2,ngrp,grpsize,glw){
  pairs.cov2=c(i1,i2)#pair of variables
  grplst=list()
  ind=1
  onegrp=rep(FALSE,ngrp)
  #loop to create a list of the index of variables in each group.
  for(g in 1:ngrp){
    grplst[[g]]=ind:cumsum(grpsize)[g]
    ind=ind+grpsize[g]
    if(pairs.cov2[1]%in%grplst[[g]]&&pairs.cov2[2]%in%grplst[[g]]){
      onegrp[g]=TRUE}
  }
  #univariate probabilities of the pair.
  prob1i1=prob1[j1,i1]; prob1i2=prob1[j2,i2]
  if(sum(onegrp)==1){
    #This is if the variables are in the same group.
    temp1=fden2[,,j1,i1];temp2=fden2[,,j2,i2]
    bifden12=glw%*%(temp1*temp2)
  }else{
    #This is if the variables are in different groups.
    temp1=glw%*%fden2[,,j1,i1];temp2=glw%*%fden2[,,j2,i2]
    bifden12=(temp1*temp2)
  }
  if(i1==i2 & j1!=j2) {
    -prob1i1*prob1i2
  }else {
    if(i1==i2 & j1==j2) {
      temp<-prob1i1
      temp*(1-temp)
    } else {
      prob2<-glw%*%as.vector(bifden12)
      prob2-prob1i1*prob1i2}
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
cov3.bifactor<-function(i1,j1,i2,j2,i3,j3,prob1,fden2,ngrp,grpsize,glw)
{
  pairs.cov3=c(i1,i2,i3)#trivariate variables.
  ngrpcov3=ngrp#it cannot be more than ngrp groups.
  grplst=list()
  ind=1
  item1=rep(FALSE,ngrpcov3);item2=rep(FALSE,ngrpcov3);item3=rep(FALSE,ngrpcov3)
  y1=rep(FALSE,ngrpcov3);y2=rep(FALSE,ngrpcov3);y3=rep(FALSE,ngrpcov3)
  for(g in 1:ngrp){#ngrp of all data
    grplst[[g]]=ind:cumsum(grpsize)[g]
    ind=ind+grpsize[g]
    if(pairs.cov3[1]%in%grplst[[g]]){item1[g]=1;y1[g]=1}else{item1[g]=FALSE;y1[g]=FALSE}
    if(pairs.cov3[2]%in%grplst[[g]]){item2[g]=1;y2[g]=1}else{item2[g]=FALSE;y2[g]=FALSE}
    if(pairs.cov3[3]%in%grplst[[g]]){item3[g]=1;y3[g]=1}else{item3[g]=FALSE;y3[g]=FALSE}
  }
  #How many groups/testlets in pairs.cov3?
  testlets=length(unique(c(which(item1==1),which(item2==1),which(item3==1))))
  #place the items according to their groups.

  grpmat=cbind((pairs.cov3)*rbind(item1,item2,item3))
  #exclude the extra unused groups.
  ngrpmat=cbind(grpmat[ , colSums(grpmat != 0) != 0])

  #What the categories used wrt to each item and for each group?
  y=cbind(c(j1,j2,j3)*rbind(y1,y2,y3))
  #exclude the extra unused groups.
  ny=cbind(y[ , colSums(y != 0) != 0])

  #The possible combinations of items
  itemgrp12=c(which(ngrpmat[1,]>0),which(ngrpmat[2,]>0))
  itemgrp13=c(which(ngrpmat[1,]>0),which(ngrpmat[3,]>0))
  itemgrp23=c(which(ngrpmat[2,]>0),which(ngrpmat[3,]>0))
  itemgrp123=c(which(ngrpmat[1,]>0),which(ngrpmat[2,]>0),which(ngrpmat[3,]>0))

  #compute corresponding mass functinos.
  if(length(unique(itemgrp12))==1){
    fden2.12=glw%*%(fden2[,,j1,i1]*fden2[,,j2,i2])
  }else{fden2.12=(glw%*%fden2[,,j1,i1])*(glw%*%fden2[,,j2,i2])}

  if(length(unique(itemgrp13))==1){
    fden2.13=glw%*%(fden2[,,j1,i1]*fden2[,,j3,i3])
  }else{fden2.13=(glw%*%fden2[,,j1,i1])*(glw%*%fden2[,,j3,i3])}

  if(length(unique(itemgrp23))==1){
    fden2.23=glw%*%(fden2[,,j2,i2]*fden2[,,j3,i3])
  }else{fden2.23=(glw%*%fden2[,,j2,i2])*(glw%*%fden2[,,j3,i3])}

  if(length(unique(itemgrp123))==1){
    #This is if items 1, 2, and 3 are in the same group
    fden2.123=glw%*%(fden2[,,j1,i1]*fden2[,,j2,i2]*fden2[,,j3,i3])
  }else if(length(unique(itemgrp123))==2){
    #This is if items 1, 2, and 3 are in two groups.
    tripairs=list()
    for(g in 1:ngrp){
      tripairs[[g]]=as.numeric(c(grpmat[,g][grpmat[,g]!=0]))
    }
    yj=c(j1,j2,j3)
    ymat=list(which(ny[,1]>0),which(ny[,2]>0))
    ymat=ymat[order(sapply(ymat,length),decreasing=F)]
    ni1= unlist(tripairs[lengths(tripairs)==1])[1]#group1
    yj1=yj[ymat[[1]]]
    ni2= unlist(tripairs[lengths(tripairs)==2])[1]#group2
    yj2=yj[ymat[[2]][1]]
    ni3= unlist(tripairs[lengths(tripairs)==2])[2]#group2
    yj3=yj[ymat[[2]][2]]
    fden2.123=((glw%*%fden2[,,yj1,ni1])*(glw%*%(fden2[,,yj2,ni2]*fden2[,,yj3,ni3])))
  }else if(length(unique(itemgrp123))==3){
    fden2.123=((glw%*%fden2[,,j1,i1])*(glw%*%fden2[,,j2,i2])*(glw%*%fden2[,,j3,i3]))
  }

  if((i1==i2 & j1!=j2) | (i1==i3 & j1!=j3) ){
    prob2<-glw%*%as.vector(fden2.23)
    -prob1[j1,i1]*prob2
  } else {
    if((i1==i2) | (i1==i3))
    { prob2<-glw%*%as.vector(fden2.23)
    prob2*(1-prob1[j1,i1])
    } else {
      prob3<-glw%*%as.vector(fden2.123)
      prob2<-glw%*%as.vector(fden2.23)
      prob1<-prob1[j1,i1]
      prob3-prob1*prob2
    }}
}

# purpose:
# to calculate the length of distinct values in (i_1,i_2,i_3,i_4)
ndistinct=function(i1,i2,i3,i4)
{ tem=unique(c(i1,i2,i3,i4))
length(tem)
}

# purpose:
# to calculate cov(Pr(y_{i1}=j1,y_{i2}=j2),Pr(y_{i3}=j3,y_{i4}=j4) for the
# 1-factor model
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
cov4.bifactor<-function(i1,j1,i2,j2,i3,j3,i4,j4,fden2,ngrp,grpsize,glw){
  nd=ndistinct(i1,i2,i3,i4)
  pairs.cov4=c(i1,i2,i3,i4)
  ngrpcov4=ngrp#it cannot be more than ngrp groups
  grplst=list()
  ind=1
  item1=rep(FALSE,ngrpcov4);item2=rep(FALSE,ngrpcov4)
  item3=rep(FALSE,ngrpcov4);item4=rep(FALSE,ngrpcov4)
  y1=rep(FALSE,ngrpcov4);y2=rep(FALSE,ngrpcov4)
  y3=rep(FALSE,ngrpcov4);y4=rep(FALSE,ngrpcov4)

  for(g in 1:ngrp){#ngrp of all data
    grplst[[g]]=ind:cumsum(grpsize)[g]
    ind=ind+grpsize[g]
    if(pairs.cov4[1]%in%grplst[[g]]){item1[g]=1;y1[g]=1}else{item1[g]=FALSE;y1[g]=FALSE}
    if(pairs.cov4[2]%in%grplst[[g]]){item2[g]=1;y2[g]=1}else{item2[g]=FALSE;y2[g]=FALSE}
    if(pairs.cov4[3]%in%grplst[[g]]){item3[g]=1;y3[g]=1}else{item3[g]=FALSE;y3[g]=FALSE}
    if(pairs.cov4[4]%in%grplst[[g]]){item4[g]=1;y4[g]=1}else{item4[g]=FALSE;y4[g]=FALSE}
  }
  ngrp.i=length(unique(c(which(item1==1),which(item2==1),which(item3==1),which(item4==1))))
  grpmat=cbind((pairs.cov4)*rbind(item1,item2,item3,item4))
  ngrpmat=cbind(grpmat[ , colSums(grpmat != 0) != 0] )

  #jmat=rbind(jind1,jind2,jind3,jind4)
  #njmatt=cbind(jmat[ , colSums(jmat != 0) != 0])
  #
  ymat=cbind(c(j1,j2,j3,j4)*rbind(y1,y2,y3,y4))
  #njmat=cbind((c(j1,j2,j3,j4))*jmat[ , colSums(jmat != 0) != 0])

  itemgrp12=c(which(ngrpmat[1,]>0),which(ngrpmat[2,]>0))
  itemgrp34=c(which(ngrpmat[3,]>0),which(ngrpmat[4,]>0))
  itemgrp123=c(which(ngrpmat[1,]>0),which(ngrpmat[2,]>0),which(ngrpmat[3,]>0))
  itemgrp124=c(which(ngrpmat[1,]>0),which(ngrpmat[2,]>0),which(ngrpmat[4,]>0))
  itemgrp1234=c(which(ngrpmat[1,]>0),which(ngrpmat[2,]>0),
                which(ngrpmat[3,]>0),which(ngrpmat[4,]>0))

  if(length(unique(itemgrp12))==1){
    fden2.12=glw%*%(fden2[,,j1,i1]*fden2[,,j2,i2])
  }else{fden2.12=(glw%*%fden2[,,j1,i1])*(glw%*%fden2[,,j2,i2])}


  if(length(unique(itemgrp34))==1){
    fden2.34=glw%*%(fden2[,,j3,i3]*fden2[,,j4,i4])
  }else{fden2.34=(glw%*%fden2[,,j3,i3])*(glw%*%fden2[,,j4,i4])}

  if(length(unique(itemgrp123))==1){
    fden2.123=glw%*%(fden2[,,j1,i1]*fden2[,,j2,i2]*fden2[,,j3,i3])
  }else if(length(unique(itemgrp123))==2){
    tripairs=list()
    jlst=list()
    for(g in 1:ngrp){
      tripairs[g]=list(c(grpmat[c(1,2,3),g][grpmat[c(1,2,3),g]!=0]))
      jlst[g]=list(c(ymat[c(1,2,3),g][ymat[c(1,2,3),g]!=0]))
    }
    jlst=jlst[lengths(jlst)>0 ];
    tripairs=tripairs[lengths(tripairs)>0 ];
    prodg1=1;prodg2=1
    for(i in 1:length(tripairs[[1]])){
      prodg1=prodg1*fden2[,,jlst[[1]][i],tripairs[[1]][i]]
    }
    for(i in 1:length(tripairs[[2]])){
      prodg2=prodg2*fden2[,,jlst[[2]][i],tripairs[[2]][i]]
    }
    fden2.123=(glw%*%prodg1)*(glw%*%prodg2)
  } else if(length(unique(itemgrp123))==3){
    fden2.123=((glw%*%fden2[,,j1,i1])*(glw%*%fden2[,,j2,i2])*(glw%*%fden2[,,j3,i3]))
  }
  if(length(unique(itemgrp124))==1){
    fden2.124=glw%*%(fden2[,,j1,i1]*fden2[,,j2,i2]*fden2[,,j4,i4])
  }else if(length(unique(itemgrp124))==2){
    tripairs=list()
    jlst=list()
    for(g in 1:ngrp){
      tripairs[[g]]=as.numeric(c(grpmat[c(1,2,4),g][grpmat[c(1,2,4),g]!=0]))
      jlst[g]=list(c(ymat[c(1,2,4),g][ymat[c(1,2,4),g]!=0]))
    }
    jlst=jlst[lengths(jlst)>0 ];
    tripairs=tripairs[lengths(tripairs)>0 ];
    prodg1=1;prodg2=1

    for(i in 1:length(tripairs[[1]])){
      prodg1=prodg1*fden2[,,jlst[[1]][i],tripairs[[1]][i]]
    }
    for(i in 1:length(tripairs[[2]])){
      prodg2=prodg2*fden2[,,jlst[[2]][i],tripairs[[2]][i]]
    }

    fden2.124=(glw%*%prodg1)*(glw%*%prodg2)

  }  else if(length(unique(itemgrp124))==3){
    fden2.124=((glw%*%fden2[,,j1,i1])*(glw%*%fden2[,,j2,i2])*(glw%*%fden2[,,j4,i4]))
  }

  if(length(unique(itemgrp1234))==4){
    fden2.1234=((glw%*%fden2[,,j1,i1])*(glw%*%fden2[,,j2,i2])*
                  (glw%*%fden2[,,j3,i3])*(glw%*%fden2[,,j4,i4]))
  }else if(length(unique(itemgrp1234))==1){
    fden2.1234=glw%*%(fden2[,,j1,i1]*fden2[,,j2,i2]*
                         fden2[,,j3,i3] * fden2[,,j4,i4])
  }else if(length(unique(itemgrp1234))==3){

    quadpairs=list()
    jlst=list()
    for(g in 1:ngrp){
      quadpairs[[g]]=as.numeric(c(grpmat[,g][grpmat[,g]!=0]))
      jlst[[g]]=as.numeric(c(ymat[,g][ymat[,g]!=0]))
    }
    jlst=jlst[lengths(jlst)>0 ];
    quadpairs=quadpairs[lengths(quadpairs)>0 ];
    prodg1=1;prodg2=1;prodg3=1
    for(i in 1:length(quadpairs[[1]])){
      prodg1=prodg1*fden2[,,jlst[[1]][i],quadpairs[[1]][i]]
    }
    for(i in 1:length(quadpairs[[2]])){
      prodg2=prodg2*fden2[,,jlst[[2]][i],quadpairs[[2]][i]]
    }
    for(i in 1:length(quadpairs[[3]])){
      prodg3=prodg3*fden2[,,jlst[[3]][i],quadpairs[[3]][i]]
    }
    fden2.1234=(glw%*%prodg1)*(glw%*%prodg2)*(glw%*%prodg3)

  }else if(length(unique(itemgrp1234))==2){
    quadpairs=list()
    jlst=list()
    for(g in 1:ngrp){
      quadpairs[[g]]=as.numeric(c(grpmat[,g][grpmat[,g]!=0]))
      jlst[[g]]=as.numeric(c(ymat[,g][ymat[,g]!=0]))
    }
    jlst=jlst[lengths(jlst)>0 ];
    quadpairs=quadpairs[lengths(quadpairs)>0 ];
    prodg1=1;prodg2=1
    for(i in 1:length(quadpairs[[1]])){
      prodg1=prodg1*fden2[,,jlst[[1]][i],quadpairs[[1]][i]]
    }
    for(i in 1:length(quadpairs[[2]])){
      prodg2=prodg2*fden2[,,jlst[[2]][i],quadpairs[[2]][i]]
    }
    fden2.1234=(glw%*%prodg1)*(glw%*%prodg2)
  }

  if(nd==4){
    prob4<-glw%*%as.vector(fden2.1234)
    prob21<-glw%*%as.vector(fden2.12)
    prob22<-glw%*%as.vector(fden2.34)
    prob4-prob21*prob22
  } else {
    if(nd==3){
      if((i1==i3 & j1==j3) | (i2==i3 & j2==j3)){
        prob3<-glw%*%as.vector(fden2.124)
        prob21<-glw%*%as.vector(fden2.12)
        prob22<-glw%*%as.vector(fden2.34)
        prob3-prob21*prob22
      } else {
        if((i1==i4 & j1==j4) | (i2==i4 & j2==j4)){
          prob3<-glw%*%as.vector(fden2.123)
          prob21<-glw%*%as.vector(fden2.12)
          prob22<-glw%*%as.vector(fden2.34)
          prob3-prob21*prob22
        } else {
          if((i1==i3 & j1!=j3) | (i1==i4 & j1!=j4) |
             (i2==i3 & j2!=j3) |(i2==i4 & j2!=j4)){
            prob21<-glw%*%as.vector(fden2.12)
            prob22<-glw%*%as.vector(fden2.34)
            -prob21*prob22
          }
        }
      }
    } else {
      if(nd==2){
        if (j1!=j3 | j2!=j4){
          prob21<-glw%*%as.vector(fden2.12)
          prob22<-glw%*%as.vector(fden2.34)
          -prob21*prob22
        } else {
          prob2<-glw%*%as.vector(fden2.12)
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
Xi2.bifactor<-function(prob1,fden2,ngrp,grpsize,glw)
{
  dims<-dim(fden2)
  d=dims[4]
  K=dims[3]
  d1=d-1
  rows<-NULL
  for(i1 in 1:d)
  { for(j1 in 2:K)
  {
    res<-NULL
    for(i2 in 1:d)
    { for(j2 in 2:K)
    { res<-c(res,cov2.bifactor(i1,j1,i2,j2,prob1,fden2,ngrp,grpsize,glw))
    }
    }
    for(i2 in 1:d1)
    { for(i3 in (i2+1):d)
    { for(j2 in 2:K)
    { for(j3 in 2:K)
    { res<-c(res,cov3.bifactor(i1,j1,i2,j2,i3,j3,prob1,fden2,ngrp,grpsize,glw))
    }}}}
    rows<-rbind(rows,res)
  }}

  for(i1 in 1:d1)
  { for(i2 in (i1+1):d)
  { for(j1 in 2:K)
  { for(j2 in 2:K)
  { res<-NULL
  for(i3 in 1:d)
  { for(j3 in 2:K)
  { res<-c(res,cov3.bifactor(i3,j3,i1,j1,i2,j2,prob1,fden2,ngrp,grpsize,glw))
  }}
  for(i3 in 1:d1)
  { for(i4 in (i3+1):d)
  { for(j3 in 2:K)
  { for(j4 in 2:K)
  { res<-c(res,cov4.bifactor(i1,j1,i2,j2,i3,j3,i4,j4,fden2,ngrp,grpsize,glw))
  }}}}
  rows<-rbind(rows,res)}}}}
  rows
}
################################################################################




# purpose:
# to calculate all the derivatives of the univariate probabilities wrt to
# (a) the cutpoints  and (b) the copula parameters for the bifactor model
# input:
# dnorma: \phi[\Phi^{-1}(cutp)]
# output: the matrix with derivatives of all univariate
# probabilities wrt to the cutpoints and copula parameters
all.der.uprob.bifactor<-function(dnorma,SpC)
{
  #number of columns in dnorma, this is the number of variables
  d=ncol(dnorma)

  #Number of parameters (assuming all have 1-param copula)
  ddot=2*d

  #numer of rows of dnorma, this is the probs for each categ excluding categ 0.
  r=nrow(dnorma)

  # negative values for dnorma, this is used for the derivatives of probs with respect to
  # the cut-points.
  neg.dnorma=-dnorma

  #number of rows for the univ. prob. with respec to the cutpoints and copula parameters
  rows.result= d*r
  # creating the matrix for derivatives of uni. prob,
  result = matrix( 0 , ncol =  d*r+ddot , nrow =  d*r)

  # The diagnoal of the matrix will always be negative and will include
  # all elements from dnorma:
  diag(result)=neg.dnorma

  # Subdiagnoal will include the other cutpoints excluding
  # the first one (I don't mean categ=0, but I mean categ=1)
  # extra element is added as zero, because there is an empty
  # cell between each variable and this is done by adding 0.
  ndnorma=dnorma[-1,]
  subdiagonal=as.numeric(rbind(ndnorma,0))[-rows.result]

  # letting the values be as subdiagonal in original matrix
  diag(result[-nrow(result),-1])=subdiagonal

  # Return the matrix with derivatives of univaraite margins
  # For two-factor copula model with BVN copulas (if SpC is TRUE) 
  # we exclude the last columnn that correspond to the nth 
  # copula that is set to independence.
  if ( is.null(SpC) ) {
    output = result
  }else {
    output = result[ , -(d*r+d+SpC)] 
  }
  return(output)
}





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
der.bivprob.cutp.bifactor<-function(i1,i2,y1,y2,i,cutp,dnorma,cden2,fden2,ngrp,grpsize,glw){

  i1i2=c(i1,i2)#pair of variables
  grplst=list()
  ind=1
  onegrp=rep(FALSE,ngrp)
  #loop to create a list of the index of variables in each group.
  for(g in 1:ngrp){
    grplst[[g]]=ind:cumsum(grpsize)[g]
    ind=ind+grpsize[g]
    if(i1i2[1]%in%grplst[[g]]&&i1i2[2]%in%grplst[[g]]){
      onegrp[g]=TRUE}
  }

  if(sum(onegrp)==1){

    if(i1==i & (y1-1)==cutp){
      inner<-dnorma[cutp+1,i]*cden2[,,cutp+1,i]*fden2[,,y2,i2]
      glw%*%as.vector(glw %*%inner)
    } else {
      if(i2==i & (y2-1)==cutp){
        inner<-dnorma[cutp+1,i]*cden2[,,cutp+1,i]*fden2[,,y1,i1]
        glw%*%as.vector(glw %*%inner)
      } else {
        if(i1==i & cutp==(y1-2)){
          inner<--dnorma[cutp+1,i]*cden2[,,cutp+1,i]*fden2[,,y2,i2]
          glw%*%as.vector(glw%*% inner)
        } else {
          if(i2==i & cutp==(y2-2)){
            inner<--dnorma[cutp+1,i]*cden2[,,cutp+1,i]*fden2[,,y1,i1]
            glw%*%as.vector(glw %*%inner)
          } else {
            0
          }
        }
      }
    }
  }else{

    if(i1==i & (y1-1)==cutp){
      innerg1<- glw%*%(dnorma[cutp+1,i]*cden2[,,cutp+1,i])
      innerg2<- glw%*%fden2[,,y2,i2]
      glw%*%as.vector(innerg1*innerg2)
    } else {
      if(i2==i & (y2-1)==cutp){
        innerg1<- glw%*%(dnorma[cutp+1,i]*cden2[,,cutp+1,i])
        innerg2<- glw%*%fden2[,,y1,i1]
        glw%*%as.vector(innerg1*innerg2)
      } else {
        if(i1==i & cutp==(y1-2)){
          innerg1<- glw%*%(-dnorma[cutp+1,i]*cden2[,,cutp+1,i])
          innerg2<- glw%*%fden2[,,y2,i2]
          glw%*%as.vector(innerg1*innerg2)
        } else {
          if(i2==i & cutp==(y2-2)){
            innerg1<- glw%*%(-dnorma[cutp+1,i]*cden2[,,cutp+1,i])
            innerg2<- glw%*%fden2[,,y1,i1]
            glw%*%as.vector(innerg1*innerg2)
          } else {
            0
          }
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
all.der.bivprob.cutp.bifactor<-function(dnorma,cden2,fden2,ngrp,grpsize,glw)
{ dims<-dim(fden2)
d=dims[4]
K=dims[3]
K2=K-2
cols<-NULL
for(i in 1:d)
{ for(cutp in 0:K2)
{ res<-NULL
for(i1 in 1:(d-1))
{ for(i2 in (i1+1):d)
{ for(y1 in 2:K)
{ for(y2 in 2:K)
  res<-c(res,der.bivprob.cutp.bifactor(i1,i2,y1,y2,i,cutp,dnorma,cden2,fden2,ngrp,grpsize,glw))}}}
cols<-cbind(cols,res)
}
}
cols
}



# purpose:
# to calculate the derivative of \Pr_{i1 i2, y1 y2} wrt group-specific theta
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
# the derivative of \Pr_{i1 i2, y1 y2} wrt thXg_j
der.bivprob.thxg.bifactor<-function(i1,i2,y1,y2,i,fden2,fdotden2,ngrp,grpsize,glw){
  i1i2=c(i1,i2)#pair of variables
  grplst=list()
  ind=1
  onegrp=rep(FALSE,ngrp)
  #loop to create a list of the index of variables in each group.
  for(g in 1:ngrp){
    grplst[[g]]=ind:cumsum(grpsize)[g]
    ind=ind+grpsize[g]
    if(i1i2[1]%in%grplst[[g]]&&i1i2[2]%in%grplst[[g]]){
      onegrp[g]=TRUE}
  }
  if(sum(onegrp)==1){
    if(i1==i){
      inner<-fden2[,,y2,i2]*fdotden2[,,y1,i]
      glw%*%as.vector(glw%*%inner)
    } else {
      if(i2==i){
        inner<-fden2[,,y1,i1]*fdotden2[,,y2,i]
        glw%*%as.vector(glw%*%inner)
      } else {
        0
      }
    }
  }else{
    if(i1==i){
      innerg1<- glw%*%fden2[,,y2,i2]
      innerg2<- glw%*%fdotden2[,,y1,i]
      glw%*%as.vector(innerg1*innerg2)
    } else {
      if(i2==i){
        innerg1<- glw%*%fden2[,,y1,i1]
        innerg2<- glw%*%fdotden2[,,y2,i]
        glw%*%as.vector(innerg1*innerg2)
      } else {
        0
      }
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
all.der.bivprob.thxg.bifactor<-function(fden2,fdotden2,ngrp,grpsize,glw,SpC)
{ dims<-dim(fden2)
d=dims[4]
K=dims[3]
cols<-NULL

#adjusting the number of columns if SpC is TRUE or FALSE.
if (is.null(SpC)){
  dF2delta = 1:d 
} else {
  dF2delta = (1:d)[-SpC]
}

for(i in dF2delta)
{res<-NULL
for(i1 in 1:(d-1))
{ for(i2 in (i1+1):d)
{ for(y1 in 2:K)
{ for(y2 in 2:K)
  res<-c(res,der.bivprob.thxg.bifactor(i1,i2,y1,y2,i,fden2,fdotden2,ngrp,grpsize,glw))}}}
cols<-cbind(cols,res)
}
cols
}

###############################################################################

# purpose:
# to calculate the derivative of \Pr_{i1 i2, y1 y2} wrt thx0
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
# the derivative of \Pr_{i1 i2, y1 y2} wrt thx0
der.bivprob.thx0.bifactor<-function(i1,i2,y1,y2,i,fden2,fbarden2,ngrp,grpsize,glw){
  i1i2=c(i1,i2)#pair of variables
  grplst=list()
  ind=1
  onegrp=rep(FALSE,ngrp)
  #loop to create a list of the index of variables in each group.
  for(g in 1:ngrp){
    grplst[[g]]=ind:cumsum(grpsize)[g]
    ind=ind+grpsize[g]
    if(i1i2[1]%in%grplst[[g]]&&i1i2[2]%in%grplst[[g]]){
      onegrp[g]=TRUE}
  }
  if(sum(onegrp)==1){
    if(i1==i){
      inner<-fden2[,,y2,i2]*fbarden2[,,y1,i]
      glw%*%as.vector(glw%*%inner)
    } else {
      if(i2==i){
        inner<-fden2[,,y1,i1]*fbarden2[,,y2,i]
        glw%*%as.vector(glw%*%inner)
      } else {
        0
      }
    }
  }else{
    if(i1==i){
      innerg1<-glw%*%fden2[,,y2,i2]
      innerg2<-glw%*%fbarden2[,,y1,i]
      glw%*%as.vector(innerg1*innerg2)
    } else {
      if(i2==i){
        innerg1<-glw%*%fden2[,,y1,i1]
        innerg2<-glw%*%fbarden2[,,y2,i]
        glw%*%as.vector(innerg1*innerg2)
      } else {
        0
      }
    }
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
all.der.bivprob.thx0.bifactor<-function(fden2,fbarden2,ngrp,grpsize,glw)
{ dims<-dim(fden2)
d=dims[4]
K=dims[3]
cols<-NULL
for(i in 1:d)
{res<-NULL
for(i1 in 1:(d-1))
{ for(i2 in (i1+1):d)
{ for(y1 in 2:K)
{ for(y2 in 2:K)
  res<-c(res,der.bivprob.thx0.bifactor(i1,i2,y1,y2,i,fden2,fbarden2,ngrp,grpsize,glw))}}}
cols<-cbind(cols,res)
}
cols
}

################################################################################
# purpose:
# to calculate  the Delta_2 matrix for the bifactor model
# input:
# dnorma: \phi[\Phi^{-1}(cutp)]
# cden2: cdens12 (see Appendix)
# fden2: fden2 (see Appendix)
# fdotden2: fdotden2 (see Appendix)
# fbarden2: fbarden2 (see Appendix)
# output:
# the Delta_2 matrix
Delta2.bifactor<-function(dnorma,cden2,fden2,fdotden2,fbarden2,ngrp,grpsize,glw,SpC)
{ #the derivatives of the univariate probabilities
  res1<-all.der.uprob.bifactor(dnorma,SpC)
  #the derivatives of the bivariate probabilities wrt cutpoints
  res2<-all.der.bivprob.cutp.bifactor(dnorma,cden2,fden2,ngrp,grpsize,glw)
  #the derivatives of the bivariate probabilities wrt copula parmaeters
  res3<-all.der.bivprob.thx0.bifactor(fden2,fbarden2,ngrp,grpsize,glw)
  res4<-all.der.bivprob.thxg.bifactor(fden2,fdotden2,ngrp,grpsize,glw,SpC)
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
M2.bifactor<-function(dat,dnorma,prob1,cden2,fden2,fdotden2,fbarden2,ngrp,grpsize,glw,SpC,iprint=FALSE)
{
  dims<-dim(fden2)
  d=dims[4]
  K=dims[3]
  n=nrow(dat)
  K1=K-1
  K1d=1:(K1*d)
  bpairs<-bivpairs(d)
  bm<-nrow(bpairs)
  V<-Xi2.bifactor(prob1,fden2,ngrp,grpsize,glw)
  if(iprint)
  { cat("\ndim(V): dimension of Xi_2 matrix\n")
    print(dim(V))
  }
  delta<-Delta2.bifactor(dnorma,cden2,fden2,fdotden2,fbarden2,ngrp,grpsize,glw,SpC)
  if(iprint)
  { cat("\ndim(delta): dimension of Delta_2 matrix\n")
    print(dim(delta))
  }
  oc.delta<-orthcomp(delta)
  inv<-solve(t(oc.delta)%*%V%*%oc.delta)
  C2<-oc.delta%*%inv%*%t(oc.delta)
  pi.r<-pihat(dat,ngrp,grpsize,prob1,fden2,glw)
  p.r<-pobs(dat)
  stat<-nrow(dat)*t(p.r-pi.r)%*%C2%*%(p.r-pi.r)
  dof<-nrow(delta)-ncol(delta)
  pvalue<-pchisq(stat,dof, lower.tail = FALSE)
  dg<-round(n*apply(matrix(abs(p.r[-K1d]-pi.r[-K1d]),bm,K1*K1,byrow=T),1,max))
  list("M_2"=stat,"df"=dof,"p-value"=pvalue,"diagnostic"=dg)
}

