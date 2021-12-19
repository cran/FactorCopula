
store.nested<-function(cutp,thX0Xg,thXgYj,
                       pcondY,pconddotY,dcopY,
                       dcopX,ddotcopX,nq,ngrp,glw,gln){
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
    pcondY.j=pcondY[[j]]
    pconddotY.j=pconddotY[[j]]
    dcopY.j=dcopY[[j]]
    for(k in 1:(K-1)){
      condcdf[,k,j]<-pcondY.j(cutp[k,j],gln,thXgYj[j])
      conddotcdf[,k,j]<-pconddotY.j(cutp[k,j],gln,thXgYj[j])
      cden[,k,j]<-dcopY.j(cutp[k,j],gln,thXgYj[j])
    }
  }

  fden<-array(NA,dim=c(nq,nq,ngrp))
  fdotden<-array(NA,dim=c(nq,nq,ngrp))
  for(g in 1:ngrp){
    dcopX.g=dcopX[[g]]
    ddotcopX.g=ddotcopX[[g]]
    for(q in 1:nq){
      fden[,q,g]=dcopX.g(gln[q],gln,thX0Xg[g])
      fdotden[,q,g]=ddotcopX.g(gln[q],gln,thX0Xg[g])
    }
  }

  pmf<-array(NA,dim=c(nq,K,d))
  pmfdot<-array(NA,dim=c(nq,K,d))

  arr0<-array(0,dim=c(nq,1,d))
  arr1<-array(1,dim=c(nq,1,d))
  condcdf<-abind(arr0,condcdf,arr1,along=2)
  conddotcdf<-abind(arr0,conddotcdf,arr0,along=2)
  for(j in 1:d)
  { for(k in 1:K)
  { pmf[,k,j]<-condcdf[,k+1,j]-condcdf[,k,j]
  pmfdot[,k,j]<-conddotcdf[,k+1,j]-conddotcdf[,k,j]
  }
  }
  list(prob1=prob1,dnorma=dnorma,pmf=pmf,pmfdot=pmfdot,
       cden=cden,fden=fden,fdotden=fdotden)
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
pihat.nested<-function(dat,ngrp,grpsize,prob1,pmf,fden,glw,gln)
{

  d=ncol(dat)
  n=nrow(dat)
  K<-length(unique(dat[,1]))

  #uninvariate probabilities
  pihat1=prob1[-1,]
  bpairs=bivpairs(d)

  res2=NULL
  for(i in 1:nrow(bpairs)){
    items=bpairs[i,]

    grplst=list()
    ind=1
    item1=rep(FALSE,ngrp);item2=rep(FALSE,ngrp)
    for(g in 1:ngrp){#ngrp of all data
      grplst[[g]]=ind:cumsum(grpsize)[g]
      ind=ind+grpsize[g]
      if(items[1]%in%grplst[[g]]){item1[g]=1}else{item1[g]=FALSE}
      if(items[2]%in%grplst[[g]]){item2[g]=1}else{item2[g]=FALSE}
    }
    #How many groups/testlets in pairs.cov3?
    testlets=length(unique(c(which(item1==1),which(item2==1))))
    #place the items according to their groups.
    grpmat=rbind(item1,item2)
    #The possible combinations of items
    itemgrp=c(which(grpmat[1,]>0),which(grpmat[2,]>0))

    if(testlets==1){
      res1<-NULL
      for(j1 in 2:K){
        for(j2 in 2:K){
          pmfprod=(pmf[,j1,items[1]]*pmf[,j2,items[2]])
          dfactor.v0gy=glw %*%(pmfprod*fden[,,itemgrp[1]])
          res1<-c(res1,glw %*% as.vector(dfactor.v0gy))
        }
      }
      res1
    }
    else{
      res1<-NULL;
      for(j1 in 2:K){
        for(j2 in 2:K){
          pmf1=pmf[,j1,items[1]]
          pmf2=pmf[,j2,items[2]]
          dfactor.v0g1y=glw%*%(pmf1*fden[,,itemgrp[1]])
          dfactor.v0g2y=glw%*%(pmf2*fden[,,itemgrp[2]])
          res1<-c(res1,glw%*%(as.vector(dfactor.v0g1y)*as.vector(dfactor.v0g2y)))
        }
      }
      res1
    }
    res2=c(res2,res1)
  }
  return(c(pihat1,res2))
}

#============================= Covariance

cov2.nested<-function(i1,j1,i2,j2,prob1,pmf,fden,ngrp,grpsize,glw,gln){
  items=c(i1,i2)
  grplst=list()
  ind=1
  item1=rep(FALSE,ngrp);item2=rep(FALSE,ngrp)
  for(g in 1:ngrp){#ngrp of all data
    grplst[[g]]=ind:cumsum(grpsize)[g]
    ind=ind+grpsize[g]
    if(items[1]%in%grplst[[g]]){item1[g]=1}else{item1[g]=FALSE}
    if(items[2]%in%grplst[[g]]){item2[g]=1}else{item2[g]=FALSE}
  }
  #How many groups/testlets in pairs.cov3?
  testlets=length(unique(c(which(item1==1),which(item2==1))))
  #place the items according to their groups.
  grpmat=rbind(item1,item2)
  #The possible combinations of items
  itemgrp=c(which(grpmat[1,]>0),which(grpmat[2,]>0))

  #univariate probabilities of the pair.
  prob1i1=prob1[j1,i1]; prob1i2=prob1[j2,i2]
  if(testlets==1){
    #This is if the variables are in the same group.
    dfactor.v0gy=fden[,,itemgrp[1]]
    temp1=pmf[,j1,i1];temp2=pmf[,j2,i2]
    bifden12=glw%*%((temp1*temp2)*dfactor.v0gy)
  }else{
    #This is if the variables are in different groups.
    dfactor.v0g1y=fden[,,itemgrp[1]]
    dfactor.v0g2y=fden[,,itemgrp[2]]
    temp1=glw%*%(pmf[,j1,i1]*dfactor.v0g1y);
    temp2=glw%*%(pmf[,j2,i2]*dfactor.v0g2y)
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


cov3.nested<-function(i1,j1,i2,j2,i3,j3,prob1,pmf,fden,ngrp,grpsize,glw,gln){
  pairs.cov3=c(i1,i2,i3)#trivariate variables.
  ngrpcov3=3#it cannot be more than 3 groups.
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
  grpmat=pairs.cov3*rbind(item1,item2,item3)
  #exclude the extra unused groups.
  ngrpmat=cbind(grpmat[ , colSums(grpmat != 0) != 0] )

  #What the categories used wrt to each item and for each group?
  y=cbind(c(j1,j2,j3)*rbind(y1,y2,y3))
  #exclude the extra unused groups.
  ny=cbind(y[ , colSums(y != 0) != 0])

  #The possible combinations of items
  itemgrp12=c(which(grpmat[1,]>0),which(grpmat[2,]>0))
  itemgrp13=c(which(grpmat[1,]>0),which(grpmat[3,]>0))
  itemgrp23=c(which(grpmat[2,]>0),which(grpmat[3,]>0))
  itemgrp123=c(which(grpmat[1,]>0),which(grpmat[2,]>0),which(grpmat[3,]>0))

  #compute corresponding mass functinos.
  if(length(unique(itemgrp12))==1){
    dfactor.v0gy=fden[,,itemgrp12[1]]
    fden2.12=glw%*%((pmf[,j1,i1]*pmf[,j2,i2])*dfactor.v0gy)
  }else{
    dfactor.v0g1y=fden[,,itemgrp12[1]]
    dfactor.v0g2y=fden[,,itemgrp12[2]]
    fden2.12=(glw%*%(dfactor.v0g1y*pmf[,j1,i1]))*( glw%*%(dfactor.v0g2y*pmf[,j2,i2]))
    }

  if(length(unique(itemgrp13))==1){
    dfactor.v0gy=fden[,,itemgrp13[1]]
    fden2.13=glw%*%(dfactor.v0gy*(pmf[,j1,i1]*pmf[,j3,i3]))
  }else{
    dfactor.v0g1y=fden[,,itemgrp13[1]]
    dfactor.v0g2y=fden[,,itemgrp13[2]]
    fden2.13=(glw%*%(dfactor.v0g1y*pmf[,j1,i1]))*(glw%*%(dfactor.v0g2y*pmf[,j3,i3]))
  }

  if(length(unique(itemgrp23))==1){
    dfactor.v0gy=fden[,,itemgrp23[1]]
    fden2.23=glw%*%(dfactor.v0gy*(pmf[,j2,i2]*pmf[,j3,i3]))
  }else{
    dfactor.v0g1y=fden[,,itemgrp23[1]]
    dfactor.v0g2y=fden[,,itemgrp23[2]]
    fden2.23=(glw%*%(dfactor.v0g1y*pmf[,j2,i2]))*(glw%*%(dfactor.v0g2y*pmf[,j3,i3]))
  }

  if(length(unique(itemgrp123))==1){
    #This is if items 1, 2, and 3 are in the same group
    dfactor.v0gy=fden[,,itemgrp123[1]]
    fden2.123=glw%*%(dfactor.v0gy*(pmf[,j1,i1]*pmf[,j2,i2]*pmf[,j3,i3]))
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

    #identifying the groups for each item
    grp23=itemgrp123[duplicated(itemgrp123)]
    grp1=itemgrp123[-which(itemgrp123==grp23)]

    dfactor.v0gy23=fden[,,grp23]
    dfactor.v0gy1=fden[,,grp1]

    fden2.123=((glw%*%(dfactor.v0gy1*pmf[,yj1,ni1]))*
                 (glw%*%(dfactor.v0gy23*(pmf[,yj2,ni2]*pmf[,yj3,ni3]))))
  }
  if(length(unique(itemgrp123))==3){
    dfactor.v0gy1=fden[,,itemgrp123[1]]
    dfactor.v0gy2=fden[,,itemgrp123[2]]
    dfactor.v0gy3=fden[,,itemgrp123[3]]

    fden2.123=((glw%*%(dfactor.v0gy1*pmf[,j1,i1]))*
                 (glw%*%(dfactor.v0gy2*pmf[,j2,i2]))*
                 (glw%*%(dfactor.v0gy3*pmf[,j3,i3])))
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
cov4.nested<-function(i1,j1,i2,j2,i3,j3,i4,j4,pmf,fden,ngrp,grpsize,glw,gln){
  nd=ndistinct(i1,i2,i3,i4)
  pairs.cov4=c(i1,i2,i3,i4)
  ngrpcov4=ngrp#it cannot be more than 3 groups
  grplst=list()
  ind=1
  item1=rep(FALSE,ngrp);item2=rep(FALSE,ngrp)
  item3=rep(FALSE,ngrp);item4=rep(FALSE,ngrp)
  y1=rep(FALSE,ngrp);y2=rep(FALSE,ngrp)
  y3=rep(FALSE,ngrp);y4=rep(FALSE,ngrp)

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

  ymat=cbind(c(j1,j2,j3,j4)*rbind(y1,y2,y3,y4))
  #exclude the extra unused groups.
  ny=cbind(ymat[ , colSums(ymat != 0) != 0])

  itemgrp12=c(which(grpmat[1,]>0),which(grpmat[2,]>0))
  itemgrp13=c(which(grpmat[1,]>0),which(grpmat[3,]>0))
  itemgrp14=c(which(grpmat[1,]>0),which(grpmat[4,]>0))
  itemgrp23=c(which(grpmat[2,]>0),which(grpmat[3,]>0))
  itemgrp24=c(which(grpmat[2,]>0),which(grpmat[4,]>0))
  itemgrp34=c(which(grpmat[3,]>0),which(grpmat[4,]>0))
  itemgrp123=c(which(grpmat[1,]>0),which(grpmat[2,]>0),which(grpmat[3,]>0))
  itemgrp124=c(which(grpmat[1,]>0),which(grpmat[2,]>0),which(grpmat[4,]>0))
  itemgrp1234=c(which(grpmat[1,]>0),which(grpmat[2,]>0),
                which(grpmat[3,]>0),which(grpmat[4,]>0))

  if(length(unique(itemgrp12))==1){
    dfactor.v0gy=fden[,,itemgrp12[1]]
    fden2.12=glw%*%(dfactor.v0gy*(pmf[,j1,i1]*pmf[,j2,i2]))
  }else{
    dfactor.v0gy1=fden[,,itemgrp12[1]]
    dfactor.v0gy2=fden[,,itemgrp12[2]]
    fden2.12=(glw%*%(dfactor.v0gy1*pmf[,j1,i1]))*(glw%*%(dfactor.v0gy2*pmf[,j2,i2]))
    }

  if(length(unique(itemgrp13))==1){
    dfactor.v0gy=fden[,,itemgrp13[1]]
    fden2.13=glw%*%(dfactor.v0gy*(pmf[,j1,i1]*pmf[,j3,i3]))
  }else{
    dfactor.v0gy1=fden[,,itemgrp13[1]]
    dfactor.v0gy2=fden[,,itemgrp13[2]]
    fden2.13=(glw%*%(dfactor.v0gy1*pmf[,j1,i1]))*(glw%*%(dfactor.v0gy2*pmf[,j3,i3]))
  }

  if(length(unique(itemgrp14))==1){
    dfactor.v0gy=fden[,,itemgrp14[1]]
    fden2.14=glw%*%(dfactor.v0gy*(pmf[,j1,i1]*pmf[,j4,i4]))
  }else{
    dfactor.v0gy1=fden[,,itemgrp14[1]]
    dfactor.v0gy2=fden[,,itemgrp14[2]]
    fden2.14=(glw%*%(dfactor.v0gy1*pmf[,j1,i1]))*(glw%*%(dfactor.v0gy2*pmf[,j4,i4]))
  }

  if(length(unique(itemgrp23))==1){
    dfactor.v0gy=fden[,,itemgrp23[1]]
    fden2.23=glw%*%(dfactor.v0gy*(pmf[,j2,i2]*pmf[,j3,i3]))
  }else{
    dfactor.v0gy1=fden[,,itemgrp23[1]]
    dfactor.v0gy2=fden[,,itemgrp23[2]]
    fden2.23=(glw%*%(dfactor.v0gy1*pmf[,j2,i2]))*(glw%*%(dfactor.v0gy2*pmf[,j3,i3]))
  }

  if(length(unique(itemgrp24))==1){
    dfactor.v0gy=fden[,,itemgrp24[1]]
    fden2.24=glw%*%(dfactor.v0gy*(pmf[,j2,i2]*pmf[,j4,i4]))
  }else{
    dfactor.v0gy1=fden[,,itemgrp24[1]]
    dfactor.v0gy2=fden[,,itemgrp24[2]]
    fden2.24=(glw%*%(dfactor.v0gy1*pmf[,j2,i2]))*(glw%*%(dfactor.v0gy2*pmf[,j4,i4]))
  }

  if(length(unique(itemgrp34))==1){
    dfactor.v0gy=fden[,,itemgrp34[1]]
    fden2.34=glw%*%(dfactor.v0gy*(pmf[,j3,i3]*pmf[,j4,i4]))
  }else{
    dfactor.v0gy1=fden[,,itemgrp34[1]]
    dfactor.v0gy2=fden[,,itemgrp34[2]]
    fden2.34=(glw%*%(dfactor.v0gy1*pmf[,j3,i3]))*(glw%*%(dfactor.v0gy2*pmf[,j4,i4]))
  }

  if(length(unique(itemgrp123))==1){
    dfactor.v0gy=fden[,,itemgrp123[1]]
    fden2.123=glw%*%(dfactor.v0gy*(pmf[,j1,i1]*pmf[,j2,i2]*pmf[,j3,i3]))

  }else if(length(unique(itemgrp123))==2){

    tripairs=list()
    jlst=list()
    for(g in 1:ngrp){
      tripairs[[g]]=as.numeric(c(grpmat[c(1,2,3),g][grpmat[c(1,2,3),g]!=0]))
      jlst[[g]]=as.numeric(c(ymat[c(1,2,3),g][ymat[c(1,2,3),g]!=0]))
    }

    jlst=jlst[lengths(jlst)>0 ];
    tripairs=tripairs[lengths(tripairs)>0 ];
    prodg1=1;prodg2=1;prodg3=1
    for(i in 1:length(tripairs[[1]])){
      prodg1=prodg1*pmf[,jlst[[1]][i],tripairs[[1]][i]]
    }
    for(i in 1:length(tripairs[[2]])){
      prodg2=prodg2*pmf[,jlst[[2]][i],tripairs[[2]][i]]
    }

    #identifying the groups for each item
    ngrp.123=which(colSums(grpmat)>0)
    dfactor.v0g1y=fden[,,ngrp.123[1]]
    dfactor.v0g2y=fden[,,ngrp.123[2]]

    fden2.123=(glw%*%(dfactor.v0g1y*prodg1))*(glw%*%(dfactor.v0g2y*prodg2))


  } else if(length(unique(itemgrp123))==3){
    dfactor.v0gy1=fden[,,itemgrp123[1]]
    dfactor.v0gy2=fden[,,itemgrp123[2]]
    dfactor.v0gy3=fden[,,itemgrp123[3]]

    fden2.123=((glw%*%(dfactor.v0gy1*pmf[,j1,i1]))*
                 (glw%*%(dfactor.v0gy2*pmf[,j2,i2]))*
                 (glw%*%(dfactor.v0gy3*pmf[,j3,i3])))
  }
  if(length(unique(itemgrp124))==1){
    dfactor.v0gy=fden[,,itemgrp124[1]]
    fden2.124=glw%*%(dfactor.v0gy*(pmf[,j1,i1]*pmf[,j2,i2]*pmf[,j4,i4]))

  }else if(length(unique(itemgrp124))==2){
    tripairs=list()
    jlst=list()
    for(g in 1:ngrp){
      tripairs[[g]]=as.numeric(c(grpmat[c(1,2,4),g][grpmat[c(1,2,4),g]!=0]))
      jlst[[g]]=as.numeric(c(ymat[c(1,2,4),g][ymat[c(1,2,4),g]!=0]))
    }

    jlst=jlst[lengths(jlst)>0 ];
    tripairs=tripairs[lengths(tripairs)>0 ];
    prodg1=1;prodg2=1;prodg3=1
    for(i in 1:length(tripairs[[1]])){
      prodg1=prodg1*pmf[,jlst[[1]][i],tripairs[[1]][i]]
    }
    for(i in 1:length(tripairs[[2]])){
      prodg2=prodg2*pmf[,jlst[[2]][i],tripairs[[2]][i]]
    }

    #identifying the groups for each item
    ngrp.124=which(colSums(grpmat)>0)
    dfactor.v0g1y=fden[,,ngrp.124[1]]
    dfactor.v0g2y=fden[,,ngrp.124[2]]

    fden2.124=(glw%*%(dfactor.v0g1y*prodg1))*(glw%*%(dfactor.v0g2y*prodg2))


  }  else if(length(unique(itemgrp124))==3){
    dfactor.v0gy1=fden[,,itemgrp124[1]]
    dfactor.v0gy2=fden[,,itemgrp124[2]]
    dfactor.v0gy4=fden[,,itemgrp124[3]]

    fden2.124=((glw%*%(dfactor.v0gy1*pmf[,j1,i1]))*
                 (glw%*%(dfactor.v0gy2*pmf[,j2,i2]))*
                 (glw%*%(dfactor.v0gy4*pmf[,j4,i4])))
  }

  if(length(unique(itemgrp1234))==4){
    dfactor.v0gy1=fden[,,itemgrp1234[1]]
    dfactor.v0gy2=fden[,,itemgrp1234[2]]
    dfactor.v0gy3=fden[,,itemgrp1234[3]]
    dfactor.v0gy4=fden[,,itemgrp1234[4]]

    fden2.1234=((glw%*%(dfactor.v0gy1*pmf[,j1,i1]))*
                  (glw%*%(dfactor.v0gy2*pmf[,j2,i2]))*
                  (glw%*%(dfactor.v0gy3*pmf[,j3,i3]))*
                  (glw%*%(dfactor.v0gy4*pmf[,j4,i4])))

  }else if(length(unique(itemgrp1234))==1){
    dfactor.v0gy=fden[,,itemgrp1234[1]]

    fden2.1234=glw%*%(dfactor.v0gy*(pmf[,j1,i1]*pmf[,j2,i2]*
                                 pmf[,j3,i3] * pmf[,j4,i4]))

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
      prodg1=prodg1*pmf[,jlst[[1]][i],quadpairs[[1]][i]]
    }
    for(i in 1:length(quadpairs[[2]])){
      prodg2=prodg2*pmf[,jlst[[2]][i],quadpairs[[2]][i]]
    }
    for(i in 1:length(quadpairs[[3]])){
      prodg3=prodg3*pmf[,jlst[[3]][i],quadpairs[[3]][i]]
    }

    #identifying the groups for each item
    ngrp.1234=which(colSums(grpmat)>0)
    dfactor.v0g1y=fden[,,ngrp.1234[1]]
    dfactor.v0g2y=fden[,,ngrp.1234[2]]
    dfactor.v0g3y=fden[,,ngrp.1234[3]]

    fden2.1234=(glw%*%(dfactor.v0g1y*prodg1))*(glw%*%(dfactor.v0g2y*prodg2))*(glw%*%(dfactor.v0g3y*prodg3))

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
      prodg1=prodg1*pmf[,jlst[[1]][i],quadpairs[[1]][i]]
    }
    for(i in 1:length(quadpairs[[2]])){
      prodg2=prodg2*pmf[,jlst[[2]][i],quadpairs[[2]][i]]
    }

    #identifying the groups for each item
    ngrp.1234=which(colSums(grpmat)>0)
    dfactor.v0g1y=fden[,,ngrp.1234[1]]
    dfactor.v0g2y=fden[,,ngrp.1234[2]]

    fden2.1234=(glw%*%(dfactor.v0g1y*prodg1))*(glw%*%(dfactor.v0g2y*prodg2))
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
# to calculate the \Xi_2 matrix for the nested model
# input:
# prob1: the (K\times d) matrix with the univariate probabilities
# fden2: fden2 (see Appendix)
# output:
# the \Xi_2 matrix
Xi2.nested<-function(prob1,pmf,fden,ngrp,grpsize,glw,gln)
{ dims<-dim(pmf)
d=dims[3]
K=dims[2]
d1=d-1
rows<-NULL
for(i1 in 1:d)
{ for(j1 in 2:K)
{
  res<-NULL
  for(i2 in 1:d)
  { for(j2 in 2:K)
  { res<-c(res,cov2.nested(i1,j1,i2,j2,prob1,pmf,fden,ngrp,grpsize,glw,gln))
  }
  }
  for(i2 in 1:d1)
  { for(i3 in (i2+1):d)
  { for(j2 in 2:K)
  { for(j3 in 2:K)
  { res<-c(res,cov3.nested(i1,j1,i2,j2,i3,j3,prob1,pmf,fden,ngrp,grpsize,glw,gln))
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
{ res<-c(res,cov3.nested(i3,j3,i1,j1,i2,j2,prob1,pmf,fden,ngrp,grpsize,glw,gln))
}}
for(i3 in 1:d1)
{ for(i4 in (i3+1):d)
{ for(j3 in 2:K)
{ for(j4 in 2:K)
{ res<-c(res,cov4.nested(i1,j1,i2,j2,i3,j3,i4,j4,pmf,fden,ngrp,grpsize,glw,gln))
}}}}
rows<-rbind(rows,res)}}}}
rows
}


################################################################################

# purpose:
# to calculate all the derivatives of the univariate probabilities wrt to
# (a) the cutpoints  and (b) the copula parameters for the nested model
# input:
# dnorma: \phi[\Phi^{-1}(cutp)]
# output: the matrix with derivatives of all univariate
# probabilities wrt to the cutpoints and copula parameters
# The number of parameters is different here from the bifactor model.
all.der.uprob.nested<-function(dnorma,ngrp,SpC)
{
  #number of columns in dnorma, this is the number of variables
  d=ncol(dnorma)

  #Number of parameters (assuming all have 1-param copula)
  ddot=d+ngrp

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
    output = result[ , -(d*r+ngrp+SpC)] 
  }
  # Return the matrix with derivatives of univaraite margins
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
der.bivprob.cutp.nested<-function(i1,i2,y1,y2,i,cutp,dnorma,cden,pmf,fden,ngrp,grpsize,glw,gln){

  i1i2=c(i1,i2)#pair of variables
  grplst=list()
  ind=1
  item1=rep(FALSE,ngrp);item2=rep(FALSE,ngrp)
  for(g in 1:ngrp){#ngrp of all data
    grplst[[g]]=ind:cumsum(grpsize)[g]
    ind=ind+grpsize[g]
    if(i1i2[1]%in%grplst[[g]]){item1[g]=1}else{item1[g]=FALSE}
    if(i1i2[2]%in%grplst[[g]]){item2[g]=1}else{item2[g]=FALSE}
  }
  ngrp.i1i2=length(unique(c(which(item1==1),which(item2==1))))
  grpmat=cbind((i1i2)*rbind(item1,item2))

  if(sum(ngrp.i1i2)==1){
    grp=which(colSums(grpmat)>0)
    if(i1==i & (y1-1)==cutp){
      inner<- ((dnorma[cutp+1,i]*cden[,cutp+1,i]*pmf[,y2,i2])*fden[,,grp])
      glw%*%as.vector(glw%*%inner)
    } else {
      if(i2==i & (y2-1)==cutp){
        inner<- ((dnorma[cutp+1,i]*cden[,cutp+1,i]*pmf[,y1,i1])*fden[,,grp])
        glw%*%as.vector(glw%*%inner)
      } else {
        if(i1==i & cutp==(y1-2)){
          inner<- ((-dnorma[cutp+1,i]*cden[,cutp+1,i]*pmf[,y2,i2])*fden[,,grp])
          glw%*%as.vector(glw%*%inner)
        } else {
          if(i2==i & cutp==(y2-2)){
            inner<- ((-dnorma[cutp+1,i]*cden[,cutp+1,i]*pmf[,y1,i1])*fden[,,grp])
            glw%*%as.vector(glw%*%inner)
          } else {
            0
          }
        }
      }
    }
  }else{

    grp=which(colSums(grpmat)>0)
    grp1=grp[1]
    grp2=grp[2]

    if(i1==i & (y1-1)==cutp){
      innerg1<- glw%*%((dnorma[cutp+1,i]*cden[,cutp+1,i])*fden[,,grp1])
      innerg2<- glw%*%(pmf[,y2,i2]*fden[,,grp2])
      glw%*%as.vector(innerg1*innerg2)
    } else {
      if(i2==i & (y2-1)==cutp){
        innerg1<- glw%*%((dnorma[cutp+1,i]*cden[,cutp+1,i])*fden[,,grp2])
        innerg2<- glw%*%(pmf[,y1,i1]*fden[,,grp1])
        glw%*%as.vector(innerg1*innerg2)
      } else {
        if(i1==i & cutp==(y1-2)){
          innerg1<- glw%*%((-dnorma[cutp+1,i]*cden[,cutp+1,i])*fden[,,grp1])
          innerg2<- glw%*%(pmf[,y2,i2]*fden[,,grp2])
          glw%*%as.vector(innerg1*innerg2)
        } else {
          if(i2==i & cutp==(y2-2)){
            innerg1<- glw%*%((-dnorma[cutp+1,i]*cden[,cutp+1,i])*fden[,,grp2])
            innerg2<- glw%*%(pmf[,y1,i1]*fden[,,grp1])
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
all.der.bivprob.cutp.nested<-function(dnorma,cden,pmf,fden,ngrp,grpsize,glw,gln)
{
dims<-dim(pmf)
d=dims[3]
K=dims[2]
K2=K-2
cols<-NULL
for(i in 1:d)
{ for(cutp in 0:K2)
{ res<-NULL
for(i1 in 1:(d-1))
{ for(i2 in (i1+1):d)
{ for(y1 in 2:K)
{ for(y2 in 2:K)
  res<-c(res,der.bivprob.cutp.nested(i1,i2,y1,y2,i,cutp,dnorma,cden,pmf,fden,ngrp,grpsize,glw,gln))}}}
cols<-cbind(cols,res)
}
}
cols
}

# purpose:
# to calculate the derivative of \Pr_{i1 i2, y1 y2} wrt group-specific theta
# for the nested model
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
der.bivprob.thxg.nested<-function(i1,i2,y1,y2,i,pmf,fden,pmfdot,ngrp,grpsize,glw,gln){
  i1i2=c(i1,i2)#pair of variables
  grplst=list()
  ind=1
  item1=rep(FALSE,ngrp);item2=rep(FALSE,ngrp)
  for(g in 1:ngrp){#ngrp of all data
    grplst[[g]]=ind:cumsum(grpsize)[g]
    ind=ind+grpsize[g]
    if(i1i2[1]%in%grplst[[g]]){item1[g]=1}else{item1[g]=FALSE}
    if(i1i2[2]%in%grplst[[g]]){item2[g]=1}else{item2[g]=FALSE}
  }
  ngrp.i1i2=length(unique(c(which(item1==1),which(item2==1))))
  grpmat=cbind((i1i2)*rbind(item1,item2))


  if(sum(ngrp.i1i2)==1){
    grp=which(colSums(grpmat)>0)
    if(i1==i){
      inner<- ((pmf[,y2,i2]*pmfdot[,y1,i])*fden[,,grp])
      glw%*%as.vector(glw%*%inner)
    } else {
      if(i2==i){
        inner<- ((pmf[,y1,i1]*pmfdot[,y2,i])*fden[,,grp])
        glw%*%as.vector(glw%*%inner)
      } else {
        0
      }
    }
  }else{

    grp=which(colSums(grpmat)>0)
    grp1=grp[1]
    grp2=grp[2]

    if(i1==i){
      innerg1<- glw%*%(pmf[,y2,i2]*fden[,,grp2])
      innerg2<- glw%*%(pmfdot[,y1,i]*fden[,,grp1])
      glw%*%as.vector(innerg1*innerg2)
    } else {
      if(i2==i){
        innerg1<- glw%*%(pmf[,y1,i1]*fden[,,grp1])
        innerg2<- glw%*%(pmfdot[,y2,i]*fden[,,grp2])
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
# for the nested model
# input:
# fden: fden (see Appendix)
# fdotden: fdotden (see Appendix)
# output:
# a matrix with the derivatives of \Pr_{i1 i2, y1 y2} wrt delta_i  for all the
# bivariate probabilities and copula parameters at the 2nd factor
all.der.bivprob.thxg.nested<-function(pmf,fden,pmfdot,ngrp,grpsize,glw,gln,SpC){
  dims<-dim(pmf)
d=dims[3]
K=dims[2]
cols<-NULL

#adjusting the number of columns if SpC is TRUE or FALSE.
if (is.null(SpC)){
  dFgdelta = 1:d 
} else {
  dFgdelta = (1:d)[-SpC]
}

for(i in dFgdelta)
{res<-NULL
for(i1 in 1:(d-1))
{ for(i2 in (i1+1):d)
{ for(y1 in 2:K)
{ for(y2 in 2:K)
  res<-c(res,der.bivprob.thxg.nested(i1,i2,y1,y2,i,pmf,fden,pmfdot,ngrp,grpsize,glw,gln))}}}
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
der.bivprob.thx0.nested<-function(i1,i2,y1,y2,i,pmf,fden,fdotden,ngrp,grpsize,glw,gln){
  i1i2=c(i1,i2)#pair of variables
  grplst=list()
  ind=1
  item1=rep(FALSE,ngrp);item2=rep(FALSE,ngrp)
  for(g in 1:ngrp){#ngrp of all data
    grplst[[g]]=ind:cumsum(grpsize)[g]
    ind=ind+grpsize[g]
    if(i1i2[1]%in%grplst[[g]]){item1[g]=1}else{item1[g]=FALSE}
    if(i1i2[2]%in%grplst[[g]]){item2[g]=1}else{item2[g]=FALSE}
  }
  ngrp.i1i2=length(unique(c(which(item1==1),which(item2==1))))
  grpmat=cbind((i1i2)*rbind(item1,item2))

  if(sum(ngrp.i1i2)==1){
    grp=which(colSums(grpmat)>0)
    if(grp==i){
      inner<- ((pmf[,y2,i2]*pmf[,y1,i1])*fdotden[,,i])
      glw%*%as.vector(glw%*%inner)
    } else {
      0
    }
  }else{

    grp=which(colSums(grpmat)>0)
    grp1=grp[1]
    grp2=grp[2]

    if(grp1==i){
      innerg1<- glw%*%(pmf[,y2,i2]*fden[,,grp2])
      innerg2<- glw%*%(pmf[,y1,i1]*fdotden[,,i])
      glw%*%as.vector(innerg1*innerg2)
    }else if(grp2==i){
        innerg1<- glw%*%(pmf[,y2,i2]*fdotden[,,i])
        innerg2<- glw%*%(pmf[,y1,i1]*fden[,,grp1])
        glw%*%as.vector(innerg1*innerg2)
    } else {
      0
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
all.der.bivprob.thx0.nested<-function(pmf,fden,fdotden,ngrp,grpsize,glw,gln)
{
dims<-dim(pmf)
d=dims[3]
K=dims[2]
cols<-NULL
for(i in 1:ngrp)
{res<-NULL
for(i1 in 1:(d-1))
{ for(i2 in (i1+1):d)
{ for(y1 in 2:K)
{ for(y2 in 2:K)
  res<-c(res,der.bivprob.thx0.nested(i1,i2,y1,y2,i,pmf,fden,fdotden,ngrp,grpsize,glw,gln))}}}
cols<-cbind(cols,res)
}
cols
}

################################################################################
# purpose:
# to calculate  the Delta_2 matrix for the 1-factor model
# input:
# dnorma: \phi[\Phi^{-1}(cutp)]
# cden: cdens (see Appendix)
# fden: fden (see Appendix)
# fdotden: fdotden (see Appendix)
# output:
# the Delta_2 matrix
Delta2.nested<-function(dnorma,cden,pmf,pmfdot,fden,fdotden,ngrp,grpsize,glw,gln,SpC)
{ #the derivatives of the univariate probabilities
  res1<-all.der.uprob.nested(dnorma,ngrp,SpC)
  #the derivatives of the bivariate probabilities wrt cutpoints
  res2<-all.der.bivprob.cutp.nested(dnorma,cden,pmf,fden,ngrp,grpsize,glw,gln)
  #the derivatives of the bivariate probabilities wrt copula parmaeters
  res3<-all.der.bivprob.thxg.nested(pmf,fden,pmfdot,ngrp,grpsize,glw,gln,SpC)
  res4=all.der.bivprob.thx0.nested(pmf,fden,fdotden,ngrp,grpsize,glw,gln)
  rbind(res1,cbind(res2,res3,res4))
}

################################################################################
# R code for orthogonal complement (used by one of my students)
orthcomp <- function(x,tol=1.e-12)
{ s <- nrow(x)
q <- ncol(x)
if (s<=q) { return('error, matrix must have more columns than rows') }
x.svd <- svd(x,s)
if (sum(abs(x.svd$d)<tol)!=0)
{ return('error, matrix must full column rank') }
return(x.svd$u[,(q+1):s])
}

# purpose:
# to calculate  the M_2 for the 1-factor model
# input:
# dat: A data matrix where the number of rows corresponds to an
# individual's response and each column represents an item
# Number of ordinal categories for each item, coded as 0,...,(ncat-1).
# Currently supported are items that have the same number of categories.
# dnorma: \phi[\Phi^{-1}(cutp)]
# prob1: the (K\times d) matrix with the univariate probabilities
# fden: fden (see Appendix)
# fdotden: fdotden (see Appendix)
# cden: cdens (see Appendix)
# iprint: debug option
# output:
# a list with
# "M_2"=M2 stat
# "df"=degrees of freedom
# "p-value"=pvalue of M2
# "diagnostic"=biavriate diagnostics
M2.nested<-function(dat,dnorma,prob1,cden,pmf,pmfdot,fden,fdotden,
                    ngrp,grpsize,glw,gln,SpC,iprint=FALSE){
dims<-dim(pmf)
d=dims[3]
K=dims[2]
n=nrow(dat)
K1=K-1
K1d=1:(K1*d)
bpairs<-bivpairs(d)
bm<-nrow(bpairs)
V<-Xi2.nested(prob1,pmf,fden,ngrp,grpsize,glw,gln)
if(iprint)
{ cat("\ndim(V): dimension of Xi_2 matrix\n")
  print(dim(V))
}
delta<-Delta2.nested(dnorma,cden,pmf,pmfdot,fden,fdotden,ngrp,grpsize,glw,gln,SpC)
if(iprint)
{ cat("\ndim(delta): dimension of Delta_2 matrix\n")
  print(dim(delta))
}
oc.delta<-orthcomp(delta)
inv<-solve(t(oc.delta)%*%V%*%oc.delta)
C2<-oc.delta%*%inv%*%t(oc.delta)
pi.r<-pihat.nested(dat,ngrp,grpsize,prob1,pmf,fden,glw,gln)
p.r<-pobs(dat)
stat<-n*t(p.r-pi.r)%*%C2%*%(p.r-pi.r)
dof<-nrow(delta)-ncol(delta)
pvalue<-pchisq(stat,dof, lower.tail = FALSE)
dg<-round(n*apply(matrix(abs(p.r[-K1d]-pi.r[-K1d]),bm,K1*K1,byrow=T),1,max))
list("M_2"=stat,"df"=dof,"p-value"=pvalue,"diagnostic"=dg)
}



