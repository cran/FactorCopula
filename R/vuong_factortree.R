#M1: base model to be compared to: e.g. BVN
#M2: Model 2

vuong.1fv = function(tau.M1,cop.M1,A.M1, tau.M2, cop.M2,A.M2, 
                     ordinal, nq, param=F){
  #Gauss legendre points and weights
  gl<-gauss.quad.prob(nq)
  gln<-gl$nodes
  glw<-gl$weights
  d<-ncol(ordinal)
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  n<-nrow(ordinal)
  cutpoints<-unifcuts(ordinal)
  #--------------------------------------------------------
  # names of copulas for Model 1
  likelihood.vuong.den=lk_f1vine(tau.M1,ydat=ordinal,
                                 A=A.M1,cutp=cutpoints,nq,gln,glw,cop.M1,param=F)
  
  # likelihood in the numerator 
  likelihood.vuong.num=lk_f1vine(tau.M2,ydat=ordinal,
                                 A=A.M2,cutp=cutpoints,nq,gln,glw,cop.M2,param=F)
  
  
  D = log(likelihood.vuong.num/likelihood.vuong.den) #log-likelihood ratio
  avrg=mean(D) #Average of the D
  
  # Z score
  z = sqrt(n) * avrg/sd(D)
  
  #P-value for  two tailed
  pvalue <- 2 * pnorm(-abs(z)) #two tailed
  
  # The intervals to be added and subtracted to create the 95% CI.
  interval =  1.96 * (1/sqrt(n)) * sd(D) 
  
  # 95% confidence interval
  CI95=c( avrg - interval, avrg + interval )
  
  # Exporting result
  result <- data.frame(z, pvalue,  CI95[1], CI95[2])
  
  names(result) <- c("z", "p.value", "CI.left", "CI.right")
  return(result)
}


#M1: base model to be compared to: e.g. BVN
#M2: Model 2

vuong.2fv = function(tau.M1,cop.M1,A.M1, tau.M2, cop.M2,A.M2, 
                     ordinal, nq, param=F){
  #Gauss legendre points and weights
  gl<-gauss.quad.prob(nq)
  gln<-gl$nodes
  glw<-gl$weights
  d<-ncol(ordinal)
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  n<-nrow(ordinal)
  cutpoints<-unifcuts(ordinal)
  #--------------------------------------------------------
  # names of copulas for Model 1
  likelihood.vuong.den=lk_f2vine(tau.M1,ydat=ordinal,
                                 A=A.M1,cutp=cutpoints,nq,gln,glw,cop.M1,param=F)
  # likelihood in the numerator 
  likelihood.vuong.num=lk_f2vine(tau.M2,ydat=ordinal,
                                 A=A.M2,cutp=cutpoints,nq,gln,glw,cop.M2,param=F)
  
  D = log(likelihood.vuong.num/likelihood.vuong.den) #log-likelihood ratio
  avrg=mean(D) #Average of the D
  
  # Z score
  z = sqrt(n) * avrg/sd(D)
  
  #P-value for  two tailed
  pvalue <- 2 * pnorm(-abs(z)) #two tailed
  
  # The intervals to be added and subtracted to create the 95% CI.
  interval =  1.96 * (1/sqrt(n)) * sd(D) 
  
  # 95% confidence interval
  CI95=c( avrg - interval, avrg + interval )
  
  # Exporting result
  result <- data.frame(z, pvalue,  CI95[1], CI95[2])
  
  names(result) <- c("z", "p.value", "CI.left", "CI.right")
  return(result)
}


#M1: 1-factor tree
#M2: 2-factor tree
vuong.1fv.2fv = function(tau.M1,cop.M1,A.M1, tau.M2, cop.M2,A.M2, 
                         ordinal, nq, param=F){
  #Gauss legendre points and weights
  gl<-gauss.quad.prob(nq)
  gln<-gl$nodes
  glw<-gl$weights
  d<-ncol(ordinal)
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  n<-nrow(ordinal)
  cutpoints<-unifcuts(ordinal)
  #--------------------------------------------------------
  # names of copulas for Model 1
  likelihood.vuong.den=lk_f1vine(tau.M1,ydat=ordinal,
                                 A=A.M1,cutp=cutpoints,nq,gln,glw,cop.M1,param=F)
  
  # likelihood in the numerator 
  likelihood.vuong.num=lk_f2vine(tau.M2,ydat=ordinal,
                                 A=A.M2,cutp=cutpoints,nq,gln,glw,cop.M2,param=F)
  
  
  D = log(likelihood.vuong.num/likelihood.vuong.den) #log-likelihood ratio
  avrg=mean(D) #Average of the D
  
  # Z score
  z = sqrt(n) * avrg/sd(D)
  
  #P-value for  two tailed
  pvalue <- 2 * pnorm(-abs(z)) #two tailed
  
  # The intervals to be added and subtracted to create the 95% CI.
  interval =  1.96 * (1/sqrt(n)) * sd(D) 
  
  # 95% confidence interval
  dim1=length(tau.M1)
  dim2=length(tau.M2)
  AIC.correction= (dim2-dim1)/n
  CI95=c( avrg - AIC.correction -interval, avrg - AIC.correction + interval )
  
  # Exporting result
  result <- data.frame(z, pvalue,  CI95[1], CI95[2])
  
  names(result) <- c("z", "p.value", "CI.left", "CI.right")
  return(result)
}


# The 1-factor likelihood to be used in vuong statistics
lk_f1= function(tau,ydat,copF1,cutp,gln,glw)
{
  d=ncol(cutp)
  # names of copulas
  cops=copulas(d1=0,d2=d,d3=0,cop_name_f1=copF1,cop_name_f2=NULL)
  #cpars
  theta=mapply(function(x,y) tau2par(x,y),x=copF1,y=tau)
  
  lk=dfcopula_f1(NULL,ydat,NULL,dcop=cops$dcop$f1,
                 pcondcop=cops$pcondcop$f1,
                 cutp,NULL,th=theta,gln,glw,param=F)
  return(lk)
} 

# The 2-factor likelihood to be used in vuong statistics
lk_f2= function(tau,ydat,copF1,copF2,cutp,gln,glw)
{
  d=ncol(cutp)
  # names of copulas
  cops=copulas(d1=0,d2=d,d3=0,cop_name_f1=copF1,cop_name_f2=copF2)
  #cpars
  theta1=mapply(function(x,y) tau2par(x,y),x=copF1,y=tau[1:d])
  theta2=mapply(function(x,y) tau2par(x,y),x=copF2,y=tau[(d+1):(d*2)])
  theta=c(theta1,theta2)
  lk=dfcopula_f2(NULL,ydat,NULL,d,dcopf1=cops$dcop$f1,dcopf2=cops$dcop$f2,
                 pcondcopf1=cops$pcondcop$f1,pcondcopf2=cops$pcondcop$f2,
                 cutp,NULL,th=theta,gln,glw,param=F)
  return(lk)
}

#M1: 1-factor
#M2: 1-factor tree
vuong.1fv.1f = function(tau.M1,cop.M1, tau.M2, cop.M2,A.M2, 
                        ordinal, nq){
  #Gauss legendre points and weights
  gl<-gauss.quad.prob(nq)
  gln<-gl$nodes
  glw<-gl$weights
  d<-ncol(ordinal)
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  n<-nrow(ordinal)
  cutpoints<-unifcuts(ordinal)
  #--------------------------------------------------------
  # names of copulas for Model 1
  likelihood.vuong.den=lk_f1(tau.M1,ydat=ordinal,cop.M1,
                             cutp=cutpoints,gln,glw)
  
  # likelihood in the numerator 
  likelihood.vuong.num=lk_f1vine(tau.M2,ydat=ordinal,
                                 A=A.M2,cutp=cutpoints,nq,gln,glw,cop.M2,param=F)
  
  
  D = log(likelihood.vuong.num/likelihood.vuong.den) #log-likelihood ratio
  avrg=mean(D) #Average of the D
  
  # Z score
  z = sqrt(n) * avrg/sd(D)
  
  #P-value for  two tailed
  pvalue <- 2 * pnorm(-abs(z)) #two tailed
  
  # The intervals to be added and subtracted to create the 95% CI.
  interval =  1.96 * (1/sqrt(n)) * sd(D) 
  
  # 95% confidence interval
  dim1=length(tau.M1)
  dim2=length(tau.M2)
  AIC.correction= (dim2-dim1)/n
  CI95=c( avrg - AIC.correction -interval, avrg - AIC.correction + interval )
  
  # Exporting result
  result <- data.frame(z, pvalue,  CI95[1], CI95[2])
  
  names(result) <- c("z", "p.value", "CI.left", "CI.right")
  return(result)
}

#M1: 2-factor
#M2: 1-factor tree
vuong.1fv.2f = function(tau.M1,cop.M1F1,cop.M1F2, tau.M2, cop.M2,A.M2, 
                        ordinal, nq){
  #Gauss legendre points and weights
  gl<-gauss.quad.prob(nq)
  gln<-gl$nodes
  glw<-gl$weights
  d<-ncol(ordinal)
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  n<-nrow(ordinal)
  cutpoints<-unifcuts(ordinal)
  #--------------------------------------------------------
  # names of copulas for Model 1
  likelihood.vuong.den=lk_f2(tau.M1,ydat=ordinal,cop.M1F1,cop.M1F2,
                             cutp=cutpoints,gln,glw)
  
  # likelihood in the numerator 
  likelihood.vuong.num=lk_f1vine(tau.M2,ydat=ordinal,
                                 A=A.M2,cutp=cutpoints,nq,gln,glw,cop.M2,param=F)
  
  
  D = log(likelihood.vuong.num/likelihood.vuong.den) #log-likelihood ratio
  avrg=mean(D) #Average of the D
  
  # Z score
  z = sqrt(n) * avrg/sd(D)
  
  #P-value for  two tailed
  pvalue <- 2 * pnorm(-abs(z)) #two tailed
  
  # The intervals to be added and subtracted to create the 95% CI.
  interval =  1.96 * (1/sqrt(n)) * sd(D) 
  
  # 95% confidence interval
  dim1=length(tau.M1)
  dim2=length(tau.M2)
  AIC.correction= (dim2-dim1)/n
  CI95=c( avrg - AIC.correction -interval, avrg - AIC.correction + interval )
  
  # Exporting result
  result <- data.frame(z, pvalue,  CI95[1], CI95[2])
  
  names(result) <- c("z", "p.value", "CI.left", "CI.right")
  return(result)
}


#M1: 1-factor
#M2: 2-factor tree
vuong.2fv.1f = function(tau.M1,cop.M1, tau.M2, cop.M2,A.M2, 
                        ordinal, nq){
  #Gauss legendre points and weights
  gl<-gauss.quad.prob(nq)
  gln<-gl$nodes
  glw<-gl$weights
  d<-ncol(ordinal)
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  n<-nrow(ordinal)
  cutpoints<-unifcuts(ordinal)
  #--------------------------------------------------------
  # names of copulas for Model 1
  likelihood.vuong.den=lk_f1(tau.M1,ydat=ordinal,cop.M1,
                             cutp=cutpoints,gln,glw)
  
  # likelihood in the numerator 
  likelihood.vuong.num=lk_f2vine(tau.M2,ydat=ordinal,
                                 A=A.M2,cutp=cutpoints,nq,gln,glw,cop.M2,param=F)
  
  
  D = log(likelihood.vuong.num/likelihood.vuong.den) #log-likelihood ratio
  avrg=mean(D) #Average of the D
  
  # Z score
  z = sqrt(n) * avrg/sd(D)
  
  #P-value for  two tailed
  pvalue <- 2 * pnorm(-abs(z)) #two tailed
  
  # The intervals to be added and subtracted to create the 95% CI.
  interval =  1.96 * (1/sqrt(n)) * sd(D) 
  
  # 95% confidence interval
  dim1=length(tau.M1)
  dim2=length(tau.M2)
  AIC.correction= (dim2-dim1)/n
  CI95=c( avrg - AIC.correction -interval, avrg - AIC.correction + interval )
  
  # Exporting result
  result <- data.frame(z, pvalue,  CI95[1], CI95[2])
  
  names(result) <- c("z", "p.value", "CI.left", "CI.right")
  return(result)
}

#M1: 2-factor
#M2: 2-factor tree
vuong.2fv.2f = function(tau.M1,cop.M1F1,cop.M1F2, tau.M2, cop.M2,A.M2, 
                        ordinal, nq){
  #Gauss legendre points and weights
  gl<-gauss.quad.prob(nq)
  gln<-gl$nodes
  glw<-gl$weights
  d<-ncol(ordinal)
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  n<-nrow(ordinal)
  cutpoints<-unifcuts(ordinal)
  #--------------------------------------------------------
  # names of copulas for Model 1
  likelihood.vuong.den=lk_f2(tau.M1,ydat=ordinal,cop.M1F1,cop.M1F2,
                             cutp=cutpoints,gln,glw)
  
  # likelihood in the numerator 
  likelihood.vuong.num=lk_f2vine(tau.M2,ydat=ordinal,
                                 A=A.M2,cutp=cutpoints,nq,gln,glw,cop.M2,param=F)
  
  
  D = log(likelihood.vuong.num/likelihood.vuong.den) #log-likelihood ratio
  avrg=mean(D) #Average of the D
  
  # Z score
  z = sqrt(n) * avrg/sd(D)
  
  #P-value for  two tailed
  pvalue <- 2 * pnorm(-abs(z)) #two tailed
  
  # The intervals to be added and subtracted to create the 95% CI.
  interval =  1.96 * (1/sqrt(n)) * sd(D) 
  
  # 95% confidence interval
  dim1=length(tau.M1)
  dim2=length(tau.M2)
  AIC.correction= (dim2-dim1)/n
  CI95=c( avrg - AIC.correction -interval, avrg - AIC.correction + interval )
  
  # Exporting result
  result <- data.frame(z, pvalue,  CI95[1], CI95[2])
  
  names(result) <- c("z", "p.value", "CI.left", "CI.right")
  return(result)
}
