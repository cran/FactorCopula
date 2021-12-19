# This functoin is the Vuong test. 
# purpose ========================
# Calculate the Vuong statistics,
# p-value and 95% confidence interval CI will be exported.
# To detect whether two models are statistically different.
# Input   ========================
# cpar.M1        : estimated parameters for the model 1
# cpar.M2        : estimated parameters for the model 2
# copF0.M1       : copula names that link item with common factor for Model 1
# copFg.M1       : copula names that link item with group factor for Model 1
# copF0.M2       : copula names that link item with common factor for Model 2
# copFg.M2       : copula names that link item with group factor for Model 2
# ordinal        : ordinal data
# ngrp           : number of groups
# grpsize        : a vector of size of each group
# Output   ========================
# Vector of z score, p-value, and 95% CI

vuong.bifactor = function(cpar.M1,copF0.M1, copFg.M1, cpar.M2, copF0.M2, copFg.M2, 
                          ordinal,ngrp, grpsize, nq, param=F){
  #Gauss legendre points and weights
  gl<-gauss.quad.prob(nq)
  gln<-gl$nodes
  glw<-gl$weights
  d=ncol(ordinal)
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  n<-nrow(ordinal)
  cutpoints<-cuts(ordinal)
  #--------------------------------------------------------
  # names of copulas for Model 1
  cops.M1=copulas_bifactor(copF0.M1, copFg.M1)
  likelihood.vuong.den=bipmf(ordinal, cpar.M1[1:d], cpar.M1[(d+1):(d*2)], cutpoints, 
                             pcondcopX0 = cops.M1$pcondcopX0, 
                             pcondcopXg = cops.M1$pcondcopXgY, nq, ngrp, grpsize, glw, 
                             gln, param)
  
  # names of copulas for Model 2
  cops.M2=copulas_bifactor(copF0.M2, copFg.M2)
  # likelihood in the numerator 
  likelihood.vuong.num=bipmf(ordinal, cpar.M2[1:d], cpar.M2[(d+1):(d*2)], cutpoints, 
                             pcondcopX0 = cops.M2$pcondcopX0, 
                             pcondcopXg = cops.M2$pcondcopXgY, nq, ngrp, grpsize, glw, 
                             gln, param)
  
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


#vuong.bifactor(cpar.M1=th.bvn,copF0.M1=rep("bvn",20), copFg.M1=rep("bvn",20),
#               cpar.M2=th,copF0.M2=copX0Xg, copFg.M2=copXgY,
#               ordinal=ydat,ngrp=ngrp, grpsize=grpsize, nq=25, param=F)


vuong.nested = function(cpar.M1,copF0.M1, copFg.M1, cpar.M2, copF0.M2, copFg.M2, 
                        ordinal,ngrp, grpsize, nq, param=F){
  #Gauss legendre points and weights
  gl<-gauss.quad.prob(nq)
  gln<-gl$nodes
  glw<-gl$weights
  d=ncol(ordinal)
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  n<-nrow(ordinal)
  cutpoints<-cuts(ordinal)
  #--------------------------------------------------------
  # names of copulas for Model 1
  cops.M1=copulas_nested(copF0.M1, copFg.M1)
  likelihood.vuong.den=nestedcop.ir(ordinal, cutpoints, dcop = cops.M1$dcopX0Xg, 
                                    pcondcop = cops.M1$pcondcopXgY, qcondcop = cops.M1$qcopX0Xg, 
                                    th.x0g = cpar.M1[1:ngrp], th.xgy = cpar.M1[(ngrp+1):(d+ngrp)],
                                    nq, ngrp, grpsize, glw, gln, param)
  
  # names of copulas for Model 2
  cops.M2=copulas_nested(copF0.M2, copFg.M2)
  # likelihood in the numerator 
  likelihood.vuong.num=nestedcop.ir(ordinal, cutpoints, dcop = cops.M2$dcopX0Xg, 
                                    pcondcop = cops.M2$pcondcopXgY, qcondcop = cops.M2$qcopX0Xg, 
                                    th.x0g = cpar.M2[1:ngrp], th.xgy = cpar.M2[(ngrp+1):(d+ngrp)],
                                    nq, ngrp, grpsize, glw, gln, param)
  
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

#vuong.nested(cpar.M1=rep(0.4,d+ngrp),copF0.M1=rep("bvn",ngrp), copFg.M1=rep("bvn",d), 
#             cpar.M2=rep(0.4,d+ngrp),copF0.M2=rep("frk",ngrp), copFg.M2=rep("frk",d), 
#             ydat,ngrp, grpsize, nq, param=F)


# M1 is nested model
vuong.bifactor.nested = function(cpar.M1,copF0.M1, copFg.M1, cpar.M2, copF0.M2, copFg.M2, 
                                 ordinal,ngrp, grpsize, nq, param=F){
  #Gauss legendre points and weights
  gl<-gauss.quad.prob(nq)
  gln<-gl$nodes
  glw<-gl$weights
  d=ncol(ordinal)
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  n<-nrow(ordinal)
  cutpoints<-cuts(ordinal)
  #--------------------------------------------------------
  # names of copulas for Model 1
  cops.M1=copulas_nested(copF0.M1, copFg.M1)
  likelihood.vuong.den=nestedcop.ir(ordinal, cutpoints, dcop = cops.M2$dcopX0Xg, 
                                    pcondcop = cops.M1$pcondcopXgY, qcondcop = cops.M1$qcopX0Xg, 
                                    th.x0g = cpar.M1[1:ngrp], th.xgy = cpar.M1[(ngrp+1):(d+ngrp)],
                                    nq, ngrp, grpsize, glw, gln, param)
  
  # names of copulas for Model 2
  cops.M2=copulas_bifactor(copF0.M2, copFg.M2)
  # likelihood in the numerator 
  likelihood.vuong.num=bipmf(ordinal, cpar.M2[1:d], cpar.M2[(d+1):(d*2)], cutpoints, 
                             pcondcopX0 = cops.M2$pcondcopX0, 
                             pcondcopXg = cops.M2$pcondcopXgY, nq, ngrp, grpsize, glw, 
                             gln, param)
  
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

# M1 is bi-factor model
vuong.nested.bifactor = function(cpar.M1,copF0.M1, copFg.M1, cpar.M2, copF0.M2, copFg.M2, 
                                 ordinal,ngrp, grpsize, nq, param=F){
  #Gauss legendre points and weights
  gl<-gauss.quad.prob(nq)
  gln<-gl$nodes
  glw<-gl$weights
  d=ncol(ordinal)
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  n<-nrow(ordinal)
  cutpoints<-cuts(ordinal)
  #--------------------------------------------------------
  # names of copulas for Model 2
  cops.M2=copulas_nested(copF0.M2, copFg.M2)
  likelihood.vuong.num=nestedcop.ir(ordinal, cutpoints, dcop = cops.M2$dcopX0Xg, 
                                    pcondcop = cops.M2$pcondcopXgY, qcondcop = cops.M2$qcopX0Xg, 
                                    th.x0g = cpar.M2[1:ngrp], th.xgy = cpar.M2[(ngrp+1):(d+ngrp)],
                                    nq, ngrp, grpsize, glw, gln, param)
  
  # names of copulas for Model 1
  cops.M1=copulas_bifactor(copF0.M1, copFg.M1)
  # likelihood in the denominator 
  likelihood.vuong.den=bipmf(ordinal, cpar.M1[1:d], cpar.M1[(d+1):(d*2)], cutpoints, 
                             pcondcopX0 = cops.M1$pcondcopX0, 
                             pcondcopXg = cops.M1$pcondcopXgY, nq, ngrp, grpsize, glw, 
                             gln, param)
  
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
