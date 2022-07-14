# This functoin is used for the Vuong test. 
#
# purpose ========================
# Calculate the Vuong statistics,
# p-value and 95% confidence interval CI will be exported.
# To detect whether two models are statistically different.
# Input   ========================
# copF1          : copula names
# cpar.bvn       : estimated parameters for the model that needs to be compared
# cpar           : estimated parameters for the model that needs to be compared

# Output   ========================
# Vector of z score, p-value, 95% CI

vuong.1f= function(cpar.bvn, cpar, copF1, continuous, ordinal, count, gl, param=FALSE){
  #Gauss legendre points and weights
  gln<-gl$nodes
  glw<-gl$weights
  
  d1=ifelse(is.null(continuous),0,ncol(continuous))
  d2=ifelse(is.null(ordinal),0,ncol(ordinal))
  d3=ifelse(is.null(count),0,ncol(count))
  d=d1+d2+d3
  #--------------------------------------------------------
  #                 Continuous variables
  #-----------------                    -------------------
  u <- NULL
  if(!is.null(continuous)){
    n<-nrow(continuous)
    u<-apply(continuous,2,rank)/(n+1)
  }
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  cutpoints<-NULL
  if(!is.null(ordinal)){
    n<-nrow(ordinal)
    cutpoints<-cuts(ordinal)
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
  }
  
  #--------------------------------------------------------
  cpar=as.list(mapply( function(x, y) if(x=="1rjoe" || x=="2rjoe" 
                    || x=="1rgum"|| x=="2rgum"){abs(y)}else{y},
                          x=copF1, y=cpar))
  #--------------------------------------------------------
  # names of copulas
  cops=copulas(d1,d2,d3,cop_name_f1=copF1,cop_name_f2=NULL)
  # likelihood in the numerator 
  likelihood.vuong.num=dfcopula_f1(u,ordinal,count,dcop=cops$dcop$f1,
                                   pcondcop=cops$pcondcop$f1,
                                   cutpoints,est_count,th=cpar,gln,glw,param)
  
  copF1.den = rep("bvn", d)
  cops.den=copulas(d1,d2,d3,cop_name_f1=copF1.den,cop_name_f2=NULL)
  likelihood.vuong.den=dfcopula_f1(u,ordinal,count,dcop=cops.den$dcop$f1,
                                   pcondcop=cops.den$pcondcop$f1,
                                   cutpoints,est_count,th=cpar.bvn,gln,glw,param)
  
  n = length(likelihood.vuong.num) #number of observation
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



vuong.2f= function(cpar.bvn, cpar, copF1, copF2, continuous, ordinal, count, gl, param=FALSE){
  #Gauss legendre points and weights
  gln<-gl$nodes
  glw<-gl$weights
  
  d1=ifelse(is.null(continuous),0,ncol(continuous))
  d2=ifelse(is.null(ordinal),0,ncol(ordinal))
  d3=ifelse(is.null(count),0,ncol(count))
  d=d1+d2+d3
  #--------------------------------------------------------
  #                 Continuous variables
  #-----------------                    -------------------
  u <- NULL
  if(!is.null(continuous)){
    n<-nrow(continuous)
    u<-apply(continuous,2,rank)/(n+1)
  }
  #--------------------------------------------------------
  #                  Ordinal variables
  #------------------                 ---------------------
  cutpoints<-NULL
  if(!is.null(ordinal)){
    n<-nrow(ordinal)
    cutpoints<-cuts(ordinal)
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
  }
  
  #----------------------------------------------------------
  #change negative to positive copula parameters
  thetaF1 = as.list(mapply( function(x, y) 
    if(x=="1rjoe" || x=="2rjoe" || x=="1rgum"|| x=="2rgum"){
      abs(y)}else{y}, x=copF1, y=cpar[[1]]))
  thetaF2 = as.list(mapply( function(x, y) 
    if(x=="1rjoe" || x=="2rjoe" || x=="1rgum"|| x=="2rgum"){
      abs(y)}else{y}, x=copF2, y=cpar[[2]]))
  
  theta = list(thetaF1, thetaF2)
  #--------------------------------------------------------
  # names of copulas
  cops=copulas(d1,d2,d3,cop_name_f1=copF1,cop_name_f2=copF2)
  # likelihood in the numerator 
  likelihood.vuong.num=dfcopula_f2(u,ordinal,count,d,dcopf1=cops$dcop$f1,dcopf2=cops$dcop$f2,
                               pcondcopf1=cops$pcondcop$f1,pcondcopf2=cops$pcondcop$f2,
                               cutpoints,est_count,th=theta,gln,glw,param)
  
  copF1.den = copF2.den = rep("bvn", d)
  cops.den=copulas(d1,d2,d3,cop_name_f1=copF1.den,cop_name_f2=copF2.den)
  likelihood.vuong.den=dfcopula_f2(u,ordinal,count,d,dcopf1=cops.den$dcop$f1,dcopf2=cops.den$dcop$f2,
                               pcondcopf1=cops.den$pcondcop$f1,pcondcopf2=cops.den$pcondcop$f2,
                               cutpoints,est_count,th=cpar.bvn,gln,glw,param)
  
  n = length(likelihood.vuong.num) #number of observation
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


