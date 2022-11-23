# Purpose: Providing the inital values for copula parameters given the name of the copula
#         
# Inputs
# copX0: copula names that link the observed variables with the 
#        common latent variable C_{Y_{jg}|X_0}
# copXg: copula names that link the univariate conditional dist. C_{Y_{jg}|X_0}
#         with the group-specific latent variable C_{Y_{jg},X_0;X_0}

# Output: initla values for given copula names
initval_tree=function(cop)
{
  d=length(cop)
  initval=rep(NA,d)
  for(j in 1:d){
    if(cop[[j]]=="bvn"){initval[[j]]=0.1}
    if(cop[[j]]=="bvt"){initval[[j]]=0.1}
    if(cop[[j]]=="bvt1"){initval[[j]]=0.1}
    if(cop[[j]]=="bvt2"){initval[[j]]=0.1}
    if(cop[[j]]=="bvt3"){initval[[j]]=0.1}
    if(cop[[j]]=="bvt4"){initval[[j]]=0.1}
    if(cop[[j]]=="bvt5"){initval[[j]]=0.1}
    if(cop[[j]]=="bvt6"){initval[[j]]=0.1}
    if(cop[[j]]=="bvt7"){initval[[j]]=0.1}
    if(cop[[j]]=="bvt8"){initval[[j]]=0.1}
    if(cop[[j]]=="bvt9"){initval[[j]]=0.1}
    
    if(cop[[j]]=="frk"){initval[[j]]=0.1}
    if(cop[[j]]=="gum"){initval[[j]]=0.1}
    if(cop[[j]]=="sgum"||cop[[j]]=="rgum"){initval[[j]]=0.1}
    if(cop[[j]]=="gum90"||cop[[j]]=="1rgum"){initval[[j]]=-0.1}
    if(cop[[j]]=="gum270"||cop[[j]]=="2rgum"){initval[[j]]=-0.1}
  }
  output=c(initval=initval)
  return(output)
}

