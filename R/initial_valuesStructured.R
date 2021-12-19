# Purpose: Providing the inital values for copula parameters given the name of the copula
#         
# Inputs
# copX0: copula names that link the observed variables with the 
#        common latent variable C_{Y_{jg}|X_0}
# copXg: copula names that link the univariate conditional dist. C_{Y_{jg}|X_0}
#         with the group-specific latent variable C_{Y_{jg},X_0;X_0}

# Output: initla values for given copula names
initval=function(copX0,copXg)
{
  copnX0Xg=c(copX0,copXg)
  d=length(c(copX0,copXg))
  initvalX0GY=rep(NA,d)
  for(j in 1:d){
    if(copnX0Xg[[j]]=="bvn"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="bvt"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="bvt1"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="bvt2"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="bvt3"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="bvt4"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="bvt5"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="bvt6"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="bvt7"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="bvt8"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="bvt9"){initvalX0GY[[j]]=0.2}
    
    if(copnX0Xg[[j]]=="frk"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="gum"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="sgum"||copnX0Xg[[j]]=="rgum"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="gum90"||copnX0Xg[[j]]=="1rgum"){initvalX0GY[[j]]=-0.2}
    if(copnX0Xg[[j]]=="gum270"||copnX0Xg[[j]]=="2rgum"){initvalX0GY[[j]]=-0.2}
    if(copnX0Xg[[j]]=="joe"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="sjoe"||copnX0Xg[[j]]=="rjoe"){initvalX0GY[[j]]=0.2}
    if(copnX0Xg[[j]]=="joe90"||copnX0Xg[[j]]=="1rjoe"){initvalX0GY[[j]]=-0.2}
    if(copnX0Xg[[j]]=="joe270"||copnX0Xg[[j]]=="2rjoe"){initvalX0GY[[j]]=-0.2}
  }
  output=c(initvalX0GY=initvalX0GY)
  return(output)
}

