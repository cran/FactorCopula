# This code is to convert the copula parameters to the (- inf, inf) space to 
# avoid convergence problems and boundary issues. Also, there is 
# a following a code to convert the estimated paramters back to their
# original code.


#Input: 
# 1
# copname: 
# 1- for Joe and Gumbel and the survival copulas: "gum" 
# 2- for BVN and t-copula: "bvn" 
# 3- for BB1, BB7, and the survival versions: "bb1","bb7"
# 2
# par: copula paramter in the original scale.
# par2: the second copula parameter for the two-parameter families.

#Output: a vector of parameter(s) after parametrisation.

linkfun=function(copname, par){
  if(copname=="gum" |copname=="rgum" |
     copname=="joe" |copname=="rjoe"){
    delta=log(par-1)
  }else if(copname=="1rgum" |copname=="2rgum" 
           |copname=="1rjoe" |copname=="2rjoe"){
    delta=log(abs(par)-1)
  }else if(copname=="bvn" | copname=="bvt1"| 
           copname=="bvt2"| copname=="bvt3"| 
           copname=="bvt4"| copname=="bvt5"| 
           copname=="bvt6"| copname=="bvt7"| 
           copname=="bvt8"| copname=="bvt9"){
      delta=log((1+par)/(1-par))
    }else if(copname=="frank" | copname=="frk"){
        delta=par
      }else if( length(par)==2){
        if(copname=="bb1" | copname=="rbb1"){
          delta=log(par[1])
          delta2=log(par[2]-1)
          delta=cbind(delta,delta2)
        }else
          if(copname=="bb8" | copname=="rbb8"){
            delta=log(par[1]-1)
            delta2=log(par[2]/(1-par[2]))
            delta=cbind(delta,delta2)
          }else
            if(copname=="bb7" | copname=="rbb7"){
              delta=log(par[1]-1)
              delta2=log(par[2])
              delta=cbind(delta,delta2)
            }  else 
              if(copname=="bb10" | copname=="rbb10" ){
                delta=log(par[1])
                delta2=log(par[2]/(1-par[2]))
                delta=cbind(delta,delta2)
              }
      }
  else{stop("number of parameters must be two for the two-parameter family")}
  return(delta)
}


#Input: 
# 1
# copname: 
# 1- for Joe and Gumbel and the survival copulas: "gum" 
# 2- for BVN and t-copula: "bvn" 
# 3- for BB1, BB7, and the survival versions: "bb1","bb7"
# 2
# delta: copula paramter after parametrisation.
# delta2: the second copula parameter for the two-parameter families after parametrisation.

#Output: a vector of parameter(s) in the original scale

linkfuninv=function(copname, delta ){
  if(copname=="gum" |copname=="rgum"
     |copname=="joe" |copname=="rjoe"){
    theta=exp(delta)+1
  }else if(copname=="1rgum" |copname=="2rgum"  
                |copname=="1rjoe" |copname=="2rjoe"){
    theta=-(exp(delta)+1)
  }else if((copname=="bvn" | copname=="bvt1"| 
        copname=="bvt2"| copname=="bvt3"| 
        copname=="bvt4"| copname=="bvt5"| 
        copname=="bvt6"| copname=="bvt7"| 
        copname=="bvt8"| copname=="bvt9")){
      theta=((exp(delta)-1)/(exp(delta)+1))
    }else
      if(copname=="frank" | copname=="frk"){
        theta=delta
      }else if( length(delta)==2){
        if(copname=="bb1" | copname=="rbb1"){
          theta1=exp(delta[1])
          theta2=exp(delta[2])+1
          theta=cbind(theta1,theta2)
        }else if(copname=="bb7" | copname=="rbb7" ){
          theta1=exp(delta[1])+1
          theta2=exp(delta[2])
          theta=cbind(theta1,theta2)
        }else if(copname=="bb8" | copname=="rbb8"){
            theta1=exp(delta[1])+1
            theta2=exp(delta[2])/(1+exp(delta[2]))
            theta=cbind(theta1,theta2)
          }else if(copname=="bb10" | copname=="rbb10" ){
                theta1=exp(delta[1])
                theta2=exp(delta[2])/(1+exp(delta[2]))
                theta=cbind(theta1,theta2)
              }
      }
  else{stop("parameters are not suitable for two-parameter family")}
  return(theta)
}




# Purpose: Convert from parameterisation to copula parameters
#          in a more convenient way, especially if we
#          use various bivariate copulas for the variables.

# Input: 
# deltas: the parameterisations of copula parameters
# copF1 : copula names for the first factor.
# copF2 : copula names for the second factor.

# Output: 
# a list of copula parameters for the first factor, and the second factor is specefied.

delta2cpar=function(deltas,copF1,copF2=NULL){
  d<-length(copF1)
  pars1=list()
  pars2=list()
  for(j in 1:(d)){
    pars1[[j]]=linkfuninv(copF1[[j]],deltas[[1]][[j]] )
    
    if(is.null(copF2)){pars2=NULL}else{
      pars2[[j]]=linkfuninv(copF2[[j]],deltas[[2]][[j]] )}
  }
  if(is.null(copF2)){output=list(f1=pars1)}else{output=list(f1=pars1,f2=pars2)}
  output
}