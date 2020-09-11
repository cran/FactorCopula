# Purpose: calculate the SE using delta method
# Input:::::::::::::::
# copname:copula names.
# delta : parameters transformed in domain R; chech "linkfun" and "linkfuninv"
# SE  : sqrt(diag(solve(model$hessian))) of the reparameterization.

#Output: a vector of standard error of the Kendall's taus.
deltamethodSE=function(copname, delta, variance, covariance){
  
  if(copname=="gum" |copname=="rgum" | copname=="1rgum" |copname=="2rgum"){
    stdrror=  exp(delta)/(exp(delta) + 1)^2 * sqrt(variance)
    
  }else if(copname=="joe" |copname=="rjoe"|copname=="1rjoe" |copname=="2rjoe"){
      stdrror= (2 * exp(delta)/(exp(delta) + 1)^2 * trigamma(2/(exp(delta) + 1) + 1) * 
          2/(2 - (exp(delta) + 1)) + (digamma(2) - digamma(2/(exp(delta) + 
          1) + 1)) * 2 * exp(delta)/(2 - (exp(delta) + 1))^2 ) * sqrt(variance)
      
    }else if(copname=="bvn"|copname=="bvt1"|copname=="bvt2"|copname=="bvt3"|
       copname=="bvt4"|copname=="bvt5"|copname=="bvt6"|copname=="bvt7"
       |copname=="bvt8"|copname=="bvt9"){
      stdrror=(2/pi * ((exp(delta)/(exp(delta) + 1) - (exp(delta) - 1) 
                     * exp(delta)/(exp(delta) + 1)^2)/sqrt(1 - ((exp(delta) - 1)/(exp(delta) + 1))^2)) )* sqrt(variance)
    }else
      if(copname=="frank" | copname=="frk"){
      stdrror= (-4/(delta^2) * (1/delta) * integrate( function(x) x/(exp(x) - 1),
                  lower = 0, upper = delta)$value + (4/delta)*
          (( (delta/(exp(delta)-1)) - (1/delta) * integrate(function(x) x/(exp(x) - 1),
                  lower = 0, upper = delta)$value )/delta) + (4/delta^2) )* sqrt(variance)
      }else if(copname=="bb1"|copname=="rbb1"){
          
        val = bb1var.delta2tau(delta,variance,covariance)
        stdrror = sqrt(val)
          
        }else if(copname=="bb7"|copname=="rbb7"){
          
          val = bb7var.delta2tau(delta, variance, covariance)
          stdrror = sqrt(val)
        }else if(copname=="bb8"||copname=="rbb8"){
          
          val = bb8var.delta2tau(delta, variance, covariance)
          stdrror = sqrt(val)
        }else if(copname=="rbb10"|copname=="bb10"){
          
          val = bb10var.delta2tau(delta, variance, covariance)
          stdrror = sqrt(val)
        }
  return(stdrror)
}


# SE using delta method.
# delta: estimated parameters using reparameteriztion
# variance: variance
# covariance: covariance
# 
bb1var.delta2tau=function(delta,variance,covariance) 
{
  delta1 = delta[1]
  delta2 = delta[2]
  
  #dFdelta1
  dFdelta1 = (2*exp(delta1))/((2+ exp(delta1))^2*(1+exp(delta2)))
  
  #dFdelta2
  dFdelta2 = (2*exp(delta2))/((2+exp(delta1))*(1+exp(delta2))^2)
  
  var.tau=dFdelta1^2 *variance[1]+dFdelta2^2 *variance[2]+2*covariance*dFdelta1*dFdelta2
  var.tau
}

# SE using delta method.
# delta: estimated parameters using reparameteriztion
# variance: variance
# covariance: covariance
# 
bb7var.delta2tau=function(delta,variance,covariance) {
  delta1 = delta[1]
  delta2 = delta[2]
  
  dFdelta1=function(s){
    ed1 = exp(delta1)
    ed2 = exp(delta2)
    ed1p1 = (1+ed1)
    ed2p1 = (1+ed2)
    s1 = 1-s
    
    tem1 = (1/(ed1p1^2))*exp(delta1-delta2) *(s1)^-ed1
    tem2 = 1-(1-(s1)^ed1p1)^ed2p1-(s1)^ed1
    tem3 = (s1)^ed1 *s
    tem4 = (1-(s1)^ed1p1)^ed2
    tem5 = (-1-ed2* (s1)^ed1p1)
    tem6 = ed1p1 *( 1+ tem4 *tem5)* log(s1)
    val = tem1 *(tem2+tem3+ tem6)
    val
  }
  
  
  dFdelta2=function(s){
    
    ed1 = exp(delta1)
    ed2 = exp(delta2)
    ed1p1 = (1+ed1)
    ed2p1 = (1+ed2)
    s1 = 1-s
    
    
    tem1 = (1/ed1p1) *exp(-delta2) *(1-s1^ed1p1)
    tem2 = tem1*s1^-ed1
    tem3 = (1-s1^ed1p1)^ed2
    tem4 = (-1+ed2 *log(1-s1^ed1p1))
    
    val =tem2 * ( 1+ tem3 *tem4)
    val
  }
  
  Tause.d1 = 4*integrate(dFdelta1, 0, 1)$v #* SE_th
  Tause.d2 = 4*integrate(dFdelta2, 0, 1)$v #* SE_th
  
  var.tau=Tause.d1^2 *variance[1]+Tause.d2^2 *variance[2]+2*covariance*Tause.d1*Tause.d2
  var.tau
}


#BB8 delta2tauSE
# SE using delta method.
# delta: estimated parameters using reparameteriztion
# variance: variance
# covariance: covariance
# 

bb8var.delta2tau=function(delta,variance,covariance) {
  
  delta1 = delta[1]
  delta2 = delta[2]
  
  dFdelta1 = function(s){
    ed1 = exp(delta1)
    ed2 = exp(delta2)
    d1s = delta1+s
    d2s = delta2+s
    exps = exp(s)
    exp_s = exp(-s)
    ed1s = exp(d1s)
    ed2s = exp(d2s)
    ed1p1 = 1+ed1
    ed2p1 = 1+ed2
    
    tem1 = -2*(1/((ed1p1)^3))
    tem2 = exp(delta1-2 *(d2s))
    tem3 = (ed2p1)* (1-(1/(ed2p1))^(ed1p1)) 
    tem4 = (1+exp_s * (-1+(1/(ed2p1))^(ed1p1)))^(-2+2/(ed1p1)) *s 
    
    tem5 = tem1 * tem2 * tem3 * tem4
    
    tem6 = (1/(ed2p1))^ed1
    tem7 = (-1-ed2+exps+ed1s+ed2s+exp(delta1+d2s)+(1/(ed2p1))^ed1)
    tem8 = log(1/(ed2p1))/(-1-ed2+exps+ed2s+(1/(ed2p1))^ed1)
    tem9 = (ed2p1-(1/(ed2p1))^ed1)
    tem10 = (ed1p1+log(1+exp_s *(-1+(1/(ed2p1))^(ed1p1))))
    val = tem5*(( tem6 *tem7 *tem8+(tem9*tem10)/(ed1p1)))
    val
  }
  
  dFdelta2 = function(s){
    ed1 = exp(delta1)
    ed2 = exp(delta2)
    d2s = delta2+s
    exps = exp(s)
    ed2s = exp(d2s)
    ed1p1 = 1+ed1
    ed2p1 = 1+ed2
    
    tem1 = 2*exp(-2 * delta2) * (ed2p1)^2 *(ed2p1-(1/(ed2p1))^ed1)
    tem2 = exp(delta1+d2s)* (1/(ed2p1))^ed1
    tem3 = ( -1+(1/(ed2p1))^ed1)
    tem4 =  (-1-ed2+exps+ed2s+(1/(ed2p1))^ed1)
    tem5 = (-1+(1/(ed2p1))^(ed1p1))
    tem6 = (1+exp(-s) *tem5 )^(2/(ed1p1)) *s
    tem7 = (-1-ed2+exps+ed2s+(1/(ed2p1))^ed1)
    tem8 = ((ed1p1)^2* tem7^3)
    val = (tem1*(tem2+tem3*tem4)*tem6 )/tem8
    val
  }
  
  Tause.d1= -4*integrate(dFdelta1, 0, 350, rel.tol = 1e-06)$v #* SE_th
  Tause.d2= -4*integrate(dFdelta2, 0, 350, rel.tol = 1e-06)$v #* SE_th
  
  var.tau=Tause.d1^2 *variance[1]+Tause.d2^2 *variance[2]+2*covariance*Tause.d1*Tause.d2
  var.tau
}



# SE using delta method.
# delta: estimated parameters using reparameteriztion
# variance: variance
# covariance: covariance
# 
bb10var.delta2tau=function(delta,variance,covariance) 
{
  delta1 = delta[1]
  delta2 = delta[2]
  
  dFdelta1=function(s){
    expd1=exp(-delta1)
    expd2=exp(delta2)
    expd21=expd2+1
    
    tem1=(1/((-expd2 + exp(s) + exp(delta2 + s))^2))
    tem2=exp(-3*delta1 + 2*s - 2*expd1*s)
    tem3=(1/expd21)^(-2 + 2*expd1)
    tem4=(1 - exp(delta2 - s)/expd21)^(-2*expd1)
    tem5=log(1/expd21) 
    tem6=log(1 - exp(delta2 - s)/expd21)
    
    dfd1=tem1*2*tem2*tem3*tem4*s*(-exp(delta1) + s - tem5+ tem6)
    dfd1
  }
  
  dFdelta2=function(s){
    expd1=exp(-delta1)
    expd2=exp(delta2)
    expd21=exp(delta2)+1
    expd2s=exp(delta2 + s)
    expd2ns=exp(delta2 - s)
    
    tem1=(1/((-expd2 + exp(s) + expd2s)^3))
    tem2=exp(-3 *delta1 + delta2 +2* s - 2 *expd1* s)
    tem3=(1/expd21)^(-1 + 2 *expd1)
    tem4=(1 + exp(delta1) + expd2 - exp(s) - expd2s)
    tem5=(1-expd2ns/expd21)^(-2* expd1)
    
    dfd2=tem1*2*tem2 *tem3* tem4* tem5* s
    dfd2
  }
  
  Tause.th= -4*integrate(dFdelta1, 0, 350, rel.tol = 1e-06)$v #* SE_th
  Tause.p= -4*integrate(dFdelta2, 0, 350, rel.tol = 1e-06)$v #* SE_th
  
  var.tau=Tause.th^2 *variance[1]+Tause.p^2 *variance[2]+2*covariance*Tause.th*Tause.p
  var.tau
}


