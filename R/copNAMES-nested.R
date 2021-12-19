# Purpose: Providing the copula functions for the 
#          observed variables.
#         
# Inputs
# d     : number of total observed variables.
# copX0: copula names that link the observed variables with the 
#        common latent variable C_{Y_{jg}|X_0}
# copXg: copula names that link the univariate conditional dist. C_{Y_{jg}|X_0}
#         with the group-specific latent variable C_{Y_{jg},X_0;X_0}

# Output: Copula functions for all variables
#         stored in lists, and to be used in the pr. mass function.
copulas_nested=function(copX0Xg,copXgY)
{
  ngrp=length(copX0Xg)
  d=length(copXgY)
  dcopX0Xg=list()
  pcondXgY=list()
  qcopX0Xg=list()
  qcondXgY=list()
  for(j in 1:d){
    if(copXgY[[j]]=="bvn"){pcondXgY[[j]]=pcond.bvn;qcondXgY[[j]]=qcond.bvn}
    if(copXgY[[j]]=="bvt"){pcondXgY[[j]]=pcond.t;qcondXgY[[j]]=qcond.t}
    if(copXgY[[j]]=="bvt1"){pcondXgY[[j]]=pcond.t1;qcondXgY[[j]]=qcond.t1}
    if(copXgY[[j]]=="bvt2"){pcondXgY[[j]]=pcond.t2;qcondXgY[[j]]=qcond.t2}
    if(copXgY[[j]]=="bvt3"){pcondXgY[[j]]=pcond.t3;qcondXgY[[j]]=qcond.t3}
    if(copXgY[[j]]=="bvt4"){pcondXgY[[j]]=pcond.t4;qcondXgY[[j]]=qcond.t4}
    if(copXgY[[j]]=="bvt5"){pcondXgY[[j]]=pcond.t5;qcondXgY[[j]]=qcond.t5}
    if(copXgY[[j]]=="bvt6"){pcondXgY[[j]]=pcond.t6;qcondXgY[[j]]=qcond.t6}
    if(copXgY[[j]]=="bvt7"){pcondXgY[[j]]=pcond.t7;qcondXgY[[j]]=qcond.t7}
    if(copXgY[[j]]=="bvt8"){pcondXgY[[j]]=pcond.t8;qcondXgY[[j]]=qcond.t8}
    if(copXgY[[j]]=="bvt9"){pcondXgY[[j]]=pcond.t9;qcondXgY[[j]]=qcond.t9}
    
    if(copXgY[[j]]=="frk"){pcondXgY[[j]]=pcond.frank;qcondXgY[[j]]=qcond.frank}
    
    if(copXgY[[j]]=="gum"){pcondXgY[[j]]=pcond.gumbel;qcondXgY[[j]]=qcond.gumbel}
    if(copXgY[[j]]=="sgum"||copXgY[[j]]=="rgum"){pcondXgY[[j]]=pcond.sgumbel;qcondXgY[[j]]=qcond.sgumbel}
    if(copXgY[[j]]=="gum90"||copXgY[[j]]=="1rgum"){pcondXgY[[j]]=pcond.gumbel.90;qcondXgY[[j]]=qcond.gumbel.90}
    if(copXgY[[j]]=="gum270"||copXgY[[j]]=="2rgum"){pcondXgY[[j]]=pcond.gumbel.270;qcondXgY[[j]]=qcond.gumbel.270}
    
    if(copXgY[[j]]=="joe"){pcondXgY[[j]]=pcond.joe;qcondXgY[[j]]=qcond.joe}
    if(copXgY[[j]]=="sjoe"||copXgY[[j]]=="rjoe"){pcondXgY[[j]]=pcond.sjoe;qcondXgY[[j]]=qcond.sjoe}
    if(copXgY[[j]]=="joe90"||copXgY[[j]]=="1rjoe"){pcondXgY[[j]]=pcond.joe.90;qcondXgY[[j]]=qcond.joe.90}
    if(copXgY[[j]]=="joe270"||copXgY[[j]]=="2rjoe"){pcondXgY[[j]]=pcond.joe.270;qcondXgY[[j]]=qcond.joe.270}
  }
  
  for(j in 1:ngrp){
    if(copX0Xg[[j]]=="bvn"){dcopX0Xg[[j]]=dbvn}
    if(copX0Xg[[j]]=="bvt"){dcopX0Xg[[j]]=dtcop}
    if(copX0Xg[[j]]=="bvt1"){dcopX0Xg[[j]]=dtcop1}
    if(copX0Xg[[j]]=="bvt2"){dcopX0Xg[[j]]=dtcop2}
    if(copX0Xg[[j]]=="bvt3"){dcopX0Xg[[j]]=dtcop3}
    if(copX0Xg[[j]]=="bvt4"){dcopX0Xg[[j]]=dtcop4}
    if(copX0Xg[[j]]=="bvt5"){dcopX0Xg[[j]]=dtcop5}
    if(copX0Xg[[j]]=="bvt6"){dcopX0Xg[[j]]=dtcop6}
    if(copX0Xg[[j]]=="bvt7"){dcopX0Xg[[j]]=dtcop7}
    if(copX0Xg[[j]]=="bvt8"){dcopX0Xg[[j]]=dtcop8}
    if(copX0Xg[[j]]=="bvt9"){dcopX0Xg[[j]]=dtcop9}
    
    if(copX0Xg[[j]]=="frk"){dcopX0Xg[[j]]=dfrank}
    
    if(copX0Xg[[j]]=="gum"){dcopX0Xg[[j]]=dgumbel}
    if(copX0Xg[[j]]=="sgum"||copX0Xg[[j]]=="rgum"){dcopX0Xg[[j]]=dsgumbel}
    if(copX0Xg[[j]]=="gum90"||copX0Xg[[j]]=="1rgum"){dcopX0Xg[[j]]=dgumbel.90}
    if(copX0Xg[[j]]=="gum270"||copX0Xg[[j]]=="2rgum"){dcopX0Xg[[j]]=dgumbel.270}
    
    if(copX0Xg[[j]]=="joe"){dcopX0Xg[[j]]=djoe}
    if(copX0Xg[[j]]=="sjoe"||copX0Xg[[j]]=="rjoe"){dcopX0Xg[[j]]=dsjoe}
    if(copX0Xg[[j]]=="joe90"||copX0Xg[[j]]=="1rjoe"){dcopX0Xg[[j]]=djoe.90}
    if(copX0Xg[[j]]=="joe270"||copX0Xg[[j]]=="2rjoe"){dcopX0Xg[[j]]=djoe.270}
  }
  
  for(j in 1:ngrp){
    if(copX0Xg[[j]]=="bvn"){qcopX0Xg[[j]]=qcond.bvn}
    if(copX0Xg[[j]]=="bvt"){qcopX0Xg[[j]]=qcond.t}
    if(copX0Xg[[j]]=="bvt1"){qcopX0Xg[[j]]=qcond.t1}
    if(copX0Xg[[j]]=="bvt2"){qcopX0Xg[[j]]=qcond.t2}
    if(copX0Xg[[j]]=="bvt3"){qcopX0Xg[[j]]=qcond.t3}
    if(copX0Xg[[j]]=="bvt4"){qcopX0Xg[[j]]=qcond.t4}
    if(copX0Xg[[j]]=="bvt5"){qcopX0Xg[[j]]=qcond.t5}
    if(copX0Xg[[j]]=="bvt6"){qcopX0Xg[[j]]=qcond.t6}
    if(copX0Xg[[j]]=="bvt7"){qcopX0Xg[[j]]=qcond.t7}
    if(copX0Xg[[j]]=="bvt8"){qcopX0Xg[[j]]=qcond.t8}
    if(copX0Xg[[j]]=="bvt9"){qcopX0Xg[[j]]=qcond.t9}
    
    if(copX0Xg[[j]]=="frk"){qcopX0Xg[[j]]=qcond.frank}
    
    if(copX0Xg[[j]]=="gum"){qcopX0Xg[[j]]=qcond.gumbel}
    if(copX0Xg[[j]]=="sgum"||copX0Xg[[j]]=="rgum"){qcopX0Xg[[j]]=qcond.sgumbel}
    if(copX0Xg[[j]]=="gum90"||copX0Xg[[j]]=="1rgum"){qcopX0Xg[[j]]=qcond.gumbel.90}
    if(copX0Xg[[j]]=="gum270"||copX0Xg[[j]]=="2rgum"){qcopX0Xg[[j]]=qcond.gumbel.270}
    
    if(copX0Xg[[j]]=="joe"){qcopX0Xg[[j]]=qcond.joe}
    if(copX0Xg[[j]]=="sjoe"||copX0Xg[[j]]=="rjoe"){qcopX0Xg[[j]]=qcond.sjoe}
    if(copX0Xg[[j]]=="joe90"||copX0Xg[[j]]=="1rjoe"){qcopX0Xg[[j]]=qcond.joe.90}
    if(copX0Xg[[j]]=="joe270"||copX0Xg[[j]]=="2rjoe"){qcopX0Xg[[j]]=qcond.joe.270}
    
  }
  
  output=list(dcopX0Xg=dcopX0Xg,pcondcopXgY=pcondXgY,qcopX0Xg=qcopX0Xg,qcondXgY=qcondXgY)
  return(output)
}
