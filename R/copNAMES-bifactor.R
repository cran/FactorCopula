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
copulas_bifactor=function(copX0,copXg)
{
  d=length(copX0)
  pcondcopX0=list()
  qcondcopX0=list()
  pcondcopXg=list()
  qcondcopXg=list()
  for(j in 1:d){
    if(copX0[[j]]=="bvn"){pcondcopX0[[j]]=pcond.bvn;qcondcopX0[[j]]=qcond.bvn}
    if(copX0[[j]]=="bvt"){pcondcopX0[[j]]=pcond.t;qcondcopX0[[j]]=qcond.t}
    if(copX0[[j]]=="bvt1"){pcondcopX0[[j]]=pcond.t1;qcondcopX0[[j]]=qcond.t1}
    if(copX0[[j]]=="bvt2"){pcondcopX0[[j]]=pcond.t2;qcondcopX0[[j]]=qcond.t2}
    if(copX0[[j]]=="bvt3"){pcondcopX0[[j]]=pcond.t3;qcondcopX0[[j]]=qcond.t3}
    if(copX0[[j]]=="bvt4"){pcondcopX0[[j]]=pcond.t4;qcondcopX0[[j]]=qcond.t4}
    if(copX0[[j]]=="bvt5"){pcondcopX0[[j]]=pcond.t5;qcondcopX0[[j]]=qcond.t5}
    if(copX0[[j]]=="bvt6"){pcondcopX0[[j]]=pcond.t6;qcondcopX0[[j]]=qcond.t6}
    if(copX0[[j]]=="bvt7"){pcondcopX0[[j]]=pcond.t7;qcondcopX0[[j]]=qcond.t7}
    if(copX0[[j]]=="bvt8"){pcondcopX0[[j]]=pcond.t8;qcondcopX0[[j]]=qcond.t8}
    if(copX0[[j]]=="bvt9"){pcondcopX0[[j]]=pcond.t9;qcondcopX0[[j]]=qcond.t9}
    
    if(copX0[[j]]=="frk"){pcondcopX0[[j]]=pcond.frank;qcondcopX0[[j]]=qcond.frank}
    
    if(copX0[[j]]=="gum"){pcondcopX0[[j]]=pcond.gumbel;qcondcopX0[[j]]=qcond.gumbel}
    if(copX0[[j]]=="sgum"||copX0[[j]]=="rgum"){pcondcopX0[[j]]=pcond.sgumbel;qcondcopX0[[j]]=qcond.sgumbel}
    if(copX0[[j]]=="gum90"||copX0[[j]]=="1rgum"){pcondcopX0[[j]]=pcond.gumbel.90;qcondcopX0[[j]]=qcond.gumbel.90}
    if(copX0[[j]]=="gum270"||copX0[[j]]=="2rgum"){pcondcopX0[[j]]=pcond.gumbel.270;qcondcopX0[[j]]=qcond.gumbel.270}
    
    if(copX0[[j]]=="joe"){pcondcopX0[[j]]=pcond.joe;qcondcopX0[[j]]=qcond.joe}
    if(copX0[[j]]=="sjoe"||copX0[[j]]=="rjoe"){pcondcopX0[[j]]=pcond.sjoe;qcondcopX0[[j]]=qcond.sjoe}
    if(copX0[[j]]=="joe90"||copX0[[j]]=="1rjoe"){pcondcopX0[[j]]=pcond.joe.90;qcondcopX0[[j]]=qcond.joe.90}
    if(copX0[[j]]=="joe270"||copX0[[j]]=="2rjoe"){pcondcopX0[[j]]=pcond.joe.270;qcondcopX0[[j]]=qcond.joe.270}
  }
  
  
  for(j in 1:d){
    if(copXg[[j]]=="bvn"){pcondcopXg[[j]]=pcond.bvn;qcondcopXg[[j]]=qcond.bvn}
    if(copXg[[j]]=="bvt"){pcondcopXg[[j]]=pcond.t;qcondcopXg[[j]]=qcond.t}
    if(copXg[[j]]=="bvt1"){pcondcopXg[[j]]=pcond.t1;qcondcopXg[[j]]=qcond.t1}    
    if(copXg[[j]]=="bvt2"){pcondcopXg[[j]]=pcond.t2;qcondcopXg[[j]]=qcond.t2}    
    if(copXg[[j]]=="bvt3"){pcondcopXg[[j]]=pcond.t3;qcondcopXg[[j]]=qcond.t3}
    if(copXg[[j]]=="bvt4"){pcondcopXg[[j]]=pcond.t4;qcondcopXg[[j]]=qcond.t4}
    if(copXg[[j]]=="bvt5"){pcondcopXg[[j]]=pcond.t5;qcondcopXg[[j]]=qcond.t5}
    if(copXg[[j]]=="bvt6"){pcondcopXg[[j]]=pcond.t6;qcondcopXg[[j]]=qcond.t6}
    if(copXg[[j]]=="bvt7"){pcondcopXg[[j]]=pcond.t7;qcondcopXg[[j]]=qcond.t7}
    if(copXg[[j]]=="bvt8"){pcondcopXg[[j]]=pcond.t8;qcondcopXg[[j]]=qcond.t8}
    if(copXg[[j]]=="bvt9"){pcondcopXg[[j]]=pcond.t9;qcondcopXg[[j]]=qcond.t9}
    
    if(copXg[[j]]=="frk"){pcondcopXg[[j]]=pcond.frank;qcondcopXg[[j]]=qcond.frank}
    
    if(copXg[[j]]=="gum"){pcondcopXg[[j]]=pcond.gumbel;qcondcopXg[[j]]=qcond.gumbel}
    if(copXg[[j]]=="sgum"||copXg[[j]]=="rgum"){pcondcopXg[[j]]=pcond.sgumbel;qcondcopXg[[j]]=qcond.sgumbel}
    if(copXg[[j]]=="gum90"||copXg[[j]]=="1rgum"){pcondcopXg[[j]]=pcond.gumbel.90;qcondcopXg[[j]]=qcond.gumbel.90}
    if(copXg[[j]]=="gum270"||copXg[[j]]=="2rgum"){pcondcopXg[[j]]=pcond.gumbel.270;qcondcopXg[[j]]=qcond.gumbel.270}
    
    if(copXg[[j]]=="joe"){pcondcopXg[[j]]=pcond.joe;qcondcopXg[[j]]=qcond.joe}
    if(copXg[[j]]=="sjoe"||copXg[[j]]=="rjoe"){pcondcopXg[[j]]=pcond.sjoe;qcondcopXg[[j]]=qcond.sjoe}
    if(copXg[[j]]=="joe90"||copXg[[j]]=="1rjoe"){pcondcopXg[[j]]=pcond.joe.90;qcondcopXg[[j]]=qcond.joe.90}
    if(copXg[[j]]=="joe270"||copXg[[j]]=="2rjoe"){pcondcopXg[[j]]=pcond.joe.270;qcondcopXg[[j]]=qcond.joe.270}
  }
  
      
output=list(pcondcopX0=pcondcopX0,qcondcopX0=qcondcopX0,pcondcopXgY=pcondcopXg,qcondcopXgY=qcondcopXg)
return(output)
}

