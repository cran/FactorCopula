# Purpose: Provide the copula functions for the all
#          variables for the M2 statistic.
#         
# Inputs
# d   : number of all variables.
# cop_name1: copula names for the first factor; starting 
#        with the continuous variables; ordinal and then the count variables.
# cop_name2: copula names for the second factor; starting 
#        with the continuous variables; ordinal and then the count variables.

# Output: Copula functions for the continuous; ordinal; and count variables(if there are)
#         stored in lists; and to be used in the density function.
M2copulas=function(d,cop_name_f1,cop_name_f2=NULL){
  
  #copula names for each factor
  cop_name<-list(first_factor=cop_name_f1,second_factor=cop_name_f2)
  cop_name<-Filter(Negate(is.null), cop_name)#negating the null.
  num_factors<-length(cop_name)#number of factors
  
  #Empty lists to save outputs
  pcondcop=list(f1=list(),f2=list())
  pconddotcop=list(f1=list(),f2=list())
  dcop=list(f1=list(),f2=list())
  for(i in 1:num_factors){
  for(j in 1:d){
    if(cop_name[[i]][[j]]=="bvn"){dcop[[i]][[j]]=dbvn;pcondcop[[i]][[j]]=pcond.bvn;pconddotcop[[i]][[j]]=pcond.2bvn}
    
    if(cop_name[[i]][[j]]=="bvt1"){dcop[[i]][[j]]=dtcop1;pcondcop[[i]][[j]]=pcond.t1;pconddotcop[[i]][[j]]=pcond.2t1}
    if(cop_name[[i]][[j]]=="bvt2"){dcop[[i]][[j]]=dtcop2;pcondcop[[i]][[j]]=pcond.t2;pconddotcop[[i]][[j]]=pcond.2t2}
    if(cop_name[[i]][[j]]=="bvt3"){dcop[[i]][[j]]=dtcop3;pcondcop[[i]][[j]]=pcond.t3;pconddotcop[[i]][[j]]=pcond.2t3}
    if(cop_name[[i]][[j]]=="bvt4"){dcop[[i]][[j]]=dtcop4;pcondcop[[i]][[j]]=pcond.t4;pconddotcop[[i]][[j]]=pcond.2t4}
    if(cop_name[[i]][[j]]=="bvt5"){dcop[[i]][[j]]=dtcop5;pcondcop[[i]][[j]]=pcond.t5;pconddotcop[[i]][[j]]=pcond.2t5}
    if(cop_name[[i]][[j]]=="bvt6"){dcop[[i]][[j]]=dtcop6;pcondcop[[i]][[j]]=pcond.t6;pconddotcop[[i]][[j]]=pcond.2t6}
    if(cop_name[[i]][[j]]=="bvt7"){dcop[[i]][[j]]=dtcop7;pcondcop[[i]][[j]]=pcond.t7;pconddotcop[[i]][[j]]=pcond.2t7}
    if(cop_name[[i]][[j]]=="bvt8"){dcop[[i]][[j]]=dtcop8;pcondcop[[i]][[j]]=pcond.t8;pconddotcop[[i]][[j]]=pcond.2t8}
    if(cop_name[[i]][[j]]=="bvt9"){dcop[[i]][[j]]=dtcop9;pcondcop[[i]][[j]]=pcond.t9;pconddotcop[[i]][[j]]=pcond.2t9}
    
    if(cop_name[[i]][[j]]=="frk"){dcop[[i]][[j]]=dfrank;pcondcop[[i]][[j]]=pcond.frank;pconddotcop[[i]][[j]]=pcond.2frank}
    if(cop_name[[i]][[j]]=="gum"){dcop[[i]][[j]]=dgumbel;pcondcop[[i]][[j]]=pcond.gumbel;pconddotcop[[i]][[j]]=pcond.2gumbel}
    if(cop_name[[i]][[j]]=="rgum"){dcop[[i]][[j]]=dsgumbel;pcondcop[[i]][[j]]=pcond.sgumbel;pconddotcop[[i]][[j]]=pcond.2sgumbel}
    if(cop_name[[i]][[j]]=="1rgum"){dcop[[i]][[j]]=dgumbel.90;pcondcop[[i]][[j]]=pcond.gumbel.90;pconddotcop[[i]][[j]]=pcond.2gumbel.90}
    if(cop_name[[i]][[j]]=="2rgum"){dcop[[i]][[j]]=dgumbel.270;pcondcop[[i]][[j]]=pcond.gumbel.270;pconddotcop[[i]][[j]]=pcond.2gumbel.270}
    if(cop_name[[i]][[j]]=="joe"){dcop[[i]][[j]]=djoe;pcondcop[[i]][[j]]=pcond.joe;pconddotcop[[i]][[j]]=pcond.2joe}
    if(cop_name[[i]][[j]]=="rjoe"){dcop[[i]][[j]]=dsjoe;pcondcop[[i]][[j]]=pcond.sjoe;pconddotcop[[i]][[j]]=pcond.2sjoe}
    if(cop_name[[i]][[j]]=="1rjoe"){dcop[[i]][[j]]=djoe.90;pcondcop[[i]][[j]]=pcond.joe.90;pconddotcop[[i]][[j]]=pcond.2joe.90}
    if(cop_name[[i]][[j]]=="2rjoe"){dcop[[i]][[j]]=djoe.270;pcondcop[[i]][[j]]=pcond.joe.270;pconddotcop[[i]][[j]]=pcond.2joe.270}
    
    if(cop_name[[i]][[j]]=="bb1"){dcop[[i]][[j]]=dbb1;pcondcop[[i]][[j]]=pcond.bb1;pconddotcop[[i]][[j]]=c(pcond.2bb1.theta,pcond.2bb1.delta) }
    if(cop_name[[i]][[j]]=="rbb1"){dcop[[i]][[j]]=dsbb1;pcondcop[[i]][[j]]=pcond.sbb1;pconddotcop[[i]][[j]]=c(pcond.2sbb1.theta,pcond.2sbb1.delta) }
    
    if(cop_name[[i]][[j]]=="bb7"){dcop[[i]][[j]]=dbb7;pcondcop[[i]][[j]]=pcond.bb7;pconddotcop[[i]][[j]]=c(pcond.2bb7.theta,pcond.2bb7.delta) }
    if(cop_name[[i]][[j]]=="rbb7"){dcop[[i]][[j]]=dsbb7;pcondcop[[i]][[j]]=pcond.sbb7;pconddotcop[[i]][[j]]=c(pcond.2sbb7.theta,pcond.2sbb7.delta) }
    
    if(cop_name[[i]][[j]]=="bb8"){dcop[[i]][[j]]=dbb8;pcondcop[[i]][[j]]=pcond.bb8;pconddotcop[[i]][[j]]=c(pcond.2bb8.theta,pcond.2bb8.delta) }
    if(cop_name[[i]][[j]]=="rbb8"){dcop[[i]][[j]]=dsbb8;pcondcop[[i]][[j]]=pcond.sbb8;pconddotcop[[i]][[j]]=c(pcond.2sbb8.theta,pcond.2sbb8.delta) }
    
    if(cop_name[[i]][[j]]=="bb10"){dcop[[i]][[j]]=dbb10;pcondcop[[i]][[j]]=pcond.bb10;pconddotcop[[i]][[j]]=c(pcond.2bb10.theta,pcond.2bb10.delta) }
    if(cop_name[[i]][[j]]=="rbb10"){dcop[[i]][[j]]=dsbb10;pcondcop[[i]][[j]]=pcond.sbb10;pconddotcop[[i]][[j]]=c(pcond.2sbb10.theta,pcond.2sbb10.delta) }
  }
  }
    output=list(dcop=dcop,pcondcop=pcondcop,
                 pconddotcop=pconddotcop)
  return(output)
} 

