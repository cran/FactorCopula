# Purpose: Providing the copula functions for the 
#          observed variables.
#         
# Inputs
# d     : number of total observed variables.
# copF1: copula names that link the observed variables with the 
#        first latent variable
# copF2: copula names that link the observed variables with the 
#        second latent variable
# copvine: copula names for the last tree.


# Output: Copula functions for all variables
#         stored in lists, and to be used in the pr. mass function.
copulas_factortree=function(copF1,copF2,copvine)
{
  d=length(copF1)
  if(any(copF2==0)){copF2=rep(0,d)}
  pcondcopF1=list()
  pcondcopF2=list()
  pcop=list()

  for(j in 1:d){
    if( copF1[[j]]=="bvn"){ pcondcopF1[[j]]=pcond.bvn}
    if( copF1[[j]]=="bvt1"){ pcondcopF1[[j]]=pcond.t1}
    if( copF1[[j]]=="bvt2"){ pcondcopF1[[j]]=pcond.t2}
    if( copF1[[j]]=="bvt3"){ pcondcopF1[[j]]=pcond.t3}
    if( copF1[[j]]=="bvt4"){ pcondcopF1[[j]]=pcond.t4}
    if( copF1[[j]]=="bvt5"){ pcondcopF1[[j]]=pcond.t5}
    if( copF1[[j]]=="bvt6"){ pcondcopF1[[j]]=pcond.t6}
    if( copF1[[j]]=="bvt7"){ pcondcopF1[[j]]=pcond.t7}
    if( copF1[[j]]=="bvt8"){ pcondcopF1[[j]]=pcond.t8}
    if( copF1[[j]]=="bvt9"){ pcondcopF1[[j]]=pcond.t9}
    if( copF1[[j]]=="frk"){ pcondcopF1[[j]]=pcond.frank}
    if( copF1[[j]]=="gum"){pcondcopF1[[j]]=pcond.gumbel}
    if( copF1[[j]]=="sgum"|| copF1[[j]]=="rgum"){ pcondcopF1[[j]]=pcond.sgumbel}
    if( copF1[[j]]=="gum90" ||copF1[[j]]=="1rgum"){ pcondcopF1[[j]]=pcond.gumbel.90}
    if( copF1[[j]]=="gum270" ||copF1[[j]]=="2rgum"){ pcondcopF1[[j]]=pcond.gumbel.270}
  }
  
  for(j in 1:d){
    if( copF2[[j]]=="bvn"){ pcondcopF2[[j]]=pcond.bvn}
    if( copF2[[j]]=="bvt1"){ pcondcopF2[[j]]=pcond.t1}
    if( copF2[[j]]=="bvt2"){ pcondcopF2[[j]]=pcond.t2}
    if( copF2[[j]]=="bvt3"){ pcondcopF2[[j]]=pcond.t3}
    if( copF2[[j]]=="bvt4"){ pcondcopF2[[j]]=pcond.t4}
    if( copF2[[j]]=="bvt5"){ pcondcopF2[[j]]=pcond.t5}
    if( copF2[[j]]=="bvt6"){ pcondcopF2[[j]]=pcond.t6}
    if( copF2[[j]]=="bvt7"){ pcondcopF2[[j]]=pcond.t7}
    if( copF2[[j]]=="bvt8"){ pcondcopF2[[j]]=pcond.t8}
    if( copF2[[j]]=="bvt9"){ pcondcopF2[[j]]=pcond.t9}
    if( copF2[[j]]=="frk"){ pcondcopF2[[j]]=pcond.frank}
    if( copF2[[j]]=="gum"){pcondcopF2[[j]]=pcond.gumbel}
    if( copF2[[j]]=="sgum"|| copF2[[j]]=="rgum"){ pcondcopF2[[j]]=pcond.sgumbel}
    if( copF2[[j]]=="gum90" ||copF2[[j]]=="1rgum"){ pcondcopF2[[j]]=pcond.gumbel.90}
    if( copF2[[j]]=="gum270" ||copF2[[j]]=="2rgum"){ pcondcopF2[[j]]=pcond.gumbel.270}
    if( copF2[[j]]==0){pcondcopF2[[j]]=NA}
  }
  
  for(j in 1:(d-1)){
    if(copvine[[j]]=="bvn"){pcop[[j]]=pbvncop}
    if(copvine[[j]]=="bvt1"){pcop[[j]]=pbvtcop1}
    if(copvine[[j]]=="bvt2"){pcop[[j]]=pbvtcop2}
    if(copvine[[j]]=="bvt3"){pcop[[j]]=pbvtcop3}
    if(copvine[[j]]=="bvt4"){pcop[[j]]=pbvtcop4}
    if(copvine[[j]]=="bvt5"){pcop[[j]]=pbvtcop5}
    if(copvine[[j]]=="bvt6"){pcop[[j]]=pbvtcop6}
    if(copvine[[j]]=="bvt7"){pcop[[j]]=pbvtcop7}
    if(copvine[[j]]=="bvt8"){pcop[[j]]=pbvtcop8}
    if(copvine[[j]]=="bvt9"){pcop[[j]]=pbvtcop9}
    if(copvine[[j]]=="frk"){pcop[[j]]=pfrank}
    if(copvine[[j]]=="gum"){pcop[[j]]=pgumbel}
    if(copvine[[j]]=="sgum"||copvine[[j]]=="rgum"){pcop[[j]]=psgumbel}
    if(copvine[[j]]=="gum90"||copvine[[j]]=="1rgum"){pcop[[j]]=pgumbel.90}
    if(copvine[[j]]=="gum270"||copvine[[j]]=="2rgum"){pcop[[j]]=pgumbel.270}
  }
  
output=list(pcondF1=pcondcopF1, pcondF2=pcondcopF2, copvine=pcop)
return(output)
}
