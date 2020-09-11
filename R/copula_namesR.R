# Purpose: Providing the copula functions for the 
#          variables (continuous, ordinal. And counts if there are any).
#         
# Inputs
# d1   : dimensions for the continuous variables.
# d2   : dimensions for the ordinal variables.
# d3   : dimensions for the count variables.
# cop_name_f1: copula names for the first factor, starting 
#        with the continuous variables, ordinal and then the count variables.
# cop_name_f2: copula names for the second factor, starting 
#        with the continuous variables, ordinal and then the count variables.

# Output: Copula functions for the continuous, ordinal, and count variables
#         stored in lists, and to be used for the likelihood function.
copulas=function(d1=0,d2=0,d3=0,cop_name_f1,cop_name_f2=NULL){
  #number of total items
  d=d1+d2+d3
  #copula names for each factor
  cop_name<-list(first_factor=cop_name_f1,second_factor=cop_name_f2)
  cop_name<-Filter(Negate(is.null), cop_name)#negating the null.
  num_factors<-length(cop_name)#number of factors
  #Empty lists to save outputs
  dcop=list(f1=list(),f2=list())
  pcondcop=list(f1=list(),f2=list())
  qcondcop=list(f1=list(),f2=list())
  for(i in 1:num_factors){
    for(j in 1:d){
      if(cop_name[[i]][[j]]=="bvn"){
        dcop[[i]][[j]]=dbvn;pcondcop[[i]][[j]]=pcond.bvn;qcondcop[[i]][[j]]=qcond.bvn
      }
      if(cop_name[[i]][[j]]=="bvt1"){
        dcop[[i]][[j]]=dtcop1;pcondcop[[i]][[j]]=pcond.t1;qcondcop[[i]][[j]]=qcond.t1
      }
      if(cop_name[[i]][[j]]=="bvt2"){
        dcop[[i]][[j]]=dtcop2;pcondcop[[i]][[j]]=pcond.t2;qcondcop[[i]][[j]]=qcond.t2
      }
      if(cop_name[[i]][[j]]=="bvt3"){
        dcop[[i]][[j]]=dtcop3;pcondcop[[i]][[j]]=pcond.t3;qcondcop[[i]][[j]]=qcond.t3
      }
      if(cop_name[[i]][[j]]=="bvt4"){
        dcop[[i]][[j]]=dtcop4;pcondcop[[i]][[j]]=pcond.t4;qcondcop[[i]][[j]]=qcond.t4
      }
      if(cop_name[[i]][[j]]=="bvt5"){
        dcop[[i]][[j]]=dtcop5;pcondcop[[i]][[j]]=pcond.t5;qcondcop[[i]][[j]]=qcond.t5
      }
      if(cop_name[[i]][[j]]=="bvt6"){
        dcop[[i]][[j]]=dtcop6;pcondcop[[i]][[j]]=pcond.t6;qcondcop[[i]][[j]]=qcond.t6
      }
      if(cop_name[[i]][[j]]=="bvt7"){
        dcop[[i]][[j]]=dtcop7;pcondcop[[i]][[j]]=pcond.t7;qcondcop[[i]][[j]]=qcond.t7
      }
      if(cop_name[[i]][[j]]=="bvt8"){
        dcop[[i]][[j]]=dtcop8;pcondcop[[i]][[j]]=pcond.t8;qcondcop[[i]][[j]]=qcond.t8
      }
      if(cop_name[[i]][[j]]=="bvt9"){
        dcop[[i]][[j]]=dtcop9;pcondcop[[i]][[j]]=pcond.t9;qcondcop[[i]][[j]]=qcond.t9
      }
      if(cop_name[[i]][[j]]=="frk"){
        dcop[[i]][[j]]=dfrank;pcondcop[[i]][[j]]=pcond.frank;qcondcop[[i]][[j]]=qcond.frank
        }
      if(cop_name[[i]][[j]]=="gum"){
        dcop[[i]][[j]]=dgumbel;pcondcop[[i]][[j]]=pcond.gumbel;qcondcop[[i]][[j]]=qcond.gumbel
        }
      if(cop_name[[i]][[j]]=="rgum"){
        dcop[[i]][[j]]=dsgumbel;pcondcop[[i]][[j]]=pcond.sgumbel;qcondcop[[i]][[j]]=qcond.sgumbel
        }
      if(cop_name[[i]][[j]]=="1rgum"){
        dcop[[i]][[j]]=dgumbel.90;pcondcop[[i]][[j]]=pcond.gumbel.90;qcondcop[[i]][[j]]=qcond.gumbel.90
        }
      if(cop_name[[i]][[j]]=="2rgum"){
        dcop[[i]][[j]]=dgumbel.270;pcondcop[[i]][[j]]=pcond.gumbel.270;qcondcop[[i]][[j]]=qcond.gumbel.270
        }
      if(cop_name[[i]][[j]]=="joe"){
        dcop[[i]][[j]]=djoe;pcondcop[[i]][[j]]=pcond.joe;qcondcop[[i]][[j]]=qcond.joe
        }
      if(cop_name[[i]][[j]]=="rjoe"){
        dcop[[i]][[j]]=dsjoe;pcondcop[[i]][[j]]=pcond.sjoe;qcondcop[[i]][[j]]=qcond.sjoe
        }
      if(cop_name[[i]][[j]]=="1rjoe"){
        dcop[[i]][[j]]=djoe.90;pcondcop[[i]][[j]]=pcond.joe.90;qcondcop[[i]][[j]]=qcond.joe.90
        }
      if(cop_name[[i]][[j]]=="2rjoe"){
        dcop[[i]][[j]]=djoe.270;pcondcop[[i]][[j]]=pcond.joe.270;qcondcop[[i]][[j]]=qcond.joe.270
      }
      
      if(cop_name[[i]][[j]]=="bb1"){
        dcop[[i]][[j]]=dbb1;pcondcop[[i]][[j]]=pcond.bb1;qcondcop[[i]][[j]]=qcond.bb1
      }
      if(cop_name[[i]][[j]]=="rbb1"){
        dcop[[i]][[j]]=dsbb1;pcondcop[[i]][[j]]=pcond.sbb1;qcondcop[[i]][[j]]=qcond.sbb1
      }
      
      if(cop_name[[i]][[j]]=="bb7"){
        dcop[[i]][[j]]=dbb7;pcondcop[[i]][[j]]=pcond.bb7;qcondcop[[i]][[j]]=qcond.bb7
      }
      if(cop_name[[i]][[j]]=="rbb7"){
        dcop[[i]][[j]]=dsbb7;pcondcop[[i]][[j]]=pcond.sbb7;qcondcop[[i]][[j]]=qcond.sbb7
      }
      
      if(cop_name[[i]][[j]]=="bb8"){
        dcop[[i]][[j]]=dbb8;pcondcop[[i]][[j]]=pcond.bb8;qcondcop[[i]][[j]]=qcond.bb8
      }
      if(cop_name[[i]][[j]]=="rbb8"){
        dcop[[i]][[j]]=dsbb8;pcondcop[[i]][[j]]=pcond.sbb8;qcondcop[[i]][[j]]=qcond.sbb8
      }
      if(cop_name[[i]][[j]]=="bb10"){
        dcop[[i]][[j]]=dbb10;pcondcop[[i]][[j]]=pcond.bb10;qcondcop[[i]][[j]]=qcond.bb10
        }
      if(cop_name[[i]][[j]]=="rbb10"){
        dcop[[i]][[j]]=dsbb10;pcondcop[[i]][[j]]=pcond.sbb10;qcondcop[[i]][[j]]=qcond.sbb10
        }
    }
  }
  return(list(dcop=dcop, pcondcop=pcondcop, qcondcop=qcondcop))
}


