# Purpose: provide initial/starting values for all variables.
# Inputs
# cop_name_f1: copula names for the first factor, starting 
#        with the continuous variables, ordinal and then the count variables.
# cop_name_f2: copula names for the second factor, starting 
#        with the continuous variables, ordinal and then the count variables.

# Output: Starting values for the continuous, ordinal, and count variables
#         stored in lists, and to be used in the optimisation.
#       
initPar=function(cop_name_f1,cop_name_f2=NULL){
  d=length(cop_name_f1)
  initF1=list()
  initF2=list()
  for(j in 1:d){
    if(cop_name_f1[[j]]=="gum" | cop_name_f1[[j]]=="rgum" | cop_name_f1[[j]]=="1rgum" | cop_name_f1[[j]]=="2rgum" | 
       cop_name_f1[[j]]=="joe"| cop_name_f1[[j]]=="rjoe" | cop_name_f1[[j]]=="1rjoe" | cop_name_f1[[j]]=="2rjoe") {
      initF1[[j]] = linkfun("gum",1.1) 
    } else if( cop_name_f1[[j]]=="bvn" | cop_name_f1[[j]]=="bvt1"| 
               cop_name_f1[[j]]=="bvt2"| cop_name_f1[[j]]=="bvt3"| 
               cop_name_f1[[j]]=="bvt4"| cop_name_f1[[j]]=="bvt5"| 
               cop_name_f1[[j]]=="bvt6"| cop_name_f1[[j]]=="bvt7"| 
               cop_name_f1[[j]]=="bvt8"| cop_name_f1[[j]]=="bvt9") { 
      initF1[[j]] = linkfun("bvn",.1) 
    } else if( cop_name_f1[[j]]=="frk"){
      initF1[[j]] = linkfun("frk",1) 
    }else if( cop_name_f1[[j]]=="bb1" | cop_name_f1[[j]]=="rbb1") {
      initF1[[j]] = linkfun("bb1", c(.2,1.2) ) 
    }else if( cop_name_f1[[j]]=="bb7"| cop_name_f1[[j]]=="rbb7") {
      initF1[[j]] = linkfun("bb7", c(1.2,0.2) ) 
    } else if( cop_name_f1[[j]]=="bb8"| cop_name_f1[[j]]=="rbb8") {
      initF1[[j]] = linkfun("bb8", c(3,0.6) ) 
    } else if( cop_name_f1[[j]]=="bb10" | cop_name_f1[[j]]=="rbb10") {
      initF1[[j]] = linkfun("bb10", c(10,0.8) ) }
    if(is.null(cop_name_f2)) {
      output=initF1 
    } else {
      if( cop_name_f2[[j]]=="gum" | cop_name_f2[[j]]=="rgum" | cop_name_f2[[j]]=="1rgum" | cop_name_f2[[j]]=="2rgum" | 
          cop_name_f2[[j]]=="joe"| cop_name_f2[[j]]=="rjoe" | cop_name_f2[[j]]=="1rjoe" | cop_name_f2[[j]]=="2rjoe") {
        initF2[[j]] = linkfun("gum",1.1) 
      } else if( cop_name_f2[[j]]=="bvn" | cop_name_f2[[j]]=="bvt1"| 
                 cop_name_f2[[j]]=="bvt2"| cop_name_f2[[j]]=="bvt3"| 
                 cop_name_f2[[j]]=="bvt4"| cop_name_f2[[j]]=="bvt5"| 
                 cop_name_f2[[j]]=="bvt6"| cop_name_f2[[j]]=="bvt7"| 
                 cop_name_f2[[j]]=="bvt8"| cop_name_f2[[j]]=="bvt9") {
        initF2[[j]] = linkfun("bvn",.1) 
      } else if( cop_name_f2[[j]]=="frk") {
        initF2[[j]] = linkfun("frk",1) 
      }else if( cop_name_f2[[j]]=="bb1" | cop_name_f2[[j]]=="rbb1") {
        initF2[[j]] = linkfun("bb1", c(.2,1.2) ) 
      }else if( cop_name_f2[[j]]=="bb7"| cop_name_f2[[j]]=="rbb7") {
        initF2[[j]] = linkfun("bb7", c(1.2,0.2) ) 
      } else if( cop_name_f2[[j]]=="bb8"| cop_name_f2[[j]]=="rbb8") {
        initF2[[j]] = linkfun("bb8", c(3,0.6) ) 
      }else if( cop_name_f2[[j]]=="bb10" | cop_name_f2[[j]]=="rbb10") {
        initF2[[j]] = linkfun("bb10", c(10,0.8) ) 
      }
      output = list(f1=initF1,f2=initF2)
    }
  }
  return(output)
}
