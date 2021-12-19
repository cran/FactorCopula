LUbound<- function(copnames) {
  d=length(copnames)
  LUB=rep(list(c(NA,NA)),d)
  
  for(j in 1:d){
    if(copnames[[j]]=="bvn"){LUB[[j]][1]= -0.95;LUB[[j]][2]=  0.95 
    }else
    if(copnames[[j]]=="bvt" |copnames[[j]]=="bvt1"| copnames[[j]]=="bvt2" | copnames[[j]]=="bvt3" |copnames[[j]]=="bvt4" |
       copnames[[j]]=="bvt5"|copnames[[j]]=="bvt6"|copnames[[j]]=="bvt7"|
       copnames[[j]]=="bvt8"|copnames[[j]]=="bvt9"){LUB[[j]][1]= -0.95;LUB[[j]][2]=  0.95 
    }else
    if(copnames[[j]]=="frk"){LUB[[j]][1]= -0.95;LUB[[j]][2]=  0.95 
    }else
    if(copnames[[j]]=="gum"|copnames[[j]]=="sgum"|copnames[[j]]=="rgum"){LUB[[j]][1]=  0.01;LUB[[j]][2]=  0.95  
    }else
    if(copnames[[j]]=="gum90"|copnames[[j]]=="gum270"||copnames[[j]]=="1rgum"|copnames[[j]]=="2rgum"){LUB[[j]][1]=-0.95; LUB[[j]][2]=  -0.01 }
    else
    if(copnames[[j]]=="joe"|copnames[[j]]=="sjoe"|copnames[[j]]=="rjoe"){LUB[[j]][1]=  0.01;LUB[[j]][2]=  0.95  }
    else
    if(copnames[[j]]=="joe90"|copnames[[j]]=="joe270"||copnames[[j]]=="1rjoe"|copnames[[j]]=="2rjoe"){LUB[[j]][1]=  -0.95;LUB[[j]][2]=  -0.01 }
  }
  
  output=c(lower.upper.B=LUB)
  return(output)
}
  
  