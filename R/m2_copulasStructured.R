M2copulas_structured=function(copX0,copXg){
  d0=length(copX0)
  dg=length(copXg)
  pcondcopX0=list()
  pcondcopXg=list()
  pcond2copX0=list()
  pcond2copXg=list()
  dcopX0=list()
  dcopXg=list()
  d2copX0=list()
  
  for(j in 1:d0){
    if(copX0[[j]]=="bvn"){pcondcopX0[[j]]=pcond.bvn;pcond2copX0[[j]]=pcond.2bvn;dcopX0[[j]]=dbvn;d2copX0[[j]]=d2bvn}
    if(copX0[[j]]=="bvt"){pcondcopX0[[j]]=pcond.t;pcond2copX0[[j]]=pcond.2t;dcopX0[[j]]=dtcop;d2copX0[[j]]=NULL}
    
    if(copX0[[j]]=="bvt1"){pcondcopX0[[j]]=pcond.t1;
    pcond2copX0[[j]]=pcond.2t1;dcopX0[[j]]=dtcop1;d2copX0[[j]]=NULL}
    
    if(copX0[[j]]=="bvt2"){pcondcopX0[[j]]=pcond.t2;
    pcond2copX0[[j]]=pcond.2t2;dcopX0[[j]]=dtcop2;d2copX0[[j]]=d2tcop2}
    
    if(copX0[[j]]=="bvt3"){pcondcopX0[[j]]=pcond.t3;
    pcond2copX0[[j]]=pcond.2t3;dcopX0[[j]]=dtcop3;d2copX0[[j]]=d2tcop3}
    if(copX0[[j]]=="bvt4"){pcondcopX0[[j]]=pcond.t4;
    pcond2copX0[[j]]=pcond.2t4;dcopX0[[j]]=dtcop4;d2copX0[[j]]=d2tcop4}
    
    if(copX0[[j]]=="bvt5"){pcondcopX0[[j]]=pcond.t5;
    pcond2copX0[[j]]=pcond.2t5;dcopX0[[j]]=dtcop5;d2copX0[[j]]=d2tcop5}
    if(copX0[[j]]=="bvt6"){pcondcopX0[[j]]=pcond.t6;
    pcond2copX0[[j]]=pcond.2t6;dcopX0[[j]]=dtcop6;d2copX0[[j]]=d2tcop6}
    
    if(copX0[[j]]=="bvt7"){pcondcopX0[[j]]=pcond.t7;
    pcond2copX0[[j]]=pcond.2t7;dcopX0[[j]]=dtcop7;d2copX0[[j]]=d2tcop7}
    if(copX0[[j]]=="bvt8"){pcondcopX0[[j]]=pcond.t8;
    pcond2copX0[[j]]=pcond.2t8;dcopX0[[j]]=dtcop8;d2copX0[[j]]=d2tcop8}
    if(copX0[[j]]=="bvt9"){pcondcopX0[[j]]=pcond.t9;
    pcond2copX0[[j]]=pcond.2t9;dcopX0[[j]]=dtcop9;d2copX0[[j]]=d2tcop9}
    
    
    if(copX0[[j]]=="frk"){pcondcopX0[[j]]=pcond.frank;pcond2copX0[[j]]=pcond.2frank;
    dcopX0[[j]]=dfrank;d2copX0[[j]]=d2frank}
    
    if(copX0[[j]]=="gum"){pcondcopX0[[j]]=pcond.gumbel;pcond2copX0[[j]]=pcond.2gumbel;
    dcopX0[[j]]=dgumbel;d2copX0[[j]]=d2gumbel}
    if(copX0[[j]]=="sgum"|copX0[[j]]=="rgum"){pcondcopX0[[j]]=pcond.sgumbel;pcond2copX0[[j]]=pcond.2sgumbel;
    dcopX0[[j]]=dsgumbel;d2copX0[[j]]=d2sgumbel}
    if(copX0[[j]]=="gum90"|copX0[[j]]=="1rgum"){pcondcopX0[[j]]=pcond.gumbel.90;pcond2copX0[[j]]=pcond.2gumbel.90;
    dcopX0[[j]]=dgumbel.90;d2copX0[[j]]=d2gumbel.90}
    if(copX0[[j]]=="gum270"|copX0[[j]]=="2rgum"){pcondcopX0[[j]]=pcond.gumbel.270;pcond2copX0[[j]]=pcond.2gumbel.270;
    dcopX0[[j]]=d2gumbel.270;d2copX0[[j]]=d2gumbel.270}
    
    if(copX0[[j]]=="joe"){pcondcopX0[[j]]=pcond.joe;pcond2copX0[[j]]=pcond.2joe;
    dcopX0[[j]]=djoe;d2copX0[[j]]=d2joe}
    if(copX0[[j]]=="sjoe"|copX0[[j]]=="rjoe"){pcondcopX0[[j]]=pcond.sjoe;pcond2copX0[[j]]=pcond.2sjoe;
    dcopX0[[j]]=dsjoe;d2copX0[[j]]=d2sjoe}
    
    if(copX0[[j]]=="joe90"|copX0[[j]]=="1rjoe"){pcondcopX0[[j]]=pcond.joe.90;pcond2copX0[[j]]=pcond.2joe.90;
    dcopX0[[j]]=djoe.90;d2copX0[[j]]=d2joe.90}
    if(copX0[[j]]=="joe270"|copX0[[j]]=="2rjoe"){pcondcopX0[[j]]=pcond.joe.270;pcond2copX0[[j]]=pcond.2joe.270;
    dcopX0[[j]]=djoe.270;d2copX0[[j]]=d2joe.270}
    
  }
  
  
  for(j in 1:dg){
    if(copXg[[j]]=="bvn"){pcondcopXg[[j]]=pcond.bvn;pcond2copXg[[j]]=pcond.2bvn;dcopXg[[j]]=dbvn}
    if(copXg[[j]]=="bvt"){pcondcopXg[[j]]=pcond.t;pcond2copXg[[j]]=pcond.2t;dcopXg[[j]]=dtcop}
    
    if(copXg[[j]]=="bvt1"){pcondcopXg[[j]]=pcond.t1;pcond2copXg[[j]]=pcond.2t1;dcopXg[[j]]=dtcop1}
    if(copXg[[j]]=="bvt2"){pcondcopXg[[j]]=pcond.t2;pcond2copXg[[j]]=pcond.2t2;dcopXg[[j]]=dtcop2}
    
    if(copXg[[j]]=="bvt3"){pcondcopXg[[j]]=pcond.t3;pcond2copXg[[j]]=pcond.2t3;dcopXg[[j]]=dtcop3}
    if(copXg[[j]]=="bvt4"){pcondcopXg[[j]]=pcond.t4;pcond2copXg[[j]]=pcond.2t4;dcopXg[[j]]=dtcop4}
    if(copXg[[j]]=="bvt5"){pcondcopXg[[j]]=pcond.t5;pcond2copXg[[j]]=pcond.2t5;dcopXg[[j]]=dtcop5}
    if(copXg[[j]]=="bvt6"){pcondcopXg[[j]]=pcond.t6;pcond2copXg[[j]]=pcond.2t6;dcopXg[[j]]=dtcop6}
    if(copXg[[j]]=="bvt7"){pcondcopXg[[j]]=pcond.t7;pcond2copXg[[j]]=pcond.2t7;dcopXg[[j]]=dtcop7}
    if(copXg[[j]]=="bvt8"){pcondcopXg[[j]]=pcond.t8;pcond2copXg[[j]]=pcond.2t8;dcopXg[[j]]=dtcop8}
    if(copXg[[j]]=="bvt9"){pcondcopXg[[j]]=pcond.t9;pcond2copXg[[j]]=pcond.2t9;dcopXg[[j]]=dtcop9}
    
    if(copXg[[j]]=="frk"){pcondcopXg[[j]]=pcond.frank;pcond2copXg[[j]]=pcond.2frank;dcopXg[[j]]=dfrank}
    
    if(copXg[[j]]=="gum"){pcondcopXg[[j]]=pcond.gumbel;pcond2copXg[[j]]=pcond.2gumbel;dcopXg[[j]]=dgumbel}
    if(copXg[[j]]=="sgum"|copXg[[j]]=="rgum"){pcondcopXg[[j]]=pcond.sgumbel;pcond2copXg[[j]]=pcond.2sgumbel;dcopXg[[j]]=dsgumbel}
    if(copXg[[j]]=="gum90"|copXg[[j]]=="1rgum"){pcondcopXg[[j]]=pcond.gumbel.90;pcond2copXg[[j]]=pcond.2gumbel.90;dcopXg[[j]]=dgumbel.90}
    if(copXg[[j]]=="gum270"|copXg[[j]]=="2rgum"){pcondcopXg[[j]]=pcond.gumbel.270;pcond2copXg[[j]]=pcond.2gumbel.270;dcopXg[[j]]=dgumbel.270}
    
    if(copXg[[j]]=="joe"){pcondcopXg[[j]]=pcond.joe;pcond2copXg[[j]]=pcond.2joe;dcopXg[[j]]=djoe}
    if(copXg[[j]]=="sjoe"|copXg[[j]]=="rjoe"){pcondcopXg[[j]]=pcond.sjoe;pcond2copXg[[j]]=pcond.2sjoe;dcopXg[[j]]=dsjoe}
    if(copXg[[j]]=="joe90"|copXg[[j]]=="1rjoe"){pcondcopXg[[j]]=pcond.joe.90;pcond2copXg[[j]]=pcond.2joe.90;dcopXg[[j]]=djoe.90}
    if(copXg[[j]]=="joe270"|copXg[[j]]=="2rjoe"){pcondcopXg[[j]]=pcond.joe.270;pcond2copXg[[j]]=pcond.2joe.270;dcopXg[[j]]=djoe.270}
    
  }
  
  output=list(pcondcopX0=pcondcopX0,pcond2copX0=pcond2copX0,dcopX0=dcopX0,d2copX0=d2copX0,
              pcondcopXg=pcondcopXg,pcond2copXg=pcond2copXg,dcopXg=dcopXg)
  return(output)
}

