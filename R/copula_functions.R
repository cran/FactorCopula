#Conditional cdf for Gaussian copula
pcond.bvn=function(v,u,theta,param=F){
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  cond.cdf<- pnorm((qnorm(v)-cpar*qnorm(u))/sqrt(1-cpar^2))
  cond.cdf[ v <= 0 | u <= 0 |  u>= 1] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}

# Gaussain density
dbvn=function(u,v,theta,param=F){     
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  
  x1=qnorm(u); x2=qnorm(v)
  qf=x1^2+x2^2-2*cpar*x1*x2
  qf=qf/(1-cpar^2)
  con=sqrt(1-cpar^2)*(2*pi)
  pdf=exp(-.5*qf)/con
  pdf=pdf/(dnorm(x1)*dnorm(x2))
  pdf
}
#=======
#Conditional cdf for Student-t copula
#Conditional cdf for Student-t copula
pcond.t = function(v,u,theta,param=F){
  
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  #  df=1
  tem = qt(cbind(u,v),df)
  
  x=tem[,1]
  y=tem[,2]
  
  term1=(df+x*x)*(1-cpar*cpar)/(df+1)
  term2 = sqrt(term1)
  quant = (y-cpar*x)/term2
  cond.cdf = pt(quant,df+1)
  cond.cdf[ v <= 0 | u <= 0 |  u>= 1] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}


pcond.t1 = function(v,u,theta,param=F){
  
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=1
  tem = qt(cbind(u,v),df)
  
  x=tem[,1]
  y=tem[,2]
  
  term1=(df+x*x)*(1-cpar*cpar)/(df+1)
  term2 = sqrt(term1)
  quant = (y-cpar*x)/term2
  cond.cdf = pt(quant,df+1)
  cond.cdf[ v <= 0 | u <= 0 |  u>= 1] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}

pcond.t2 = function(v,u,theta,param=F){
  
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=2
  tem = qt(cbind(u,v),df)
  
  x=tem[,1]
  y=tem[,2]
  
  term1=(df+x*x)*(1-cpar*cpar)/(df+1)
  term2 = sqrt(term1)
  quant = (y-cpar*x)/term2
  cond.cdf = pt(quant,df+1)
  cond.cdf[ v <= 0 | u <= 0 |  u>= 1] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}

pcond.t3 = function(v,u,theta,param=F){
  
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=3
  tem = qt(cbind(u,v),df)
  
  x=tem[,1]
  y=tem[,2]
  
  term1=(df+x*x)*(1-cpar*cpar)/(df+1)
  term2 = sqrt(term1)
  quant = (y-cpar*x)/term2
  cond.cdf = pt(quant,df+1)
  cond.cdf[ v <= 0 | u <= 0 |  u>= 1] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}

pcond.t4 = function(v,u,theta,param=F){
  
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=4
  tem = qt(cbind(u,v),df)
  
  x=tem[,1]
  y=tem[,2]
  
  term1=(df+x*x)*(1-cpar*cpar)/(df+1)
  term2 = sqrt(term1)
  quant = (y-cpar*x)/term2
  cond.cdf = pt(quant,df+1)
  cond.cdf[ v <= 0 | u <= 0 |  u>= 1] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}

pcond.t5 = function(v,u,theta,param=F){
  
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=5
  tem = qt(cbind(u,v),df)
  
  x=tem[,1]
  y=tem[,2]
  
  term1=(df+x*x)*(1-cpar*cpar)/(df+1)
  term2 = sqrt(term1)
  quant = (y-cpar*x)/term2
  cond.cdf = pt(quant,df+1)
  cond.cdf[ v <= 0 | u <= 0 |  u>= 1] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}

pcond.t6 = function(v,u,theta,param=F){
  
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=6
  tem = qt(cbind(u,v),df)
  
  x=tem[,1]
  y=tem[,2]
  
  term1=(df+x*x)*(1-cpar*cpar)/(df+1)
  term2 = sqrt(term1)
  quant = (y-cpar*x)/term2
  cond.cdf = pt(quant,df+1)
  cond.cdf[ v <= 0 | u <= 0 |  u>= 1] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}


pcond.t7 = function(v,u,theta,param=F){
  
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=7
  tem = qt(cbind(u,v),df)
  
  x=tem[,1]
  y=tem[,2]
  
  term1=(df+x*x)*(1-cpar*cpar)/(df+1)
  term2 = sqrt(term1)
  quant = (y-cpar*x)/term2
  cond.cdf = pt(quant,df+1)
  cond.cdf[ v <= 0 | u <= 0 |  u>= 1] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}

pcond.t8 = function(v,u,theta,param=F){
  
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=8
  tem = qt(cbind(u,v),df)
  
  x=tem[,1]
  y=tem[,2]
  
  term1=(df+x*x)*(1-cpar*cpar)/(df+1)
  term2 = sqrt(term1)
  quant = (y-cpar*x)/term2
  cond.cdf = pt(quant,df+1)
  cond.cdf[ v <= 0 | u <= 0 |  u>= 1] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}

pcond.t9 = function(v,u,theta,param=F){
  
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=9
  tem = qt(cbind(u,v),df)
  
  x=tem[,1]
  y=tem[,2]
  
  term1=(df+x*x)*(1-cpar*cpar)/(df+1)
  term2 = sqrt(term1)
  quant = (y-cpar*x)/term2
  cond.cdf = pt(quant,df+1)
  cond.cdf[ v <= 0 | u <= 0 |  u>= 1] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}


# the t copula density
dtcop=function(u,v,theta,param=F){    
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  #df=1
  tem1=cbind(u,v)
  tem2=qt(tem1,df)
  tem3=dt(tem2,df)
  x<-tem2[,1]
  y<-tem2[,2]
  z<-tem3[,1]
  w<-tem3[,2]
  x2=x*x
  y2=y*y
  cpar2=cpar*cpar
  val<-df/2*(1+(x2+y2-2*cpar*x*y)/(df*(1-cpar2)))^(-(df+2)/2)/(df*pi*z*w*sqrt(1-cpar2))
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}


dtcop1=function(u,v,theta,param=F){    
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=1
  tem1=cbind(u,v)
  tem2=qt(tem1,df)
  tem3=dt(tem2,df)
  x<-tem2[,1]
  y<-tem2[,2]
  z<-tem3[,1]
  w<-tem3[,2]
  x2=x*x
  y2=y*y
  cpar2=cpar*cpar
  val<-df/2*(1+(x2+y2-2*cpar*x*y)/(df*(1-cpar2)))^(-(df+2)/2)/(df*pi*z*w*sqrt(1-cpar2))
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

dtcop2=function(u,v,theta,param=F){    
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=2
  tem1=cbind(u,v)
  tem2=qt(tem1,df)
  tem3=dt(tem2,df)
  x<-tem2[,1]
  y<-tem2[,2]
  z<-tem3[,1]
  w<-tem3[,2]
  x2=x*x
  y2=y*y
  cpar2=cpar*cpar
  val<-df/2*(1+(x2+y2-2*cpar*x*y)/(df*(1-cpar2)))^(-(df+2)/2)/(df*pi*z*w*sqrt(1-cpar2))
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

dtcop3=function(u,v,theta,param=F){    
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=3
  tem1=cbind(u,v)
  tem2=qt(tem1,df)
  tem3=dt(tem2,df)
  x<-tem2[,1]
  y<-tem2[,2]
  z<-tem3[,1]
  w<-tem3[,2]
  x2=x*x
  y2=y*y
  cpar2=cpar*cpar
  val<-df/2*(1+(x2+y2-2*cpar*x*y)/(df*(1-cpar2)))^(-(df+2)/2)/(df*pi*z*w*sqrt(1-cpar2))
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

dtcop4=function(u,v,theta,param=F){    
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=4
  tem1=cbind(u,v)
  tem2=qt(tem1,df)
  tem3=dt(tem2,df)
  x<-tem2[,1]
  y<-tem2[,2]
  z<-tem3[,1]
  w<-tem3[,2]
  x2=x*x
  y2=y*y
  cpar2=cpar*cpar
  val<-df/2*(1+(x2+y2-2*cpar*x*y)/(df*(1-cpar2)))^(-(df+2)/2)/(df*pi*z*w*sqrt(1-cpar2))
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

dtcop5=function(u,v,theta,param=F){    
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=5
  tem1=cbind(u,v)
  tem2=qt(tem1,df)
  tem3=dt(tem2,df)
  x<-tem2[,1]
  y<-tem2[,2]
  z<-tem3[,1]
  w<-tem3[,2]
  x2=x*x
  y2=y*y
  cpar2=cpar*cpar
  val<-df/2*(1+(x2+y2-2*cpar*x*y)/(df*(1-cpar2)))^(-(df+2)/2)/(df*pi*z*w*sqrt(1-cpar2))
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

dtcop6=function(u,v,theta,param=F){    
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=6
  tem1=cbind(u,v)
  tem2=qt(tem1,df)
  tem3=dt(tem2,df)
  x<-tem2[,1]
  y<-tem2[,2]
  z<-tem3[,1]
  w<-tem3[,2]
  x2=x*x
  y2=y*y
  cpar2=cpar*cpar
  val<-df/2*(1+(x2+y2-2*cpar*x*y)/(df*(1-cpar2)))^(-(df+2)/2)/(df*pi*z*w*sqrt(1-cpar2))
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}


dtcop7=function(u,v,theta,param=F){    
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=7
  tem1=cbind(u,v)
  tem2=qt(tem1,df)
  tem3=dt(tem2,df)
  x<-tem2[,1]
  y<-tem2[,2]
  z<-tem3[,1]
  w<-tem3[,2]
  x2=x*x
  y2=y*y
  cpar2=cpar*cpar
  val<-df/2*(1+(x2+y2-2*cpar*x*y)/(df*(1-cpar2)))^(-(df+2)/2)/(df*pi*z*w*sqrt(1-cpar2))
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

dtcop8=function(u,v,theta,param=F){    
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=8
  tem1=cbind(u,v)
  tem2=qt(tem1,df)
  tem3=dt(tem2,df)
  x<-tem2[,1]
  y<-tem2[,2]
  z<-tem3[,1]
  w<-tem3[,2]
  x2=x*x
  y2=y*y
  cpar2=cpar*cpar
  val<-df/2*(1+(x2+y2-2*cpar*x*y)/(df*(1-cpar2)))^(-(df+2)/2)/(df*pi*z*w*sqrt(1-cpar2))
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

dtcop9=function(u,v,theta,param=F){    
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if(param){cpar=((exp(theta)-1)/(exp(theta)+1))}else{cpar=theta}
  df=9
  tem1=cbind(u,v)
  tem2=qt(tem1,df)
  tem3=dt(tem2,df)
  x<-tem2[,1]
  y<-tem2[,2]
  z<-tem3[,1]
  w<-tem3[,2]
  x2=x*x
  y2=y*y
  cpar2=cpar*cpar
  val<-df/2*(1+(x2+y2-2*cpar*x*y)/(df*(1-cpar2)))^(-(df+2)/2)/(df*pi*z*w*sqrt(1-cpar2))
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

#=======

#Conditional cdf for Gumbel copula
pcond.gumbel=function ( v,u, theta , param=F){
  
  u[u>=1]=.9999999999999
  v[v>=1]=.9999999999999
  u[u<=0]=.0000000000001
  v[v<=0]=.0000000000001
  
  if(param){cpar=exp(theta)+1}else{cpar=theta}
  
  b = - log(u)
  z = - log(v)
  
  term1 = b^cpar
  term2 = z^cpar
  
  summation = term1 + term2
  
  inner.part = summation^(1/cpar)
  cond.cdf = exp(- inner.part)
  cond.cdf = cond.cdf * (1 + term2/term1)^( (1/cpar) - 1)
  cond.cdf = cond.cdf/u
  cond.cdf
}

# Gumbel density
dgumbel=function(u,v,theta,param=F){
  u[u>=1]=.9999999999999
  v[v>=1]=.9999999999999
  u[u<=0]=.0000000000001
  v[v<=0]=.0000000000001
  
  if(param){cpar=exp(theta)+1}else{cpar=theta}
  
  l1= -log(u); l2= -log(v);
  tem1=l1^cpar; tem2=l2^cpar; sm=tem1+tem2; tem=sm^(1./cpar);
  cdf=exp(-tem);
  pdf=cdf*tem*tem1*tem2*(tem+cpar-1.);
  pdf=pdf/(sm*sm*l1*l2*u*v);
  pdf
  
}
#=======
pcond.frank=function ( v,u, cpar,param=F){
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  cpar[cpar == 0] = 1e-10
  cpar1 = 1 - exp(-cpar)
  tem = 1 - exp(-cpar * u)
  ccdf = (1 - tem)/(cpar1/(1 - exp(-cpar * v)) - tem)
  ccdf
}

dfrank=function(u,v,cpar,param=F){     
  u[u==1]=.9999999
  v[v==1]=.9999999
  u[u==0]=.0000001
  v[v==0]=.0000001
  
  t1=1.-exp(-cpar);
  tem1=exp(-cpar*u); tem2=exp(-cpar*v);
  pdf=cpar*tem1*tem2*t1;
  tem=t1-(1.-tem1)*(1.-tem2);
  pdf=pdf/(tem*tem);
  pdf
}
#========

pcond.joe=function (v, u, theta,param=F){    
  u[u==1]=.9999999
  v[v==1]=.9999999
  u[u==0]=.0000001
  v[v==0]=.0000001
  
  if(param){cpar=exp(theta)+1}else{cpar=theta}
  
  temv = (1 - v)^cpar
  temu = (1 - u)^cpar
  ccdf = 1 + temv/temu - temv
  ccdf = ccdf^(-1 + 1/cpar)
  ccdf = ccdf * (1 - temv)
  ccdf
}


djoe=function (u, v, theta,param=F){ 
  u[u==1]=.9999999
  v[v==1]=.9999999
  u[u==0]=.0000001
  v[v==0]=.0000001
  
  if(param){cpar=exp(theta)+1}else{cpar=theta}
  
  f1 = 1 - u
  f2 = 1 - v
  tem1 = f1^cpar
  tem2 = f2^cpar
  sm = tem1 + tem2 - tem1 * tem2
  tem = sm^(1/cpar)
  pdf = tem * ((cpar - 1) * tem1 * tem2 + tem1 * tem1 * tem2 +
                 tem1 * tem2 * tem2 - tem1 * tem1 * tem2 * tem2)
  pdf = pdf/(sm * sm)
  pdf = pdf/(f1 * f2)
  pdf
}
#===========================================================
#===========================================================
#======== Rotated versions by 180 (Survival)
#===========================================================
#===========================================================
# Survival Gumbel
dsgumbel=function(u,v,theta,param=F){
  u=1-u; v=1-v;
  
  if(param){cpar=exp(theta)+1}else{cpar=theta}
  
  l1= -log(u); l2= -log(v);
  tem1=l1^cpar; tem2=l2^cpar; sm=tem1+tem2; tem=sm^(1/cpar);
  cdf=exp(-tem);
  pdf=cdf*tem*tem1*tem2*(tem+cpar-1.);
  pdf=pdf/(sm*sm*l1*l2*u*v);
  pdf[ u <= 0 | v <= 0 ] <- 0
  pdf[ u >= 1 | v >= 1] <- 0
  pdf
}

pcond.sgumbel=function(v,u,theta,param=F){
  u[u==1]=.9999999
  v[v==1]=.9999999
  u[u==0]=.0000001
  v[v==0]=.0000001
  u1=1-u
  v1=1-v
  cond.cdf=1-pcond.gumbel(v1,u1,theta,param)
  cond.cdf[ v <= 0| u <= 0 | u >= 1 ] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}

# Survival Joe
dsjoe=function (u, v, theta,param=F){ 
  if(param){cpar=exp(theta)+1}else{cpar=theta}
  u=1-u
  v=1-v
  f1 = 1 - u
  f2 = 1 - v
  tem1 = f1^cpar
  tem2 = f2^cpar
  sm = tem1 + tem2 - tem1 * tem2
  tem = sm^(1/cpar)
  pdf = tem * ((cpar - 1) * tem1 * tem2 + tem1 * tem1 * tem2 +
                 tem1 * tem2 * tem2 - tem1 * tem1 * tem2 * tem2)
  pdf = pdf/(sm * sm)
  pdf = pdf/(f1 * f2)
  pdf[ u <= 0 | v <= 0 ] <- 0
  pdf[ u >= 1 | v >= 1] <- 0
  pdf
}

pcond.sjoe=function (v, u, theta,param=F){   
  u[u==1]=.9999999
  v[v==1]=.9999999
  u[u==0]=.0000001
  v[v==0]=.0000001
  
  if(param){cpar=exp(theta)+1}else{cpar=theta}
  
  u=1-u
  v=1-v
  temv = (1 - v)^cpar
  temu = (1 - u)^cpar
  ccdf = 1 + temv/temu - temv
  ccdf = ccdf^(-1 + 1/cpar)
  ccdf = ccdf * (1 - temv)
  
  ccdf=1-ccdf
  ccdf[ v <= 0 | u <= 0 | u >= 1] <- 0
  ccdf[ v == 1 ] <- 1
  ccdf
}
#===========================================================
#===========================================================
#======== Rotated versions by 90 and 270
#===========================================================
#===========================================================

#======== Gumbel.90 Gumbel.270
pcond.gumbel.90=function(v,u,theta,param=F){
  u[u==1]=.9999999
  v[v==1]=.9999999
  u[u==0]=.0000001
  v[v==0]=.0000001
  v1=1-v
  cond.cdf= 1-pcond.gumbel(v1, u, theta, param)
  cond.cdf[ v <= 0 | u <= 0 | u >= 1] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}


pcond.gumbel.270=function(v,u,theta,param=F){
  u[u==1]=.9999999
  v[v==1]=.9999999
  u[u==0]=.0000001
  v[v==0]=.0000001
  u1=1-u
  cond.cdf=pcond.gumbel(v, u1, theta, param)
  cond.cdf[ v <= 0 | u <= 0 | u >= 1] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}

dgumbel.90=function(u,v,theta,param=F)
{ 
  u = 1 - u
  val<-  dgumbel(u,v,theta,param)
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

dgumbel.270=function(u,v,theta,param=F){
  
  v = 1 - v
  val<- dgumbel(u,v,theta,param)
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

#======== Joe.90 Joe.270
pcond.joe.90=function (v,u, theta,param=F){  
  u[u==1]=.9999999
  v[v==1]=.9999999
  u[u==0]=.0000001
  v[v==0]=.0000001
  v=1-v
  if(param){cpar=exp(theta)+1}else{cpar=theta}
  
  temv = (1 - v)^cpar
  temu = (1 - u)^cpar
  ccdf = 1 + temv/temu - temv
  ccdf = ccdf^(-1 + 1/cpar)
  ccdf = ccdf * (1 - temv)
  ccdf[ v <= 0 | u <= 0 | u >= 1] <- 0
  ccdf[ v == 1 ] <- 1
  1-ccdf
}

pcond.joe.270=function (v,u, theta,param=F){  
  u[u==1]=.9999999
  v[v==1]=.9999999
  u[u==0]=.0000001
  v[v==0]=.0000001
  u=1-u
  if(param){cpar=exp(theta)+1}else{cpar=theta}
  
  temv = (1 - v)^cpar
  temu = (1 - u)^cpar
  ccdf = 1 + temv/temu - temv
  ccdf = ccdf^(-1 + 1/cpar)
  ccdf = ccdf * (1 - temv)
  ccdf[ v <= 0 | u <= 0 | u >= 1] <- 0
  ccdf[ v == 1 ] <- 1
  ccdf
}

djoe.90=function (u, v, theta,param=F){     
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  u=1-u
  f1 = 1 - u
  f2 = 1 - v
  
  if(param){cpar=exp(theta)+1}else{cpar=theta}
  
  tem1 = f1^cpar
  tem2 = f2^cpar
  sm = tem1 + tem2 - tem1 * tem2
  tem = sm^(1/cpar)
  pdf = tem * ((cpar - 1) * tem1 * tem2 + tem1 * tem1 * tem2 +
                 tem1 * tem2 * tem2 - tem1 * tem1 * tem2 * tem2)
  pdf = pdf/(sm * sm)
  pdf = pdf/(f1 * f2)
  pdf[ u <= 0 | v <= 0 ] <- 0
  pdf[ u >= 1 | v >= 1] <- 0
  pdf
}

djoe.270=function (u, v, theta,param=F){     
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  v=1-v
  f1 = 1 - u
  f2 = 1 - v
  
  if(param){cpar=exp(theta)+1}else{cpar=theta}
  
  tem1 = f1^cpar
  tem2 = f2^cpar
  sm = tem1 + tem2 - tem1 * tem2
  tem = sm^(1/cpar)
  pdf = tem * ((cpar - 1) * tem1 * tem2 + tem1 * tem1 * tem2 +
                 tem1 * tem2 * tem2 - tem1 * tem1 * tem2 * tem2)
  pdf = pdf/(sm * sm)
  pdf = pdf/(f1 * f2)
  pdf[ u <= 0 | v <= 0 ] <- 0
  pdf[ u >= 1 | v >= 1] <- 0
  pdf
}

#BB1
#Conditional cdf for bb1 copula
pcond.bb1 =function (v,u, cpar,param=F)
{    
  u[u==1]=.9999999
  v[v==1]=.9999999
  u[u==0]=.0000001
  v[v==0]=.0000001
  
  if (is.matrix(cpar)) {
    theta = cpar[, 1]
    delta = cpar[, 2]
  }
  else {
    theta = cpar[1]
    delta = cpar[2]
  }
  
  if(param){th=exp(theta);de=exp(delta)+1}else{th=theta;de=delta}
  
  de1 = 1/de
  th1 = 1/th
  ut = (u^(-th) - 1)
  vt = (v^(-th) - 1)
  x = ut^de
  y = vt^de
  sm = x + y
  smd = sm^(de1)
  tem = (1 + smd)^(-th1 - 1)
  ccdf = tem * smd * x * (ut + 1)/sm/ut/u
  ccdf[ v <= 0 | u <= 0 | u >= 1 ] <- 0
  ccdf[ v == 1 ] <- 1
  ccdf
}


pcond.sbb1=function(v,u,theta,param=F){
  u[u==1]=.9999999
  v[v==1]=.9999999
  u[u==0]=.0000001
  v[v==0]=.0000001
  u1=1-u
  v1=1-v
  cond.cdf=1-pcond.bb1(v1,u1,theta,param)
  cond.cdf[ v <= 0| u <= 0 | u >= 1 ] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}

dbb1=function (u, v, cpar,param=F){
  if (is.matrix(cpar)) {
    theta = cpar[, 1]
    delta = cpar[, 2]
  }
  else {
    theta = cpar[1]
    delta = cpar[2]
  }
  
  if(param){th=exp(theta);de=exp(delta)+1}else{th=theta;de=delta}
  
  de1 = 1/de
  th1 = 1/th
  ut = (u^(-th) - 1)
  vt = (v^(-th) - 1)
  x = ut^de
  y = vt^de
  sm = x + y
  smd = sm^(de1)
  tem = (1 + smd)^(-th1 - 2) * (th * (de - 1) + (th * de + 1) * smd)
  pdf = tem * smd * x * y * (ut + 1) * (vt + 1)/sm/sm/ut/vt/u/v
  pdf[ u <= 0 | v <= 0 ] <- 0
  pdf[ u == 1 | v >= 1] <- 0
  pdf
}

dsbb1=function(u, v, cpar,param=F){
  u1=1-u
  v1=1-v
  val=dbb1(u1, v1, cpar,param)
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
} 

#Conditional cdf for bb7 copula
pcond.bb7=function (v, u, cpar,param=F){
  if (is.matrix(cpar)) {
    theta = cpar[, 1]
    delta = cpar[, 2]
  }
  else {
    theta = cpar[1]
    delta = cpar[2]
  }
  
  if(param){th=exp(theta)+1;de=exp(delta)}else{th=theta;de=delta}
  
  de1 = 1/de
  th1 = 1/th
  ut = 1 - (1 - u)^th
  vt = 1 - (1 - v)^th
  x = ut^(-de) - 1
  y = vt^(-de) - 1
  sm = x + y + 1
  smd = sm^(-de1)
  tem = (1 - smd)^(th1 - 1)
  ccdf = tem * smd * (x + 1) * (1 - ut)/sm/ut/(1 - u)
  ccdf[ v <= 0 | u <= 0 | u >= 1 ] <- 0
  ccdf[ v == 1 ] <- 1
  ccdf
  
}

pcond.sbb7=function(v,u,theta,param=F){
  u[u==1]=.9999999
  v[v==1]=.9999999
  u[u==0]=.0000001
  v[v==0]=.0000001
  u1=1-u
  v1=1-v
  cond.cdf=1-pcond.bb7(v1,u1,theta,param)
  cond.cdf[ v <= 0| u <= 0 | u >= 1 ] <- 0
  cond.cdf[ v == 1 ] <- 1
  cond.cdf
}

dbb7=function (u, v, cpar,param=F){
  if (is.matrix(cpar)) {
    theta = cpar[, 1]
    delta = cpar[, 2]
  }
  else {
    theta = cpar[1]
    delta = cpar[2]
  }
  
  if(param){th=exp(theta)+1;de=exp(delta)}else{th=theta;de=delta}
  
  de1 = 1/de
  th1 = 1/th
  ut = 1 - (1 - u)^th
  vt = 1 - (1 - v)^th
  x = ut^(-de) - 1
  y = vt^(-de) - 1
  sm = x + y + 1
  smd = sm^(-de1)
  tem = (1 - smd)^(th1 - 2)
  tem = tem * (th * (de + 1) - (th * de + 1) * smd)
  pdf = tem * smd * (x + 1) * (y + 1) * (1 - ut) * (1 - vt)/sm/sm/ut/vt/(1 -u)/(1 - v)
  pdf[ u <= 0 | v <= 0 ] <- 0
  pdf[ u >= 1 | v >= 1] <- 0
  pdf
}

dsbb7=function(u, v, cpar,param=F){
  u1=1-u
  v1=1-v
  val=dbb7(u1, v1, cpar,param)
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
} 
# BB8
pcond.bb8=function (v, u, cpar,param=F)
{
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if (is.matrix(cpar)) {
    theta = cpar[, 1]
    delta = cpar[, 2]
  }
  else {
    theta = cpar[1]
    delta = cpar[2]
  }
  
  if(param){th=exp(theta)+1;de=(exp(delta))/(exp(delta)+1)}else{th=theta;de=delta}
  
  ut = (1 - de * u)^th
  x = 1 - ut
  y = 1 - (1 - de * v)^th
  eta1 = 1/(1 - (1 - de)^th)
  tem = (1 - eta1 * x * y)^(1/th - 1)
  den = (1 - de * u)/ut
  ccdf = eta1 * y * tem/den
  ccdf
}





dbb8=function (u, v, cpar,param=F)
{
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  if (is.matrix(cpar)) {
    theta = cpar[, 1]
    delta = cpar[, 2]
  }
  else {
    theta = cpar[1]
    delta = cpar[2]
  }
  
  if(param){th=exp(theta)+1;de=(exp(delta))/(exp(delta)+1)}else{th=theta;de=delta}
  
  ut = (1 - de * u)^th
  vt = (1 - de * v)^th
  x = 1 - ut
  y = 1 - vt
  eta1 = 1/(1 - (1 - de)^th)
  tem = (1 - eta1 * x * y)^(1/th - 2)
  pdf = eta1 * de * tem * (th - eta1 * x * y) * ut * vt/(1 - 
                                                           de * u)/(1 - de * v)
  pdf
}


pcond.sbb8=function (v, u, cpar,param=F){
  
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  u=1-u
  v=1-v
  
  if (is.matrix(cpar)) {
    theta = cpar[, 1]
    delta = cpar[, 2]
  }
  else {
    theta = cpar[1]
    delta = cpar[2]
  }
  
  if(param){th=exp(theta)+1;de=(exp(delta))/(exp(delta)+1)}else{th=theta;de=delta}
  
  ut = (1 - de * u)^th
  x = 1 - ut
  y = 1 - (1 - de * v)^th
  eta1 = 1/(1 - (1 - de)^th)
  tem = (1 - eta1 * x * y)^(1/th - 1)
  den = (1 - de * u)/ut
  ccdf = eta1 * y * tem/den
  1-ccdf
}





dsbb8=function (u, v, cpar,param=F)
{
  u[u==1]=.9999999999999
  v[v==1]=.9999999999999
  u[u==0]=.0000000000001
  v[v==0]=.0000000000001
  
  
  u=1-u
  v=1-v
  
  if (is.matrix(cpar)) {
    theta = cpar[, 1]
    delta = cpar[, 2]
  }
  else {
    theta = cpar[1]
    delta = cpar[2]
  }
  
  if(param){th=exp(theta)+1;de=(exp(delta))/(exp(delta)+1)}else{th=theta;de=delta}
  
  ut = (1 - de * u)^th
  vt = (1 - de * v)^th
  x = 1 - ut
  y = 1 - vt
  eta1 = 1/(1 - (1 - de)^th)
  tem = (1 - eta1 * x * y)^(1/th - 2)
  pdf = eta1 * de * tem * (th - eta1 * x * y) * ut * vt/(1 - 
                                                           de * u)/(1 - de * v)
  pdf
}

#BB10====================================
dbb10=function (u, v, cpar,param=F){
  if (is.matrix(cpar)) {
    th = cpar[, 1]
    ppi = cpar[, 2]
  }
  else {
    th = cpar[1]
    ppi = cpar[2]
  }
  
  if(param){th=exp(th);ppi=((exp(ppi))/(exp(ppi)+1))}else{th=th;ppi=ppi}
  ut = u^th
  vt = v^th
  tem = 1 - ppi * (1 - ut) * (1 - vt)
  ttem = tem^(-1/th)
  pdf = ttem/tem/tem
  pdf = pdf * (1 - ppi + ppi * (1 + th) * ut * vt - ppi * (1 - 
                                                             ppi) * (1 - ut) * (1 - vt))
  pdf[ u <= 0 | v <= 0 ] <- 0
  pdf[ u >= 1 | v >= 1] <- 0
  pdf
}

dsbb10=function(u, v, cpar,param=F){
  u1=1-u
  v1=1-v
  val=dbb10(u1, v1, cpar,param)
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
} 

pcond.bb10=function (v, u, cpar,param=F){
  if (is.matrix(cpar)) {
    th = cpar[, 1]
    ppi = cpar[, 2]
  }
  else {
    th = cpar[1]
    ppi = cpar[2]
  }
  if(param){th=exp(th);ppi=((exp(ppi))/(exp(ppi)+1))}else{th=th;ppi=ppi}
  ut = u^th
  vt = v^th
  tem = 1 - ppi * (1 - ut) * (1 - vt)
  ttem = tem^(-1/th)
  ccdf = ttem/tem * v * (1 - ppi + ppi * vt)
  ccdf[ v <= 0 | u <= 0 |  u>= 1] <- 0
  ccdf[ v == 1 ] <- 1
  ccdf
}

pcond.sbb10=function(v, u, cpar,param=F) 
{
  v1=1-v
  u1=1-u
  val=1-pcond.bb10(v1, u1, cpar,param) 
  #val[ v <= 0 | u <= 0 |  u>= 1] <- 0
  #val[ v == 1 ] <- 1
  val
}

#================================================================================
#================================================================================
#================================================================================
#================================================================================
# derivative of the h function with respect to the copula parameter
#================================================================================
pcond.2frank<-function(u,v,a)
{ expa<-exp(a)
num<-(exp(a + a*v)*(exp(2*a*u)*(-1 + v) + expa*v -
                      exp(a*u)*(-1 + u - expa*u + v + expa*v)))
den<-(exp(a*(u + v)) - expa*(-1 + exp(a*u) + exp(a*v)))^2
val<-num/den
val[ u <= 0 | v <= 0 ] <- 0
val[ u >= 1 | v >= 1] <- 0
val
}
#===========================================
pcond.2bvn<-function(u,v,r)
{ qnormv<-qnorm(v)
tem1<-qnorm(u)-r*qnormv
tem<-1-r^2
tem2<-sqrt(tem)
tem3<-(-qnormv*tem2 + tem1*r/tem2)/tem
val<-dnorm(tem1/tem2)*tem3
val[ u <= 0 | v <= 0 ] <- 0
val[ u >= 1 | v >= 1] <- 0
val
}
#===========================================
pcond.2t<-function(u,v,r){
  #df=1
  qtv<-qt(v,df)
  tem1<-qt(u,df)-r*qtv
  tem<-(1-r^2)*(df+qtv^2)/(df+1)
  tem2<-sqrt(tem)
  tem3<--r*(df+qtv^2)/(df+1)/tem2
  tem4<-(-qtv*tem2 - tem1*tem3)/tem
  val<-dt(tem1/tem2,df+1)*tem4
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}


pcond.2t1<-function(u,v,r){
  df=1
  qtv<-qt(v,df)
  tem1<-qt(u,df)-r*qtv
  tem<-(1-r^2)*(df+qtv^2)/(df+1)
  tem2<-sqrt(tem)
  tem3<--r*(df+qtv^2)/(df+1)/tem2
  tem4<-(-qtv*tem2 - tem1*tem3)/tem
  val<-dt(tem1/tem2,df+1)*tem4
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

pcond.2t2<-function(u,v,r){
  df=2
  qtv<-qt(v,df)
  tem1<-qt(u,df)-r*qtv
  tem<-(1-r^2)*(df+qtv^2)/(df+1)
  tem2<-sqrt(tem)
  tem3<--r*(df+qtv^2)/(df+1)/tem2
  tem4<-(-qtv*tem2 - tem1*tem3)/tem
  val<-dt(tem1/tem2,df+1)*tem4
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

pcond.2t3<-function(u,v,r){
  df=3
  qtv<-qt(v,df)
  tem1<-qt(u,df)-r*qtv
  tem<-(1-r^2)*(df+qtv^2)/(df+1)
  tem2<-sqrt(tem)
  tem3<--r*(df+qtv^2)/(df+1)/tem2
  tem4<-(-qtv*tem2 - tem1*tem3)/tem
  val<-dt(tem1/tem2,df+1)*tem4
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

pcond.2t4<-function(u,v,r){
  df=4
  qtv<-qt(v,df)
  tem1<-qt(u,df)-r*qtv
  tem<-(1-r^2)*(df+qtv^2)/(df+1)
  tem2<-sqrt(tem)
  tem3<--r*(df+qtv^2)/(df+1)/tem2
  tem4<-(-qtv*tem2 - tem1*tem3)/tem
  val<-dt(tem1/tem2,df+1)*tem4
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}


pcond.2t5<-function(u,v,r){
  df=5
  qtv<-qt(v,df)
  tem1<-qt(u,df)-r*qtv
  tem<-(1-r^2)*(df+qtv^2)/(df+1)
  tem2<-sqrt(tem)
  tem3<--r*(df+qtv^2)/(df+1)/tem2
  tem4<-(-qtv*tem2 - tem1*tem3)/tem
  val<-dt(tem1/tem2,df+1)*tem4
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

pcond.2t6<-function(u,v,r){
  df=6
  qtv<-qt(v,df)
  tem1<-qt(u,df)-r*qtv
  tem<-(1-r^2)*(df+qtv^2)/(df+1)
  tem2<-sqrt(tem)
  tem3<--r*(df+qtv^2)/(df+1)/tem2
  tem4<-(-qtv*tem2 - tem1*tem3)/tem
  val<-dt(tem1/tem2,df+1)*tem4
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}


pcond.2t7<-function(u,v,r){
  df=7
  qtv<-qt(v,df)
  tem1<-qt(u,df)-r*qtv
  tem<-(1-r^2)*(df+qtv^2)/(df+1)
  tem2<-sqrt(tem)
  tem3<--r*(df+qtv^2)/(df+1)/tem2
  tem4<-(-qtv*tem2 - tem1*tem3)/tem
  val<-dt(tem1/tem2,df+1)*tem4
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

pcond.2t8<-function(u,v,r){
  df=8
  qtv<-qt(v,df)
  tem1<-qt(u,df)-r*qtv
  tem<-(1-r^2)*(df+qtv^2)/(df+1)
  tem2<-sqrt(tem)
  tem3<--r*(df+qtv^2)/(df+1)/tem2
  tem4<-(-qtv*tem2 - tem1*tem3)/tem
  val<-dt(tem1/tem2,df+1)*tem4
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

pcond.2t9<-function(u,v,r){
  df=9
  qtv<-qt(v,df)
  tem1<-qt(u,df)-r*qtv
  tem<-(1-r^2)*(df+qtv^2)/(df+1)
  tem2<-sqrt(tem)
  tem3<--r*(df+qtv^2)/(df+1)/tem2
  tem4<-(-qtv*tem2 - tem1*tem3)/tem
  val<-dt(tem1/tem2,df+1)*tem4
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

#===========================================
pcond.2gumbel=function (v, u, cpar){
  x = -log(u)
  y = -log(v)
  tx = log(x)
  ty = log(y)
  xd = x^cpar
  yd = y^cpar
  sm = xd + yd
  tem = sm^(1/cpar)
  deriv12 = tem * xd * yd * (tem + cpar - 1)/(sm * sm * x *
                                                y)
  lcdf = -tem
  lpdf = lcdf + x + y + log(deriv12)
  logs = log(sm)
  dlsq = cpar * cpar
  dlcu = dlsq * cpar
  sder1 = xd * tx + yd * ty
  mder1 = tem * sder1/(sm * cpar) - tem * logs/dlsq
  logm = log(tem)
  msq = tem * tem
  usq = u * u
  lccdf = x - tem + (1 - cpar) * (logm - tx)
  mu = -tem * xd/(u * sm * x)
  lcder1 = -mder1 + (1 - cpar) * mder1/tem - logm + tx
  lcder1u = -(cpar - 1) * mu/tem - mu - (cpar - 1)/(x * u) -
    1/u
  ccdf = exp(lccdf)
  ccdfdu = ccdf * lcder1u
  ccdfdpar = ccdf * lcder1
  
  ccdfdpar[ v <= 0 | u <= 0 ] <- 0
  ccdfdpar[ v >= 1 | u >= 1] <- 0
  ccdfdpar
}
#===========================================
pcond.2sgumbel<-function(v,u,a)
{ u=1-u
v=1-v
val<--pcond.2gumbel(v,u,a)
val[ v <= 0 | u <= 0 ] <- 0
val[ v >= 1 | u >= 1] <- 0
val
}
#===========================================
pcond.2gumbel.90=function(v,u,a)
{ 
  v=1-v
  val<--pcond.2gumbel(v,u,a)
  val[ v <= 0 | u <= 0 ] <- 0
  val[ v >= 1 | u >= 1] <- 0
  val
}
#===========================================
pcond.2gumbel.270=function(v,u,a){ 
  u=1-u
  val<-pcond.2gumbel(v,u,a)
  val[ v <= 0 | u <= 0 ] <- 0
  val[ v >= 1 | u >= 1] <- 0
  val
}
#===========================================
pcond.2joe=function(v,u,cpar){
  nv=1 - v
  nu=1 - u
  ncpar=1/cpar
  nvcpar=(nv)^cpar
  nucpar=(nu)^cpar
  nncpar=-1 + ncpar
  nv.nu=nvcpar/nucpar
  lnv=log(nv)
  nvnc.cpar=1 + nv.nu - nvcpar
  nvnc.cpar2=nvnc.cpar^(nncpar)
  lnucpar=(nucpar * log(nu))
  lnvcpar=log((1 + nvcpar/(nu)^cpar - nvcpar))
  temp1=nncpar * (nvcpar * lnv/nucpar - nvcpar * lnucpar/(nucpar)^2 - nvcpar * lnv)
  temp2=nvnc.cpar^(nncpar - 1) * temp1 - nvnc.cpar2 * (lnvcpar * (ncpar^2))
  
  cdfdrv=temp2 * (1 - nvcpar) - (nvnc.cpar2) * (nvcpar * lnv)
  cdfdrv[ v <= 0 | u <= 0 ] <- 0
  cdfdrv[ v >= 1 | u >= 1] <- 0
  return(cdfdrv)
}
#===========================================
pcond.2sjoe<-function(v,u,a)
{   u=1-u
v=1-v
val<--pcond.2joe(v,u,a)
val[ v <= 0 | u <= 0 ] <- 0
val[ v >= 1 | u >= 1] <- 0
val
}
#===========================================
pcond.2joe.90 = function(v,u,cpar)
{ 
  v=1-v
  val<--pcond.2joe(v,u,cpar)
  val[ v <= 0 | u <= 0 ] <- 0
  val[ v >= 1 | u >= 1] <- 0
  val
}
#===========================================
pcond.2joe.270=function(v,u,cpar){
  u=1-u
  val<-pcond.2joe(v,u,cpar)
  val[ v <= 0 | u <= 0 ] <- 0
  val[ v >= 1 | u >= 1] <- 0
  val
}

#=====================================
# Derivatives wrt cpar for bbx copulas

# bb1 derivation
bb1deriv= function(u1,u2,theta,delta)
{ t1=u1^(-theta); t2=u2^(-theta)
t1a=t1-1.; t2a=t2-1.
tu1=-log(u1); tu2=-log(u2)
ttu1=log(t1a); ttu2=log(t2a)
td01=(t1a)^delta; td02=(t2a)^delta
td11=td01/t1a; td12=td02/t2a;
td21=td11/t1a; td22=td12/t2a;
s=td01+td02
sder1th= delta*(td11*t1*tu1+td12*t2*tu2)
sder1dl= td01*ttu1+td02*ttu2
m=s^(1./delta); m1=m/s
ts=log(s)
dlsq=delta*delta; dlcu=delta*dlsq
mder1th=m1*sder1th/delta
mder1dl=m1*sder1dl/delta - m*ts/dlsq
m1der1dl= mder1dl/s - m*sder1dl/s^2
list(sm=s,mexp=m,mder1th=mder1th,mder1dl=mder1dl)
}

pcond.2bb1.theta=function ( v,u, cpar)
{
  theta = cpar[1]
  delta = cpar[2]
  mder = bb1deriv(u, v, theta, delta)
  m = mder$mexp
  mder1th = mder$mder1th
  mder1dl = mder$mder1dl
  sm = mder$sm
  mp1 = 1 + m
  msq = m * m
  thsq = theta * theta
  lmp1 = log(mp1)
  logm = log(m)
  t2 = u^(-theta)
  tu2 = -log(u)
  t2a = t2 - 1
  lt2a = log(t2a)
  cf1 = 1 + 1/theta
  dl1n = 1 - delta
  lcdf = -cf1 * lmp1 + (dl1n) * (logm - lt2a) + (theta + 1) *
    tu2
  lcder1th = lmp1/thsq - cf1 * mder1th/mp1
  lcder1th = lcder1th + (dl1n) * (mder1th/m - t2 * tu2/t2a) +
    tu2
  lcder1dl = -cf1 * mder1dl/mp1 + (dl1n) * mder1dl/m - logm +
    lt2a
  dmdx = m/sm/delta
  x = t2a^delta
  dxdu = -delta * theta * x/t2a * t2/u
  lcder1u = -cf1/mp1 * dmdx * dxdu + dl1n/m * dmdx * dxdu -
    dl1n/delta/x * dxdu - (theta + 1)/u
  ccdf = exp(lcdf)
  ccdfdth = ccdf * lcder1th
  ccdfddl = ccdf * lcder1dl
  ccdfdu = ccdf * lcder1u
  ccdfdth[ v <= 0 | u <= 0 ] <- 0
  ccdfdth[ v >= 1 | u >= 1] <- 0
  ccdfdth
}

pcond.2bb1.delta=function ( v,u, cpar)
{
  theta = cpar[1]
  delta = cpar[2]
  mder = bb1deriv(u, v, theta, delta)
  m = mder$mexp
  mder1th = mder$mder1th
  mder1dl = mder$mder1dl
  sm = mder$sm
  mp1 = 1 + m
  msq = m * m
  thsq = theta * theta
  lmp1 = log(mp1)
  logm = log(m)
  t2 = u^(-theta)
  tu2 = -log(u)
  t2a = t2 - 1
  lt2a = log(t2a)
  cf1 = 1 + 1/theta
  dl1n = 1 - delta
  lcdf = -cf1 * lmp1 + (dl1n) * (logm - lt2a) + (theta + 1) *
    tu2
  lcder1th = lmp1/thsq - cf1 * mder1th/mp1
  lcder1th = lcder1th + (dl1n) * (mder1th/m - t2 * tu2/t2a) +
    tu2
  lcder1dl = -cf1 * mder1dl/mp1 + (dl1n) * mder1dl/m - logm +
    lt2a
  dmdx = m/sm/delta
  x = t2a^delta
  dxdu = -delta * theta * x/t2a * t2/u
  lcder1u = -cf1/mp1 * dmdx * dxdu + dl1n/m * dmdx * dxdu -
    dl1n/delta/x * dxdu - (theta + 1)/u
  ccdf = exp(lcdf)
  ccdfdth = ccdf * lcder1th
  ccdfddl = ccdf * lcder1dl
  ccdfdu = ccdf * lcder1u
  
  ccdfddl[ v <= 0 | u <= 0 ] <- 0
  ccdfddl[ v >= 1 | u >= 1] <- 0
  ccdfddl
}


pcond.2sbb1.theta=function(v,u,cpar)
{ u=1-u
v=1-v
val<--pcond.2sbb1.theta(v,u,cpar)
val[ v <= 0 | u <= 0 ] <- 0
val[ v >= 1 | u >= 1] <- 0
val
}

pcond.2sbb1.delta=function(v,u,cpar)
{ u=1-u
v=1-v
val<--pcond.2sbb1.delta(v,u,cpar)
val[ v <= 0 | u <= 0 ] <- 0
val[ v >= 1 | u >= 1] <- 0
val
}




pcond.2bb7.theta=function(v,u,cpar){
  th=cpar[1]
  de=cpar[2]
  
  val=((((((((1 - (((1 - (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 
                     1 + 1)^(-(1/de)))^(((1/th) - 1) - 1) * (((1/th) - 1) * ((((1 - 
                                                                                  (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 1 
                                                                              + 1)^((-(1/de)) - 1) * ((-(1/de)) * ((1 - (1 - v)^th)^((-de) - 1) * ((-de) * 
                                                                                                                                                     ((1 - v)^th * log((1 - v)))) + (1 - (1 - u)^th)^((-de) - 
                                                                                                                                                                                                        1) * ((-de) * ((1 - u)^th * log((1 - u)))))))) - (1 - (((1 - 
                                                                                                                                                                                                                                                                   (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 1 + 1)^(-(1/de)))^((1/th) - 
                                                                                                                                                                                                                                                                                                                                          1) * (log((1 - (((1 - (1 - u)^th)^(-de) - 1) + (1 - (1 - 
                                                                                                                                                                                                                                                                                                                                                                                                 v)^th)^(-de) - 1 + 1)^(-(1/de)))) * (1/th^2))) * ((((1 - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                        (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 1 + 1)^(-(1/de))) - 
              ((1 - (((1 - (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 
                       1 + 1)^(-(1/de)))^((1/th) - 1)) * ((((1 - (1 - u)^th)^(-de) - 
                                                              1) + (1 - (1 - v)^th)^(-de) - 1 + 1)^((-(1/de)) - 1) * 
                                                            ((-(1/de)) * ((1 - (1 - v)^th)^((-de) - 1) * ((-de) * 
                                                                                                            ((1 - v)^th * log((1 - v)))) + (1 - (1 - u)^th)^((-de) - 
                                                                                                                                                               1) * ((-de) * ((1 - u)^th * log((1 - u)))))))) * 
             ((1 - (1 - u)^th)^(-de) - 1 + 1) - ((1 - (((1 - (1 - u)^th)^(-de) - 
                                                          1) + (1 - (1 - v)^th)^(-de) - 1 + 1)^(-(1/de)))^((1/th) - 
                                                                                                             1)) * ((((1 - (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 
                                                                                                                       1 + 1)^(-(1/de))) * ((1 - (1 - u)^th)^((-de) - 1) * ((-de) * 
                                                                                                                                                                              ((1 - u)^th * log((1 - u)))))) * (1 - (1 - (1 - u)^th)) + 
            ((1 - (((1 - (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 
                     1 + 1)^(-(1/de)))^((1/th) - 1)) * ((((1 - (1 - u)^th)^(-de) - 
                                                            1) + (1 - (1 - v)^th)^(-de) - 1 + 1)^(-(1/de))) * ((1 - 
                                                                                                                  (1 - u)^th)^(-de) - 1 + 1) * ((1 - u)^th * log((1 - u))))/(((1 - 
                                                                                                                                                                                 (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 1 + 1) + 
           ((1 - (((1 - (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 
                    1 + 1)^(-(1/de)))^((1/th) - 1)) * ((((1 - (1 - u)^th)^(-de) - 
                                                           1) + (1 - (1 - v)^th)^(-de) - 1 + 1)^(-(1/de))) * ((1 - 
                                                                                                                 (1 - u)^th)^(-de) - 1 + 1) * (1 - (1 - (1 - u)^th)) * 
           ((1 - (1 - v)^th)^((-de) - 1) * ((-de) * ((1 - v)^th * 
                                                       log((1 - v)))) + (1 - (1 - u)^th)^((-de) - 1) * ((-de) * 
                                                                                                          ((1 - u)^th * log((1 - u)))))/(((1 - (1 - u)^th)^(-de) - 
                                                                                                                                            1) + (1 - (1 - v)^th)^(-de) - 1 + 1)^2)/(1 - (1 - u)^th) + 
          ((1 - (((1 - (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 
                   1 + 1)^(-(1/de)))^((1/th) - 1)) * ((((1 - (1 - u)^th)^(-de) - 
                                                          1) + (1 - (1 - v)^th)^(-de) - 1 + 1)^(-(1/de))) * ((1 - 
                                                                                                                (1 - u)^th)^(-de) - 1 + 1) * (1 - (1 - (1 - u)^th))/(((1 - 
                                                                                                                                                                         (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 1 + 
                                                                                                                                                                       1) * ((1 - u)^th * log((1 - u)))/(1 - (1 - u)^th)^2)/(1 - 
                                                                                                                                                                                                                               u))
  val[ v <= 0 | u <= 0 ] <- 0
  val[ v >= 1 | u >= 1] <- 0
  val
}



pcond.2bb7.delta=function(v,u,cpar){
  th=cpar[1]
  de=cpar[2]
  
  val=((((((1 - (((1 - (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 
                   1 + 1)^(-(1/de)))^((1/th) - 1)) * ((((1 - (1 - u)^th)^(-de) - 
                                                          1) + (1 - (1 - v)^th)^(-de) - 1 + 1)^(-(1/de)) * (log((((1 - 
                                                                                                                     (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 1 + 1)) * 
                                                                                                              (1/de^2)) - (((1 - (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 
                                                                                                                             1 + 1)^((-(1/de)) - 1) * ((-(1/de)) * ((1 - (1 - v)^th)^(-de) * 
                                                                                                                                                                      log((1 - (1 - v)^th)) + (1 - (1 - u)^th)^(-de) * log((1 - 
                                                                                                                                                                                                                              (1 - u)^th))))) - (1 - (((1 - (1 - u)^th)^(-de) - 1) + (1 - 
                                                                                                                                                                                                                                                                                        (1 - v)^th)^(-de) - 1 + 1)^(-(1/de)))^(((1/th) - 1) - 1) * 
            (((1/th) - 1) * ((((1 - (1 - u)^th)^(-de) - 1) + (1 - (1 - 
                                                                     v)^th)^(-de) - 1 + 1)^(-(1/de)) * (log((((1 - (1 - u)^th)^(-de) - 
                                                                                                                1) + (1 - (1 - v)^th)^(-de) - 1 + 1)) * (1/de^2)) - (((1 - 
                                                                                                                                                                         (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 1 + 
                                                                                                                                                                       1)^((-(1/de)) - 1) * ((-(1/de)) * ((1 - (1 - v)^th)^(-de) * 
                                                                                                                                                                                                            log((1 - (1 - v)^th)) + (1 - (1 - u)^th)^(-de) * log((1 - 
                                                                                                                                                                                                                                                                    (1 - u)^th)))))) * ((((1 - (1 - u)^th)^(-de) - 1) + (1 - 
                                                                                                                                                                                                                                                                                                                           (1 - v)^th)^(-de) - 1 + 1)^(-(1/de)))) * ((1 - (1 - u)^th)^(-de) - 
                                                                                                                                                                                                                                                                                                                                                                       1 + 1) - ((1 - (((1 - (1 - u)^th)^(-de) - 1) + (1 - (1 - 
                                                                                                                                                                                                                                                                                                                                                                                                                              v)^th)^(-de) - 1 + 1)^(-(1/de)))^((1/th) - 1)) * ((((1 - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 1 + 1)^(-(1/de))) * 
           ((1 - (1 - u)^th)^(-de) * log((1 - (1 - u)^th)))) * (1 - 
                                                                  (1 - (1 - u)^th))/(((1 - (1 - u)^th)^(-de) - 1) + (1 - (1 - 
                                                                                                                            v)^th)^(-de) - 1 + 1) + ((1 - (((1 - (1 - u)^th)^(-de) - 
                                                                                                                                                              1) + (1 - (1 - v)^th)^(-de) - 1 + 1)^(-(1/de)))^((1/th) - 
                                                                                                                                                                                                                 1)) * ((((1 - (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 
                                                                                                                                                                                                                           1 + 1)^(-(1/de))) * ((1 - (1 - u)^th)^(-de) - 1 + 1) * (1 - 
                                                                                                                                                                                                                                                                                     (1 - (1 - u)^th)) * ((1 - (1 - v)^th)^(-de) * log((1 - (1 - 
                                                                                                                                                                                                                                                                                                                                               v)^th)) + (1 - (1 - u)^th)^(-de) * log((1 - (1 - u)^th)))/(((1 - 
                                                                                                                                                                                                                                                                                                                                                                                                              (1 - u)^th)^(-de) - 1) + (1 - (1 - v)^th)^(-de) - 1 + 1)^2)/(1 - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                             (1 - u)^th)/(1 - u))
  val[ v <= 0 | u <= 0 ] <- 0
  val[ v >= 1 | u >= 1] <- 0
  val
}



pcond.2sbb7.theta=function(v,u,cpar)
{ u=1-u
v=1-v
val<--pcond.2bb7.theta(v,u,cpar)
val[ v <= 0 | u <= 0 ] <- 0
val[ v >= 1 | u >= 1] <- 0
val
}

pcond.2sbb7.delta=function(v,u,cpar)
{ u=1-u
v=1-v
val<--pcond.2bb7.delta(v,u,cpar)
val[ v <= 0 | u <= 0 ] <- 0
val[ v >= 1 | u >= 1] <- 0
val
}




pcond.2bb8.theta=function(v,u,cpar){
  th=cpar[1]
  de=cpar[2]
  deu = de * u
  dev = de * v
  deu1 = 1 - deu
  dev1 = 1 - dev
  de1 = 1 - de
  
  val=((((de1)^th * log((de1))/(1 - (de1)^th)^2 * (1 - (1 - 
                                                          dev)^th) - (1/(1 - (de1)^th)) * ((dev1)^th * 
                                                                                             log((dev1)))) * ((1 - (1/(1 - (de1)^th)) * (1 - 
                                                                                                                                           (deu1)^th) * (1 - (dev1)^th))^(1/th - 1)) - (1/(1 - 
                                                                                                                                                                                             (de1)^th)) * (1 - (dev1)^th) * ((1 - (1/(1 - (1 - 
                                                                                                                                                                                                                                             de)^th)) * (1 - (deu1)^th) * (1 - (dev1)^th))^(1/th - 
                                                                                                                                                                                                                                                                                              1) * (log((1 - (1/(1 - (de1)^th)) * (1 - (deu1)^th) * 
                                                                                                                                                                                                                                                                                                           (1 - (dev1)^th))) * (1/th^2)) + (1 - (1/(1 - (de1)^th)) * 
                                                                                                                                                                                                                                                                                                                                              (1 - (deu1)^th) * (1 - (dev1)^th))^((1/th - 1) - 
                                                                                                                                                                                                                                                                                                                                                                                    1) * ((1/th - 1) * (((de1)^th * log((de1))/(1 - (1 - 
                                                                                                                                                                                                                                                                                                                                                                                                                                       de)^th)^2 * (1 - (deu1)^th) - (1/(1 - (de1)^th)) * 
                                                                                                                                                                                                                                                                                                                                                                                                           ((deu1)^th * log((deu1)))) * (1 - (dev1)^th) - 
                                                                                                                                                                                                                                                                                                                                                                                                          (1/(1 - (de1)^th)) * (1 - (deu1)^th) * ((1 - de * 
                                                                                                                                                                                                                                                                                                                                                                                                                                                     v)^th * log((dev1)))))))/((deu1)/((1 - de * 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          u)^th)) + (1/(1 - (1 - de)^th)) * (1 - (dev1)^th) * 
         ((1 - (1/(1 - (de1)^th)) * (1 - (deu1)^th) * (1 - 
                                                         (dev1)^th))^(1/th - 1)) * ((deu1) * ((1 - 
                                                                                                 deu)^th * log((deu1)))/((deu1)^th)^2)/((1 - 
                                                                                                                                           deu)/((deu1)^th))^2)
  
  val[ v <= 0 | u <= 0 ] <- 0
  val[ v >= 1 | u >= 1] <- 0
  val
  
}

pcond.2bb8.delta=function(v,u,cpar){
  th=cpar[1]
  de=cpar[2]
  deu = de * u
  dev = de * v
  deu1 = 1 - deu
  dev1 = 1 - dev
  
  val=((((1/(1 - (1 - de)^th)) * ((dev1)^(th - 1) * (th * v)) - 
           (1 - de)^(th - 1) * th/(1 - (1 - de)^th)^2 * (1 - (1 - de * 
                                                                v)^th)) * ((1 - (1/(1 - (1 - de)^th)) * (1 - (1 - de * 
                                                                                                                u)^th) * (1 - (dev1)^th))^(1/th - 1)) - (1/(1 - (1 - 
                                                                                                                                                                   de)^th)) * (1 - (dev1)^th) * ((1 - (1/(1 - (1 - de)^th)) * 
                                                                                                                                                                                                    (1 - (deu1)^th) * (1 - (dev1)^th))^((1/th - 1) - 
                                                                                                                                                                                                                                          1) * ((1/th - 1) * (((1/(1 - (1 - de)^th)) * ((deu1)^(th - 
                                                                                                                                                                                                                                                                                                  1) * (th * u)) - (1 - de)^(th - 1) * th/(1 - (1 - de)^th)^2 * 
                                                                                                                                                                                                                                                                 (1 - (deu1)^th)) * (1 - (dev1)^th) + (1/(1 - 
                                                                                                                                                                                                                                                                                                            (1 - de)^th)) * (1 - (deu1)^th) * ((dev1)^(th - 
                                                                                                                                                                                                                                                                                                                                                         1) * (th * v))))))/((deu1)/((deu1)^th)) + (1/(1 - 
                                                                                                                                                                                                                                                                                                                                                                                                         (1 - de)^th)) * (1 - (dev1)^th) * ((1 - (1/(1 - (1 - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                            de)^th)) * (1 - (deu1)^th) * (1 - (dev1)^th))^(1/th - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             1)) * (u/((deu1)^th) - (deu1) * ((deu1)^(th - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        1) * (th * u))/((deu1)^th)^2)/((deu1)/((1 - de * 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  u)^th))^2)
  
  val[ v <= 0 | u <= 0 ] <- 0
  val[ v >= 1 | u >= 1] <- 0
  val
  
}



pcond.2sbb8.theta=function(v,u,cpar)
{ u=1-u
v=1-v
val<--pcond.2bb8.theta(v,u,cpar)
val[ v <= 0 | u <= 0 ] <- 0
val[ v >= 1 | u >= 1] <- 0
val
}

pcond.2sbb8.delta=function(v,u,cpar)
{ u=1-u
v=1-v
val<--pcond.2bb8.delta(v,u,cpar)
val[ v <= 0 | u <= 0 ] <- 0
val[ v >= 1 | u >= 1] <- 0
val
}



pcond.2bb10.theta=function(v,u,cpar){
  th=cpar[1]
  delta=cpar[2]
  ut=u^th
  vt=v^th
  logu=log(u)
  logv=log(v)
  utn=-1+ut
  vtn=-1+vt
  deltauv=delta *utn* vtn
  
  num=(delta *th *(1+th) *(ut *logu+vt *(logv-ut *(logu+logv))))
  denm=(-1+delta *utn* vtn)
  exp1=(1+delta *vtn)*(num/denm+log(1-deltauv))
  val=v*(1-deltauv)^(-((1+th)/th))*(delta* vt *logv+exp1/(th^2))
  val[ v <= 0 | u <= 0 ] <- 0
  val[ v >= 1 | u >= 1] <- 0
  val
}

pcond.2bb10.delta=function(v,u,cpar){
  th=cpar[1]
  ppi=cpar[2]
  ut=u^th
  vt=v^th
  utn=-1+ut
  vtn=-1+vt
  
  exp1=(1-ppi *utn *vtn)^(-2-1/th)
  dccdf=v*vtn*exp1*(-1+ ((1+th)*ut) + (ppi*utn*vtn))
  
  dccdf.ppi=dccdf/th 
  dccdf.ppi[ v <= 0 | u <= 0 ] <- 0
  dccdf.ppi[ v >= 1 | u >= 1] <- 0
  dccdf.ppi
}


pcond.2sbb10.delta=function(v,u,cpar){
  th=cpar[1]
  ppi=cpar[2]
  v1=1-v
  u1=1-u
  val=-pcond.2bb10.delta(v1,u1,cpar)
  val[ v <= 0 | u <= 0 ] <- 0
  val[ v >= 1 | u >= 1] <- 0
  val
}

pcond.2sbb10.theta=function(v,u,cpar){
  v1=1-v
  u1=1-u
  val=-pcond.2bb10.theta(v1, u1, cpar)
  val[ v <= 0 | u <= 0 ] <- 0
  val[ v >= 1 | u >= 1] <- 0
  val
}


#derivative of copula densities wrt to parameters
d2gumbel=function(u,v,th){
  lu=-log(u)
  lv=-log(v)
  lvth=(lv)^th
  luth=(lu)^th
  lulv=luth+lvth
  ith=1/th
  lulvth=(lulv)^(ith)
  llu=log(lu)
  llv=log(lv)
  lluth=luth *llu
  llvth=lvth *llv
  term2=(lu)^(-1+th)*(lulv)^(-2+ith)
  term3=(-1+th+lulvth)
  term4=(-1+th+lulvth)
  term5=1/(u*v)*exp(-lulvth)*term2*(lv)^(-1+th)
  term6=(lulv)^(-1+ith)
  term7=(th*lluth-(lulv)*log(lulv)+th*llvth)
  term8=(-1+th+lulvth)
  term9=((-2+ith) *(lluth+llvth))
  term10=-(log(lulv)/th^2)
  term11=(-1+th+lulvth) 
  term12=(lulv)^(-1+ith)
  term13=(th *lluth-(lulv)* log(lulv)+th *llvth)
  term14=1+term3*llu+term4*llv+(1/(th^2))*term6 *term7
  term15=(1/(th^2))*term8* term12 * term13
  val=term5*(term14-term15+term11*(term10+term9/lulv))
  val
}

d2frank=function(u,v,cpar){
  uv=u + v  
  nuv=- u + v
  uvneg= u - v
  cparuv=cpar*(uv)
  denom=(-exp(cpar) + exp(cpar + cpar*u) - 
           exp(cparuv) + exp(cpar + cpar*v))^3
  tem1=cpar*(1+nuv)
  tem2=exp(cpar*(2 + v))
  tem3=exp(cpar*(2 + u))
  tem4=cpar*(-1 + uv)
  tem5=(tem2*(1 + cpar*(uvneg)) + 
          exp(cpar + cpar*u)*(-1 + cpar*(1 + uvneg)) + 
          tem3*(1 + cpar*(nuv)) + 
          exp(cpar + cpar*v)*(-1 + tem1) + 
          exp(cpar*(1 + uv))*(-1 + cpar*(-2 + uv)) + 
          exp(cparuv)*(1 - tem4) + 
          exp(cpar)*(1 + tem4) - 
          exp(2*cpar)*(1 + cparuv))
  
  val=(exp(cpar*(1 + uv))*tem5)/denom
  val  
}

d2sgumbel<-function(u,v,a)
{ u=1-u
v=1-v
val<-d2gumbel(u,v,a)
val[ u <= 0 | v <= 0 ] <- 0
val[ u >= 1 | v >= 1] <- 0
val
}

d2gumbel.90<-function(u,v,a){
  u=1-u
  val<-d2gumbel(u,v,a)
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

d2gumbel.270<-function(u,v,a){
  v=1-v
  val<-d2gumbel(u,v,a)
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}


#library(VineCopula)
#d2tcop=function(u1,u2,rho){
#  cop <- BiCop(family = 2, par = rho, par2 = df)
#  val=NULL
#  for(i in 1:length(u2)){
#    val=c(val,BiCopDeriv(u1=u1, u2=u2[i], cop, deriv = "par"))
#  }
#  val
#}

#d2tcop1=function(u1,u2,rho){
#  df=1
#  cop <- BiCop(family = 2, par = rho, par2 = df)
#  val=NULL
#  for(i in 1:length(u2)){
#    val=c(val,BiCopDeriv(u1=u1, u2=u2[i], cop, deriv = "par"))
#  }
#  val
#}

d2tcop2=function(u1,u2,rho){
  df=2.000000000001
  cop <- BiCop(family = 2, par = rho, par2 = df)
  val=NULL
  for(i in 1:length(u2)){
    val=c(val,BiCopDeriv(u1=u1, u2=u2[i], cop, deriv = "par"))
  }
  val
}

d2tcop3=function(u1,u2,rho){
  df=3
  cop <- BiCop(family = 2, par = rho, par2 = df)
  val=NULL
  for(i in 1:length(u2)){
    val=c(val,BiCopDeriv(u1=u1, u2=u2[i], cop, deriv = "par"))
  }
  val
}

d2tcop4=function(u1,u2,rho){
  df=4
  cop <- BiCop(family = 2, par = rho, par2 = df)
  val=NULL
  for(i in 1:length(u2)){
    val=c(val,BiCopDeriv(u1=u1, u2=u2[i], cop, deriv = "par"))
  }
  val
}

d2tcop5=function(u1,u2,rho){
  df=5
  cop <- BiCop(family = 2, par = rho, par2 = df)
  val=NULL
  for(i in 1:length(u2)){
    val=c(val,BiCopDeriv(u1=u1, u2=u2[i], cop, deriv = "par"))
  }
  val
}

d2tcop6=function(u1,u2,rho){
  df=6
  cop <- BiCop(family = 2, par = rho, par2 = df)
  val=NULL
  for(i in 1:length(u2)){
    val=c(val,BiCopDeriv(u1=u1, u2=u2[i], cop, deriv = "par"))
  }
  val
}

d2tcop7=function(u1,u2,rho){
  df=7
  cop <- BiCop(family = 2, par = rho, par2 = df)
  val=NULL
  for(i in 1:length(u2)){
    val=c(val,BiCopDeriv(u1=u1, u2=u2[i], cop, deriv = "par"))
  }
  val
}

d2tcop8=function(u1,u2,rho){
  df=8
  cop <- BiCop(family = 2, par = rho, par2 = df)
  val=NULL
  for(i in 1:length(u2)){
    val=c(val,BiCopDeriv(u1=u1, u2=u2[i], cop, deriv = "par"))
  }
  val
}

d2tcop9=function(u1,u2,rho){
  df=9
  cop <- BiCop(family = 2, par = rho, par2 = df)
  val=NULL
  for(i in 1:length(u2)){
    val=c(val,BiCopDeriv(u1=u1, u2=u2[i], cop, deriv = "par"))
  }
  val
}

#library(VineCopula)
d2bvn=function(u1,u2,rho){
  cop <- BiCop(family = 1, par = rho)
  val=NULL
  for(i in 1:length(u2)){
    val=c(val,BiCopDeriv(u1=u1, u2=u2[i], cop, deriv = "par"))
  }
  val
}


#library(VineCopula)
d2joe=function(u1,u2,rho){
  cop <- BiCop(family = 6, par = rho)
  val=NULL
  for(i in 1:length(u2)){
    val=c(val,BiCopDeriv(u1=u1, u2=u2[i], cop, deriv = "par"))
  }
  val
}

d2sjoe=function(u1,u2,rho){
  cop <- BiCop(family = 16, par = rho)
  val=NULL
  for(i in 1:length(u2)){
    val=c(val,BiCopDeriv(u1=u1, u2=u2[i], cop, deriv = "par"))
  }
  val
}

d2joe.90<-function(u,v,a){
  u=1-u
  val<-d2joe(u,v,a)
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}

d2joe.270<-function(u,v,a){
  v=1-v
  val<-d2joe(u,v,a)
  val[ u <= 0 | v <= 0 ] <- 0
  val[ u >= 1 | v >= 1] <- 0
  val
}


#===========================================
#===========================================
# Quantile of conditional copulas
#===========================================
#===========================================
#conditional quantiel for Frank copula
qcond.frank = function (p, u, cpar)
{
  cpar0 = exp(-cpar)
  cpar1 = 1 - cpar0
  etem = exp(-cpar * u + log(1/p - 1))
  tem = 1 - cpar1/(etem + 1)
  v = (-log(tem))/cpar
  isinf = is.infinite(v)
  v[isinf] = (-log(cpar0 + etem[isinf]))/cpar
  v
}

#conditional quantiel for Gaussian copula
qcond.bvn = function (p, v, cpar)
{
  x = qnorm(p)
  y = qnorm(v)
  tem = x * sqrt(1 - cpar * cpar) + cpar * y
  pnorm(tem)
}

#conditional quantiel for Student-t copula
qcond.t = function (p, u, cpar){
  #df=1
  rho = cpar
  x2 = qt(p, df + 1)
  x1 = qt(u, df)
  tem = x2 * sqrt((df + x1 * x1) * (1 - rho * rho)/(df + 1)) +
    rho * x1
  pt(tem, df)
}

qcond.t1 = function (p, u, cpar){
  df=1
  rho = cpar
  x2 = qt(p, df + 1)
  x1 = qt(u, df)
  tem = x2 * sqrt((df + x1 * x1) * (1 - rho * rho)/(df + 1)) +
    rho * x1
  pt(tem, df)
}

qcond.t2 = function (p, u, cpar){
  df=2
  rho = cpar
  x2 = qt(p, df + 1)
  x1 = qt(u, df)
  tem = x2 * sqrt((df + x1 * x1) * (1 - rho * rho)/(df + 1)) +
    rho * x1
  pt(tem, df)
}

qcond.t3 = function (p, u, cpar){
  df=3
  rho = cpar
  x2 = qt(p, df + 1)
  x1 = qt(u, df)
  tem = x2 * sqrt((df + x1 * x1) * (1 - rho * rho)/(df + 1)) +
    rho * x1
  pt(tem, df)
}

qcond.t4 = function (p, u, cpar){
  df=4
  rho = cpar
  x2 = qt(p, df + 1)
  x1 = qt(u, df)
  tem = x2 * sqrt((df + x1 * x1) * (1 - rho * rho)/(df + 1)) +
    rho * x1
  pt(tem, df)
}

qcond.t5 = function (p, u, cpar){
  df=5
  rho = cpar
  x2 = qt(p, df + 1)
  x1 = qt(u, df)
  tem = x2 * sqrt((df + x1 * x1) * (1 - rho * rho)/(df + 1)) +
    rho * x1
  pt(tem, df)
}

qcond.t6 = function (p, u, cpar){
  df=6
  rho = cpar
  x2 = qt(p, df + 1)
  x1 = qt(u, df)
  tem = x2 * sqrt((df + x1 * x1) * (1 - rho * rho)/(df + 1)) +
    rho * x1
  pt(tem, df)
}

qcond.t7 = function (p, u, cpar){
  df=7
  rho = cpar
  x2 = qt(p, df + 1)
  x1 = qt(u, df)
  tem = x2 * sqrt((df + x1 * x1) * (1 - rho * rho)/(df + 1)) +
    rho * x1
  pt(tem, df)
}

qcond.t8 = function (p, u, cpar){
  df=8
  rho = cpar
  x2 = qt(p, df + 1)
  x1 = qt(u, df)
  tem = x2 * sqrt((df + x1 * x1) * (1 - rho * rho)/(df + 1)) +
    rho * x1
  pt(tem, df)
}

qcond.t9 = function (p, u, cpar){
  df=9
  rho = cpar
  x2 = qt(p, df + 1)
  x1 = qt(u, df)
  tem = x2 * sqrt((df + x1 * x1) * (1 - rho * rho)/(df + 1)) +
    rho * x1
  pt(tem, df)
}
#conditional quantiel for Gumbel copula
qcond.gumbel=function (p, u, cpar, eps = 1e-06, mxiter = 30)
{
  x = -log(u)
  cpar1 = cpar - 1
  con = log(p) - x - cpar1 * log(x)
  z = x * (2^(1/cpar))
  mxdif = 1
  iter = 0
  while (mxdif > eps & iter < mxiter) {
    g = z + cpar1 * log(z) + con
    gp = 1 + cpar1/z
    diff = g/gp
    z = z - diff
    iter = iter + 1
    while (z <= x) {
      diff = diff/2
      z = z + diff
    }
    mxdif = abs(diff)
  }
  y = ((z^cpar) - (x^cpar))^(1/cpar)
  return(exp(-y))
}
#conditional quantiel for survival Gumbel copula
qcond.sgumbel = function(p, u, cpar, eps = 1e-06, mxiter = 30)
{
  v1 = qcond.gumbel(1 - p, 1 - u, cpar, eps, mxiter)
  1 - v1
}


qcond.gumbel.90=function(p, u, cpar, eps = 1e-06, mxiter = 30)
{
  v1 = qcond.gumbel(1 - p,  u, cpar, eps, mxiter)
  1 - v1
}


qcond.gumbel.270=function(p, u, cpar, eps = 1e-06, mxiter = 30)
{
  v1 = qcond.gumbel( p, 1- u, cpar, eps, mxiter)
  v1
}


#conditional quantiel for Joe copula
qcond.joe=function (p, u, cpar, eps = 1e-06, mxiter = 30)
{
  cpar1 = cpar - 1
  cpartem = -cpar1/(1 + cpar1)
  cpar1inv = -1/cpar1
  ubar = 1 - u
  ud = ubar^cpar
  cpari = 1/cpar
  tem = (1 - p)^cpartem - 1
  tem = tem * (1 - u)^(-cpar1) + 1
  v = tem^cpar1inv
  v = 1 - v
  diff = 1
  iter = 0
  while (max(abs(diff)) > eps & iter < mxiter) {
    vbar = 1 - v
    vd = vbar^cpar
    sm = ud + vd - ud * vd
    smd = sm^cpari
    c21 = 1 + vd/ud - vd
    c21 = c21^(-1 + cpari)
    c21 = c21 * (1 - vd)
    pdf = smd * ud * vd * (cpar1 + sm)/sm/sm/ubar/vbar
    iter = iter + 1
    if (any(is.nan(pdf)) | any(is.nan(c21))) {
      diff = diff/(-2)
    }
    else diff = (c21 - p)/pdf
    v = v - diff
    while (min(v) <= 0 | max(v) >= 1 | max(abs(diff)) > 0.25) {
      diff = diff/2
      v = v + diff
    }
  }
  v
}

#conditional quantiel for survival Gumbel copula
qcond.sjoe = function(p, u, cpar, eps = 1e-06, mxiter = 30)
{
  v1 = qcond.joe(1 - p, 1 - u, cpar, eps, mxiter)
  1 - v1
}

qcond.joe.90 = function(p, u, cpar, eps = 1e-06, mxiter = 30)
{
  v1 = qcond.joe( 1-p, u, cpar, eps, mxiter)
  1-v1
}

qcond.joe.270 = function(p, u, cpar, eps = 1e-06, mxiter = 30)
{
  v1 = qcond.joe( p, 1-u, cpar, eps, mxiter)
  v1
}

#conditional quantiel for BB1 copula
qcond.bb1=function (p, u, cpar, eps = 1e-06, mxiter = 30, iprint = F)
{
  if (is.matrix(cpar)) {
    th = cpar[, 1]
    de = cpar[, 2]
  }
  else {
    th = cpar[1]
    de = cpar[2]
  }
  
  de1 = 1/de
  th1 = 1/th
  ut = (u^(-th) - 1)
  x = ut^de
  den = (1 + ut)^(-th1 - 1) * ut/x
  pden = p * den
  y = pden^(-1/(th1 * de1 + 1)) - x
  if (x < 1e-05) {
    y = min(x * (p^(-de/(de - 1)) - 1), 1e-05)
    if (iprint)
      cat("\\nsmall x case ", x, y, "\\n")
  }
  if (x > 1e+05) {
    r = p^(-de * th/(1 + de * th)) - 1
    y = r * x
    if (iprint)
      cat("\\nlarge x case ", x, y, "\\n")
    eps = eps * 1e-04
  }
  if (de <= 1.1) {
    thr = -th/(1 + th)
    tem = (p^thr) - 1
    tem = tem * (u^(-th)) + 1
    y = (tem - 1)^de
  }
  else if (th < 0.2) {
    v = qcond.gumbel(p, u, de)
    y = ((v^(-th)) - 1)^de
  }
  diff = 1
  iter = 0
  while (abs(diff/y) > eps & iter < mxiter) {
    sm = x + y
    smd = sm^de1
    G21 = (1 + smd)^(-th1 - 1) * smd/sm
    gpdf = -G21
    gpdf = gpdf/(1 + smd)/sm/de/th
    gpdf = gpdf * (th * (de - 1) + (th * de + 1) * smd)
    iter = iter + 1
    diff = (G21 - pden)/gpdf
    y = y - diff
    if (iprint)
      cat(iter, y, diff, "\\n")
    while (y <= 0) {
      diff = diff/2
      y = y + diff
    }
  }
  v = (y^de1 + 1)^(-th1)
  if (iprint)
    cat(v, "ended at iter. ", iter, "\\n")
  return(v)
}

qcond.sbb1 = function(p, u, cpar, eps = 1e-06, mxiter = 30){
  p1 = 1 - p
  u1 = 1 - u
  v1 = qcond.bb1(p1, u1, cpar, eps, mxiter)
  1 - v1
}

qcond.bb7=function (p, u, cpar, eps = 1e-06, mxiter = 30, iprint = F) 
{
  if (is.matrix(cpar)) {
    th = cpar[, 1]
    de = cpar[, 2]
  }
  else {
    th = cpar[1]
    de = cpar[2]
  }
  de1 = 1/de
  th1 = 1/th
  ut = 1 - (1 - u)^th
  x = ut^(-de) - 1
  if (x <= 0) {
    v = 0.999999
    if (u > v) 
      v = u
    if (iprint) 
      cat("\n**** x below 1 ", p, u, x, "\n")
    return(v)
  }
  den = (1 - ut)^(th1 - 1) * ut/(x + 1)
  pden = p * den
  rhs = pden * (1 - (2 * x + 1)^(-de1))^(1 - th1)
  y = rhs^(-de/(de + 1)) - 1 - x
  if (y <= 0) 
    y = 0.1
  if (th > 3 & (u > 0.8 | p > 0.9)) {
    diff = 1
    v = u
    v1 = 0
    v2 = 1
    vold = v
    if (iprint) 
      cat("\n (p,u)=", p, u, "\n")
    while (diff > eps) {
      ccdf = pcond.bb7(v, u, cpar)
      di = ccdf - p
      if (di < 0) {
        v1 = v
      }
      else {
        v2 = v
      }
      diff = v2 - v1
      vold = v
      v = (v1 + v2)/2
      if (iprint) 
        cat(vold, ccdf, di, "\n")
    }
    return(v)
  }
  if (x < 1e-05) {
    epsx = de * (1 - ut)
    tem = p * (1 - (1 + de1) * epsx)
    tem = tem^(-th/(th - 1)) - 1
    epsy = tem * epsx
    if (epsy > 1e-05) 
      epsy = 1e-05
    y = epsy
    if (iprint) 
      cat("\nsmall x case ", x, y, "\n")
  }
  if (th < 1.01) {
    thr = -th/(1 + th)
    tem = (p^thr) - 1
    tem = tem * (u^(-th)) + 1
    y = (tem - 1)^de
  }
  else if (de < 0.1) {
    v = 1 - qcond.joe(p, u, de)
    y = v^th
    y = (1 - y)^(-de) - 1
  }
  diff = 1
  iter = 0
  while (max(abs(diff/y)) > eps & iter < mxiter) {
    sm = x + y + 1
    smd = sm^(-de1)
    G21 = (1 - smd)^(th1 - 1) * smd/sm
    gpdf = -G21
    gpdf = gpdf/(1 - smd)/sm/de/th
    gpdf = gpdf * (th * (de + 1) - (th * de + 1) * smd)
    iter = iter + 1
    diff = (G21 - pden)/gpdf
    y = y - diff
    while (min(y) <= 0 | max(abs(diff)) > 5) {
      diff = diff/2
      y = y + diff
    }
    if (iprint) 
      cat(iter, y, diff, "\n")
  }
  v = 1 - (1 - (y + 1)^(-de1))^th1
  if (iprint) 
    cat(v, "ended at iter. ", iter, "\n")
  v
}

qcond.sbb7 = function(p, u, cpar, eps = 1e-06, mxiter = 30){
  p1 = 1 - p
  u1 = 1 - u
  v1 = qcond.bb7(p1, u1, cpar, eps, mxiter)
  1 - v1
}

#------- bb8
pbb8=function (u, v, cpar) 
{
  if (is.matrix(cpar)) {
    vth = cpar[, 1]
    de = cpar[, 2]
  }
  else {
    vth = cpar[1]
    de = cpar[2]
  }
  x = 1 - (1 - de * u)^vth
  y = 1 - (1 - de * v)^vth
  eta1 = 1/(1 - (1 - de)^vth)
  tem = (1 - eta1 * x * y)^(1/vth)
  cdf = (1 - tem)/de
  cdf
}
qcond.bb8=function (p, u, cpar, eps = 1e-06, mxiter = 30, iprint = F) 
{
  vth = cpar[1]
  de = cpar[2]
  vth1 = 1/vth
  eta = (1 - (1 - de)^vth)
  eta1 = 1/eta
  be = 4 * pbb8(0.5, 0.5, cpar) - 1
  ut = (1 - de * u)^vth
  x = 1 - ut
  con = log(eta1) + (1 - vth1) * log(1 - x) - log(p)
  iter = 0
  diff = 1
  y = 0.5 * x
  if (be >= 0.8) 
    y = x
  while (iter < mxiter & abs(diff) > eps) {
    tem = 1 - eta1 * x * y
    h = log(y) + (vth1 - 1) * log(tem) + con
    hp = 1/y + (1 - vth1) * eta1 * x/tem
    diff = h/hp
    y = y - diff
    if (iprint) 
      cat(iter, diff, y, "\\n")
    if (any(is.nan(y))) {
      if (iprint) 
        cat("***", p, u, cpar, x, "\\n")
      return(u)
    }
    while (min(y) <= 0 | max(y) >= eta) {
      diff = diff/2
      y = y + diff
    }
    iter = iter + 1
  }
  if (iprint & iter >= mxiter) 
    cat("***did not converge\\n")
  v = (1 - y)^vth1
  v = (1 - v)/de
  v
}

qcond.sbb8 = function(p, u, cpar, eps = 1e-06, mxiter = 30){
  p1 = 1 - p
  u1 = 1 - u
  v1 = qcond.bb8(p1, u1, cpar, eps, mxiter)
  1 - v1
}

#===BB10
pbb10=function (u, v, cpar) 
{
  if (is.matrix(cpar)) {
    th = cpar[, 1]
    ppi = cpar[, 2]
  }
  else {
    th = cpar[1]
    ppi = cpar[2]
  }
  ut = u^th
  vt = v^th
  tem = 1 - ppi * (1 - ut) * (1 - vt)
  cdf = u * v * tem^(-1/th)
  cdf
}
#conditional quantiel for BB10 copula
qcond.bb10=function (p, u, cpar, eps = 1e-06, mxiter = 30, iprint = F) 
{
  th = cpar[1]
  ppi = cpar[2]
  th1 = 1/th
  be = 4 * pbb10(0.5, 0.5, cpar) - 1
  ut = u^th
  mxdif = 1
  iter = 0
  diff = 1
  v = 0.7 * u
  if (be >= 0.8) 
    v = u
  while (mxdif > eps & iter < mxiter) {
    vt = v^th
    tem = 1 - ppi * (1 - ut) * (1 - vt)
    ttem = tem^(-th1)
    ccdf = ttem/tem * v * (1 - ppi + ppi * vt)
    pdf = ttem/tem/tem
    pdf = pdf * (1 - ppi + ppi * (1 + th) * ut * vt - ppi * 
                   (1 - ppi) * (1 - ut) * (1 - vt))
    h = ccdf - p
    hp = pdf
    diff = h/hp
    v = v - diff
    iter = iter + 1
    while (min(v) <= 0 | max(v) >= 1) {
      diff = diff/2
      v = v + diff
    }
    if (iprint) 
      cat(iter, diff, v, "\\n")
    mxdif = abs(diff)
  }
  if (iprint & iter >= mxiter) {
    cat("***did not converge\\n")
    cat("p=", p, " u=", u, " theta=", th, " pi=", ppi, " lastv=", 
        v, "\\n")
  }
  v
}
#conditional quantiel for survival BB10 copula
qcond.sbb10=function (p, u, cpar, eps = 1e-06, mxiter = 30, iprint = F) 
{
  p=1-p
  u=1-u
  th = cpar[1]
  ppi = cpar[2]
  th1 = 1/th
  be = 4 * pbb10(0.5, 0.5, cpar) - 1
  ut = u^th
  mxdif = 1
  iter = 0
  diff = 1
  v = 0.7 * u
  if (be >= 0.8) 
    v = u
  while (mxdif > eps & iter < mxiter) {
    vt = v^th
    tem = 1 - ppi * (1 - ut) * (1 - vt)
    ttem = tem^(-th1)
    ccdf = ttem/tem * v * (1 - ppi + ppi * vt)
    pdf = ttem/tem/tem
    pdf = pdf * (1 - ppi + ppi * (1 + th) * ut * vt - ppi * 
                   (1 - ppi) * (1 - ut) * (1 - vt))
    h = ccdf - p
    hp = pdf
    diff = h/hp
    v = v - diff
    iter = iter + 1
    while (min(v) <= 0 | max(v) >= 1) {
      diff = diff/2
      v = v + diff
    }
    if (iprint) 
      cat(iter, diff, v, "\\n")
    mxdif = abs(diff)
  }
  if (iprint & iter >= mxiter) {
    cat("***did not converge\\n")
    cat("p=", p, " u=", u, " theta=", th, " pi=", ppi, " lastv=", 
        v, "\\n")
  }
  1-v
}



