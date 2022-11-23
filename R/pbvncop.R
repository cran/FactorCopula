# bivariate t copula cdf
# cpar = copula parameter (rho,nu) with -1<rho<1, nu>0
#   nu is an integer here because input to pbvt is integer
pbvtcop=function(u,v,cpar)
{ rho=cpar[1]; nu=cpar[2]
if(nu<1) nu=1
nu=floor(nu) # input to pbvt is an integer >= 1.
# need correction for boundary, check for better fix later (pbvn is OK)
u[u>=1]=.9999999
v[v>=1]=.9999999
u[u<=0]=.0000001
v[v<=0]=.0000001
xt=qt(u,nu); yt=qt(v,nu); 
# note that pbvt returns -1 for nu<1.
out=pbvt(xt,yt,cpar) # function in this library
out
}

# bivariate normal copula
# cpar = copula parameter with -1<cpar<1
pbvncop=function(u,v,cpar)
{ # endpoint corrections to prevent NaN
  u[1-u<1.e-9]=1-1.e-9
  v[1-v<1.e-9]=1-1.e-9
  u[u<1.e-9]=1.e-9
  v[v<1.e-9]=1.e-9
  x=qnorm(u)
  y=qnorm(v)
  cdf=pbnorm(x,y,cpar)  # function in this library
  cdf
}