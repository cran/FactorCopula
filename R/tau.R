### Kendall's tau conversion for copula parameters
### The choice of copulas:
### Gaussian="bvn"
### t-copula="bvt",
### Gumbel="gum"
### Joe="joe"
### Frank="frk"

par2tau = function(copulaname,cpar){
    par=cpar[1]
    if(length(cpar)==2){par2=cpar[2]}
    #=========================================================================================
    #========================================GAUSSIAN=========================================
    if(copulaname=="bvn" )
    {
        par2tau.bvn = function(par)
        {
            if( any(abs(par) >= 1) )
            stop("The parameter of the Gaussian copula has to be in the interval (-1,1).")
            2/pi * asin(par)
        }
        tau=par2tau.bvn(par)
    }
    #=========================================================================================
    #========================================STUDENT-T========================================
    if(copulaname=="bvt1"||copulaname=="bvt2" ||copulaname=="bvt3"||
        copulaname=="bvt4"||copulaname=="bvt5"||copulaname=="bvt6"||
        copulaname=="bvt7"||copulaname=="bvt8"||copulaname=="bvt9")
    {
        par2tau.bvt = function(par)
        {
            if( any(abs(par) >= 1))
            stop("The parameter of the Student-t copula has to be in the interval (-1,1).")
            2/pi * asin(par)
        }
        tau=par2tau.bvt(par)
    }
    #=========================================================================================
    #==========================================GUMBEL=========================================
    if(copulaname=="gum" || copulaname=="rgum")
    {
        par2tau.gumbel = function(par)
        {
            if( any(par < 1))
            stop("The parameter of the Gumbel copula has to be in the interval [1,inf).")
            1 - 1/par
        }
        tau=par2tau.gumbel(par)
    }
    #=========================================================================================
    #==========================================GUMBEL=========================================
    if(copulaname=="1rgum" || copulaname=="2rgum")
    {
      par2tau.rgumbel = function(par)
      {
        if( any(par > 1))
          stop("The parameter of the Gumbel copula has to be in the interval (inf,-1].")
        -1 - 1/par
      }
      tau=par2tau.rgumbel(par)
    }
    #=========================================================================================
    #===========================================FRANK=========================================
    if(copulaname=="frk" || copulaname=="frank")
    {
        par2tau.frank = function(par)
        {
            if (any(par == 0) )
            stop("The parameter of the Frank copula has to be unequal to 0.")
            f = function(x)
            {
                x/(exp(x) - 1)
            }
            d= length(par)
            tau=rep(NA,d)
            for(j in 1:d){
                tau[j] = 1 - 4/par[j] + 4/par[j]^2 * integrate(f, lower = 0,upper = par[j])$value
            }
            tau
        }
        tau=par2tau.frank(par)
    }
    
    #=========================================================================================
    #===========================================JOE=========================================
    if(copulaname=="joe" || copulaname=="rjoe")
    {
        par2tau.joe = function(par)
        {
        if ( any(par < 1))
        stop("The parameter of the Joe copula has to be in the interval (1,inf).")
         param1 <- 2/par + 1
        tem <- digamma(2) - digamma(param1)
        tau <- 1 + tem * 2/(2 - par)
        tau[par == 2] <- 1 - trigamma(2)
        tau
    }
    tau=par2tau.joe(par)
    }
    #=========================================================================================
    #===========================================JOE=========================================
    if(copulaname=="1rjoe" |copulaname=="2rjoe")
    {
      par2tau.rjoe = function(par)
      {
        if ( any(par > -1))
          stop("The parameter of the Joe copula has to be in the interval (-inf,1).")
        par <- abs(par)
        param1 <- 2/par + 1
        tem <- digamma(2) - digamma(param1)
        tau <- 1 + tem * 2/(2 - par)
        tau[par == 2] <- 1 - trigamma(2)
        return(-tau)
      }
      tau=par2tau.rjoe(par)
    }
    #=========================================================================================
    #===========================================BB1=========================================
    if(copulaname=="bb1"|copulaname=="rbb1")
    {
         par2tau.bb1 = function(par,par2)
         {
             if(par2 < 1)
             stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
              1 - 2/(par2 * (par + 2))
    }
    tau=par2tau.bb1(par,par2)
    }
  #=========================================================================================
  #===========================================BB7=========================================
    if(copulaname=="bb7"|copulaname=="rbb7")
    {
      par2tau.bb7 = function(par,par2)
      {
        if(par < 1)
          stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
        
        kt <- function(t) {
          ((1 - (1 - t)^par)^-par2 - 1)/(-par * par2 * (1 -
          t)^(par - 1) * (1 - (1 - t)^par)^(-par2 - 1))
        }
        1 + 4 * integrate(kt, 0, 1)$value
      }
     tau = par2tau.bb7(par,par2)
    }
    #=========================================================================================
    #===========================================BB8=========================================
    if(copulaname=="bb8"|copulaname=="rbb8"){
      par2tau.bb8 = function(par,par2)
      {
        if(par < 1)stop("The 1st par of the BB8/sBB8 copula has to be in the interval [1,inf).")
        if(par2 > 1 || par2 <= 0 )stop("The 2nd par of the BB8/sBB8 copula has to be in the interval (0,1].")
        vth = par
        de = par2
        eta = 1 - (1 - de)^vth
        psider = function(s) {
          s1 = exp(-s)
          psid = (eta/vth/de) * ((1 - eta * s1)^(1/vth - 1)) * 
            s1
          s * psid^2
        }
        tem = integrate(psider, 0, Inf, rel.tol = 1e-06)
        tau = 1 - 4 * tem$value
        tau
      }
      tau=par2tau.bb8(par,par2)
    }
  #=========================================================================================
  #===========================================BB10=========================================
  if(copulaname=="bb10"|copulaname=="rbb10"){
    par2tau.bb10 = function(par,par2)
    {
      th = par
      p = par2
      if(p>0.99999999){p=0.9999998}
      
      qq = (1 - p)^(1/th)
      psider = function(s) {
        s1 = exp(-s)
        sth = exp(-s/th)
        psid = (qq/th) * ((1 - p * s1)^(-1/th - 1)) * sth
        s * psid^2
      }
      tem = integrate(psider, 0, Inf, rel.tol = 1e-06)
      tau = 1 - 4 * tem$value
      tau
    }
    tau=par2tau.bb10(par,par2)
  }
    #=========================================================================================
    #=========================================================================================
    return(tau)
}


############################################################################################
############################################################################################
# Converting tau to copula parameters
# the copula families are Guassian, Gumbel,Survival Gumbel and Frank

tau2par = function(copulaname,tau)
{
    #=========================================================================================
    #========================================GAUSSIAN/t=======================================
  if(copulaname=="bvn"||copulaname=="bvt1"||copulaname=="bvt2" ||copulaname=="bvt3"||
     copulaname=="bvt4"||copulaname=="bvt5"||copulaname=="bvt6"||
     copulaname=="bvt7"||copulaname=="bvt8"||copulaname=="bvt9")
    {
        tau2par.bvn=function(tau){
            sin(pi * tau/2)
        }
        par = tau2par.bvn(tau)
    }
    #=========================================================================================
    #==========================================GUMBEL=========================================
    if(copulaname=="gum" | copulaname=="rgum")
    {
        tau2par.gumbel = function(tau){
          if (any(tau< 0))return(1.00)
            1/(1 - tau)
        }
        par = tau2par.gumbel(tau)
    }
  #==========================================Rotated GUMBEL============================
  if(copulaname=="1rgum" |copulaname=="2rgum")
  {
    tau2par.rgumbel = function(tau)
    {
      if (any(tau> 0))return(1.00)
      tau=abs(tau)
      1/(1 - tau)
    }
    par = -tau2par.rgumbel(tau)
  }
    #=========================================================================================
    #===========================================FRANK=========================================
    if(copulaname=="frk" | copulaname=="frank")
    {
      tau2par.frank = function(tau)
      {
        d=length(tau)
        
        a=rep(NA,d)
        for(j in 1:d){
          
          if ((tau[j]) > 0.99999) return(Inf)
          a[j] <- 1
          
          if (tau[j] < 0) {
            a[j] <- -1
            tau[j] <- -tau[j]
          }
        }
        
        
        frac1<-function(t)
        {
          t/(exp(t)-1)
        }
        
        debye1<-function(theta)
        {
          (1/theta)*integrate(frac1, 0, theta)$value
        }
        
        d<-length(tau)
        v<-rep(NA,d)
        for(j in 1:d){
          
          v[j] <- uniroot(function(x) tau[j] - (1 - 4/x + 4/x * debye1(x)),
                          lower = 0 + .Machine$double.eps^0.5, upper = 5e5,
                          tol = .Machine$double.eps^0.5)$root
        }
        return(a*v)
        
      }
        par = tau2par.frank(tau)
    }
  #=======================================================================================
  #===========================================JOE=========================================
  if(copulaname=="joe" | copulaname=="rjoe")
  {  
    tau2par.joe <- function(tau) {
    if ( any((tau) < 0)) {
      return(1.00)
    } else {
      
      tauF <- function(a) {
        1 + 4/a^2 * integrate(function(x) log(x) * x * (1 - x)^(2 * (1 - a)/a), 0, 1)$value
      }
       d = length(tau)
      v=rep(NA,d)
      for(j in 1:d){
      v[j] <- uniroot(function(x) tau[j] - tauF(x), lower = 1, upper = 500, tol = .Machine$double.eps^.5)$root
      }
      v
    }
  }
    par = tau2par.joe(tau)
  }
  #=========================================================================================
  #===========================================Rotated Joe=============================
  if(copulaname=="1rjoe" |copulaname=="2rjoe"){  
    tau2par.rjoe <- function(tau) {
      if (tau > 0) {
        return(1.00)
      } else {
        tau=abs(tau)
        tauF <- function(a) {
          1 + 4/a^2 * integrate(function(x) log(x) * x * (1 - x)^(2 * (1 - a)/a), 0, 1)$value
        }
        d = length(tau)
        v=rep(NA,d)
        for(j in 1:d){
          v[j] <- uniroot(function(x) tau[j] - tauF(x), lower = 1, upper = 500, tol = .Machine$double.eps^.5)$root
        }
        v
      }
    }
    par = -tau2par.rjoe(tau)
  }
    #=========================================================================================
    #=========================================================================================
    return(par)
}

