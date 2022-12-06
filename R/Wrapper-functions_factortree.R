#Wrapper function for factor tree copula
# (1) MLE 1-factor tree copula model.
# (2) MLE 2-factor tree copula model.
# (4) Simulation code for 1-2-factor tree.
# (5) Model selection.

#A function returns qcondcopula that matches the name
#of the inputs "copnames"
qcond_cop=function(copnames)
{ 
  cop_names=c("bvn","bvt1","bvt2","bvt3","bvt4","bvt5","bvt6","bvt7","bvt8","bvt9",
    "frk","gum","rgum","1rgum","2rgum")
  qcop=c(qcond.bvn,qcond.t1,qcond.t2,qcond.t3,qcond.t4,qcond.t5,qcond.t6
         ,qcond.t7,qcond.t8,qcond.t9,qcond.frank,qcond.gumbel,qcond.sgumbel,
         qcond.gumbel.90,qcond.gumbel.270)
  output=qcop[which(copnames==cop_names)][[1]]
  return(output)
}
  

#Simulation from 1-factor tree copula
#n: simulation times
#A: vine array
#theta1: copula parameter for first factor
#delta: copula parameter for first vine tree
#copname1: a single copula name 
#(e.g., "bvn", "gum", etc.) for first factor
#copname2: a single copula name 
#(e.g., "bvn", "gum", etc.) for second factor
#copnametree: a single copula name 
#(e.g., "bvn", "gum", etc.) for vine tree
r1factortree=function(n, d, A, copname1, copnametree, theta1, delta,K){
  
  qcondfactor=qcond_cop(copname1)
  qcondvine=qcond_cop(copnametree)
  u = simfacvine1(nsim=n, A, parfactor=theta1, 
                          parvine=delta, qcondfactor,qcondvine)
  y = continuous2ordinal(u,K)
}

r2factortree=function(n, d, A, copname1, copname2, copnametree,
                      theta1, theta2, delta,K){
  
  qcondfactor1=qcond_cop(copname1)
  qcondfactor2=qcond_cop(copname2)
  qcondvine=qcond_cop(copnametree)
  u = simfacvine2(nsim=n, A, parfactor1=theta1,parfactor2=theta2, 
                              parvine=delta, qcondfactor1,qcondfactor2,qcondvine)
  y = continuous2ordinal(u,K)
}

#> A
#[,1] [,2] [,3] [,4] [,5]
#[1,]    1    1    2    3    4
#[2,]    0    2    1    2    3
#[3,]    0    0    3    1    2
#[4,]    0    0    0    4    1
#[5,]    0    0    0    0    5
#Dvine
#d=5
#A=matrix(0,d,d)
#A[1,]=c(1,c(1:(d-1)))
#diag(A)=1:d
#only the first row and diagonal are used for 1-truncated vine
#set.seed(1)
#ydat1=r1factortree(n=500, d, A, copname1="frk", copnametree="frk", 
#             theta=cpar1, delta=cpar2, K=5)

#set.seed(1)
#ydat2=r2factortree(n=500, d, A, copname1="gum", copname2="gum", 
#                   copnametree="gum",theta1=cpar1,theta2=cpar1,
#                   delta=cpar2, K=5)
  

#------------------------------------------------------------
# MLE Estimation of 1-factor tree copula for item response
# Input:
#y: data in ordinal scale
#A: vine array
#cop: copula names
#E.g., 
#copnmF = rep("bvn",d)
#copnmV = rep("bvn",d-1)
#cop=c(copnmF,copnmV)
# gl: Gauss legendre quadrature points and weights.

# Outputs:
# 1) "cutpoints": cutpoints for ordinal data.
# 2) "loglik": log-likelihood value.
# 3) "taus": estimated Kendall taus.
# 4) "SEs": Standard errors of the estimated taus.

mle1FactorTree=function(y, A, cop, gl, hessian = F, print.level = 0) 
{
  gln <- gl$nodes
  glw <- gl$weights
  nq <- length(glw)
  d <- ncol(y)
  n <- nrow(y)
  cutpoints <- unifcuts(y)
  initpar = initval_tree(cop)
  mle <- nlm(loglk_f1vine,p=initpar,ydat=y, A, 
             cutp=unifcuts(y),nq=nq,gln,glw,
             copnm=cop,print.level=print.level,hessian=hessian)
  if (hessian == TRUE) {
    SEtaus = rep(NA, 2*d-1 )
    try(SEtaus <- sqrt(diag(solve(mle$h))), T)
  }
  else {
    SEtaus = NULL
  }
  list(cutpoints = cutpoints, taus = mle$e, SEs = SEtaus, loglik = -mle$m)
}


# MLE Estimation of 2-factor tree copula for item response
# Input:
#y: data in ordinal scale
#A: vine array
#cop: copula names
#E.g., 
#copnmF1 = rep("bvn",d)
#copnmF2 = rep("bvn",d)
#copnmV = rep("bvn",d-1)
#cop=c(copnmF1,copnmF2,copnmV)
# gl: Gauss legendre quadrature points and weights.

# Outputs:
# 1) "cutpoints": cutpoints for ordinal data.
# 2) "loglik": log-likelihood value.
# 3) "taus": estimated Kendall taus.
# 4) "SEs": Standard errors of the estimated taus.
mle2FactorTree=function(y, A, cop, gl, hessian = F, print.level = 0) 
{
  gln <- gl$nodes
  glw <- gl$weights
  nq <- length(glw)
  d <- ncol(y)
  n <- nrow(y)
  cutpoints <- unifcuts(as.matrix(y))
  initpar = initval_tree(cop)
  mle <- nlm(loglk_f2vine,p=initpar,ydat=y, A, 
             cutp=unifcuts(y),nq=nq,gln,glw,
             copnm=cop,print.level=print.level,hessian=hessian)
  if (hessian == TRUE) {
    SEtaus = rep(NA, 3*d-1 )
    try(SEtaus <- sqrt(diag(solve(mle$h))), T)
  }
  else {
    SEtaus = NULL
  }
  list(cutpoints = cutpoints, taus = mle$e, SEs = SEtaus, loglik = -mle$m)
}


#mle2=mleF2V(y, A, cop=c(rep("gum",3*d-1)), gl=gauss.quad.prob(15), 
#           hessian = F, print.level = 2)

#------------------------------------------------------------
#Vine tree selection
#Input---
#y: n \times d matrix  with ordinal data 
#rmat: polychoric correlation matrix
#alg: select from the following 
# 1:  1-factor tree partial algorithm.
# 2:  2-factor tree partial algorithm.
# 3:  factor tree polychoric algorithm.
# Note that if selected 3, then any other correlation matrix
# would work as well.

#output---
#F1treeA: matrix filled only in 1st row and diagonal
#        value in each row are connected with its corresponding 
#        diagnoal value.(vine array)
#Pcor=pmat_f1: Partial correlation matrix (Y_j1Y_j2|V_1V_2)
selectFactorTrVine=function(y,rmat,alg){
  if(alg==1){
    output=select1FTr.partial(y,rmat)
    }else if(alg==2){
      output=select2FTr.partial(y,rmat)
      }else if(alg==3){
        output=selectVineTree(y,rmat)
        }else{
          output="select a viable number for model selection"
        }
  return(output)
}


copselect1FactorTree=function(y, A, copnames, gl){
  d=ncol(y)
  copula_f1=copselect1Factor(y, copnames, gl, d)
  f1copulas=unlist(copula_f1$"1st factor")
  vinetreecopulas=copselect1Ftree(y, A, f1copnames=f1copulas, vinecopnames=copnames, gl, d)
  
  return(list("1Factor copulas" = f1copulas,
              "Vine tree copulas" = vinetreecopulas$"Vine tree copulas", 
              "estimated taus" = vinetreecopulas$"estimated taus",
              "Log-likelihood"=vinetreecopulas$"Log-likelihood"))
}

copselect2FactorTree=function(y, A, copnames, gl){
  d=ncol(y)
  copula_f2=copselect2Factor(y, copnames, gl, d)
  f2copulas1=copula_f2$"1st factor"
  f2copulas2=copula_f2$"2nd factor"
  vinetreecopulas=copselect2Ftree(y, A, f1copnames=f2copulas1, f2copnames=f2copulas2, vinecopnames=copnames, gl, d)
  
  return(list("1Factor copulas" = f2copulas1,"2Factor copulas" = f2copulas2,
              "Vine tree copulas" = vinetreecopulas$"Vine tree copulas", 
              "estimated taus" = vinetreecopulas$"estimated taus",
              "Log-likelihood"=vinetreecopulas$"Log-likelihood"))
  
}



#---------------------------------------- Vuong's test
# purpose ========================
# Calculate the Vuong statistics,
# p-value and 95% confidence interval CI will be exported.
# To detect whether two models are statistically different.
# Input   ========================
#models: choose a number in [1,...,4]
#         1: M1:1-factor tree  , M2:1-factor tree
#         2: M1:2-factor tree  , M2:2-factor tree
#         3: M1:1-factor tree  , M2:2-factor tree
#         4: M1:1-factor       , M2:1-factor tree
#         5: M1:1-factor       , M2:2-factor tree
#         6: M1:2-factor       , M2:1-factor tree
#         7: M1:2-factor       , M2:2-factor tree
#tau.m1: vector of copula paramters in Kendall's tau for model 1.
#tau.m2: vector of copula paramters in Kendall's tau for model 2.
# copnames.m1: vector of copula families names for model 1.
# copnames.m2: vector of copula families names for model 2.
#y: item response data
#A.m1: vine array for model 1, if it is 1-factor tree, or 2-factor tree.
#      otherwise keep as NULL
#A.m2: vine array for model 2, if it is 1-factor tree, or 2-factor tree.
#      
# Output   ========================
# Vector of z score, p-value, 95% CI
vuong_FactorTree = function(models, y, A.m1=NULL, tau.m1, copnames.m1, 
                            tau.m2, copnames.m2,A.m2, 
                            nq){
  d=ncol(y)
  
  if(models==1){
    out=vuong.1fv(tau.M1=tau.m1,cop.M1=copnames.m1,A.M1=A.m1, 
                  tau.M2=tau.m2,cop.M2=copnames.m2,A.M2=A.m2, 
                             ordinal=y, nq, param=F)
  }else if(models==2){
    out=vuong.2fv(tau.M1=tau.m1,cop.M1=copnames.m1,A.M1=A.m1, 
                  tau.M2=tau.m2,cop.M2=copnames.m2,A.M2=A.m2, 
                  ordinal=y, nq, param=F)
    
  }else if(models==3){
    out=vuong.1fv.2fv(tau.M1=tau.m1,cop.M1=copnames.m1,A.M1=A.m1, 
                  tau.M2=tau.m2,cop.M2=copnames.m2,A.M2=A.m2, 
                  ordinal=y, nq, param=F)
    
  }else if(models==4){
    out=vuong.1fv.1f(tau.M1=tau.m1,cop.M1=copnames.m1,
                 tau.M2=tau.m2,cop.M2=copnames.m2,A.M2=A.m2, 
                 ordinal=y, nq)
  }else if(models==5){
    cop.m1f1=copnames.m1[1:d]
    cop.m1f2=copnames.m1[(d+1):(2*d)]
    out=vuong.1fv.2f(tau.M1=tau.m1,cop.M1F1=cop.m1f1,cop.M1F2=cop.m1f2,
                 tau.M2=tau.m1,cop.M2=copnames.m2,A.M2=A.m2, 
                 ordinal=y, nq)
  }else if(models==6){
    out=vuong.2fv.1f(tau.M1=tau.m1,cop.M1=copnames.m1,
                 tau.M2=tau.m2,cop.M2=copnames.m2,A.M2=A.m2, 
                 ordinal=y, nq)
  }else if(models==7){
    cop.m1f1=copnames.m1[1:d]
    cop.m1f2=copnames.m1[(d+1):(2*d)]
    out=vuong.2fv.2f(tau.M1=tau.m1,cop.M1F1=cop.m1f1,cop.M1F2=cop.m1f2,
                 tau.M2=tau.m1,cop.M2=copnames.m2,A.M2=A.m2,
                 ordinal=y, nq)
  }else{
    out="select a viable number of models in Vuong's test"
  }
  return(out)
}
