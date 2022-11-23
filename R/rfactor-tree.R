#nsim: simulation times
#A: vine array
#parfactor: copula parameter for first factor
#parvine: copula parameter for first vine tree
#qcondfactor:qconditional copula for first factor 
#qcondvine:qconditional copula for vine tree

simfacvine1 = function (nsim, A, parfactor, 
                        parvine, qcondfactor,qcondvine)
{
  d = ncol(A)
  p = matrix(runif(nsim * d), nsim, d)
  qq = array(0, c(nsim, d, d))
  u = matrix(0, nsim, d)
  
  u[, 1] = p[, 1]
  qq[, 1, 1] = p[, 1]
  qq[, 2, 2] = p[, 2]
  
  for(i in 1:nsim){
    u[i, 2] = qcondvine(p[i, 2], p[i, 1], parvine[1])
    qq[i, 1, 2] = u[i, 2]
  }
  for (j in 3:d) {
    for(i in 1:nsim){
      qq[i, 2, j] = p[i, j]
      qq[i, 1, j] = qcondvine(qq[i, 2, j], u[i, A[1, j]], parvine[j-1] )
      u[i, j] = qq[i, 1, j]
    }
  }
  
  # latent variable ~ U(0,1)
  v = runif(nsim)
  y = matrix(0, nsim, d)
  for(j in 1:d){
    th1 = parfactor[j]
    for(i in 1:nsim){
      y[i,j] = qcondfactor(u[i,j] , v[i]  , th1)
    }
  }
  y
}

#nsim: simulation times
#A: vine array
#parfactor1: copula parameter for first factor
#parfactor2: copula parameter for second factor
#parvine: copula parameter for first vine tree
#qcondfactor1:qconditional copula for first factor 
#qcondfactor2:qconditional copula for second factor 
#qcondvine:qconditional copula for vine tree
simfacvine2 = function (nsim, A, parfactor1,parfactor2, 
                        parvine, qcondfactor1,qcondfactor2,qcondvine)
{
  d = ncol(A)
  p = matrix(runif(nsim * d), nsim, d)
  qq = array(0, c(nsim, d, d))
  u = matrix(0, nsim, d)
  
  u[, 1] = p[, 1]
  qq[, 1, 1] = p[, 1]
  qq[, 2, 2] = p[, 2]
  
  for(i in 1:nsim){
    u[i, 2] = qcondvine(p[i, 2], p[i, 1], parvine[1])
    qq[i, 1, 2] = u[i, 2]
  }
  
  for (j in 3:d) {
    for(i in 1:nsim){
      qq[i, 2, j] = p[i, j]
      qq[i, 1, j] = qcondvine(qq[i, 2, j], u[i, A[1, j]], parvine[j-1] )
      u[i, j] = qq[i, 1, j]
    }
  }
  
  # latent variables ~ U(0,1)
  v1 = runif(nsim)
  v2 = runif(nsim)
  y = matrix(0, nsim, d)
  for(j in 1:d){
    th_f1 = parfactor1[j]
    th_f2 = parfactor2[j]
    for(i in 1:nsim){
      tem2 = qcondfactor2(u[i,j] , v2[i] , th_f2)
      y[i,j] = qcondfactor1(tem2 , v1[i]  , th_f1)
    }
  }
  return(y)
}

