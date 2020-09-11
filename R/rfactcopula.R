# Simulation for one-factor and two-factor copula models
# These models can be constructed from different
# bivaraite linking copulas.
#-------------------------------------------
# simulations for one-factor copula
#-------------------------------------------
r1Fcopula=function(n,d,theta,qcondcop){
  if (is.array(theta)) {
    theta =  as.vector(theta)
  }
  if (length(theta) != d) {
    return("The number of parameters are not
           equal to the dimensions (d)")
  }
  v = runif(n)
  u = matrix(NA,n,d)
  for (j in 1:d) {
    qcondcopula=qcondcop[[j]]
    th_j = theta[[j]]
    for (i in 1:n) {
      V= v[i] 
      X= runif(1)
      u[i,j] = qcondcopula(X,V,th_j)
    }
  }
  return(u)
}

#-------------------------------------------
# simulations for two-factor copula
#-------------------------------------------
r2Fcopula=function (n, d, theta, delta, qcondfact1, qcondfact2){
  if (is.array(theta)) {
    theta =  as.vector(theta)
  }
  if (is.array(delta)) {
    delta =  as.vector(delta)
  }
  h1 = runif(n)
  h2 = runif(n)
  u = matrix(NA, n, d)
  for (j in 1:d) {
    qcondcopula2=qcondfact2[[j]]
    qcondcopula1=qcondfact1[[j]]
    th_j=theta[[j]]
    del_j=delta[[j]]
    for (i in 1:n) {
      v1 = h1[i]
      v2 = h2[i]
      X = runif(1)
      q2 = qcondcopula2(X, v2, del_j)
      u[i, j] = qcondcopula1(q2, v1, th_j)
    }
  }
  u
}