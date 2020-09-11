#transform continuous to ordinal
continuous2ordinal = function(continuous, categ){
  if( is.matrix(continuous) == FALSE ){continuous=matrix(continuous)}
  n = nrow(continuous) 
  u = apply(continuous, 2, rank)/(n+1)
  d= ncol(u)
  cutpts=matrix(rep((seq(1:(categ-1))/(categ)),d),ncol=d)
  ydat = matrix(NA, ncol = d, nrow = n)
  for (j in 1:d) {
    y = cut(u[, j], c(0, cutpts[, j], 1))
    ydat[, j] = as.integer(y) - 1
  }
  ydat
}

#transform count to ordinal
count2ordinal=function(count, categ){
  if( is.matrix(count) == FALSE ){count=matrix(count)}
  d= ncol(count)
  n = nrow(count) 
  ydat = matrix(NA, ncol = d, nrow = n)
  for (j in 1:d) {
    cutpts=seq(from = min(count[,j]),to = max(count[,j]), 
               length.out = categ + 1)
    y = cut(count[,j], breaks=cutpts,include.lowest = T)
    ydat[, j] = as.integer(y) - 1
  }
  ydat
}

#transform continuous to ordinal for simulation
continuous2ordinal.sim = function(u, categ){
  if( is.matrix(u) == FALSE ){u=matrix(u)}
  d= ncol(u)
  cutpts=matrix(rep((seq(1:(categ-1))/(categ)),d),ncol=d)
  ydat = matrix(NA, ncol = d, nrow = nrow(u))
  for (j in 1:d) {
    y = cut(u[, j], c(0, cutpts[, j], 1))
    ydat[, j] = as.integer(y) - 1
  }
  ydat
}