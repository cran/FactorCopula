simulation.bifactor=function(nn , grsize , qcond0,qcondg , param){
  d = sum(grsize)
  mgrp = length(grsize)
  th00 = param[1:d]
  thgg = param[(d + 1):(2 * d)]
  
  zdata = matrix(0, nrow = nn, ncol = d)
  z0 = runif(nn)
  z = matrix(runif(nn * mgrp), ncol = mgrp)
  
  ind = 0
  for (jg in 1:mgrp) {
    ind1 = ind + 1
    ind2 = ind + grsize[jg]
    ind = ind + grsize[jg]
    for (ij in ind1:ind2){
      qcond00=qcond0[[ij]]
      qcondgg=qcondg[[ij]]
      
      for (i in 1:nn) {
        q1=qcondgg(runif(1),z[i,jg],thgg[ij])
        zdata[i,ij]=qcond00(q1, z0[i],th00[ij])
        
      }
    }
  }
  zdata
}

