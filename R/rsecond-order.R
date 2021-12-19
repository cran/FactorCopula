simulation.nested=function(nn,grsize,qcond0,qcondg,param)
{
  d = sum(grsize)
  mgrp = length(grsize)
  
  z0 = runif(nn)
  z = matrix(0, nrow = nn, ncol = mgrp)
  zdata = matrix(0, nrow = nn, ncol = d)

    for (jg in 1:mgrp) {
      qcond.jg=qcond0[[jg]]
      for (i in 1:nn) {
        z[i, jg] = qcond.jg(runif(1), z0[i], param[jg])
      }
    }
  
  ind = 0
  for (jg in 1:mgrp) {
    ind1 = ind + 1
    ind2 = ind + grsize[jg]
    ind = ind + grsize[jg]
    for (ij in ind1:ind2) {
      ijm = ij + mgrp
      qcondg.ij=qcondg[[ij]]
        for (i in 1:nn) {
          zdata[i, ij] = qcondg.ij(runif(1), z[i, jg], 
                                  param[ijm])
        }
    }
  }
  zdata
}
