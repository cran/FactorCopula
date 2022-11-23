#Vine tree selection techniques using Prim's algorithm.

#library(optrees)
#library(FactorCopula) for loadings of BVN factor model
#library(statmod)



#Input---
#y: n \times d matrix  with ordinal data 
#rmat: polychoric correlation matrix

#output---
#F1treeA: matrix filled only in 1st row and diagonal
#        value in each row are connected with its corresponding 
#        diagnoal value.(vine array)
#Pcor=pmat_f1: Partial correlation matrix (Y_j1Y_j2|V_1V_2)

select1FTr.partial = function(y,rmat){
  d=ncol(y)
  n=nrow(y)
  #rmat=polychoric0(y)$p
  nq=15
  gl=gauss.quad.prob(nq)
  #----
  fa1 <- mle1factor(NULL, y, NULL, rep("bvn",d), gl)
  fa1_load <- unlist(fa1$cpar$f1)
  #partial corr
  #Y_j1Y_j2|V_0
  ilength = (d*(d-1)/2)
  pmat_f1=matrix(NA, nrow = d, ncol = d)
  for(i in 1:ilength){
    cmb=combn(1:d,2)
    cmb=t(cmb)
    j1=cmb[i,1]
    j2=cmb[i,2]
    
    pj1j2=rmat[j1,j2]
    pj1V = fa1_load[j1]
    pj2V = fa1_load[j2]
    pmat_f1[j2,j1] = (pj1j2 - (pj1V*pj2V))/(sqrt((1-pj1V^2)*(1-pj2V^2)))
    pmat_f1[j1,j2] = pmat_f1[j2,j1]
    diag(pmat_f1)=1
  }
  
  #get values in the lower tri
  pmat.valf1=pmat_f1[lower.tri(pmat_f1)]

  #weights of edges
  wghtf1 = log(1-pmat.valf1^2)
  arcsf1 = cbind(t(combn(1:d,2)) , wghtf1)
  
  # Minimum cost spanning tree
  outmstf1=getMinimumSpanningTree(1:d, arcsf1, algorithm = "Prim",show.graph = F)

  A1 = matrix(NA,d,d)
  diag(A1) = outmstf1$tree.nodes
  A1[1,2:d]= outmstf1$tree.arcs[,1]
  A1[is.na(A1)] = 0
  
  return(list(F1treeA=A1,Pcor=pmat_f1))
}


