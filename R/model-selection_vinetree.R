# Vine tree selection using Prim's algorithm.
#library(optrees)

#Input---
#y: n \times d matrix  with ordinal data 
#rmat: polychoric correlation matrix

#output---
#VineTreeA: vine array,
#matrix filled only in 1st row and diagonal
#        value in each row are connected with its corresponding 
#        diagnoal value

selectVineTree = function(y,rmat){
  d=ncol(y)
  n=nrow(y)
  #rmat=polychoric0(y)$p
  #get values in the lower tri
  pmat.valf1=rmat[lower.tri(rmat)]
  
  #weights of edges
  wghtf1 = log(1-pmat.valf1^2)
  arcsf1=cbind(t(combn(1:d,2)) , wghtf1)
  
  # Minimum cost spanning tree
  outmstf1=getMinimumSpanningTree(1:d, arcsf1, algorithm = "Prim",show.graph = F)
  
  A1 = matrix(NA,d,d)
  diag(A1) = outmstf1$tree.nodes
  A1[1,2:d]= outmstf1$tree.arcs[,1]
  A1[is.na(A1)] = 0
  
  return(list(VineTreeA=A1))
}

