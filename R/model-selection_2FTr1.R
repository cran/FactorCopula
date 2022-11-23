#Vine Tree selection for 2-factor tree
#library(optrees) for getMinimumSpanningTree

#Input:
#y:ordinal data
#rmat:polychoric correlation

#Outut:
#F2treeA=A2: vine array,
#matrix filled only in 1st row and diagonal
#        value in each row are connected with its corresponding 
#        diagnoal value
#Pcor=pmat_f2: Partial correlation matrix (Y_j1Y_j2|V_1V_2)
select2FTr.partial = function(y,rmat){
  d=ncol(y)
  n=nrow(y)
  #rmat=polychoric0(y)$p
  #fa1 <- factanal(covmat = rmat, factors = 1, n)
  #fa2 <- factanal(covmat = rmat, factors = 2, n, rotation = "varimax")
  nq=15
  gl=gauss.quad.prob(nq)
  #----
  fa2 <- mle2factor.bvn( NULL, y, NULL, rep("bvn",d), rep("bvn",d), gl,SpC=d,print.level = 0)
  fa2_load1 <-unlist(fa2$cpar$f1)
  fa2_pc2 <- unlist(fa2$cpar$f2)
  load2=fa2_pc2*sqrt(1-fa2_load1^2)
  #partial corr
  #Y_j1Y_j2|V_0
  ilength = (d*(d-1)/2)
  pmat_f2=matrix(NA, nrow = d, ncol = d)
  for(i in 1:ilength){
    cmb=combn(1:d,2)
    cmb=t(cmb)
    j1=cmb[i,1]
    j2=cmb[i,2]
    
    pj1j2=rmat[j1,j2]
    
    #2-factor 
    pj1V_2f = fa2_load1[j1]
    pj2V_2f = fa2_load1[j2]
    pj1pj2_f2 = (pj1j2 - (pj1V_2f*pj2V_2f))/sqrt( (1-pj1V_2f^2)*(1-pj2V_2f^2) )
    
    partial_f2_j1 = fa2_pc2[j1]#/sqrt(1-fa2_load1[j1]^2)
    partial_f2_j2 = fa2_pc2[j2]#/sqrt(1-fa2_load1[j2]^2)
    pmat_f2[j2,j1] = (pj1pj2_f2 - (partial_f2_j1 * partial_f2_j2)
    )/( sqrt((1-partial_f2_j1^2)*(1-partial_f2_j2^2)))
    pmat_f2[j1,j2] = pmat_f2[j2,j1]
    diag(pmat_f2)=1
  }
  
  #get values in the lower tri
  pmat.valf2=pmat_f2[lower.tri(pmat_f2)]
  wghtf2 = log(1-pmat.valf2^2)
  arcsf2=cbind(t(combn(1:d,2)) , wghtf2)
  # Minimum cost spanning tree
  outmstf2=getMinimumSpanningTree(1:d, arcsf2, algorithm = "Prim",show.graph = F)
  
  A2 = matrix(NA,d,d)
  diag(A2) = outmstf2$tree.nodes
  A2[1,2:d]= outmstf2$tree.arcs[,1]
  A2[is.na(A2)] = 0
  
  return(list(F2treeA=A2,Pcor=pmat_f2))
}
