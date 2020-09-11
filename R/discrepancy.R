# This code is to calculate the discrepancy measures

discrepancy = function(cormat, n, f3 = F){
fa1 <- factanal(covmat = cormat, factors = 1, n.obs = n, rotation = "varimax")
dif1=round(cormat - (fa1$loadings %*% t(fa1$loadings) + diag(fa1$uniquenesses)), 3)
max.ad.f1=max(abs(dif1))
avg.ad.f1=mean(abs(dif1))

fa2 <- factanal(covmat = cormat, factors = 2, n.obs = n, rotation = "varimax")
dif2=round(cormat - (fa2$loadings %*% t(fa2$loadings) + diag(fa2$uniquenesses)), 3)
max.ad.f2=max(abs(dif2))
avg.ad.f2=mean(abs(dif2))

if(f3){
fa3 <- factanal(covmat = cormat, factors = 3, n.obs = n, rotation = "varimax")
dif3=round(cormat - (fa3$loadings %*% t(fa3$loadings) + diag(fa3$uniquenesses)), 3)
max.ad.f3=max(abs(dif3))
avg.ad.f3=mean(abs(dif3))
}else{
  max.ad.f3=NA
  avg.ad.f3=NA
}
#====================================================================  
#
d = ncol(cormat)
modMat1=(fa1$loadings %*% t(fa1$loadings) + diag(fa1$uniquenesses))
deltahat1=log( det(modMat1) ) - log( det(cormat) ) +
  sum(diag((solve(modMat1)%*%cormat))) - d

modMat2=(fa2$loadings %*% t(fa2$loadings) + diag(fa2$uniquenesses))
deltahat2=log( det(modMat2) ) - log( det(cormat) ) +
  sum(diag((solve(modMat2)%*%cormat))) - d

if(f3){
modMat3=(fa3$loadings %*% t(fa3$loadings) + diag(fa3$uniquenesses))
deltahat3=log( det(modMat3) ) - log( det(cormat) ) +
  sum(diag((solve(modMat3)%*%cormat))) - d
}else{
  deltahat3 = NA
}
#====================================================================   
discrepancy.f1=c(max.ad.f1,avg.ad.f1,deltahat1)
discrepancy.f2=c(max.ad.f2,avg.ad.f2,deltahat2)
discrepancy.f3=c(max.ad.f3,avg.ad.f3,deltahat3)
Dmat = rbind(discrepancy.f1, discrepancy.f2, discrepancy.f3)
colnames(Dmat) = c("D1", "D2", "D3")
rownames(Dmat) = c("f1", "f2", "f3")
return(Dmat)
}




