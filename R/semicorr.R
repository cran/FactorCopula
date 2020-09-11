#(semi-)correlations

semicorr=function(dat,type){
  d=ncol(dat)
  n=nrow(dat)
  out=NULL
  for(j1 in 1:(d-1)){
    for(j2 in (j1+1):d){
      #print(c(j1,j2))
      bivdat=apply(dat[,c(j1,j2)],2,as.numeric)
      if(type[j1]==1 & type[j2]==1){
        u=apply(bivdat,2,rank)/(n+1)
        bivdat=qnorm(u)
        corr=polycor::hetcor(bivdat[,1:2])$c[1,2]
        if(corr>0){
          i1 = (bivdat[, 1] < 0 & bivdat[, 2] < 0)
          i2 = (bivdat[, 1] > 0 & bivdat[, 2] > 0)
          lcorr=polycor::hetcor(bivdat[i1,1:2])$c[1,2]
          ucorr=polycor::hetcor(bivdat[i2,1:2])$c[1,2]
        } else {
          i1 = (bivdat[, 1] > 0 & bivdat[, 2] < 0)
          i2 = (bivdat[, 1] < 0 & bivdat[, 2] > 0)
          lcorr=polycor::hetcor(bivdat[i1,1:2])$c[1,2]
          ucorr=polycor::hetcor(bivdat[i2,1:2])$c[1,2]
        }
      } else {
        if(type[j1]==1 & type[j2]==2){
          u=rank(bivdat[,1])/(n+1)
          zdat=qnorm(u)
          bivdat=cbind(zdat,bivdat[,2])
          thres=median(unique(as.numeric(bivdat[,2])))
          corr=polycor::polyserial(bivdat[,1], bivdat[,2])
          if(corr>0)
          { i1 = (bivdat[,1]< 0 & bivdat[,2]<= thres)
          i2 = (bivdat[,1]> 0 & bivdat[,2]>= thres)
          lcorr=polycor::polyserial(bivdat[i1,1], bivdat[i1,2])
          ucorr=polycor::polyserial(bivdat[i2,1], bivdat[i2,2])
          }
          else {
            i1 = (bivdat[,1]> 0 & bivdat[,2]<= thres)
            i2 = (bivdat[,1]< 0 & bivdat[,2]>= thres)
            lcorr=polycor::polyserial(bivdat[i1,1], bivdat[i1,2])
            ucorr=polycor::polyserial(bivdat[i2,1], bivdat[i2,2])
            
          }} else {
            if(type[j1]==2 & type[j2]==2)
            { thres=c(median(unique(bivdat[,1])),median(unique(bivdat[,2])))
            #corr=mixedCor(as.data.frame(bivdat),p=1:2, correct=0, global=FALSE)$rho[1,2]
            corr=polycor::polychor(bivdat[,1], bivdat[,2])
            
            if(corr>0)
            { i1 = (bivdat[,1]<= thres[1] & bivdat[,2]<= thres[2])
            i2 =   (bivdat[,1]>= thres[1] & bivdat[,2]>= thres[2])
            lcorr=polycor::polychor(bivdat[i1,1], bivdat[i1,2])
            ucorr=polycor::polychor(bivdat[i2,1], bivdat[i2,2])
            
            }
            else {
              i1 = (bivdat[,1]>= thres[1] & bivdat[,2]<= thres[2])
              i2 = (bivdat[,1]<= thres[1] & bivdat[,2]>= thres[2])
              lcorr=polycor::polychor(bivdat[i1,1], bivdat[i1,2])
              ucorr=polycor::polychor(bivdat[i2,1], bivdat[i2,2])
            }}
          }}
      out = rbind(out,c(j1,j2,corr,lcorr,ucorr))
  }}
  colnames(out) = c(" ", " ", "rho", "lrho", "urho")
  out
}
