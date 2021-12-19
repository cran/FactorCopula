#selection algorithm for the second-order/nested factor copuls

Select.nested = function(ydat,cutp, ngrp,grpsize, allcop,
                         gl,gln, glw,nq,n,d, SpC = NULL){

    #====== The possible copulas for all variables
    length.allcop = length(allcop)
    # creating a list for common factor
    copMatF0=matrix(allcop,nrow = length.allcop,ncol = ngrp)

    grplst=list()
    ind=1
    onegrp=rep(FALSE,ngrp)
    for(g in 1:ngrp){
      grplst[[g]]=ind:cumsum(grpsize)[g]
      ind=ind+grpsize[g]
    }

    coplngthF0=coplngthFg=length.allcop

    #starting copulas
    copF0= rep("bvn",ngrp)
    copFg= rep("bvn",d)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    #========================================================================
    inner0<-matrix(NA, nrow = coplngthF0, ncol = ngrp+1)
    for(i in 1:coplngthF0)
    {
      print(c("iF0",i))
      copF0<-copMatF0[i,]
      initpar<-initval(copF0,copFg)
      modF0=list(estimate=rep(NA,d+ngrp),minimum=NA)
      try(modF0<-nlm(nestedllk,initpar,ydat=ydat,cutp=cutp,
                     copX0=copF0,copXgY=copFg,
                     nq=nq,ngrp=ngrp,grpsize=grpsize,glw, gln,SpC,param=F,
                     hessian=F),silent=TRUE)
      #AIC
      aic.F0 = 2*modF0$minimum+2*length(modF0$e)
      inner0[i, ] = c(aic.F0, copF0)
    }
    index0<-which.min(cbind(as.numeric(noquote(rbind(inner0))[,1])))
    copF0<-inner0[index0,-1]
    for(j in 1:ngrp){
      group_items=grplst[[j]]
      print(c("Group",j))
      
      innerg = matrix(NA, nrow = coplngthFg, ncol = 2*d+ngrp+1)
      for( ii in 1:coplngthFg)
      {
        print(c("iFg",ii))
        copFg[group_items]<-rep(allcop[ii],grpsize[j])
        initpar<-initval(copF0,copFg)
        modFg=list(estimate=rep(NA,d+ngrp),minimum=NA)
        try(modFg<-nlm(nestedllk,initpar,ydat=ydat,cutp=cutp,
                       copX0=copF0,copXgY=copFg,
                       nq=nq,ngrp=ngrp,grpsize=grpsize,glw, gln,SpC,param=F,
                       hessian=F),silent=TRUE)
        #AIC
        aic.Fg= 2*modFg$minimum+2*length(modFg$e)
        innerg[ii, ]=c(aic.Fg, copFg, modFg$e)
      }
      indexg=which.min(cbind(as.numeric(noquote(rbind(innerg))[,1])))
      copFg[group_items]<- innerg[indexg,c(group_items+1)]

      AIC.nested<- innerg[indexg,1]
      d.innerg=ncol(innerg)
      estimated.taus<- as.numeric(innerg[indexg,(2+d):d.innerg])
    }
    return(list("common factor" = copF0, "group-specific factor" = copFg,
                "AIC" = AIC.nested, "estimated taus" = estimated.taus))
  }


