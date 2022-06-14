
RF=function(fit){
  UseMethod("RF",fit)
}


RF.Pella=function(fit){

  r=fit[[1]][1]
  K=fit[[1]][2]
  p=fit[[1]][3]
  K<-as.numeric(K)
  Bmsy<-as.numeric(K*(1/(p+1))^(1/p))
  MSY<-as.numeric(r*K*(1/(p+1))^(1+1/p))
  x<-as.numeric(MSY/K)
  Fmsy=as.numeric((r/p)*(1-(1/(p+1))))

  fit$K=K
  fit$B_MSY=Bmsy
  fit$F_MSY=Fmsy
  fit$MSY=MSY
  fit$MSYoverK=x
  fit
}

RF.Schaefer=function(fit){

  r=fit[[1]][1]
  K=fit[[1]][2]

  K<-as.numeric(K)
  Bmsy<-as.numeric(K/2)
  MSY<-as.numeric((1/4)*r*K)
  x<-as.numeric(MSY/K)
  Fmsy=as.numeric(r/2)

  fit$K=K
  fit$B_MSY=Bmsy
  fit$F_MSY=Fmsy
  fit$MSY=MSY
  fit$MSYoverK=x
  fit
}

RF.Pella_Mult = RF.Pella_Add = RF.Pella_Mult_m = RF.Pella_Add_m = function(fit){

  r=fit[1]
  K=fit[2]
  p=fit[length(fit)]
  K<-as.numeric(K)
  Bmsy<-as.numeric(K*(1/(p+1))^(1/p))
  MSY<-as.numeric(r*K*(1/(p+1))^(1+1/p))
  x<-as.numeric(MSY/K)
  Fmsy=as.numeric(r/(p+1))

  ref_pts=data.frame(K,Bmsy,Fmsy,MSY,MSYoverK=x)
  ref_pts
}

RF.Schaefer_Mult = RF.Schaefer_Add = RF.Schaefer_Mult_m = RF.Schaefer_Add_m = function(fit){

  r=fit[1]
  K=fit[2]

  K<-as.numeric(K)
  Bmsy<-as.numeric(K/2)
  MSY<-as.numeric((1/4)*r*K)
  x<-as.numeric(MSY/K)
  Fmsy=as.numeric(r/2)

  ref_pts=data.frame(K,Bmsy,Fmsy,MSY,MSYoverK=x)
  ref_pts
}

