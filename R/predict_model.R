
predict_model=function(fit){
  UseMethod("predict_model",fit)
}


predict_model.Pella=function(fit){
  x=fit$data[[2]]
  r=fit[[1]][1]
  K=fit[[1]][2]
  p=fit[[1]][3]
  cut=K
  av <- seq(0, cut, length.out = 3*length(x))
  (r/p)*av*(1-(av/K)^(p))
}

predict_model.Schaefer=function(fit){
  x=fit$data[[2]]
  r=fit[[1]][1]
  K=fit[[1]][2]
  cut=K
  av <- seq(0, cut, length.out = 3*length(x))
  (r)*av*(1-(av/K))
}


predict_model.Schaefer_Year=function(fit){
  x=fit$data[[2]]
  r=fit$Parameter_estimates[1]
  K=fit$Parameter_estimates[2]

  as.numeric((r)*x*(1-(x/K)))
}


predict_model.Pella_Year=function(fit){
  x=fit$data[[2]]
  r=fit$Parameter_estimates[1]
  K=fit$Parameter_estimates[2]
  p=fit$Parameter_estimates[3]
  as.numeric((r/p)*x*(1-(x/K)^(p)))
}

predict_model.Pella_Mult=function(fit){
  x=fit$data[[2]]
  data=fit$data
  env=data$env
  r=fit[[1]][1]
  K=fit[[1]][2]
  c=fit[[1]][3]
  p=fit[[1]][4]
  as.numeric(exp(1)^{c*env}*((r/p)*x*(1-(x/K)^(p))))
}

predict_model.Pella_Add=function(fit){
  x=fit$data[[2]]
  data=fit$data
  env=data$env
  r=fit[[1]][1]
  K=fit[[1]][2]
  c=fit[[1]][3]
  p=fit[[1]][4]
  as.numeric((r/p)*x*(1-(x/K)^(p))+c*env*x)
}

predict_model.Schaefer_Mult=function(fit){
  x=fit$data[[2]]
  data=fit$data
  env=data$env
  r=fit[[1]][1]
  K=fit[[1]][2]
  c=fit[[1]][3]
  as.numeric(exp(1)^{c*env}*((r)*x*(1-(x/K))))
}

predict_model.Schaefer_Add=function(fit){
  x=fit$data[[2]]
  data=fit$data
  env=data$env
  r=fit[[1]][1]
  K=fit[[1]][2]
  c=fit[[1]][3]
  as.numeric((r)*x*(1-(x/K))+c*env*x)
}





predict_model.Pella_Mult_m=function(fit){
  x=fit$data[[2]]
  data=fit$data
  env=as.matrix(data$env)
  c_length=ncol(fit$data$env)
  r=fit[[1]][1]
  K=fit[[1]][2]
  c=fit[[1]][3:(2+c_length)]
  p=fit[[1]][length(fit[[1]])]
  as.numeric(exp(1)^{env %*% c}*((r/p)*x*(1-(x/K)^(p))))
}

predict_model.Pella_Add_m=function(fit){
  x=fit$data[[2]]
  data=fit$data
  env=as.matrix(data$env)
  c_length=ncol(fit$data$env)
  r=fit[[1]][1]
  K=fit[[1]][2]
  c=fit[[1]][3:(2+c_length)]
  p=fit[[1]][length(fit[[1]])]
  as.numeric((r/p)*x*(1-(x/K)^(p)) + (env %*% c)*x)
}

predict_model.Schaefer_Mult_m=function(fit){
  x=fit$data[[2]]
  data=fit$data
  env=as.matrix(data$env)
  r=fit[[1]][1]
  K=fit[[1]][2]
  c=fit[[1]][3:length(fit[[1]])]
  as.numeric(exp(1)^{env %*% c}*((r)*x*(1-(x/K))))
}

predict_model.Schaefer_Add_m=function(fit){
  x=fit$data[[2]]
  data=fit$data
  env=as.matrix(data$env)
  r=fit[[1]][1]
  K=fit[[1]][2]
  c=fit[[1]][3:length(fit[[1]])]
  as.numeric((r)*x*(1-(x/K)) + (env %*% c)*x)
}





