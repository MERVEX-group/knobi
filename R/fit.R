fitting<-function(Data){
  UseMethod("fitting",Data)
}


fitting.Pella<-function(Data){

  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_p<-Data$start_p

  my_model <- function(x, r, K, p) {
    (r/p)*x*(1-(x/K)^(p))
  }

  params <- c(r=start_r ,K=start_K ,p=start_p)

  fun <- function(p,x,y) {
    sum(
      (y-my_model(x, p["r"], p["K"], p["p"]))^2

    )
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",upper=c(2,Inf,3.5),lower=c(0.05,0,0.25))
  return(out)
}


fitting.Schaefer<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K

  my_model <- function(x, r, K) {
    (r*x*(1-(x/K)))
  }

  params <- c(r=start_r ,K=start_K)

  fun <- function(p,x,y) {
    sum(
      (y-my_model(x, p["r"], p["K"]))^2

    )
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",upper=c(2,Inf),lower=c(0.05,0))
  return(out)
}

fitting.Schaefer_Mult<-function(Data){

  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c


  my_model <- function(x,env, r, K,c) {
    exp(1)^{c*env}*((r*x*(1-(x/K))))
  }

  params <- c(r=start_r ,K=start_K,c=start_c)
  # minimize total sum of squares of residuals
  fun <- function(p,x,y,env) {
    sum(
      (y-my_model(x,env, p["r"], p["K"],p["c"]))^2

    )
  }

  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,20),lower=c(0.05,0,-20),env=as.numeric(Data$env))
  return(out$par)

}

fitting.Schaefer_Add<-function(Data){

  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c


  my_model <- function(x,env, r, K,c) {
    ((r*x*(1-(x/K)))+c*env*x)
  }

  params <- c(r=start_r ,K=start_K,c=start_c)
  # minimize total sum of squares of residuals
  fun <- function(p,x,y,env) {
    sum(
      (y-my_model(x,env, p["r"], p["K"],p["c"]))^2

    )
  }

  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,1000000),lower=c(0.05,0,-1000000),env=as.numeric(Data$env))
  return(out$par)
}


fitting.Pella_Mult<-function(Data){

  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c
  start_p<-Data$start_p


  my_model <- function(x,env, r, K,c,p) {
    exp(1)^{c*env}*((r/p)*x*(1-(x/K)^(p)))
  }

  params <- c(r=start_r ,K=start_K,c=start_c,p=start_p)
  # minimize total sum of squares of residuals
  fun <- function(p,x,y,env) {
    sum(
      (y-my_model(x,env, p["r"], p["K"],p["c"],p["p"]))^2

    )
  }

  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,20,3.5),lower=c(0.05,0,-20,0.25),env=as.numeric(Data$env))
  return(out$par)
}

fitting.Pella_Add<-function(Data){

  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c
  start_p<-Data$start_p


  my_model <- function(x,env, r, K,c,p) {
    ((r/p)*x*(1-(x/K)^(p))+c*env*x)
  }

  params <- c(r=start_r ,K=start_K,c=start_c,p=start_p)
  # minimize total sum of squares of residuals
  fun <- function(p,x,y,env) {
    sum(
      (y-my_model(x,env, p["r"], p["K"],p["c"],p["p"]))^2

    )
  }

  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,1000000,3.5),lower=c(0.05,0,-1000000,0.25),env=as.numeric(Data$env))
  return(out$par)
}


fitting.Pella_Mult_2<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_p<-Data$start_p
  start_c<-Data$start_c

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]

  my_model <- function(x,env1,env2, r, K,c1,c2,p) {
    exp(1)^{c1*env1+c2*env2}*((r/p)*x*(1-(x/K)^(p)))
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c,p=start_p)

  fun <- function(p,x,y,env1,env2) {
    sum(
      (y-my_model(x,env1,env2, p["r"], p["K"],p["c1"],p["c2"],p["p"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(20,ncol(env)),3.5),lower=c(0.05,0,rep(-20,ncol(env)),0.25),
                        env1=env1,env2=env2)

  return(out$par)
}


fitting.Pella_Mult_3<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_p<-Data$start_p
  start_c<-Data$start_c

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]; env3<-env[,3]

  my_model <- function(x,env1,env2, env3, r, K,c1,c2,c3,p) {
    exp(1)^{c1*env1+c2*env2+c3*env3}*((r/p)*x*(1-(x/K)^(p)))
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c,c3=start_c,p=start_p)

  fun <- function(p,x,y,env1,env2,env3) {
    sum(
      (y-my_model(x,env1,env2,env3, p["r"], p["K"],p["c1"],p["c2"],p["c3"],p["p"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(20,ncol(env)),3.5),lower=c(0.05,0,rep(-20,ncol(env)),0.25),
                        env1=env1,env2=env2,env3=env3)
  return(out$par)
}


fitting.Pella_Mult_4<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_p<-Data$start_p
  start_c<-Data$start_c

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]; env3<-env[,3]; env4<-env[,4]

  my_model <- function(x,env1,env2, env3, env4,r, K,c1,c2,c3,c4,p) {
    exp(1)^{c1*env1+c2*env2+c3*env3+c4*env4}*((r/p)*x*(1-(x/K)^(p)))
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c,c3=start_c,c4=start_c,p=start_p)

  fun <- function(p,x,y,env1,env2,env3,env4) {
    sum(
      (y-my_model(x,env1,env2,env3,env4, p["r"], p["K"],p["c1"],p["c2"],p["c3"],p["c4"],p["p"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(20,ncol(env)),3.5),lower=c(0.05,0,rep(-20,ncol(env)),0.25),
                        env1=env1,env2=env2,env3=env3,env4=env4)
  return(out$par)
}


fitting.Pella_Mult_5<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_p<-Data$start_p
  start_c<-Data$start_c

  env<-Data$env


  env1<-env[,1]; env2<-env[,2]; env3<-env[,3]; env4<-env[,4]; env5<-env[,5]

  my_model <- function(x,env1,env2, env3, env4, env5, r, K,c1,c2,c3,c4,c5,p) {
    exp(1)^{c1*env1+c2*env2+c3*env3+c4*env4+c5*env5}*((r/p)*x*(1-(x/K)^(p)))
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c,c3=start_c,c4=start_c,c5=start_c,p=start_p)

  fun <- function(p,x,y,env1,env2,env3,env4,env5) {
    sum(
      (y-my_model(x,env1,env2,env3,env4,env5, p["r"], p["K"],p["c1"],p["c2"],p["c3"],p["c4"],p["c5"],p["p"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(20,ncol(env)),3.5),lower=c(0.05,0,rep(-20,ncol(env)),0.25),
                        env1=env1,env2=env2,env3=env3,env4=env4,env5=env5)

  return(out$par)
}


fitting.Schaefer_Mult_2<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]

  my_model <- function(x,env1,env2, r, K,c1,c2) {
    exp(1)^{c1*env1+c2*env2}*((r*x*(1-(x/K))))
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c)

  fun <- function(p,x,y,env1,env2) {
    sum(
      (y-my_model(x,env1,env2, p["r"], p["K"],p["c1"],p["c2"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(20,ncol(env))),lower=c(0.05,0,rep(-20,ncol(env))),
                        env1=env1,env2=env2)
  return(out$par)
}



fitting.Schaefer_Mult_3<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]; env3<-env[,3]

  my_model <- function(x,env1,env2, env3, r, K,c1,c2,c3) {
    exp(1)^{c1*env1+c2*env2+c3*env3}*((r*x*(1-(x/K))))
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c,c3=start_c)

  fun <- function(p,x,y,env1,env2,env3) {
    sum(
      (y-my_model(x,env1,env2,env3, p["r"], p["K"],p["c1"],p["c2"],p["c3"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(20,ncol(env))),lower=c(0.05,0,rep(-20,ncol(env))),
                        env1=env1,env2=env2,env3=env3)
  return(out$par)
}


fitting.Schaefer_Mult_4<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]; env3<-env[,3]; env4<-env[,4]

  my_model <- function(x,env1,env2, env3, env4,r, K,c1,c2,c3,c4) {
    exp(1)^{c1*env1+c2*env2+c3*env3+c4*env4}*((r*x*(1-(x/K))))
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c,c3=start_c,c4=start_c)

  fun <- function(p,x,y,env1,env2,env3,env4) {
    sum(
      (y-my_model(x,env1,env2,env3,env4, p["r"], p["K"],p["c1"],p["c2"],p["c3"],p["c4"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(20,ncol(env))),lower=c(0.05,0,rep(-20,ncol(env))),
                        env1=env1,env2=env2,env3=env3,env4=env4)
  return(out$par)
}


fitting.Schaefer_Mult_5<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]; env3<-env[,3]; env4<-env[,4]; env5<-env[,5]

  my_model <- function(x,env1,env2, env3, env4, env5, r, K,c1,c2,c3,c4,c5) {
    exp(1)^{c1*env1+c2*env2+c3*env3+c4*env4+c5*env5}*((r*x*(1-(x/K))))
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c,c3=start_c,c4=start_c,c5=start_c)

  fun <- function(p,x,y,env1,env2,env3,env4,env5) {
    sum(
      (y-my_model(x,env1,env2,env3,env4,env5, p["r"], p["K"],p["c1"],p["c2"],p["c3"],p["c4"],p["c5"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(20,ncol(env))),lower=c(0.05,0,rep(-20,ncol(env))),
                        env1=env1,env2=env2,env3=env3,env4=env4,env5=env5)
  return(out$par)
}


fitting.Pella_Add_2<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c
  start_p<-Data$start_p

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]

  my_model <- function(x,env1,env2, r, K,c1,c2,p) {
    ((r/p)*x*(1-(x/K)^(p))+(c1*env1+c2*env2)*x)
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c,p=start_p)

  fun <- function(p,x,y,env1,env2) {
    sum(
      (y-my_model(x,env1,env2, p["r"], p["K"],p["c1"],p["c2"],p["p"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(1000000,ncol(env)),3.5),lower=c(0.05,0,rep(-1000000,ncol(env)),0.25),
                        env1=env1,env2=env2)
  return(out$par)
}


fitting.Pella_Add_3<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c
  start_p<-Data$start_p

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]; env3<-env[,3]

  my_model <- function(x,env1,env2, env3, r, K,c1,c2,c3,p) {
    ((r/p)*x*(1-(x/K)^(p))+(c1*env1+c2*env2+c3*env3)*x)
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c,c3=start_c,p=start_p)

  fun <- function(p,x,y,env1,env2,env3) {
    sum(
      (y-my_model(x,env1,env2,env3, p["r"], p["K"],p["c1"],p["c2"],p["c3"],p["p"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(1000000,ncol(env)),3.5),lower=c(0.05,0,rep(-1000000,ncol(env)),0.25),
                        env1=env1,env2=env2,env3=env3)
  return(out$par)
}


fitting.Pella_Add_4<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c
  start_p<-Data$start_p

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]; env3<-env[,3]; env4<-env[,4]

  my_model <- function(x,env1,env2, env3, env4,r, K,c1,c2,c3,c4,p) {
    ((r/p)*x*(1-(x/K)^(p))+(c1*env1+c2*env2+c3*env3+c4*env4)*x)
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c,c3=start_c,c4=start_c,p=start_p)

  fun <- function(p,x,y,env1,env2,env3,env4) {
    sum(
      (y-my_model(x,env1,env2,env3,env4, p["r"], p["K"],p["c1"],p["c2"],p["c3"],p["c4"],p["p"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(1000000,ncol(env)),3.5),lower=c(0.05,0,rep(-1000000,ncol(env)),0.25),
                        env1=env1,env2=env2,env3=env3,env4=env4)
  return(out$par)
}

fitting.Pella_Add_5<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c
  start_p<-Data$start_p

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]; env3<-env[,3]; env4<-env[,4]; env5<-env[,5]

  my_model <- function(x,env1,env2, env3, env4, env5, r, K,c1,c2,c3,c4,c5,p) {
    ((r/p)*x*(1-(x/K)^(p))+(c1*env1+c2*env2+c3*env3+c4*env4+c5*env5)*x)
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c,c3=start_c,c4=start_c,c5=start_c,p=start_p)

  fun <- function(p,x,y,env1,env2,env3,env4,env5) {
    sum(
      (y-my_model(x,env1,env2,env3,env4,env5, p["r"], p["K"],p["c1"],p["c2"],p["c3"],p["c4"],p["c5"],p["p"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(1000000,ncol(env)),3.5),lower=c(0.05,0,rep(-1000000,ncol(env)),0.25),
                        env1=env1,env2=env2,env3=env3,env4=env4,env5=env5)
  return(out$par)
}



fitting.Schaefer_Add_2<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]

  my_model <- function(x,env1,env2, r, K,c1,c2) {
    ((r*x*(1-(x/K)))+(c1*env1+c2*env2))*x
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c)

  fun <- function(p,x,y,env1,env2) {
    sum(
      (y-my_model(x,env1,env2, p["r"], p["K"],p["c1"],p["c2"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(1000000,ncol(env))),lower=c(0.05,0,rep(-1000000,ncol(env))),
                        env1=env1,env2=env2)
  return(out$par)
}


fitting.Schaefer_Add_3<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]; env3<-env[,3]

  my_model <- function(x,env1,env2, env3, r, K,c1,c2,c3) {
    ((r*x*(1-(x/K)))+(c1*env1+c2*env2+c3*env3)*x)
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c,c3=start_c)

  fun <- function(p,x,y,env1,env2,env3) {
    sum(
      (y-my_model(x,env1,env2,env3, p["r"], p["K"],p["c1"],p["c2"],p["c3"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(1000000,ncol(env))),lower=c(0.05,0,rep(-1000000,ncol(env))),
                        env1=env1,env2=env2,env3=env3)
  return(out$par)
}


fitting.Schaefer_Add_4<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]; env3<-env[,3]; env4<-env[,4]

  my_model <- function(x,env1,env2, env3, env4,r, K,c1,c2,c3,c4) {
    ((r*x*(1-(x/K)))+(c1*env1+c2*env2+c3*env3+c4*env4)*x)
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c,c3=start_c,c4=start_c)

  fun <- function(p,x,y,env1,env2,env3,env4) {
    sum(
      (y-my_model(x,env1,env2,env3,env4, p["r"], p["K"],p["c1"],p["c2"],p["c3"],p["c4"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(1000000,ncol(env))),lower=c(0.05,0,rep(-1000000,ncol(env))),
                        env1=env1,env2=env2,env3=env3,env4=env4)
  return(out$par)
}


fitting.Schaefer_Add_5<-function(Data){
  data<-Data$data
  start_r<-Data$start_r
  start_K<-Data$start_K
  start_c<-Data$start_c

  env<-Data$env

  env1<-env[,1]; env2<-env[,2]; env3<-env[,3]; env4<-env[,4]; env5<-env[,5]

  my_model <- function(x,env1,env2, env3, env4, env5, r, K,c1,c2,c3,c4,c5) {
    ((r*x*(1-(x/K)))+(c1*env1+c2*env2+c3*env3+c4*env4+c5*env5)*x)
  }
  params <- c(r=start_r ,K=start_K,c1=start_c,c2=start_c,c3=start_c,c4=start_c,c5=start_c)

  fun <- function(p,x,y,env1,env2,env3,env4,env5) {
    sum(
      (y-my_model(x,env1,env2,env3,env4,env5, p["r"], p["K"],p["c1"],p["c2"],p["c3"],p["c4"],p["c5"]))^2)
  }
  out <- optimx::optimr(params,fun,x=data$x,y=data$y,method="L-BFGS-B",
                        upper=c(2,Inf,rep(1000000,ncol(env))),lower=c(0.05,0,rep(-1000000,ncol(env))),
                        env1=env1,env2=env2,env3=env3,env4=env4,env5=env5)
  return(out$par)
}
