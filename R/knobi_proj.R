#' @title KBPM projections
#'
#' @description This function projects the stock biomass  (or stock spawning biomass) time series and then the surplus production  based on the selected catch or fishing mortality values for the projected years.
#'
#' @param knobi_results The output object of \code{\link{knobi_fit}} function (main function).
#' @param env_results Optional. The output object of \code{\link{knobi_env}} function.
#' @param n_y Optional. The number of years for projections. If not provided, enter end_y in the next argument.
#' @param end_y Optional. The end year of our projections. If not provided, enter n_y in the argument above. If both are provided, the function uses n_y.
#' @param Ct Optional. Vector, data frame or matrix containing the values of the catch for the projected. Different catch scenarios are allowed. The values for each should be provided in each of the columns. The length vector or row number should be equal to the number of years. Projections can be based on selected catch or fishing mortality values, then only one of the arguments, Ct or f, must be introduced.
#' @param f Optional. Vector, data frame or matrix containing the values of the fishing mortality for the projected. Different catch scenarios are allowed. The values for each should be provided in each of the columns. The length vector or row number should be equal to the number of years. Projections can be based on selected catch or fishing mortality values, then only one of the arguments, Ct or f, must be introduced.
#' @param env Optional. If the multicovar argument of \code{\link{knobi_env}} is FALSE, a vector, data frame or matrix containing the values of the environmental covariates (unstandardized) for the projection years (rows) and the different catch or fishing mortality settings (columns).  On the other hand, if the multicovar argument of \code{\link{knobi_env}} is TRUE, the current argument must be a list, and each argument must be a data frame corresponding to each catch or fishing mortality setting containing the values of the environmental covariates for that scenario.
#' @param plot_out Logical. TRUE means that a file with the  environmental fit plots is created. By default this argument is FALSE.
#' @param plot_dir Optional directory for creating the folder and save the plots. Required when plot_out=TRUE. The default value is the input of this argument in the knobi_fit function.
#' @param plot_filename Optional name of the folder that will contain the plots. Required when plot_out=TRUE. The default value is the input of this argument in the knobi_fit function.
#'
#' @return A list containing the projection results is provided. \itemize{
#' \item base_model: Three-dimensional array containing projection information for each catch or fishing mortality scenario that is projected using the KBPM base model (see \code{\link{knobi_fit}} details). Rows correspond to years, columns to the derived quantities (biomass or SSB, surplus production, F and catch) and third dimension to the projection scenarios.
#' \item additive_model: Four-dimensional array containing projections information for each catch or fishing mortality scenario and for each environmental scenario considering the KBPM additive model (see \code{\link{knobi_env}} details). The rows correspond to the years, the columns to the derived quantities (biomass or SSB, surplus production, F and catch), the third dimension to the environmental projection scenarios and the fourth one to the catch scenarios. It is only returned if the arguments env_results and env are provided.
#' \item multiplicative_model: Four-dimensional array containing projections information for each catch or fishing mortality scenario and for each environmental scenario considering the KBPM multiplicative model (see \code{\link{knobi_env}} details). The rows correspond to the years, the columns to the derived quantities (biomass or SSB, surplus production, F and catch), the third dimension to the environmental projection scenarios and the fourth one to the catch scenarios. It is only returned if the arguments env_results and env are provided.
#' \item biomass: Data frame containing historical biomass and their projections for all scenarios.
#' \item catch: Data frame containing historical catch and their projection values for all scenarios.
#' \item f: Data frame containing the historical fishing mortality and their projection values for all scenarios.
#' \item SP: Data frame containing the historical surplus production and their projection values for all scenarios.}
#' The plots results are displayed in the plot window and are also saved (if plot_out="TRUE") in the  provided directory or in the same directory as \code{link{knobi_fit}}.
#' If environmental information is not provided, four plots are presented in a panel reporting the biomass, surplus production, catch and fishing mortality projections for each catch or fishing mortality scenario.
#' If environmental information is provided, the plots are presented in a panel for each catch or fishing mortality scenario reporting   biomass, surplus production, catch and fishing mortality projections in each of the environmental scenarios.
#'
#' @author
#' \itemize{
#' \item{Anxo Paz}
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' \item{M. Grazia Pennino}
#' }
#'
#' @examples
#'
#' \dontrun{
#'
#' # First, run the example of knobi_fit and knobi_env function
#'
#' # Then, create the data frame containing the selected catch for the projected years.
#' # In this illustration, catch values for all the projected years are constant and
#' # equal to the catch value in the last year multiplied by 1, 1.2 and 0.8.
#'
#' catch=rep(knobi_results$data$Catch[length(knobi_results$data$Catch)],5)
#'
#' Ct=data.frame(catch=catch,
#'               catch08=0.8*catch,
#'               catch12=1.2*catch)
#'
#'
#' # Then, create the data frame containing the environmental covariable values
#' # for the projected years.
#'
#' env=data.frame(AMO1=c(0.2,0.3,0.4,0.5,0.6),
#'                AMO2=c(0.05,0.15,0.25,0.35,0.45),
#'                AMO3=c(-0.1,-0.2,-0.3,-0.4,-0.5))
#'
#'
#' # Based on the previous objects we can apply the projection function.
#'
#' knobi_proj(knobi_results, knobi_environmental, Ct=Ct, n_y=5, env=env)
#'
#' # An example without environmental information for the projections is given below.
#' # The number of projection years is also extended.
#'
#' Ct=rbind(Ct,Ct[1:3,])
#' knobi_proj(knobi_results, Ct=Ct, end_y=2027)
#'
#'  # Alternatively, projections can be based on fishing mortality.
#'  # The scenarios presented below have been created from the estimated F_msy of
#'  # knobi_fit analysis.
#'
#' fmsy=knobi_results$fit$RP$F_MSY
#' ff=rep(fmsy,5)
#' f=data.frame(f=ff,f12=ff*1.2,f08=ff*0.8)
#'
#' knobi_proj(knobi_results, f=f, n_y=5, env_results=env_results, env=env)
#'
#' # The example without environmental information for the projections is given below.
#' # The number of projection years is also extended.
#'
#' f=rbind(f,f[1:3,])
#' knobi_proj(knobi_results, f=f, end_y=2027)
#'
#'
#' # In case of multicovar=T in knobi_env, a list is required which each
#' # item is a data frame for each environmental scenario
#'
#' env=list(climate_1=data.frame(AMO=c(0.2,0.2,0.3,0.3,0.4),
#'                               Tmax_Vigo=c(19,19,20,20,21)),
#'          climate_2=data.frame(AMO=c(0.2,0.3,0.4,0.5,0.6),
#'                               Tmax_Vigo=c(19,20,21,22,23)))
#'
#' knobi_proj(knobi_results, knobi_environmental2, Ct=Ct[1:5,], n_y=5, env=env)
#' }
#'
#' @export

knobi_proj<-function(knobi_results, env_results=NULL, Ct=NULL, f=NULL, env=NULL,
                     end_y=NULL, n_y=NULL,
                     plot_out=F, plot_filename=NULL, plot_dir=NULL){


  years=knobi_results$data$years
  lastyear=years[length(years)]

  if(plot_out==T){
    old_dir=getwd()
    if (is.null(plot_dir)) {plot_dir=knobi_results$control$plot_settings$plot_dir}
    setwd(plot_dir)
    if (is.null(plot_filename)){plot_filename=knobi_results$control$plot_settings$plot_filename}
    if (plot_filename %in% list.dirs(full.names=F)){
      setwd(paste0(plot_dir,"/",plot_filename))} else {
        dir.create(plot_filename)
        setwd(paste0(plot_dir,"/",plot_filename))}
  }

  if(is.null(n_y)==T){
    newyears=c(lastyear:end_y)
  } else {
    newyears=c(lastyear:c(lastyear+n_y))
  }

  ly=length(newyears)


  if(is.null(Ct)==T){

    if(is.null(f)==T){stop('You must provide catch or f time series')}

    if(ly-1!=(nrow(as.matrix(f)))){stop('Length of f time series is different than the length of years vector')}

  } else {

    if(is.null(f)==F){stop('You must provide only catch or f time series, not both')}

    Ct=Ct

    if(ly-1!=(nrow(as.matrix(Ct)))){stop('Length of catch time series is different than the length of years vector')}

  }


  if(knobi_results$control$method=="Biomass"){
    B0=knobi_results$data$Biomass[length(knobi_results$data$Biomass)]
    Bt_ini=knobi_results$data$Biomass
  } else {
    B0=knobi_results$data$Spawning_Biomass[length(knobi_results$data$Spawning_Biomass)]
    Bt_ini=knobi_results$data$Spawning_Biomass
  }


  params=knobi_results$fit$Parameter_estimates
  r=params[1]
  K=params[2]
  if(knobi_results$control$pella){
    p=params[3]
  } else {p=1}



  if(is.null(env_results)==F){

    additive=env_results$model_env_Additive
    if(knobi_results$control$pella){
      p_a=additive[length(additive)]
      n_p=3
    } else {
      p_a=1
      n_p=2}
    r_a=additive[1]
    K_a=additive[2]
    cs_length=length(additive)-n_p
    c_a=additive[3:c(2+cs_length)]

    multiplicative=env_results$model_env_Multiplicative
    if(knobi_results$control$pella){
      p_m=multiplicative[length(multiplicative)]
      n_p=3
    } else {
      p_m=1
      n_p=2}
    r_m=multiplicative[1]
    K_m=multiplicative[2]
    cs_length=length(multiplicative)-n_p
    c_m=multiplicative[3:c(2+cs_length)]
  }



  if(is.null(f)==T){

    model <- function(Bt1,Bt,Ct,K,r,p) {
      Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))-Ct-Bt1
    }

    model_a <- function(Bt1,Bt,Ct,K,r,p,c,Xt) {
      Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))+Xt%*%c*Bt-Ct-Bt1
    }

    model_m <- function(Bt1,Bt,Ct,K,r,p,c,Xt) {
      Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))*exp(Xt%*%c)-Ct-Bt1
    }

  } else {

    model <- function(Bt1,Bt,ef,K,r,p) {
      Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))-ef*((Bt1+Bt)/2)-Bt1
    }

    model_a <- function(Bt1,Bt,ef,K,r,p,c,Xt) {
      Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))+Xt%*%c*Bt-ef*((Bt1+Bt)/2)-Bt1
    }

    model_m <- function(Bt1,Bt,ef,K,r,p,c,Xt) {
      Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))*exp(Xt%*%c)-ef*((Bt1+Bt)/2)-Bt1
    }
  }



  n_esc=max(ncol(f),ncol(Ct),ifelse(is.vector(Ct)==T,1,NA),ifelse(is.vector(f)==T,1,NA),na.rm=T)


  if(n_esc!=1){
    if(is.null(f)==T){
      if(is.data.frame(Ct)){sc_names=names(Ct)} else {
        sc_names=NULL
        for(i in 1:n_esc){
          sc_names=c(sc_names,paste0(i,"_projection"))}
      }
      colnames(Ct)=sc_names
    } else {
      if(is.data.frame(f)){sc_names=names(f)} else {
        sc_names=NULL
        for(i in 1:n_esc){
          sc_names=c(sc_names,paste0(i,"_projection"))}
      }
      colnames(f)=sc_names
    }
  } else {
    if(is.null(f)==T){
      Ct=matrix(Ct,ncol=1)
      colnames(Ct)="Projection"
      sc_names=colnames(Ct)
    } else {
      f=matrix(f,ncol=1)
      colnames(f)="Projection"
      sc_names=colnames(f)
    }
  }


  base_Bt=array(c(B0,rep(0,ly-1)),c(ly,n_esc)); colnames(base_Bt)=sc_names
  base_SP=array(rep(0,ly-1),c(ly-1,n_esc)); colnames(base_SP)=sc_names
  base_Baver=array(rep(0,ly-1),c(ly-1,n_esc)); colnames(base_Baver)=sc_names

  for(i in 1:(ly-1)){

    for(j in sc_names){

      Bi=base_Bt[i,j]

      if(is.null(f)==T){
        Ci=Ct[i,j]
        v=stats::uniroot(model,c(0,K),Bt=Bi,Ct=Ci,K=K,r=r,p=p)
      } else {
        v=stats::uniroot(model,c(0,K),Bt=Bi,ef=f[i,j],K=K,r=r,p=p)
        Ci=f[i,j]*((Bi+v$root)/2)
      }

      Bt1=v$root

      base_Bt[c(i+1),j]=Bt1

      base_SP[i,j]=as.numeric(Bt1-Bi+Ci)
      base_Baver[i,j]=(Bt1+Bi)/2
    }
  }

  if(is.null(f)==T){
    base_Ct=Ct
    base_f=Ct/base_Baver
  } else {
    base_f=f
    base_Ct=f*base_Baver
  }



  if(is.null(env_results)==F & is.null(env)==T) {stop('Environmental data is required')}
  if(is.null(env_results)==T & is.null(env)==F) {stop('Environmental fit results are required')}


  if(is.null(env_results)==F){

    env=env

    if(is.data.frame(env)==F & is.matrix(env)==F){

      n_env_esc=length(env)

      add_Bt=list(); mult_Bt=list()
      add_SP=list(); mult_SP=list()
      add_Baver=list(); mult_Baver=list()
      add_Ct=list(); mult_Ct=list()
      add_f=list(); mult_f=list()

      for(n in 1:n_env_esc){

        Xt=env[[n]]

        if(cs_length!=(ncol(Xt))){
          stop('Number of environmental covariables is different than ist number in knobi_env')}

        for(j in ncol(Xt)){
          Xt[,j]=(Xt[,j]-attr(env_results$environmental_variables[,j],
                              "scaled:center"))/attr(env_results$environmental_variables[,j],
                                                     "scaled:scale")
        }


        add_Bt[[n]]=array(c(B0,rep(0,ly-1)),c(ly,n_esc)); colnames(add_Bt[[n]])=sc_names
        add_SP[[n]]=array(rep(0,ly-1),c(ly-1,n_esc)); colnames(add_SP[[n]])=sc_names
        add_Baver[[n]]=array(rep(0,ly-1),c(ly-1,n_esc)); colnames(add_Baver[[n]])=sc_names

        mult_Bt[[n]]=array(c(B0,rep(0,ly-1)),c(ly,n_esc)); colnames(mult_Bt[[n]])=sc_names
        mult_SP[[n]]=array(rep(0,ly-1),c(ly-1,n_esc)); colnames(mult_SP[[n]])=sc_names
        mult_Baver[[n]]=array(rep(0,ly-1),c(ly-1,n_esc)); colnames(mult_Baver[[n]])=sc_names

        for(i in 1:(ly-1)){

          for(j in sc_names){

            Bi_a=add_Bt[[n]][i,j]

            if(is.null(f)==T){
              Ci=Ct[i,j]
              v=stats::uniroot(model_a,c(0,K_a),Bt=Bi_a,Ct=Ci,K=K_a,r=r_a,p=p_a,c=c_a,Xt=as.matrix(Xt[i,]))
            } else {
              v=stats::uniroot(model_a,c(0,K_a),Bt=Bi_a,ef=f[i,j],K=K_a,r=r_a,p=p_a,c=c_a,Xt=as.matrix(Xt[i,]))
              Ci=f[i,j]*((Bi_a+v$root)/2)
            }

            Bt1=v$root

            add_Bt[[n]][c(i+1),j]=Bt1
            add_SP[[n]][i,j]=as.numeric(Bt1-Bi+Ci)
            add_Baver[[n]][i,j]=(Bt1+Bi)/2

            Bi_m=mult_Bt[[n]][i,j]

            if(is.null(f)==T){
              Ci=Ct[i,j]
              v=stats::uniroot(model_m,c(0,K_m),Bt=Bi_m,Ct=Ci,K=K_m,r=r_m,p=p_m,c=c_m,Xt=as.matrix(Xt[i,]))
            } else {
              v=stats::uniroot(model_m,c(0,K_m),Bt=Bi_m,ef=f[i,j],K=K_m,r=r_m,p=p_m,c=c_m,Xt=as.matrix(Xt[i,]))
              Ci=f[i,j]*((Bi_m+v$root)/2)
            }

            Bt1=v$root

            mult_Bt[[n]][c(i+1),j]=Bt1
            mult_SP[[n]][i,j]=as.numeric(Bt1-Bi+Ci)
            mult_Baver[[n]][i,j]=(Bt1+Bi)/2

          }
        }

        if(is.null(f)==T){

          add_Ct[[n]]=Ct
          add_f[[n]]=Ct/add_Baver[[n]]

          mult_Ct[[n]]=Ct
          mult_f[[n]]=Ct/mult_Baver[[n]]

        } else {

          add_f[[n]]=f
          add_Ct[[n]]=f*add_Baver[[n]]

          mult_f[[n]]=f
          mult_Ct[[n]]=f*mult_Baver[[n]]

        }


      }

      if(is.null(names(env))==F){
        names(add_Baver)=names(add_Bt)=names(add_Ct)=names(add_f)=names(add_SP)=
          names(mult_Baver)=names(mult_Bt)=names(mult_Ct)=names(mult_f)=names(mult_SP)=names(env)
      } else {
        for(n in 1:n_env_esc){
          names(env)[[n]]=paste0("env_scenario",n)
        }
        names(add_Baver)=names(add_Bt)=names(add_Ct)=names(add_f)=names(add_SP)=
          names(mult_Baver)=names(mult_Bt)=names(mult_Ct)=names(mult_f)=names(mult_SP)=names(env)
      }

    }  else {

      if(is.vector(env)==T){
        env=matrix(env,ncol=1)
        colnames(env)=colnames(env_results$scaled_environmental_var)
      } else {
        env=as.matrix(env)
        if(is.null(colnames(env))==T){
          for(i in 1:ncol(env)){
            colnames(env)[i]=paste0(colnames(env_results$scaled_environmental_var),"_",i)
          }
        }}

      n_env_esc=ncol(env)

      add_Bt=list(); mult_Bt=list()
      add_SP=list(); mult_SP=list()
      add_Baver=list(); mult_Baver=list()
      add_Ct=list(); mult_Ct=list()
      add_f=list(); mult_f=list()

      for(n in 1:n_env_esc){

        Xt=env[,n]

        Xt=(Xt-attr(env_results$scaled_environmental_var,
                    "scaled:center"))/attr(env_results$scaled_environmental_var,"scaled:scale")

        add_Bt[[n]]=array(c(B0,rep(0,ly-1)),c(ly,n_esc)); colnames(add_Bt[[n]])=sc_names
        add_SP[[n]]=array(rep(0,ly-1),c(ly-1,n_esc)); colnames(add_SP[[n]])=sc_names
        add_Baver[[n]]=array(rep(0,ly-1),c(ly-1,n_esc)); colnames(add_Baver[[n]])=sc_names

        mult_Bt[[n]]=array(c(B0,rep(0,ly-1)),c(ly,n_esc)); colnames(mult_Bt[[n]])=sc_names
        mult_SP[[n]]=array(rep(0,ly-1),c(ly-1,n_esc)); colnames(mult_SP[[n]])=sc_names
        mult_Baver[[n]]=array(rep(0,ly-1),c(ly-1,n_esc)); colnames(mult_Baver[[n]])=sc_names

        for(i in 1:(ly-1)){

          for(j in sc_names){

            Bi_a=add_Bt[[n]][i,j]

            if(is.null(f)==T){
              Ci=Ct[i,j]
              v=stats::uniroot(model_a,c(0,K_a),Bt=Bi_a,Ct=Ci,K=K_a,r=r_a,p=p_a,c=c_a,Xt=Xt[i])
            } else {
              v=stats::uniroot(model_a,c(0,K_a),Bt=Bi_a,ef=f[i,j],K=K_a,r=r_a,p=p_a,c=c_a,Xt=Xt[i])
              Ci=f[i,j]*((Bi_a+v$root)/2)
            }

            Bt1=v$root

            add_Bt[[n]][c(i+1),j]=Bt1
            add_SP[[n]][i,j]=as.numeric(Bt1-Bi+Ci)
            add_Baver[[n]][i,j]=(Bt1+Bi)/2

            Bi_m=mult_Bt[[n]][i,j]

            if(is.null(f)==T){
              Ci=Ct[i,j]
              v=stats::uniroot(model_m,c(0,K_m),Bt=Bi_m,Ct=Ct[i,j],K=K_m,r=r_m,p=p_m,c=c_m,Xt=Xt[i])
            } else {
              v=stats::uniroot(model_m,c(0,K_m),Bt=Bi_m,ef=f[i,j],K=K_m,r=r_m,p=p_m,c=c_m,Xt=Xt[i])
              Ci=f[i,j]*((Bi_m+v$root)/2)
            }

            Bt1=v$root

            mult_Bt[[n]][c(i+1),j]=Bt1
            mult_SP[[n]][i,j]=as.numeric(Bt1-Bi+Ci)
            mult_Baver[[n]][i,j]=(Bt1+Bi)/2

          }
        }

        if(is.null(f)==T){

          add_Ct[[n]]=Ct
          add_f[[n]]=Ct/add_Baver[[n]]

          mult_Ct[[n]]=Ct
          mult_f[[n]]=Ct/mult_Baver[[n]]

        } else {

          add_f[[n]]=f
          add_Ct[[n]]=f*add_Baver[[n]]

          mult_f[[n]]=f
          mult_Ct[[n]]=f*mult_Baver[[n]]

        }


      }

      names(add_Baver)=names(add_Bt)=names(add_Ct)=names(add_f)=names(add_SP)=
        names(mult_Baver)=names(mult_Bt)=names(mult_Ct)=names(mult_f)=names(mult_SP)=colnames(env)
    }

  }



  proj_years=newyears[-1]
  total_years=c(years,proj_years,proj_years[length(proj_years)]+1)


  base_model=array(NA,dim=c(length(total_years),4,n_esc),
                   dimnames=list(total_years,c("Biomass","Catch","F","SP"),sc_names))


  for(i in 1:n_esc){
    base_model[,1,i]=c(Bt_ini,base_Bt[-1,i])
    base_model[,2,i]=c(knobi_results$data$Catch,base_Ct[,i],NA)
    base_model[,3,i]=c(knobi_results$data$F_output,base_f[,i],NA)
    base_model[,4,i]=c(knobi_results$data$SP,base_SP[,i],NA)
  }


  if(is.null(env_results)==F){

    additive_model=array(NA,dim=c(length(total_years),4,n_env_esc,n_esc),
                         dimnames=list(total_years,c("Biomass","Catch","F","SP"),colnames(env),sc_names))

    for(i in 1:n_esc){
      for(j in 1:n_env_esc){
        additive_model[,1,j,i]=c(Bt_ini,add_Bt[[j]][-1,i])
        additive_model[,2,j,i]=c(knobi_results$data$Catch,add_Ct[[j]][,i],NA)
        additive_model[,3,j,i]=c(knobi_results$data$F_output,add_f[[j]][,i],NA)
        additive_model[,4,j,i]=c(knobi_results$data$SP,add_SP[[j]][,i],NA)
      }
    }


    multiplicative_model=array(NA,dim=c(length(total_years),4,n_env_esc,n_esc),
                               dimnames=list(total_years,c("Biomass","Catch","F","SP"),colnames(env),sc_names))

    for(i in 1:n_esc){
      for(j in 1:n_env_esc){
        multiplicative_model[,1,j,i]=c(Bt_ini,mult_Bt[[j]][-1,i])
        multiplicative_model[,2,j,i]=c(knobi_results$data$Catch,mult_Ct[[j]][,i],NA)
        multiplicative_model[,3,j,i]=c(knobi_results$data$F_output,mult_f[[j]][,i],NA)
        multiplicative_model[,4,j,i]=c(knobi_results$data$SP,mult_SP[[j]][,i],NA)
      }
    }
  }


  base_plot=data.frame(years=c(years[(length(years)-9):length(years)],proj_years[1]),
                       biomass=Bt_ini[(length(years)-9):(length(years)+1)],
                       catch=c(knobi_results$data$Catch[(length(years)-9):length(years)],NA),
                       f=c(knobi_results$data$F_output[(length(years)-9):length(years)],NA),
                       SP=c(knobi_results$data$SP[(length(years)-9):length(years)],NA))


  biomass=data.frame(years=c(years,proj_years[1]),
                     biomass=Bt_ini,
                     type=rep("data",length(Bt_ini)),
                     scenario=rep("input",length(Bt_ini)),
                     env_scenario=rep("input",length(Bt_ini)),
                     model=rep("input",length(Bt_ini)))

  for(i in 1:n_esc){
    scenario=sc_names[i]
    new_b=data.frame(years=c(proj_years[-1],proj_years[length(proj_years)]+1),
                     biomass=base_Bt[-1,i],
                     type=rep("forecast",ly-1),
                     scenario=rep(scenario,ly-1),
                     env_scenario=rep("without_env",ly-1),
                     model=rep("base_model",ly-1))
    biomass=rbind(biomass,new_b)
  }


  if(is.null(env_results)==T){

    biomass=biomass[,-c(5,6)]

    biomass_plot=biomass[(length(years)-10):nrow(biomass),-3]
    biomass_plot=rbind(biomass_plot,data.frame(years=rep(biomass$years[length(years)+1],n_esc),
                                               biomass=rep(biomass$biomass[length(years)+1],n_esc),
                                               scenario=sc_names))

    vec=min(biomass_plot$years)
    vec1=max(biomass_plot$years)
    vec2=min(biomass_plot$biomass)
    vec3=max(biomass_plot$biomass)

    biomass_plots = ggplot2::ggplot(data=biomass_plot,ggplot2::aes(x=years,y=biomass)) +
      ggplot2::xlim(vec,vec1) + ggplot2::ylim(vec2,vec3) +
      ggplot2::geom_line(lwd=1.03,
                         ggplot2::aes(color=scenario)) +
      ggplot2::theme_bw() +
      ggplot2::scale_linetype_manual(values = c(1,2,3)) +
      ggplot2::geom_vline(xintercept = biomass$years[length(years)+1], linetype = "longdash") +
      ggplot2::labs(title="Biomass projections",
                    x = "Year", y = "Biomass",
                    color="Projection Scenario:",linetype="Model") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::theme(legend.position = c(.15,.85), legend.background = ggplot2::element_rect(fill = "transparent"),
                     plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

  } else {

    for(i in 1:n_esc){
      for(j in 1:n_env_esc){
        i_name=sc_names[i]
        if(is.list(env)==T){
          j_name=names(env)[j]
        } else {
        j_name=colnames(env)[j]}
        scenario=rep(i_name,2*(ly-1))
        env_scenario=rep(j_name,2*(ly-1))
        factor=rep("forecast",2*(ly-1))
        model=c(rep("additive",ly-1),rep("multiplicative",ly-1))
        ij_biomass=c(add_Bt[[j]][-1,i],mult_Bt[[j]][-1,i])
        new_b=data.frame(years=rep(c(proj_years[-1],proj_years[length(proj_years)]+1),2),
                         biomass=ij_biomass,
                         type=factor,
                         scenario=scenario,
                         env_scenario=env_scenario,
                         model=model)
        biomass=rbind(biomass,new_b)
      }}

    biomass_plot=list()
    b_years=c(proj_years,proj_years[length(proj_years)]+1)

    for(i in sc_names){
      i_biomass=base_Bt[,i]
      model=rep("Base KBPM",length(base_Bt[,i]))
      biomass_plot[[i]]=data.frame(years=b_years,biomass=i_biomass,model=model)
    }

    for(i in sc_names){
      env_scenario=rep("Without environmental effect",length(base_Bt[,i]))
      biomass_plot[[i]]$env_scenario=env_scenario
    }


    for(i in sc_names){
      for(j in 1:n_env_esc){
        ij_biomass=c(add_Bt[[j]][,i],mult_Bt[[j]][,i])
        model=c(rep("Environmental additive",ly),rep("Environmental multiplicative",ly))
        if(is.list(env)==T){
          j_name=names(env)[j]
        } else {
          j_name=colnames(env)[j]}
        env_scenario=rep(j_name,2*(ly))
        new_b=data.frame(years=b_years,biomass=ij_biomass,model=model,env_scenario=env_scenario)
        biomass_plot[[i]]=rbind(biomass_plot[[i]],new_b)
      }}

    biomass_plots=list()

    for(i in sc_names){

      i_b_plot=biomass_plot[[i]]
      i_b_plot=rbind(data.frame(years=base_plot$years,
                                biomass=base_plot$biomass,
                                model="input",env_scenario="input"),
                     i_b_plot)


      vec=min(i_b_plot$years)
      vec1=max(i_b_plot$years)
      vec2=min(i_b_plot$biomass)
      vec3=max(i_b_plot$biomass)

      biomass_plots[[i]] = ggplot2::ggplot(data=i_b_plot,ggplot2::aes(x=years,y=biomass)) +
        ggplot2::xlim(vec,vec1) + ggplot2::ylim(vec2,vec3) +
        ggplot2::geom_path(data=subset(i_b_plot,model!="input"),lwd=1.03,
                           ggplot2::aes(linetype=model,color=env_scenario)) +
        ggplot2::theme_bw() +
        ggplot2::geom_vline(xintercept = biomass$years[length(years)+1],
                            linetype = "longdash") +
        ggplot2::geom_line(data=subset(i_b_plot,model=="input"),lwd=1.03,
                           ggplot2::aes(x=years,y=biomass)) +
        ggplot2::scale_linetype_manual(values = c(1,2,3)) +
        ggplot2::labs(title="Biomass projections",
                      subtitle=i,x = "Year", y = "Biomass",
                      color="Environmental scenario",linetype="Model") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::guides(linetype=ggplot2::guide_legend(keywidth = 2.1, keyheight = 1))
    }
  }


  Ct_ini=knobi_results$data$Catch

  catch=data.frame(years=c(years),
                   catch=Ct_ini,
                   type=rep("data",length(Ct_ini)),
                   scenario=rep("input",length(Ct_ini)),
                   env_scenario=rep("input",length(Ct_ini)),
                   model=rep("input",length(Ct_ini)))

  for(i in 1:n_esc){
    scenario=sc_names[i]
    new_c=data.frame(years=proj_years,
                     catch=base_Ct[,i],
                     type=rep("forecast",ly-1),
                     scenario=rep(scenario,ly-1),
                     env_scenario=rep("without_env",ly-1),
                     model=rep("base_model",ly-1))
    catch=rbind(catch,new_c)
  }


  if(is.null(env_results)==T){

    catch=catch[,-c(5,6)]

    catch_plot=catch[(length(years)-10):nrow(catch),-3]
    catch_plot=rbind(catch_plot,data.frame(years=rep(catch$years[length(years)],n_esc),
                                               catch=rep(catch$catch[length(years)],n_esc),
                                               scenario=sc_names))

    vec=min(catch_plot$years)
    vec1=max(catch_plot$years)
    vec2=min(catch_plot$catch)
    vec3=max(catch_plot$catch)

    catch_plots = ggplot2::ggplot(data=catch_plot,ggplot2::aes(x=years,y=catch)) +
      ggplot2::xlim(vec,vec1) + ggplot2::ylim(vec2,vec3) +
      ggplot2::geom_line(lwd=1.03,
                         ggplot2::aes(color=scenario)) +
      ggplot2::theme_bw() +
      ggplot2::scale_linetype_manual(values = c(1,2,3)) +
      ggplot2::geom_vline(xintercept = catch$years[length(years)], linetype = "longdash") +
      ggplot2::labs(title="Catch projections",
                    x = "Year", y = "Catch",
                    color= "Projection Scenario:", linetype= "Model") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::theme(legend.position = c(.15,.85), legend.background = ggplot2::element_rect(fill = "transparent"),
                     plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

  } else {

    for(i in 1:n_esc){
      for(j in 1:n_env_esc){
        i_name=sc_names[i]
        if(is.list(env)==T){
          j_name=names(env)[j]
        } else {
          j_name=colnames(env)[j]}
        scenario=rep(i_name,2*(ly-1))
        env_scenario=rep(j_name,2*(ly-1))
        factor=rep("forecast",2*(ly-1))
        model=c(rep("additive",ly-1),rep("multiplicative",ly-1))
        ij_c=c(add_Ct[[j]][,i],mult_Ct[[j]][,i])
        new_c=data.frame(years=rep(proj_years,2),
                         catch=ij_c,
                         type=factor,
                         scenario=scenario,
                         env_scenario=env_scenario,
                         model=model)
        catch=rbind(catch,new_c)
      }}


    catch_plot=list()
    c_years=c(proj_years[1]-1,proj_years)

    for(i in sc_names){
      i_c=c(Ct_ini[length(Ct_ini)],base_Ct[,i])
      model=rep("Base KBPM",length(base_Ct[,i])+1)
      catch_plot[[i]]=data.frame(years=c_years,catch=i_c,model=model)
    }

    for(i in sc_names){
      env_scenario=rep("Without environmental effect",length(base_Ct[,i])+1)
      catch_plot[[i]]$env_scenario=env_scenario
    }


    for(i in sc_names){
      for(j in 1:n_env_esc){
        ij_c=c(Ct_ini[length(Ct_ini)],add_Ct[[j]][,i],Ct_ini[length(Ct_ini)],mult_Ct[[j]][,i])
        model=c(rep("Environmental additive",ly),rep("Environmental multiplicative",ly))
        if(is.list(env)==T){
          j_name=names(env)[j]
        } else {
          j_name=colnames(env)[j]}
        env_scenario=rep(j_name,2*(ly))
        new_c=data.frame(years=c_years,catch=ij_c,model=model,env_scenario=env_scenario)
        catch_plot[[i]]=rbind(catch_plot[[i]],new_c)
      }}



    catch_plots=list()

    for(i in sc_names){

      i_c_plot=catch_plot[[i]]
      i_c_plot=rbind(data.frame(years=base_plot$years[-length(base_plot$years)],
                                catch=base_plot$catch[-length(base_plot$years)],
                                model="input",env_scenario="input"),
                     i_c_plot)


      vec=min(i_c_plot$years)
      vec1=max(i_c_plot$years)
      vec2=min(i_c_plot$catch)
      vec3=max(i_c_plot$catch)

      catch_plots[[i]] = ggplot2::ggplot(data=i_c_plot,ggplot2::aes(x=years,y=catch)) +
        ggplot2::xlim(vec,vec1) + ggplot2::ylim(vec2,vec3) +
        ggplot2::geom_path(data=subset(i_c_plot,model!="input"),lwd=1.03,
                           ggplot2::aes(linetype=model,color=env_scenario)) +
        ggplot2::theme_bw() +
        ggplot2::geom_line(data=subset(i_c_plot,model=="input"),lwd=1.03,
                           ggplot2::aes(x=years,y=catch)) +
        ggplot2::geom_vline(xintercept = catch$years[length(years)], linetype = "longdash") +
        ggplot2::scale_linetype_manual(values = c(1,2,3)) +
        ggplot2::labs(title="Catch projections",
                      subtitle=i, x = "Year", y = "Catch",
                      color= "Environmental scenario", linetype= "Model") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::guides(linetype=ggplot2::guide_legend(keywidth = 2.1, keyheight = 1))
    }
  }



  f_ini=knobi_results$data$F_output

  f=data.frame(years=c(years),
               f=f_ini,
               type=rep("data",length(f_ini)),
               scenario=rep("input",length(f_ini)),
               env_scenario=rep("input",length(f_ini)),
               model=rep("input",length(f_ini)))

  for(i in 1:n_esc){
    scenario=sc_names[i]
    new_f=data.frame(years=proj_years,
                     f=base_f[,i],
                     type=rep("forecast",ly-1),
                     scenario=rep(scenario,ly-1),
                     env_scenario=rep("without_env",ly-1),
                     model=rep("base_model",ly-1))
    f=rbind(f,new_f)
  }


  if(is.null(env_results)==T){

    f=f[,-c(5,6)]

    f_plot=f[(length(years)-10):nrow(f),-3]
    f_plot=rbind(f_plot,data.frame(years=rep(f$years[length(years)],n_esc),
                                   f=rep(f$f[length(years)],n_esc),
                                   scenario=sc_names))

    vec=min(f_plot$years)
    vec1=max(f_plot$years)
    vec2=min(f_plot$f)
    vec3=max(f_plot$f)

    f_plots = ggplot2::ggplot(data=f_plot,ggplot2::aes(x=years,y=f)) +
      ggplot2::xlim(vec,vec1) + ggplot2::ylim(vec2,vec3) +
      ggplot2::geom_line(lwd=1.03,
                         ggplot2::aes(color=scenario)) +
      ggplot2::theme_bw() +
      ggplot2::scale_linetype_manual(values = c(1,2,3)) +
      ggplot2::geom_vline(xintercept = f$years[length(years)], linetype = "longdash") +
      ggplot2::labs(title="Fishing mortality projections",
                    x = "Year", y = "Fishing mortality",
                    color= "Projection Scenario:", linetype= "Model") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::theme(legend.position = c(.15,.85), legend.background = ggplot2::element_rect(fill = "transparent"),
                     plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

  } else {

    for(i in 1:n_esc){
      for(j in 1:n_env_esc){
        i_name=sc_names[i]
        if(is.list(env)==T){
          j_name=names(env)[j]
        } else {
          j_name=colnames(env)[j]}
        scenario=rep(i_name,2*(ly-1))
        env_scenario=rep(j_name,2*(ly-1))
        factor=rep("forecast",2*(ly-1))
        model=c(rep("additive",ly-1),rep("multiplicative",ly-1))
        ij_f=c(add_f[[j]][,i],mult_f[[j]][,i])
        new_f=data.frame(years=rep(proj_years,2),
                         f=ij_f,
                         type=factor,
                         scenario=scenario,
                         env_scenario=env_scenario,
                         model=model)
        f=rbind(f,new_f)
      }}


    f_plot=list()
    f_years=c_years

    for(i in sc_names){
      i_f=c(f_ini[length(f_ini)],base_f[,i])
      model=rep("Base KBPM",length(base_f[,i])+1)
      f_plot[[i]]=data.frame(years=f_years,f=i_f,model=model)
    }

    for(i in sc_names){
      env_scenario=rep("Without environmental effect",length(base_f[,i])+1)
      f_plot[[i]]$env_scenario=env_scenario
    }


    for(i in sc_names){
      for(j in 1:n_env_esc){
        ij_f=c(f_ini[length(f_ini)],add_f[[j]][,i],f_ini[length(f_ini)],mult_f[[j]][,i])
        model=c(rep("Environmental additive",ly),rep("Environmental multiplicative",ly))
        if(is.list(env)==T){
          j_name=names(env)[j]
        } else {
          j_name=colnames(env)[j]}
        env_scenario=rep(j_name,2*(ly))
        new_f=data.frame(years=f_years,f=ij_f,model=model,env_scenario=env_scenario)
        f_plot[[i]]=rbind(f_plot[[i]],new_f)
      }}


    f_plots=list()

    for(i in sc_names){

      i_f_plot=f_plot[[i]]
      i_f_plot=rbind(data.frame(years=base_plot$years[-length(base_plot$years)],
                                f=base_plot$f[-length(base_plot$years)],
                                model="input",env_scenario="input"),
                     i_f_plot)


      vec=min(i_f_plot$years)
      vec1=max(i_f_plot$years)
      vec2=min(i_f_plot$f)
      vec3=max(i_f_plot$f)

      f_plots[[i]] = ggplot2::ggplot(data=i_f_plot,ggplot2::aes(x=years,y=f)) +
        ggplot2::xlim(vec,vec1) + ggplot2::ylim(vec2,vec3) +
        ggplot2::geom_path(data=subset(i_f_plot,model!="input"),lwd=1.03,
                           ggplot2::aes(linetype=model,color=env_scenario)) +
        ggplot2::theme_bw() +
        ggplot2::geom_vline(xintercept = f$years[length(years)], linetype = "longdash") +
        ggplot2::geom_line(data=subset(i_f_plot,model=="input"),lwd=1.03,
                           ggplot2::aes(x=years,y=f)) +
        ggplot2::scale_linetype_manual(values = c(1,2,3)) +
        ggplot2::labs(title="Fishing mortality projections",
                      subtitle=i, x = "Year", y = "Fishing mortality",
                      color= "Environmental scenario", linetype= "Model") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::guides(linetype=ggplot2::guide_legend(keywidth = 2.1, keyheight = 1))
    }
  }


  SP_ini=knobi_results$data$SP

  SP=data.frame(years=c(years),
                SP=SP_ini,
                type=rep("data",length(SP_ini)),
                scenario=rep("input",length(SP_ini)),
                env_scenario=rep("input",length(SP_ini)),
                model=rep("input",length(SP_ini)))

  for(i in 1:n_esc){
    scenario=sc_names[i]
    new_sp=data.frame(years=proj_years,
                      SP=base_SP[,i],
                      type=rep("forecast",ly-1),
                      scenario=rep(scenario,ly-1),
                      env_scenario=rep("without_env",ly-1),
                      model=rep("base_model",ly-1))
    SP=rbind(SP,new_sp)
  }


  if(is.null(env_results)==T){

    SP=SP[,-c(5,6)]

    SP_plot=SP[(length(years)-10):nrow(SP),-3]
    SP_plot=rbind(SP_plot,data.frame(years=rep(SP$years[length(years)],n_esc),
                                     SP=rep(SP$SP[length(years)],n_esc),
                                     scenario=sc_names))

    vec=min(SP_plot$years)
    vec1=max(SP_plot$years)
    vec2=min(SP_plot$SP)
    vec3=max(SP_plot$SP)

    SP_plots = ggplot2::ggplot(data=SP_plot,ggplot2::aes(x=years,y=SP)) +
      ggplot2::xlim(vec,vec1) + ggplot2::ylim(vec2,vec3) +
      ggplot2::geom_line(lwd=1.03,
                         ggplot2::aes(color=scenario)) +
      ggplot2::theme_bw() +
      ggplot2::scale_linetype_manual(values = c(1,2,3)) +
      ggplot2::geom_vline(xintercept = SP$years[length(years)], linetype = "longdash") +
      ggplot2::labs(title="Surplus Production projections",
                    x ="Year", y = "Surplus Production",
                    color="Projection Scenario:",linetype="Model") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::theme(legend.position = c(.15,.85), legend.background = ggplot2::element_rect(fill = "transparent"),
                     plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

  } else {

    for(i in 1:n_esc){
      for(j in 1:n_env_esc){
        i_name=sc_names[i]
        if(is.list(env)==T){
          j_name=names(env)[j]
        } else {
          j_name=colnames(env)[j]}
        scenario=rep(i_name,2*(ly-1))
        env_scenario=rep(j_name,2*(ly-1))
        factor=rep("forecast",2*(ly-1))
        model=c(rep("additive",ly-1),rep("multiplicative",ly-1))
        ij_sp=c(add_SP[[j]][,i],mult_SP[[j]][,i])
        new_sp=data.frame(years=rep(proj_years,2),
                          SP=ij_sp,
                          type=factor,
                          scenario=scenario,
                          env_scenario=env_scenario,
                          model=model)
        SP=rbind(SP,new_sp)
      }}


    SP_plot=list()
    sp_years=c_years

    for(i in sc_names){
      i_sp=c(SP_ini[length(SP_ini)],base_SP[,i])
      model=rep("Base KBPM",length(base_SP[,i])+1)
      SP_plot[[i]]=data.frame(years=sp_years,SP=i_sp,model=model)
    }

    for(i in sc_names){
      env_scenario=rep("Without environmental effect",length(base_SP[,i])+1)
      SP_plot[[i]]$env_scenario=env_scenario
    }

    for(i in sc_names){
      for(j in 1:n_env_esc){
        ij_sp=c(SP_ini[length(SP_ini)],add_SP[[j]][,i],SP_ini[length(SP_ini)],mult_SP[[j]][,i])
        model=c(rep("Environmental additive",ly),rep("Environmental multiplicative",ly))
        if(is.list(env)==T){
          j_name=names(env)[j]
        } else {
          j_name=colnames(env)[j]}
        env_scenario=rep(j_name,2*(ly))
        new_sp=data.frame(years=sp_years,SP=ij_sp,model=model,env_scenario=env_scenario)
        SP_plot[[i]]=rbind(SP_plot[[i]],new_sp)
      }}


    SP_plots=list()

    for(i in sc_names){

      i_sp_plot=SP_plot[[i]]
      i_sp_plot=rbind(data.frame(years=base_plot$years[-length(base_plot$years)],
                                 SP=base_plot$SP[-length(base_plot$years)],
                                 model="input",env_scenario="input"),
                      i_sp_plot)


      vec=min(i_sp_plot$years)
      vec1=max(i_sp_plot$years)
      vec2=min(i_sp_plot$SP)
      vec3=max(i_sp_plot$SP)

      SP_plots[[i]] = ggplot2::ggplot(data=i_sp_plot,ggplot2::aes(x=years,y=SP)) +
        ggplot2::xlim(vec,vec1) + ggplot2::ylim(vec2,vec3) +
        ggplot2::geom_path(data=subset(i_sp_plot,model!="input"),lwd=1.03,
                           ggplot2::aes(linetype=model,color=env_scenario)) +
        ggplot2::theme_bw() +
        ggplot2::geom_vline(xintercept = SP$years[length(years)], linetype = "longdash") +
        ggplot2::geom_line(data=subset(i_sp_plot,model=="input"),lwd=1.03,
                           ggplot2::aes(x=years,y=SP)) +
        ggplot2::scale_linetype_manual(values = c(1,2,3)) +
        ggplot2::labs(title="Surplus Production projections",
                      subtitle=i, x = "Year", y = "Surplus Production",
                      color= "Environmental scenario", linetype= "Model") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::guides(linetype=ggplot2::guide_legend(keywidth = 2.1, keyheight = 1))
    }
  }


if(is.null(env_results)==T){

  forecast=list(base_model=base_model,biomass=biomass,catch=catch,f=f,SP=SP)

  if(plot_out==T){
    grDevices::jpeg("projections.jpeg",width=2500, height=2000,res=300)
    ggpubr::ggarrange(biomass_plots, SP_plots,catch_plots,f_plots, nrow = 2, ncol=2,
                      common.legend = TRUE, legend="bottom")
    grDevices::dev.off()
    cat(paste0("\n Plot successfully saved in '",getwd(),"'"),". \n")
    setwd(old_dir)
  }

  f1=ggpubr::ggarrange(biomass_plots, SP_plots,catch_plots, f_plots, nrow = 2, ncol=2,
                       common.legend = TRUE, legend="bottom")

  print(f1)

} else {

  forecast=list(base_model=base_model,additive_model=additive_model,
                multiplicative_model=multiplicative_model,
                biomass=biomass,catch=catch,f=f,SP=SP)

  for(i in sc_names){

    if(plot_out==T){
      grDevices::jpeg(paste0("projections_",i,".jpeg"),width=2500, height=2000,res=300)
      ggpubr::ggarrange(biomass_plots[[i]], SP_plots[[i]],catch_plots[[i]],f_plots[[i]],
                        nrow = 2, ncol=2, common.legend = TRUE, legend="bottom")
      grDevices::dev.off()
    }

    f1=ggpubr::ggarrange(biomass_plots[[i]], SP_plots[[i]],catch_plots[[i]],f_plots[[i]],
                         nrow = 2, ncol=2, common.legend = TRUE, legend="right")

    print(f1)
  }
}

class(forecast)="knobi"

if(plot_out==T){
cat(paste0("\n Plots successfully saved in '",getwd(),"'"),". \n")
setwd(old_dir)
}

return(forecast)

}
