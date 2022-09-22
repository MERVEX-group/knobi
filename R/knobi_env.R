
#' @title KBPM environmental analysis
#'
#' @description Analyse and model the relationships between surplus production and environmental covariables to test whether productivity changes in response to environmental fluctuations.  Environmental effects are included as additive and multiplicative effects in the general KBPM formulation (see details).
#'
#' @param knobi_results The output object of \code{\link{knobi_fit}} function (main function).
#' @param environmental A list containing the following data and settings: \itemize{
#' \item data: data frame containing the values of each one of the environmental variable(s) in one column. Each row represents a year.
#' \item years: time series of years corresponding to the environmental variable(s).
#' \item lag: optional numerical vector providing the used lag value(s) in the relation among the base KBPM surplus production (SP) residuals and the environmental variable(s). This means that the residuals(SP_{t}) is related to X_{t-lag} being X the environmental variable. The length of this argument must be equal to the number of environmental variables included.
#' \item nlag: if lag value is not provided, this argument is used to test all the lags smaller or equal to nlag (numerical vector) through cor.test function. The lag corresponding to the highest pearson correlation among the base KBPM SP residuals and the lagged environmental covariable values is considered in the environmental model.
#' \item start_c: optional numerical vector providing the start values of the environmental c parameter for the optimization of the additive and multiplicative models, respectively. By default, start_c=c(1,1). See details.
#' \item selected_var: optional character. By default, the fit is done using the environmental values according to the lag derived from the previous arguments. However, if this argument is equal to the name of the environmental variable no lag is applied to its values.
#' \item multicovar: optional logical. TRUE if you want to fit the environmental model including all the input environmental covariables, up to a maximum of 5. By default this argument is FALSE, which means that only the environmental covariable reporting the highest pearson correlation is included (after lagging it if corresponds).}
#' @param plot_out Logical. TRUE means that a file with the  environmental fit plots is created. By default this argument is FALSE.
#' @param plot_dir Optional directory for creating the folder and save the plots. Required when plot_out=TRUE. The default value is the input of this argument in the knobi_fit function.
#' @param plot_filename Optional name of the folder that will contain the plots. Required when plot_out=TRUE. The default value is the input of this argument in the knobi_fit function.
#'
#' @details
#' Additive environmental model adds the following term on the right hand of equation (1) or (2) detailed in \code{\link{knobi_fit}} function: \eqn{cX_{t}B_{t}}, being \eqn{X_{t}}the environmental variable and \eqn{B_{t}} the biomass or SSB at time \eqn{t}.
#' Multiplicative environmental model multiplies the right hand of equation (1) or (2) by \eqn{exp(cX_{t})}.
#' The subscript t denotes the time (years).
#'
#' @return A list containing the environmental analysis is provided. \itemize{
#' \item selected_lag: Estimated lag corresponding to the one reporting the highest correlation between the environmental variable and the surplus production. Derived if lag is not fixed.
#' \item fixed_lag: Input value of 'lag' argument.
#' \item lag_cor: Correlation between the environmental variable(s) value and the base KBPM SP residuals (after lagging the environmental one if corresponds).
#' \item selected_var: Environmental variable used in the fit, chosen by the user or the one derived from the highest pearson correlation procedure. In case that argument 'multicovar' is omitted, 'NULL' or equal to 'FALSE'.
#' \item model_env_Multiplicative: Estimates of the multiplicative model parameters.
#' \item model_env_Additive: Estimates of the additive model parameters.
#' \item ref_pts: Reference points (RPs) estimates for each model assuming X_t=0 (see vignettes).
#' \item scaled_environmental_var: Standardized variable used in the fit, with the 'scale' and 'center' attributes.
#' \item environmental_variables: Standardized covariables used in the fit (if 'multicovar=TRUE'), with the 'scale' and 'center' attributes.
#' \item plots3D: List with the 3D plots objects.
#' \item error: List of performance and accuracy: \itemize{
#' \item residuals: Pearson's residuals from the fit calculated as (observations-estimates)/sd(observations) for each model (base KBPM, additive model and multiplicative model).
#' \item error_table: Array of performance and accuracy (observed vs. estimated) measures for each model: Standard error of the regression (SER), coefficient of determination (R-squared), adjusted coefficient of determination (adj-R-squared), Akaike information criterion (AIC), root-mean-squared error (RMSE), mean absolute percentage error (MAPE) and the value of the F statistic corresponding to the comparison of each environmental model respect to the base model (F-value) and its corresponding p-value (Pr(>F)).}}
#' Result plots are shown in the plot window and also saved (if plot_out="TRUE") on the provided directory or in the same directory as knobi_fit.
#' The first plot reports the correlation analysis between the environmental variable(s) and the KBPM data. The second one reports the fitted values of the base model (no environmental information) and of the environmental ones.
#' If multicovar=FALSE, 3D plots reporting the surplus production curve conditioned to a grid of environmental values are also reported.
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
#' # First, run the example of knobi_fit function
#'
#' # Then, provide environmental data series
#'
#' Env=data.frame(years=seq(1973,2020),AMO=c(-0.236,-0.441,-0.32,-0.385,-0.21,-0.201,-0.132,
#' -0.041,-0.098,-0.235,-0.093,-0.23,-0.29,-0.297,0.044,-0.028,-0.106,-0.061,-0.155,-0.242,
#' -0.234,-0.2,0.112,-0.082,0.028,0.349,0.094,0.004,0.095,0.041,0.207,0.182,0.268,0.242,
#' 0.123,0.114,0.015,0.325,0.078,0.189,0.142,0.077,0.09,0.318,0.291,0.045,0.15,0.279),
#' TMax_Vigo=c(16.5,16.7,17.1,16.3,16.4,16.7,17.3,17.3,17.4,18,17.5,17.5,17.6,18,17,18.4,
#' 17.8,19.6,19.1,18,17.9,17.7,17.7,19.4,18.2,19.7,19.2,18.6,18,18.3,18.5,18.8,18.7,18.7,
#' 18.5,17.9,17.4,19.2,19,19.7,18,19.1,19.4,20,19.5,20.2,18.8,18.6))
#'
#' # The environmental data series must start in the first year of the KBPM fit data
#' # minus the provided nlag or lag
#' years=knobi_results$data$years # See knobi_fit example to obtain the knobi_results object
#' ind=which(Env[,1]==years[1])
#' ind1=which(Env[,1]==years[length(years)])
#' nlag=5
#' Env=Env[(ind-nlag):ind1,]
#'
#' # Now we create the environmental list
#' environmental=list()
#' environmental$data=data.frame(AMO=Env$AMO,Tmax=Env$TMax_Vigo)
#' environmental$years=Env$years
#' environmental$nlag=c(nlag,nlag)
#'
#' knobi_environmental<-knobi_env(knobi_results,environmental)
#' knobi_environmental
#' knobi_environmental$plots3D$additive_plot
#'
#' environmental$multicovar=T
#' knobi_env(knobi_results,environmental,plot_out=T)
#' }
#'
#' @export



knobi_env<-function(knobi_results,environmental,plot_out=F,plot_filename=NULL,plot_dir=NULL){

  env0=as.data.frame(environmental$data)
  env_names=names(env0)
  y_env=environmental$years
  if(is.null(environmental$start_c)){start_c=c(1,1)} else {
    start_c=environmental$start_c}

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

  df=data.frame(KBPM_residuals=knobi_results$fit$error$residuals,
                x = knobi_results$data$Average_Biomass,
                y = knobi_results$data$SP, Year=knobi_results$data$years)


  f_year=knobi_results$data$years[1]
  l_year=knobi_results$data$years[length(knobi_results$data$years)]

  env=list()
  res_env=list()
  df_env=df

  if (is.null(environmental$lag)){

    nlag=environmental$nlag
    data_env=list()

    for(j in env_names){
      data_env[[j]]=df
      if(!is.na(knobi_results$data$Recruitment[1])){
        data_env[[j]]=data.frame(data_env[[j]],knobi_results$data$Recruitment)}

      ind=which(y_env==f_year)
      ind1=which(y_env==l_year)
      if(length(ind)>0){data_env[[j]]$env0=env0[ind:ind1,j]} else {
        warning('The length of the environmental variable is not enough to use the number of input lags')}

      nlag_ind=which(env_names==j)
      nlag_j=nlag[nlag_ind]

      for (i in 1:nlag_j){
        ind=which(y_env==(f_year-i))
        ind1=which(y_env==(l_year-i))
        if(!is.na(knobi_results$data$Recruitment[1])){
          if(length(ind)>0){data_env[[j]][,6+i]=env0[[j]][ind:ind1]} else {
            warning('The length of the environmental variable is not enough to use the number of input lags')}
        } else {
          if(length(ind)>0){data_env[[j]][,5+i]=env0[[j]][ind:ind1]} else {
            warning('The length of the environmental variable is not enough to use the number of input lags')}
        }
      }

      vec_env=1:nlag_j
      for (i in 1:nlag_j){
        vec_env[i]=paste0(j,"_lag",i)
      }

      if(knobi_results$control$method=="Biomass"){
        if(!is.na(knobi_results$data$Recruitment[1])){
          colnames(data_env[[j]])=c("KBPM_residuals","SP","B","years","R",paste0(j,"_lag0"),vec_env)
        } else {colnames(data_env[[j]])=c("KBPM_residuals","SP","B","years",paste0(j,"_lag0"),vec_env)}
      } else {
        if(!is.na(knobi_results$data$Recruitment[1])){
          colnames(data_env[[j]])=c("KBPM_residuals","SP","SSB","years","R",paste0(j,"_lag0"),vec_env)
        } else {colnames(data_env[[j]])=c("KBPM_residuals","SP","SSB","years",paste0(j,"_lag0"),vec_env)}
      }

      env[[j]]<-data_env[[j]][,c("KBPM_residuals",paste0(j,"_lag0"),vec_env)]
      p.mat <- cor.mtest(env[[j]])
      corrplot::corrplot(round(cor(as.matrix(env[[j]]),use="na.or.complete"),2),title=j,mar=c(0,0,2,0),method="number",type="lower",diag=T,p.mat = p.mat, sig.level = 0.05)
      if(plot_out==T){
        plotname=paste0("corplot_",j,".jpeg")
        grDevices::jpeg(plotname,width=2500, height=2500,res=300)
        corrplot::corrplot(round(cor(as.matrix(env[[j]]),use="na.or.complete"),2),title=j,mar=c(0,0,2,0),method="number",type="lower",diag=T,p.mat = p.mat, sig.level = 0.05)
        grDevices::dev.off()
      }

      cor=round(cor(as.matrix(env[[j]]),use="na.or.complete"),4)
      cor=cor[,1]; cor=cor[-1]

      res_env$selected_lag[[j]]=which(max(abs(cor))==abs(cor))[1]-1
      res_env$lag_cor[[j]]=cor[[which(max(abs(cor))==abs(cor))[1]]]
      df_env[,j]=scale(env[[j]][,res_env$selected_lag[[j]]+2])

    }
  } else {

    res_env=list()
    df_env=df

    for(j in env_names){
      lag_ind=which(env_names==j)
      lag_j=environmental$lag[lag_ind]

      ind=which(y_env==(f_year-lag_j))
      ind1=which(y_env==(l_year-lag_j))
      if(length(ind)>0){df_env[,j]=scale(env0[[j]][ind:ind1])} else{
        warning('The length of the environmental variable is not enough to use the number of input lags')}
      res_env$fixed_lag[[j]]=lag_j
      cor=round(cor(as.matrix(df_env[,c("KBPM_residuals",j)]),use="na.or.complete"),4)
      cor=cor[1,2]
      res_env$lag_cor[[j]]=cor}
  }

  if(is.null(environmental$multicovar)){
    multicovar=F
  } else {
    multicovar=environmental$multicovar
  }

  if(multicovar==F){

    if(is.null(environmental$selected_var)){
      select<-which.max(abs(as.data.frame(res_env$lag_cor))[,1])
      selected_var=names(res_env$lag_cor[select])
    } else {
      selected_var=environmental$selected_var
    }

    res_env$selected_var=selected_var

    Data=list(env=df_env[,selected_var],data=df, start_r=as.numeric(knobi_results$fit$Parameter_estimates[[1]]),
              start_K=as.numeric(knobi_results$fit$Parameter_estimates[[2]]),
              start_c=start_c[1], start_p=1)

    if (knobi_results$control$pella){

      # Pella Mult

      class(Data)="Pella_Mult"
      model_Pella_Mult <-fitting(Data)
      res_env$model_env_Multiplicative = model_Pella_Mult
      class(model_Pella_Mult) = class(Data)
      ref_pts_mult = RF(model_Pella_Mult)

      # Pella Additive

      class(Data)="Pella_Add"
      model_Pella_Add <-fitting(Data)
      res_env$model_env_Additive = model_Pella_Add
      class(model_Pella_Add) = class(Data)
      ref_pts_add = RF(model_Pella_Add)

      res_env$ref_pts = list(ref_pts_mult=ref_pts_mult,ref_pts_add=ref_pts_add)

      model_Pella_Mult=list(model_Pella_Mult)
      model_Pella_Mult$data$SP=df$y
      if(knobi_results$control$method=="Biomass"){model_Pella_Mult$data$B=df$x} else {
        model_Pella_Mult$data$SSB=df$x}
      model_Pella_Mult$data$env=as.numeric(df_env[,selected_var])
      class(model_Pella_Mult)="Pella_Mult"

      model_Pella_Add=list(model_Pella_Add)
      model_Pella_Add$data$SP=df$y
      if(knobi_results$control$method=="Biomass"){model_Pella_Add$data$B=df$x} else {
        model_Pella_Add$data$SSB=df$x}
      model_Pella_Add$data$env=as.numeric(df_env[,selected_var])
      class(model_Pella_Add)="Pella_Add"

      class(knobi_results$fit)="Pella_Year"
      bv <- predict_model(knobi_results$fit)
      bv1 <- predict_model(model_Pella_Mult)
      bv2 <- predict_model(model_Pella_Add)

    } else {

      # Schaefer Mult

      class(Data)="Schaefer_Mult"
      model_Schaefer_Mult <-fitting(Data)
      res_env$model_env_Multiplicative = model_Schaefer_Mult
      class(model_Schaefer_Mult) = class(Data)
      ref_pts_mult = RF(model_Schaefer_Mult)

      # Schaefer Additive

      class(Data)="Schaefer_Add"
      model_Schaefer_Add <-fitting(Data)
      res_env$model_env_Additive = model_Schaefer_Add
      class(model_Schaefer_Add) = class(Data)
      ref_pts_add = RF(model_Schaefer_Add)

      res_env$ref_pts = list(ref_pts_mult=ref_pts_mult,ref_pts_add=ref_pts_add)

      model_Schaefer_Mult=list(model_Schaefer_Mult)
      model_Schaefer_Mult$data$SP=df$y
      if(knobi_results$control$method=="Biomass"){ model_Schaefer_Mult$data$B=df$x} else {
        model_Schaefer_Mult$data$SSB=df$x}
      model_Schaefer_Mult$data$env=as.numeric(df_env[,selected_var])
      class(model_Schaefer_Mult)="Schaefer_Mult"


      model_Schaefer_Add=list(model_Schaefer_Add)
      model_Schaefer_Add$data$SP=df$y
      if(knobi_results$control$method=="Biomass"){model_Schaefer_Add$data$B=df$x} else {
        model_Schaefer_Add$data$SSB=df$x}
      model_Schaefer_Add$data$env=as.numeric(df_env[,selected_var])
      class(model_Schaefer_Add)="Schaefer_Add"

      class(knobi_results$fit)="Schaefer_Year"
      bv <- predict_model(knobi_results$fit)
      bv1 <- predict_model(model_Schaefer_Mult)
      bv2 <- predict_model(model_Schaefer_Add)
    }

    res_env$scaled_environmental_var=df_env[,selected_var]
    colnames(res_env$scaled_environmental_var)=selected_var
    rownames(res_env$scaled_environmental_var)=df_env$Year-as.numeric(res_env$selected_lag[selected_var])

    x <- knobi_results$data$Average_Biomass
    y <- df_env[,selected_var]
    z <- knobi_results$data$SP

    r_a<-res_env$model_env_Additive[1]
    K_a<-res_env$model_env_Additive[2]
    c_a<-res_env$model_env_Additive[3]
    if(knobi_results$control$pella==T){
      p_a<-res_env$model_env_Additive[4]
    } else {p_a=1}

    r_m<-res_env$model_env_Multiplicative[1]
    K_m<-res_env$model_env_Multiplicative[2]
    c_m<-res_env$model_env_Multiplicative[3]
    if(knobi_results$control$pella==T){
      p_m<-res_env$model_env_Multiplicative[4]
    } else {p_m=1}


    grid.lines = 400
    cut_a=max(K_a+c_a*K_a*max(y),K_a+c_a*K_a*min(y))
    x.pred_a <- seq(0, cut_a, length.out = grid.lines)
    x.pred_m <- seq(0, K_m, length.out = grid.lines)

    y.pred <- c(seq(min(y),0,length.out = grid.lines/2),seq(0, max(y), length.out = grid.lines/2))
    xy_a <- expand.grid(x = x.pred_a, y = y.pred)
    xy_m <- expand.grid(x = x.pred_m, y = y.pred)

    z.pred_a <- (r_a/p_a)*xy_a$x*(1-(xy_a$x/K_a)^(p_a))+c_a*xy_a$y*xy_a$x
    z.pred_m <- exp(1)^{c_m*xy_m$y}*((r_m/p_m)*xy_m$x*(1-(xy_m$x/K_m)^p_m))
    z.pred_a <- matrix(z.pred_a,nrow=grid.lines,ncol=grid.lines)
    z.pred_m <- matrix(z.pred_m,nrow=grid.lines,ncol=grid.lines)

    y2 = y * attr(y, 'scaled:scale') + attr(y, 'scaled:center')
    y.pred2 = y.pred * attr(y, 'scaled:scale') + attr(y, 'scaled:center')

    cat("\n This can take a while... \n")

    plot3D::scatter3D(x, y2, z, bty="b2", pch = 19, cex = 1, cex.axis = 0.6,
                      colkey = list(length = 0.5, width = 0.5, cex.clab = 0.8, cex.axis = 0.8),
                      xlim = c(0,max(x.pred_a)), zlim = c(0,max(z.pred_a)*1.5),
                      theta = 28, phi = 20, ticktype = "detailed", clim = c(0,max(z.pred_a,z)),
                      xlab = "SSB", ylab = selected_var, zlab = "SP", clab = "Surplus production",
                      surf = list(x = x.pred_a, y = y.pred2, z = z.pred_a, facets = NA, fit = z),
                      main = "Additive model: Production curve", sub = knobi_results$data$Stock)
    plot3d_add<-grDevices::recordPlot()

    if (plot_out==T){
      grDevices::jpeg("additive_model.jpeg",width=2500, height=2500,res=300)
      grDevices::replayPlot(plot3d_add)
      grDevices::dev.off()
    }

    cat("\n ... Just a little more. May the 4th be with you ... \n")

    plot3D::scatter3D(x, y2, z, bty="b2", pch = 19, cex = 1, cex.axis = 0.6,
                      colkey = list(length = 0.5, width = 0.5, cex.clab = 0.8, cex.axis = 0.8),
                      xlim = c(0, max(x.pred_m)), zlim = c(0,max(z.pred_m)*1.5),
                      theta = 28, phi = 18, ticktype = "detailed", clim=c(0,max(z.pred_m,z)),
                      xlab = "SSB", ylab = selected_var, zlab = "SP", clab = "Surplus production",
                      surf = list(x = x.pred_m, y = y.pred2, z = z.pred_m, facets = NA, fit = z),
                      main = "Multiplicative model: Production curve", sub = knobi_results$data$Stock)
    plot3d_mult<-grDevices::recordPlot()

    if (plot_out==T){
      grDevices::jpeg("multiplicative_model.jpeg",width=2500, height=2500,res=300)
      grDevices::replayPlot(plot3d_mult)
      grDevices::dev.off()
    }

    res_env$plots3D=list(additive_plot=plot3d_add,multiplicative_plot=plot3d_mult)


  } else {

    envs=df_env[,-c(1:4)]

    cor=round(cor(as.matrix(envs),use="na.or.complete"),4)
    for(i in 1:ncol(cor)){
      cor[i,i]=0
    }
    cor=max(abs(cor))
    if(cor>=0.35){
      p.mat <- cor.mtest(envs)
      corrplot::corrplot(round(cor(as.matrix(envs),use="na.or.complete"),2),
                         title="Covariables correlation",mar=c(0,0,2,0),method="number",
                         type="lower",diag=T,p.mat = p.mat, sig.level = 0.05)
      print(round(cor(as.matrix(envs),use="na.or.complete"),4))
      cat("\n")
      warning("The covariables are highly correlated (see corrplot). To avoid estimation issues, consider removing variables or setting multicovar=FALSE.")
    }

    if(ncol(envs)==1){stop("The number of covariables is too short. Almost two environmental variables are required.")}
    if(ncol(envs)>=6){stop("The number of covariables is too big. Consider reducing the amount of environmental variables.")}

    Data=list(env=envs,data=df, start_r=as.numeric(knobi_results$fit$Parameter_estimates[[1]]),
              start_K=as.numeric(knobi_results$fit$Parameter_estimates[[2]]),
              start_c=start_c[1], start_p=1)

    if (knobi_results$control$pella){

      class(Data)=paste0("Pella_Mult_",ncol(envs))
      model_Pella_Mult <-fitting(Data)

      class(Data)=paste0("Pella_Add_",ncol(envs))
      model_Pella_Add <-fitting(Data)

      res_env$model_env_Multiplicative = model_Pella_Mult
      res_env$model_env_Additive = model_Pella_Add

      class(model_Pella_Mult) = "Pella_Mult_m"
      ref_pts_mult = RF(model_Pella_Mult)
      class(model_Pella_Add) = "Pella_Add_m"
      ref_pts_add = RF(model_Pella_Add)
      res_env$ref_pts = list(ref_pts_mult=ref_pts_mult,ref_pts_add=ref_pts_add)

      model_Pella_Mult=list(model_Pella_Mult)
      model_Pella_Mult$data$SP=df$y
      if(knobi_results$control$method=="Biomass"){model_Pella_Mult$data$B=df$x} else {
        model_Pella_Mult$data$SSB=df$x}
      model_Pella_Mult$data$env=envs
      class(model_Pella_Mult)="Pella_Mult_m"

      model_Pella_Add=list(model_Pella_Add)
      model_Pella_Add$data$SP=df$y
      if(knobi_results$control$method=="Biomass"){model_Pella_Add$data$B=df$x} else {
        model_Pella_Add$data$SSB=df$x}
      model_Pella_Add$data$env=envs
      class(model_Pella_Add)="Pella_Add_m"

      class(knobi_results$fit)="Pella_Year"
      bv <- predict_model(knobi_results$fit)
      bv1 <- predict_model(model_Pella_Mult)
      bv2 <- predict_model(model_Pella_Add)

    } else {

      class(Data)=paste0("Schaefer_Mult_",ncol(envs))
      model_Schaefer_Mult <-fitting(Data)

      class(Data)=paste0("Schaefer_Add_",ncol(envs))
      model_Schaefer_Add <-fitting(Data)

      res_env$model_env_Multiplicative = model_Schaefer_Mult
      res_env$model_env_Additive = model_Schaefer_Add

      class(model_Schaefer_Mult) = "Schaefer_Mult_m"
      ref_pts_mult = RF(model_Schaefer_Mult)
      class(model_Schaefer_Add) = "Schaefer_Add_m"
      ref_pts_add = RF(model_Schaefer_Add)
      res_env$ref_pts = list(ref_pts_mult=ref_pts_mult,ref_pts_add=ref_pts_add)

      model_Schaefer_Mult=list(model_Schaefer_Mult)
      model_Schaefer_Mult$data$SP=df$y
      if(knobi_results$control$method=="Biomass"){ model_Schaefer_Mult$data$B=df$x} else {
        model_Schaefer_Mult$data$SSB=df$x}
      model_Schaefer_Mult$data$env=envs
      class(model_Schaefer_Mult)="Schaefer_Mult_m"

      model_Schaefer_Add=list(model_Schaefer_Add)
      model_Schaefer_Add$data$SP=df$y
      if(knobi_results$control$method=="Biomass"){model_Schaefer_Add$data$B=df$x} else {
        model_Schaefer_Add$data$SSB=df$x}
      model_Schaefer_Add$data$env=envs
      class(model_Schaefer_Add)="Schaefer_Add_m"

      class(knobi_results$fit)="Schaefer_Year"
      bv <- predict_model(knobi_results$fit)
      bv1 <- predict_model(model_Schaefer_Mult)
      bv2 <- predict_model(model_Schaefer_Add)
    }
    res_env$environmental_variables=envs
  }

  Environmental=res_env

  envplot_df=data.frame(SP=c(df$y,bv,bv2,bv1),Year=rep(df$Year,4),
                        factor=c(rep("Observed",length(bv)),rep("Base model",length(bv)),
                                 rep("Environmental Additive",length(bv)),
                                 rep("Environmental Multiplicative",length(bv))))

  env_plot=ggplot2::ggplot(data=envplot_df,ggplot2::aes(x=Year,y=SP,color=factor)) + ggplot2::theme_bw() +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::ylim(min(envplot_df$SP),max(envplot_df$SP)) +
    ggplot2::labs(title="Environmental fits", subtitle=knobi_results$data$Stock,
                  y="Surplus Production") +
    ggplot2::theme(legend.position = c(0.15,0.85), plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),legend.title=ggplot2::element_blank(),
                   legend.background = ggplot2::element_rect(fill = "transparent"))

  print(env_plot)

  if (plot_out==T){
    p <- grDevices::recordPlot()
    grDevices::jpeg("fits_env.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()}

  Environmental$error <- error(knobi_results, Environmental, plot_out)

  if (plot_out==T){
    cat(paste0("\n ... Done! :) \n \n Plots successfully saved in '",getwd(),"'"),". \n")
    setwd(old_dir)
  } else {cat("\n ... Done! :) \n")}

  class(Environmental)="knobi"

  return(Environmental)

}
