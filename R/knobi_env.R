
#' @title KBPM environmental analysis
#'
#' @description Analyze and model the relationships between surplus production (SP) and environmental covariable(s) to test whether productivity changes in response to environmental fluctuations in three steps: (1) the analysis of the correlation between the environmental variable(s) at different delays (lags) and the KBPM residuals through the pearson's correlation or through the fits of autoregressive models; (2) the choice of the environmental variable and its time lag in the correlation with surplus production that will be included in the models; (3) the KBPM environmental fit, where environmental effects are included as additive and multiplicative effects in the KBPM formulation (see details).
#'
#' @param knobi_results The output object of \code{\link{knobi_fit}} function (main function).
#' @param data A list containing the following data: \itemize{
#' \item env: data frame containing the values of each one of the environmental variable(s) in one column. Rows represent years.
#' \item years: years in which the environmental variable(s) are reported.}
#' @param control Optional. List containing the following settings: \itemize{
#' \item nlag: this argument is used to test, in the correlation analysis, the lags smaller or equal to 'nlag' (natural number). The lag corresponding to the highest correlation among the KBPM residuals and the corresponding time lagged environmental covariable is considered in the environmental model (unless otherwise specified in 'selected_var'). This means that the KBPM residuals_{t} are related to X_{t-lag} being X the environmental variable and lag the selected value from sequence {0,1,...,nlag}. By default, 'nlag=3'. See details.
#' \item lag: optional numerical vector providing the lag value(s) to consider in the relation between the KBPM surplus production and the environmental variable(s). The length of this argument must be equal to the number of environmental variables included. Applies only if 'nlag' argument is not provided.
#' \item start_c: optional. Numerical vector providing the starting values of the environmental parameter 'c' for the additive and multiplicative models, respectively. By default, start_c=c(1,1). See details.
#' \item ar_cor: optional. Logical. By default this argument is FALSE, meaning that the correlation between the KBPM residuals and the environmental variable(s) is analyzed through a pearson correlation test, as described above. If this argument is "TRUE", the relationship  between the KBPM residuals and the environmental variables is analyzed by fitting autoregressive (AR) models including each one of the environmental time lagged variables as explanatory covariables. The environmental variable whose model reports the lowest Akaike information criterion (AIC) is selected to be included in the environmental KBPM fit. See details.
#' \item plot3d: optional. Logical. If this argument is TRUE, 3D plots reporting the surplus production curve conditioned to a grid of environmental values are provided. FALSE by default.
#' \item selected_var: optional. Character. By default, the fit is carried out  including the highest correlated environmental variable at the corresponding time lag. However, if this argument is equal to the name of one of the environmental variables, this variable is used in the environmental fit considering the lag derived from the correlation analysis.
#' \item multicovar: optional. Logical. TRUE  means that the environmental model includes all the input environmental covariables at the same time, up to a maximum of 5. By default this argument is FALSE, which means that only the environmental variable reporting the highest correlation is included (after lagging it if corresponds).}
#' @param plot_out Logical. TRUE means that a file with the plot of the retrospective fits is created. The default value is the argument input in the \code{\link{knobi_fit}} function.
#' @param plot_dir Optional. Directory to create the folder and save the plots. Required when plot_out=TRUE. The default value is the argument input in the \code{\link{knobi_fit}} function.
#' @param plot_filename Optional. Name of the folder that will contain the plots. Required when plot_out=TRUE. The default value is the argument input in the \code{\link{knobi_fit}} function.
#'
#' @details
#' Additive environmental model adds the following term on the right hand of equation (1) or (2) described in \code{\link{knobi_fit}} function: \eqn{cX_{t}B_{t}}, being \eqn{X_{t}}the environmental variable and \eqn{B_{t}} the biomass or SSB at time \eqn{t}.
#' Multiplicative environmental model multiplies the right hand of equation (1) or (2) by \eqn{exp(cX_{t})}.
#' The subscript t denotes the time (years).
#'
#' If ar_cor argument is "TRUE", for the correlation analysis between the KBPM residuals and the environmental variable(s) is carried out as follows.
#' First, an AR model is fitted to the residuals.
#' \deqn{r_t=\sum_{i=1}^{p}\beta_{i}r_{t-i}+\epsilon_{t}}
#' being \eqn{r_t} the KBPM base residual for year \eqn{t} and \eqn{p} the AR model order, estimated as the maximum time lag at which the absolute value of the partial autocorrelation value  is large than \eqn{qnorm(0.975)/\sqrt(length(r))}.
#' Then, AR models are fitted considering each one of the lagged environmental variable(s),
#' \deqn{r_{t,lag}=\sum_{i=1}^{p}\beta_{i}r_{t-i}+X_{t,lag}+\epsilon_{t}}, for \eqn{lag=0,1,...,nlag}.
#' being \eqn{X_{t,lag}} the lagged environmental variable at year {t} and time lag. Then, we have an autoregressive model for each of the lagged environmental variables.
#' The AIC values of the above models are compared, and the lagged environmental variable whose model reports the lowest AIC is used in the KBPM fit, except if the argument 'lag' is used.
#'
#' @return A list containing the results of the three-step environmental analysis is provided. \itemize{
#' \item selected_lag: data frame providing the time lag of the environmental variable(s) in the KBPM fit.
#' \item lag_cor: correlation between the environmental variable(s) and the KBPM residuals for each one of the time lags.
#' \item env_aic: if 'ar_cor=TRUE', AIC values of each one of the autoregressive models (see details).
#' \item selected_var: environmental variable used in the fit.
#' \item model_env_Multiplicative: estimates of the multiplicative model parameters.
#' \item model_env_Additive: estimates of the additive model parameters.
#' \item ref_pts: reference points (RPs) estimates for each model assuming X_t=0 (see vignettes).
#' \item scaled_environmental_var: standardized environmental variable(s) used in the fit, that means that the variable(s) are transformed subtracting its mean and dividing by its standard deviation (sd) with the 'scale' and 'center' attributes, which are the sd and the mean of the variable respectively.
#' \item plots3D: list with the 3D plots objects (if 'plot3d=TRUE').
#' \item goodness_of_fit: list of goodness-of-fit measures of the fits: \itemize{
#' \item residuals: pearson residuals from the fit for each model (base KBPM, additive model and multiplicative model).
#' \item error_table: array of performance and accuracy measures for each model: Standard error of the regression (SER), coefficient of determination (R-squared), adjusted coefficient of determination (adj-R-squared), Akaike information criterion (AIC), root-mean-squared error (RMSE) and mean absolute percentage error (MAPE) of observed and estimated surplus production values; and the statistic and p-value of the F-test whose null hypothesis establishes that the AR model without environmental covariable(s) is statistically better than the environmental AR model.}}
#' Results are also reported through plots which are displayed in the plot window and also saved (if 'plot_out=TRUE') in the provided directory or in the same directory as knobi_fit's plots.
#' The first plot reports the correlation analysis between the environmental variable(s) and the KBPM residuals. The second one reports the fitted SP values derived from the base model (no environmental information) and from the environmental ones.
#' If 'multicovar=FALSE' and 'plot3d=TRUE', 3D plots reporting the surplus production curve conditioned to a grid of environmental values are also reported.
#' If 'multicovar=TRUE', a plot with the pearson correlation between the environmental variables is reported too.
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
#' Env<-data.frame(years=seq(1973,2020),AMO=c(-0.236,-0.441,-0.32,-0.385,-0.21,-0.201,-0.132,
#' -0.041,-0.098,-0.235,-0.093,-0.23,-0.29,-0.297,0.044,-0.028,-0.106,-0.061,-0.155,-0.242,
#' -0.234,-0.2,0.112,-0.082,0.028,0.349,0.094,0.004,0.095,0.041,0.207,0.182,0.268,0.242,
#' 0.123,0.114,0.015,0.325,0.078,0.189,0.142,0.077,0.09,0.318,0.291,0.045,0.15,0.279),
#' TMax_Vigo=c(16.5,16.7,17.1,16.3,16.4,16.7,17.3,17.3,17.4,18,17.5,17.5,17.6,18,17,18.4,
#' 17.8,19.6,19.1,18,17.9,17.7,17.7,19.4,18.2,19.7,19.2,18.6,18,18.3,18.5,18.8,18.7,18.7,
#' 18.5,17.9,17.4,19.2,19,19.7,18,19.1,19.4,20,19.5,20.2,18.8,18.6))
#'
#' # The environmental data series must start in the first year of the KBPM fit data
#' # minus the provided nlag or lag
#' years<-knobi_results$data$years # See knobi_fit example to obtain the knobi_results object
#' ind<-which(Env[,1]==years[1])
#' ind1<-which(Env[,1]==years[length(years)])
#' nlag<-5
#' Env<-Env[(ind-nlag):ind1,]
#'
#' # Now we create the environmental list
#' data<-list(env=data.frame(AMO=Env$AMO,Tmax=Env$TMax_Vigo),
#'            years=Env$years)
#' control<-list(nlag=nlag)
#'
#' knobi_environmental<-knobi_env(knobi_results,data,control)
#' knobi_environmental
#' }
#'
#' @export



knobi_env<-function(knobi_results,data,control=NULL,plot_out=FALSE,plot_filename=NULL,plot_dir=NULL){

  if(is.null(control$ar_cor)==TRUE) {control$ar_cor=FALSE}

  env0<-as.data.frame(data$env)
  env_names<-names(env0)
  y_env<-data$years
  if(is.null(control$start_c)){start_c<-c(1,1)} else {
    start_c<-control$start_c}

  if(plot_out==T){
    old_dir<-getwd()
    if (is.null(plot_dir)) {plot_dir<-knobi_results$control$plot_settings$plot_dir}
    setwd(plot_dir)
    if (is.null(plot_filename)){plot_filename<-knobi_results$control$plot_settings$plot_filename}
    if (plot_filename %in% list.dirs(full.names=FALSE)){
      setwd(paste0(plot_dir,"/",plot_filename))} else {
        dir.create(plot_filename)
        setwd(paste0(plot_dir,"/",plot_filename))}
  }

  df<-data.frame(KBPM_residuals=knobi_results$fit$goodness_of_fit$residuals,
                 x = knobi_results$data$Average_Biomass,
                 y = knobi_results$data$SP, Year=knobi_results$data$years)


  f_year<-knobi_results$data$years[1]
  l_year<-knobi_results$data$years[length(knobi_results$data$years)]

  env<-list()
  res_env<-list()
  df_env<-df

  res_env$selected_lag<-array(NA,dim=c(length(env_names),2))
  colnames(res_env$selected_lag)<-c("lag","correlation")
  rownames(res_env$selected_lag)<-env_names

  if (is.null(control$lag)){

    if(is.null(control$nlag)){
      lag=3
    } else {
      lag<-control$nlag
    }

  } else {

    lag<-max(control$lag)

    res_env$selected_lag[,1]<-control$lag

  }

  data_env<-list()

  res_env$lag_cor<-array(NA,dim=c(length(env_names),lag+1))

  vec_env<-"lag_0"
  for (i in 1:lag){
    vec_env<-c(vec_env,paste0("lag_",i))
  }

  colnames(res_env$lag_cor)<-vec_env
  rownames(res_env$lag_cor)<-env_names


  for(j in env_names){

    data_env[[j]]<-df

    if(!is.na(knobi_results$data$Recruitment[1])){
      data_env[[j]]<-data.frame(data_env[[j]],knobi_results$data$Recruitment)}

    ind<-which(y_env==f_year)
    ind1<-which(y_env==l_year)
    if(length(ind)>0){data_env[[j]]$env0<-env0[ind:ind1,j]} else {
      warning('The length of the environmental variable is not enough to use the number of input lags')}

    for (i in 1:lag){
      ind<-which(y_env==(f_year-i))
      ind1<-which(y_env==(l_year-i))
      if(!is.na(knobi_results$data$Recruitment[1])){
        if(length(ind)>0){data_env[[j]][,6+i]<-env0[[j]][ind:ind1]} else {
          warning('The length of the environmental variable is not enough to use the number of input lags')}
      } else {
        if(length(ind)>0){data_env[[j]][,5+i]<-env0[[j]][ind:ind1]} else {
          warning('The length of the environmental variable is not enough to use the number of input lags')}
      }
    }

    if(knobi_results$control$method=="Biomass"){
      if(!is.na(knobi_results$data$Recruitment[1])){
        colnames(data_env[[j]])<-c("KBPM_residuals","SP","B","years","R",vec_env)
      } else {colnames(data_env[[j]])<-c("KBPM_residuals","SP","B","years",vec_env)}
    } else {
      if(!is.na(knobi_results$data$Recruitment[1])){
        colnames(data_env[[j]])<-c("KBPM_residuals","SP","SSB","years","R",vec_env)
      } else {colnames(data_env[[j]])<-c("KBPM_residuals","SP","SSB","years",vec_env)}
    }

    env[[j]]<-data_env[[j]][,c("KBPM_residuals",vec_env)]

    cor<-round(cor(as.matrix(env[[j]]),use="na.or.complete"),4)
    cor<-cor[,1]; cor<-cor[-1]

    if (is.null(control$lag)){

      res_env$selected_lag[j,2]<-cor[[which(max(abs(cor))==abs(cor))[1]]]
      res_env$selected_lag[j,1]<-which(max(abs(cor))==abs(cor))[1]-1

    } else {

      res_env$selected_lag[j,2]<-cor[res_env$selected_lag[j,1]+1]

    }

    df_env[,j]<-scale(env[[j]][,res_env$selected_lag[j,1]+2])

    res_env$lag_cor[j,]<-cor

  }

  if(control$ar_cor==TRUE){

    KBPM_residuals<-knobi_results$fit$goodness_of_fit$residuals
    pacf_res<-stats::pacf(KBPM_residuals,plot=FALSE)$acf[,1,1]

    ref<-stats::qnorm(0.975)/sqrt(length(KBPM_residuals))

    auto<-max(which(abs(pacf_res)>=ref),0)
    if(auto==0){
      warning("KBPM base residuals are not autocorrelated. An AR(0) model is fitted for SP residuals.")
    }

    fit_base <- stats::arima0(KBPM_residuals, order = c(auto, 0, 0))

    prenv_aic<-array(NA,dim=c(length(env_names),lag+1))
    colnames(prenv_aic)<-vec_env; rownames(prenv_aic)<-env_names

    for(j in env_names){
      for(i in vec_env){
        env_fit <- stats::arima0(KBPM_residuals, order = c(auto, 0, 0), xreg = env[[j]][,i])
        prenv_aic[j,i]<-env_fit$aic
      }
    }

    env_aic<-cbind(base=rep(fit_base$aic,length(env_names)),prenv_aic)
    res_env$env_aic<-env_aic

    env_aic_c<-env_aic-fit_base$aic+2
    min_aic<-env_aic_c[which.min(env_aic_c)]
    if(min_aic>=0){
      warning("AR models considering environmental variable(s) do not really improve AR model considering only the residuals")
    }

    colnames(res_env$selected_lag)[2]<-"aic"

    if (is.null(control$lag)){

      for(j in env_names){
        res_env$selected_lag[j,1]<-which(env_aic == min(env_aic[j,-1]), arr.ind=TRUE)[2]-2
        res_env$selected_lag[j,2]<-min(env_aic[j,-1])
      }

    } else {

      for(j in env_names){
        res_env$selected_lag[j,2]<-env_aic[j,(res_env$selected_lag[j,1]+2)]
      }
    }

    for(j in env_names){
      df_env[,j]<-scale(env[[j]][,res_env$selected_lag[j,1]+2])
    }


  }

  if(control$ar_cor==FALSE){

    lagf<-NULL
    corlist<-NULL

    for(i in vec_env){
      lagf<-c(lagf,rep(i,length(env_names)))
      corlist<-c(corlist,res_env$lag_cor[,i])
    }

    envcorplot_df<-data.frame(correlation=corlist,lag=lagf,factor=rep(env_names,length(vec_env)))

    envcorplot<-ggplot2::ggplot(data=envcorplot_df,ggplot2::aes_(x=~lag,y=~correlation,group=~factor,color=~factor)) +
      ggplot2::theme_bw() + ggplot2::geom_point() + ggplot2::geom_line(linetype = "dashed") + ggplot2::ylim(-1,1) +
      ggplot2::labs(title="Environmental correlation with base KBPM SP residuals", subtitle=knobi_results$data$Stock,
                    y="Correlation",x="") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5),legend.title=ggplot2::element_blank(),
                     legend.background = ggplot2::element_rect(fill = "transparent"))

    print(envcorplot)

  } else {

    lagf<-NULL
    corlist<-NULL

    for(i in vec_env){
      lagf<-c(lagf,rep(i,length(env_names)))
      corlist<-c(corlist,res_env$env_aic[,i])
    }

    maximo<-max(env_aic)
    minimo<-min(env_aic)

    envcorplot_df<-data.frame(AIC=corlist,lag=lagf,factor=rep(env_names,length(vec_env)))

    envcorplot<-ggplot2::ggplot(data=envcorplot_df,ggplot2::aes_(x=~lag,y=~AIC,group=~factor,color=~factor)) +
      ggplot2::theme_bw() + ggplot2::geom_point() + ggplot2::geom_line(linetype = "dashed") + ggplot2::ylim(minimo,maximo) +
      ggplot2::labs(title="AIC comparison", subtitle=knobi_results$data$Stock,
                    y="AIC",x="") +
      ggplot2::geom_hline(yintercept = fit_base$aic) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5),legend.title=ggplot2::element_blank(),
                     legend.background = ggplot2::element_rect(fill = "transparent"))

    print(envcorplot)

  }


  if (plot_out==TRUE){
    p <- grDevices::recordPlot()
    grDevices::jpeg("corplot.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()}

  if(is.null(control$multicovar)){
    multicovar<-FALSE
  } else {
    multicovar<-control$multicovar
  }

  if(multicovar==FALSE){

    if(is.null(control$selected_var)){

      if(control$ar_cor==FALSE){
        selected_var<-env_names[which.max(abs(res_env$selected_lag[,2]))]
      } else {
        selected_var<-env_names[which.min(abs(res_env$selected_lag[,2]))]
      }

    } else {

      selected_var<-control$selected_var

    }

    res_env$selected_var<-selected_var

    Data<-list(env=df_env[,selected_var],data=df, start_r=as.numeric(knobi_results$fit$Parameter_estimates[[1]]),
               start_K=as.numeric(knobi_results$fit$Parameter_estimates[[2]]),
               start_c=start_c[1], start_p=1)

    if (knobi_results$control$pella){

      # Pella Mult

      class(Data)<-"Pella_Mult"
      model_Pella_Mult <-fitting(Data)
      res_env$model_env_Multiplicative <- model_Pella_Mult
      class(model_Pella_Mult) <- class(Data)
      ref_pts_mult <- RF(model_Pella_Mult)

      # Pella Additive

      class(Data)<-"Pella_Add"
      model_Pella_Add <-fitting(Data)
      res_env$model_env_Additive <- model_Pella_Add
      class(model_Pella_Add) <- class(Data)
      ref_pts_add <- RF(model_Pella_Add)

      res_env$ref_pts <- list(ref_pts_mult=ref_pts_mult,ref_pts_add=ref_pts_add)

      model_Pella_Mult<-list(model_Pella_Mult)
      model_Pella_Mult$data$SP<-df$y
      if(knobi_results$control$method=="Biomass"){model_Pella_Mult$data$B<-df$x} else {
        model_Pella_Mult$data$SSB<-df$x}
      model_Pella_Mult$data$env<-as.numeric(df_env[,selected_var])
      class(model_Pella_Mult)<-"Pella_Mult"

      model_Pella_Add<-list(model_Pella_Add)
      model_Pella_Add$data$SP<-df$y
      if(knobi_results$control$method=="Biomass"){model_Pella_Add$data$B<-df$x} else {
        model_Pella_Add$data$SSB<-df$x}
      model_Pella_Add$data$env<-as.numeric(df_env[,selected_var])
      class(model_Pella_Add)<-"Pella_Add"

      class(knobi_results$fit)<-"Pella_Year"
      bv <- predict_model(knobi_results$fit)
      bv1 <- predict_model(model_Pella_Mult)
      bv2 <- predict_model(model_Pella_Add)

    } else {

      # Schaefer Mult

      class(Data)<-"Schaefer_Mult"
      model_Schaefer_Mult <-fitting(Data)
      res_env$model_env_Multiplicative <- model_Schaefer_Mult
      class(model_Schaefer_Mult) <- class(Data)
      ref_pts_mult <- RF(model_Schaefer_Mult)

      # Schaefer Additive

      class(Data)<-"Schaefer_Add"
      model_Schaefer_Add <-fitting(Data)
      res_env$model_env_Additive <- model_Schaefer_Add
      class(model_Schaefer_Add) <- class(Data)
      ref_pts_add <- RF(model_Schaefer_Add)

      res_env$ref_pts <- list(ref_pts_mult=ref_pts_mult,ref_pts_add=ref_pts_add)

      model_Schaefer_Mult<-list(model_Schaefer_Mult)
      model_Schaefer_Mult$data$SP<-df$y
      if(knobi_results$control$method=="Biomass"){ model_Schaefer_Mult$data$B<-df$x} else {
        model_Schaefer_Mult$data$SSB<-df$x}
      model_Schaefer_Mult$data$env<-as.numeric(df_env[,selected_var])
      class(model_Schaefer_Mult)<-"Schaefer_Mult"


      model_Schaefer_Add<-list(model_Schaefer_Add)
      model_Schaefer_Add$data$SP<-df$y
      if(knobi_results$control$method=="Biomass"){model_Schaefer_Add$data$B<-df$x} else {
        model_Schaefer_Add$data$SSB<-df$x}
      model_Schaefer_Add$data$env<-as.numeric(df_env[,selected_var])
      class(model_Schaefer_Add)<-"Schaefer_Add"

      class(knobi_results$fit)<-"Schaefer_Year"
      bv <- predict_model(knobi_results$fit)
      bv1 <- predict_model(model_Schaefer_Mult)
      bv2 <- predict_model(model_Schaefer_Add)
    }

    res_env$scaled_environmental_var<-df_env[,selected_var]
    colnames(res_env$scaled_environmental_var)<-selected_var
    rownames(res_env$scaled_environmental_var)<-df_env$Year-as.numeric(res_env$selected_lag[selected_var,1])


    if(is.null(control$plot3d)){plots3d<-FALSE} else {plots3d<-control$plot3d}

    if(plots3d==TRUE){

      x <- knobi_results$data$Average_Biomass
      y <- df_env[,selected_var]
      z <- knobi_results$data$SP

      r_a<-res_env$model_env_Additive[1]
      K_a<-res_env$model_env_Additive[2]
      c_a<-res_env$model_env_Additive[3]
      if(knobi_results$control$pella==TRUE){
        p_a<-res_env$model_env_Additive[4]
      } else {p_a<-1}

      r_m<-res_env$model_env_Multiplicative[1]
      K_m<-res_env$model_env_Multiplicative[2]
      c_m<-res_env$model_env_Multiplicative[3]
      if(knobi_results$control$pella==TRUE){
        p_m<-res_env$model_env_Multiplicative[4]
      } else {p_m<-1}


      grid.lines <- 400
      cut_a<-max(K_a+c_a*K_a*max(y),K_a+c_a*K_a*min(y))
      x.pred_a <- seq(0, cut_a, length.out = grid.lines)
      x.pred_m <- seq(0, K_m, length.out = grid.lines)

      y.pred <- c(seq(min(y),0,length.out = grid.lines/2),seq(0, max(y), length.out = grid.lines/2))
      xy_a <- expand.grid(x = x.pred_a, y = y.pred)
      xy_m <- expand.grid(x = x.pred_m, y = y.pred)

      z.pred_a <- (r_a/p_a)*xy_a$x*(1-(xy_a$x/K_a)^(p_a))+c_a*xy_a$y*xy_a$x
      z.pred_m <- exp(1)^{c_m*xy_m$y}*((r_m/p_m)*xy_m$x*(1-(xy_m$x/K_m)^p_m))
      z.pred_a <- matrix(z.pred_a,nrow=grid.lines,ncol=grid.lines)
      z.pred_m <- matrix(z.pred_m,nrow=grid.lines,ncol=grid.lines)

      y2 <- y * attr(y, 'scaled:scale') + attr(y, 'scaled:center')
      y.pred2 <- y.pred * attr(y, 'scaled:scale') + attr(y, 'scaled:center')

      cat("\n This can take a while... \n")

      plot3D::scatter3D(x, y2, z, bty="b2", pch = 19, cex = 1, cex.axis = 0.6,
                        colkey = list(length = 0.5, width = 0.5, cex.clab = 0.8, cex.axis = 0.8),
                        xlim = c(0,max(x.pred_a)), zlim = c(0,max(z.pred_a)*1.5),
                        theta = 28, phi = 20, ticktype = "detailed", clim = c(0,max(z.pred_a,z)),
                        xlab = "SSB", ylab = selected_var, zlab = "SP", clab = "Surplus production",
                        surf = list(x = x.pred_a, y = y.pred2, z = z.pred_a, facets = NA, fit = z),
                        main = "Additive model: Production curve", sub = knobi_results$data$Stock)
      plot3d_add<-grDevices::recordPlot()

      if (plot_out==TRUE){
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

      if (plot_out==TRUE){
        grDevices::jpeg("multiplicative_model.jpeg",width=2500, height=2500,res=300)
        grDevices::replayPlot(plot3d_mult)
        grDevices::dev.off()
      }

      cat(paste0("\n ... Done! :) \n"))


      res_env$plots3D<-list(additive_plot=plot3d_add,multiplicative_plot=plot3d_mult)
    }

  } else {

    envs<-df_env[,-c(1:4)]

    cor<-round(cor(as.matrix(envs),use="na.or.complete"),4)
    for(i in 1:ncol(cor)){
      cor[i,i]<-0
    }
    cor<-max(abs(cor))
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

    Data<-list(env=envs,data=df, start_r=as.numeric(knobi_results$fit$Parameter_estimates[[1]]),
               start_K=as.numeric(knobi_results$fit$Parameter_estimates[[2]]),
               start_c=start_c[1], start_p=1)

    if (knobi_results$control$pella){

      class(Data)<-paste0("Pella_Mult_",ncol(envs))
      model_Pella_Mult <-fitting(Data)

      class(Data)<-paste0("Pella_Add_",ncol(envs))
      model_Pella_Add <-fitting(Data)

      res_env$model_env_Multiplicative <- model_Pella_Mult
      res_env$model_env_Additive <- model_Pella_Add

      class(model_Pella_Mult) <- "Pella_Mult_m"
      ref_pts_mult <- RF(model_Pella_Mult)
      class(model_Pella_Add) <- "Pella_Add_m"
      ref_pts_add <- RF(model_Pella_Add)
      res_env$ref_pts <- list(ref_pts_mult=ref_pts_mult,ref_pts_add=ref_pts_add)

      model_Pella_Mult<-list(model_Pella_Mult)
      model_Pella_Mult$data$SP<-df$y
      if(knobi_results$control$method=="Biomass"){model_Pella_Mult$data$B<-df$x} else {
        model_Pella_Mult$data$SSB<-df$x}
      model_Pella_Mult$data$env<-envs
      class(model_Pella_Mult)<-"Pella_Mult_m"

      model_Pella_Add<-list(model_Pella_Add)
      model_Pella_Add$data$SP<-df$y
      if(knobi_results$control$method=="Biomass"){model_Pella_Add$data$B<-df$x} else {
        model_Pella_Add$data$SSB<-df$x}
      model_Pella_Add$data$env<-envs
      class(model_Pella_Add)<-"Pella_Add_m"

      class(knobi_results$fit)<-"Pella_Year"
      bv <- predict_model(knobi_results$fit)
      bv1 <- predict_model(model_Pella_Mult)
      bv2 <- predict_model(model_Pella_Add)

    } else {

      class(Data)<-paste0("Schaefer_Mult_",ncol(envs))
      model_Schaefer_Mult <-fitting(Data)

      class(Data)<-paste0("Schaefer_Add_",ncol(envs))
      model_Schaefer_Add <-fitting(Data)

      res_env$model_env_Multiplicative <- model_Schaefer_Mult
      res_env$model_env_Additive <- model_Schaefer_Add

      class(model_Schaefer_Mult) <- "Schaefer_Mult_m"
      ref_pts_mult <- RF(model_Schaefer_Mult)
      class(model_Schaefer_Add) <- "Schaefer_Add_m"
      ref_pts_add <- RF(model_Schaefer_Add)
      res_env$ref_pts <- list(ref_pts_mult=ref_pts_mult,ref_pts_add=ref_pts_add)

      model_Schaefer_Mult<-list(model_Schaefer_Mult)
      model_Schaefer_Mult$data$SP<-df$y
      if(knobi_results$control$method=="Biomass"){ model_Schaefer_Mult$data$B<-df$x} else {
        model_Schaefer_Mult$data$SSB<-df$x}
      model_Schaefer_Mult$data$env<-envs
      class(model_Schaefer_Mult)<-"Schaefer_Mult_m"

      model_Schaefer_Add<-list(model_Schaefer_Add)
      model_Schaefer_Add$data$SP<-df$y
      if(knobi_results$control$method=="Biomass"){model_Schaefer_Add$data$B<-df$x} else {
        model_Schaefer_Add$data$SSB<-df$x}
      model_Schaefer_Add$data$env<-envs
      class(model_Schaefer_Add)<-"Schaefer_Add_m"

      class(knobi_results$fit)<-"Schaefer_Year"
      bv <- predict_model(knobi_results$fit)
      bv1 <- predict_model(model_Schaefer_Mult)
      bv2 <- predict_model(model_Schaefer_Add)
    }
    res_env$scaled_environmental_var<-envs
  }

  Environmental<-res_env

  envplot_df<-data.frame(SP=c(df$y,bv,bv2,bv1),Year=rep(df$Year,4),
                         factor=c(rep("Observed",length(bv)),rep("Base model",length(bv)),
                                  rep("Environmental Additive",length(bv)),
                                  rep("Environmental Multiplicative",length(bv))))

  env_plot<-ggplot2::ggplot(data=envplot_df,ggplot2::aes_(x=~Year,y=~SP,color=~factor)) + ggplot2::theme_bw() +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::ylim(min(envplot_df$SP),max(envplot_df$SP)) +
    ggplot2::labs(title="Environmental fits", subtitle=knobi_results$data$Stock,
                  y="Surplus Production") +
    ggplot2::theme(legend.position = c(0.15,0.85), plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),legend.title=ggplot2::element_blank(),
                   legend.background = ggplot2::element_rect(fill = "transparent"))

  print(env_plot)

  if (plot_out==TRUE){
    p <- grDevices::recordPlot()
    grDevices::jpeg("fits_env.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()}

  Environmental$goodness_of_fit <- knobi_error(knobi_results, Environmental, plot_out)

  if (plot_out==TRUE){
    cat(paste0("\n Plots successfully saved in '",getwd(),"'"),". \n")
    setwd(old_dir)
  }

  class(Environmental)<-"knobi"

  return(Environmental)

}
