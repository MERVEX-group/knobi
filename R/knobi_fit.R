#' @title Known biomass Production Model (KBPM) fit
#'
#' @description This function, that is the main function of the knobi package, fits a type of surplus production models named known-biomass production models (KBPM) (MacCall, 2002). The surplus production curve is fitting using the the catch time series and the biomass or SSB (Spawning Stock Biomass) derived from the fit of other stock assessment model.
#'
#' @param data A list containing the data. \itemize{
#' \item years: time series of years corresponding to the catch time series.
#' \item Catch: time series of catch estimates from a stock assessment model.
#' \item Biomass: time series of biomass estimates from a stock assessment model. If it is available, other case introduce SSB in the next argument.
#' \item Spawning_Biomass: time series of SSB estimates from a stock assessment model. If it is available, other case introduce Biomass in the previous argument.
#' \item Stock: optional. Character string with the stock name for the plot subtitles.
#' \item Recruitment: optional. Time series of recruitment from a stock assessment model. See details.
#' \item F_input: optional. Time series of F estimates from a stock assessment model. See details.
#' \item classF_input: optional. Character indicating the type of F_input estimate. See details.
#' \item RP: optional. Values for any of the following biological reference points derived from a stock assessment model (see details): \itemize{
#' \item F_MSY: estimate of the fishing mortality at Maximum Sustainable Yield (MSY).
#' \item B_MSY: estimate of biomass at MSY (or SSB depending on method argument, see control list).
#' \item MSY: estimate of MSY.
#' \item K: estimate of the virgin biomass.}}
#' @param control A list containing the control parameters. \itemize{
#' \item pella: Logical. TRUE means that Pella-Tomlinson model is used. FALSE means that Schaefer model is fitted (by default). See details.
#' \item start_r: optional start value of the model parameter r (growth rate parameter, i.e., intrinsic rate of natural increase). See details.
#' \item start_K: optional start value of the model parameter K (the maximum population size for growth to be positive, i.e., the virgin biomass concept related to the carrying capacity). See details.
#' \item start_p: optional start value of the model parameter p (shape parameter) when Pella-Tomlinson model is used. See details.
#' \item method: establishes if the fits is carried using "SSB" or "Biomass". The argument is only required if both time series, Spawning_Biomass and Biomass, are provided to the function.}
#' @param plot_out Logical. TRUE means that files with the plots of the input time series and the fit results are created. By default FALSE.
#' @param plot_dir Optional directory for creating the folder and save the plots. Required when plot_out=TRUE.
#' @param plot_filename Optional name of the folder that will contain the plots. By default, "knobi_results". Required when plot_out=TRUE.
#'
#'
#' @details
#'
#' @return An output list updated with the fitting information. The control is the input one updated with the plot settings information. The data is also updated including the annual average biomass (mean of two consecutive years), in $data$Average_Biomass, the surplus production, in $data$SP, and the F estimates derive from knobi_fit, in $data$F_output. The fit results of the KBPM are in the list $fit, which contains: \itemize{
#' \item Parameter_estimates: Estimates of the model parameters.
#' \item data: The data used for the model.
#' \item RP: KBPM estimates of the biological reference points. \itemize{
#' \item K:  KBPM estimate of virgin biomass.
#' \item B_MSY: KBPM estimate of biomass at maximum sustainable yield.
#' \item F_MSY: KBPM estimate of fishing mortality at maximum sustainable yield.
#' \item MSY: KBPM estimate of maximum sustainable yield.
#' \item MSYoverK: ratio of MSY and K. }
#' \item optimr: A list of some results provided by \code{\link[optimr]{optimr}}: \itemize{
#' \item value: The value of the function corresponding to the parameter estimation.
#' \item convergence: An integer code. ‘0’ indicates successful completion in the optimization.
#' \item message: A character string giving any additional information returned by the optimizer, or NULL}
#' \item rror: List of performance and accuracy: \itemize{
#' \item residuals: Pearson's residuals from the fit calculated as (observations-estimates)/sqrt(estimates).
#' \item error_table: Data frame measures of estimates accuracy (error measures comparing observed an estimated values) and model performance: \itemize{
#' \item SER: Standard error of the regression, calculated as the root of the rate between the residual sum of squares and the degrees of freedom of the regression.
#' \item R-squared: Coefficient of determination.
#' \item adj-R-squared: Adjusted coefficient of determination.
#' \item AIC: Akaike information criterion.
#' \item RMSE: Root mean squared error.
#' \item MAPE: Mean absolute percentage error.}}}
#' The plots are shown in the plot window and also saved (if plot_out=TRUE) on the provided directory or in the current directory (if the directory is not provided). The following input quantities are plotted: time series of fishing mortality, SSB, surplus production and catch time series. Also plots of catch over fishing mortality, fishing mortality over SSB, and catch over SSB time series are available with a smooth line from a "loess" regression. Plot of input and output time series of fishing mortality with horizontal lines at fishing mortality at MSY ( one line if input F_MSY is NULL) is provided. The analogous SSB plot is also reported. On the other hand, the fitted surplus production curve is plotted twice with the SSB and SP observations (first one) and with the catch and SP observations (second one).
#'
#' @author
#' \itemize{
#' \item{Anxo Paz}
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' \item{M. Grazia Pennino}
#' }
#' @details The KBPMs implemented in the current package are explained below.
#' Schaefer model (1):
#' \deqn{SP_{t} = r B_{t} (1-(B_{t}/K))}
#' where \eqn{SP_{t}} is the surplus production, \eqn{B_{t}} is the biomass or \eqn{SSB_{t}} averaged (mean of two consecutive years), \eqn{r} is the population growth rate parameter, and \eqn{K} is the virgin biomass. The subscript \eqn{t} denotes the time (years).
#' Pella and Tomlinson model (2):
#' \deqn{SP_{t} = (r/p) B_{t} (1-(B_{t}/K)^{p})}
#' where \eqn{SP_{t}} is the surplus production, \eqn{B_{t}} is the biomass or \eqn{SSB_{t}} averaged,  \eqn{r} is the population growth rate parameter, \eqn{K} is the virgin biomass and \eqn{p} is the asymmetry parameter. The subscript \eqn{t} denotes the time (years).
#'
#' The recruitment values, the F_input estimates and the type of F estimation are included to have an overview of the stock status but are not required for the KBPM fitting. Similarly, the input reference points are only used for comparison.
#' On the other hand, KBPMs have also proven their usefulness for the multispecies management objectives. KBPM approach can be applied to analyze the dynamic of the total aggregated biomass and catch of all targeted fish species in a community through knobi_fit function.
#'
#' @examples
#'
#' library(knobi)
#'
#' # First step, getting the data from the ICES package
#' # install.packages("icesSAG")
#' library(icesSAG)
#'
#' summary_data <- getSAG(stock = "Hake", year = 2021)
#' Database <- subset(summary_data, summary_data[,17] == "hke.27.3a46-8abd")
#' Database <- Database[-nrow(Database),]
#'
#'
#' # Then, create the data object.
#'
#' data<-list()
#' data$Spawning_Biomass=Database$SSB # We take the SSB in our Database.
#' data$Catch=Database$catches # We take the catch in our Database.
#' data$F_input=Database$F # We take the F in our Database.
#' # Reference points estimates from ICES stock assessment model:
#' # ICES. 2021. Working Group for the Bay of Biscay and the Iberian Waters Ecoregion
#' # (WGBIE). ICES Scientific Reports. 3:48.1101 pp.
#' data$RP=list(F_MSY=0.259, B_MSY=207398, MSY=75052, B_0=NA)
#' # In this case, B_MSY is equal to SSB_MSY, since control$method="SSB"
#' # (see control list below).
#' data$classF_input="average" # Character indicating the type of F.
#' data$years=Database$Year    # Years corresponding to the catch values
#' # (can be different than the years corresponding to SSB or biomass).
#'
#'
#' # Now we define the control.
#'
#' control=list()
#' control$pella="TRUE" # Logical. TRUE means that Pella-Tomlinson model is used.
#'                          # FALSE means that Schaefer model is employed.
#' control$method="SSB" # Information for the fit: "SSB" or "Biomass".
#'
#'
#' # Finally, we can fit the model
#' knobi_results<-knobi_fit(data,control,plot_out=TRUE,plot_filename="results")
#' knobi_results$fit
#'
#'
#' # Fitting multispecific KBPM
#'
#' # Below, a multistock approximation aggregating the
#' # northern and southern stocks of sardine is performed.
#'
#' # Firstly, read southern stock data
#' sardine1 <- getSAG(stock = "pil.27.8c9a", year = 2021)
#' sardine1 <- sardine1[1:44,]
#' sardine1 <- sardine1[,c(1,6,8,12)]
#'
#' # Secondly, read northern stock data
#' sardine2 <- getSAG(stock = "pil.27.8abd", year = 2021)
#' sardine2 <- sardine2[-nrow(sardine2),]
#' sardine2 <- sardine2[,c(1,6,8,12)]
#'
#' # Extract common years of data in both stocks
#' index <- which(sardine1$Year %in% sardine2$Year)
#' sardine1 <- sardine1[index,]
#'
#' # Create a data.frame where the SSB and the catch are
#' # the sum of such data in the two stocks
#' years<-sardine1$Year
#' sardine <- data.frame(years=years,SSB=sardine1$SSB+sardine2$SSB,
#'                       catch=sardine1$catches+sardine2$catches)
#'
#' # Once the total SSB and catch are available
#' # we follow previous KBPM illustration
#' data<-list()
#' data$Spawning_Biomass=sardine$SSB
#' data$Catch=sardine$catch
#' data$years=sardine$years
#'
#' control=list()
#' control$pella="TRUE"
#' control$method="SSB"
#'
#' knobi_results2<-knobi_fit(data,control)
#' knobi_results2$fit
#'
#'
#'
#' @references
#' MacCall, A. (2002). Use of Known-Biomass Production Models to Determine Productivity of West Coast Groundsh Stocks. North American Journal of Fisheries Management, 22, 272-279.
#' @export

knobi_fit<-function(data,control=NULL,plot_out=F,plot_filename=NULL,plot_dir=NULL){

  if(plot_out==T){
    old_dir=getwd()
    if (is.null(plot_dir)) {plot_dir=old_dir}
    setwd(plot_dir)
    if (is.null(plot_filename)){plot_filename="knobi_results"}
    if (plot_filename %in% list.dirs(full.names=F)){
      setwd(paste0(plot_dir,"/",plot_filename))} else {
        dir.create(plot_filename)
        setwd(paste0(plot_dir,"/",plot_filename))}
    plot_settings=list(plot_filename=plot_filename,plot_dir=plot_dir)
    control$plot_settings=plot_settings
  }

  if(is.null(control)){control=list()}
  set.seed(123)

  # Check input data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  C=data$Catch
  years=data$years

  if(is.null(data$Biomass)){
    data$Biomass=NA
    control$method="SSB"
  }

  if(is.null(data$Spawning_Biomass)){
    data$Spawning_Biomass=NA
    control$method="Biomass"
  }

  if(is.na(data$Biomass[1])==T){control$method="SSB"}
  if(is.na(data$Spawning_Biomass[1])==T){control$method="Biomass"}

  if(length(years)!=(length(C))){stop('Length of catch time series is different than the length of years vector')}

  if(control$method=="Biomass"){
    B=data$Biomass
    if(length(B)!=(length(C)+1)){warning('The length of the catch time series is reduced according to biomass time series length')}
  }

  if(control$method=="SSB"){
    B=data$Spawning_Biomass
    if(length(B)!=(length(C)+1)){warning('The length of the catch time series is reduced according to spawning biomass time series length')}
  }

  if(length(B)!=(length(C)) & length(B)!=(length(C)+1)){
    stop('The biomass MUST be provided for the same years as the catch or for such years and the next one')}


  # Compute SP values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SP=B[-1]; l=length(B); B_aver=B[-1]
  for (i in 1:(l-1)){
    SP[i]=as.numeric(B[i+1]-B[i]+C[i])
    B_aver[i]=(B[i+1]+B[i])/2
  }
  B=B_aver

  # Update correctly dimension ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (length(C)!=length(B)){
    lb=length(B)
    C=C[1:lb]
    years=years[1:lb]
    data$Catch=C
    data$years=years
  }

  if(is.null(data$F_input)){data$F_input=NA} else {
    data$F_input=data$F_input[1:length(data$Catch)]}

  if(is.null(data$Recruitment)){data$Recruitment=NA} else {
    data$Recruitment=data$Recruitment[1:length(data$Catch)]}

  # Save SP, F and average biomass ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  data$Average_Biomass=B_aver
  data$SP=SP
  F_out=C/B_aver

  # Care with 0 values
  ind=which(B_aver==0)
  F_out[ind]=NA
  data$F_output=F_out

  # Plots input data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Plot trends
  data=data
  class(data)=gsub(" ", "", paste("KBPM_", control$method))
  class(control)=gsub(" ", "", paste("KBPM_", control$method))

  plotInput(data,plot_out)


  # Adjust the model  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Define new data.frame with the required info
  df <- data.frame(x = data$Average_Biomass, y = data$SP, Year=years)
  Year=years

  # Fit the model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(is.null(control$pella)){control$pella="FALSE"}

  if(is.null(control$start_r)){start_r=0.5000000} else {start_r=control$start_r}
  if(is.null(control$start_K)){start_K=max(df$x)} else {start_K=control$start_K}
  if(is.null(control$start_p)){start_p=1.0000000} else {start_p=control$start_p}

  Data=list(data=df,start_r=start_r,start_K=start_K,start_p=start_p)
  if (control$pella){
    class(Data)="Pella"
  } else{class(Data)="Schaefer"}
  model=fitting(Data)

  fit=list(model$par)
  names(fit)=c("Parameter_estimates")

  fit=fit
  fit$optimr=list()
  fit$optimr$value=model$value
  fit$optimr$convergence=model$convergence
  fit$optimr$message=model$message
  fit$data$SP=df$y
  if(control$method=="Biomass"){fit$data$B=df$x} else {fit$data$SSB=df$x}
  class(fit)=class(Data)


  # Output plots  fit ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  x=fit$data[[2]]
  r=fit[[1]][1]
  K=fit[[1]][2]
  if (control$pella){
    n=fit[[1]][3]
  } else{n=2}
  cut=K
  av <- seq(0, cut, length.out = 3*length(x))
  bv <- predict_model(fit); fit_base=fit


  if(!is.null(data$Stock)){subtitle=data$Stock} else {subtitle=NULL}

  df_aux=data.frame(av,bv)
  if(control$method=="Biomass"){
    xtit="SP curve and observed Biomass and SP"
    xaxis="Biomass"
    xleg="observed biomass"
  } else {
    xtit="SP curve and observed SSB and SP"
    xaxis="Spawn Biomass (SSB)"
    xleg="observed SSB"
  }

  vec=min(df_aux$bv,df$y)
  vec1=max(df_aux$bv,df$y)

  fit_plot=ggplot2::ggplot(data=df_aux,ggplot2::aes(x=av,y=bv)) + ggplot2::theme_bw() +
    ggplot2::geom_line(ggplot2::aes(size=1.5)) + ggplot2::ylim(vec,vec1) +
    ggplot2::geom_point(data=df[c(1,nrow(df)),],ggplot2::aes(x=x,y=y,size=3,color=Year)) +
    ggplot2::geom_text(data=df[c(1,nrow(df)),],ggplot2::aes(x=x,y=y,size=2,color=Year,
                                                            label=Year,vjust=-1), show.legend = FALSE) +
    ggplot2::geom_point(data=df,ggplot2::aes(x=x,y=y,color=Year)) +
    ggplot2::geom_path(data=df,ggplot2::aes(x=x,y=y,color=Year)) +
    ggplot2::labs(title=xtit,x =xaxis, y = "Surplus Production (SP)") +
    ggplot2::guides(size="none",col=ggplot2::guide_legend(title=xleg)) +
    ggplot2::scale_color_gradient(breaks=c(Year[1],Year[length(Year)])) +
    ggplot2::theme(legend.position = c(.9,.85), legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   axis.line=ggplot2::element_line())
  if(!is.null(subtitle)){
    fit_plot=fit_plot+ggplot2::labs(subtitle=subtitle)
  }
  print(fit_plot)

  if(plot_out==T){
    p <- grDevices::recordPlot()
    grDevices::jpeg("fit.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()
  }



  df$C=C

  vec=min(df_aux$bv,df$C)
  vec1=max(df_aux$bv,df$C)

  fitc_plot=ggplot2::ggplot(data=df_aux,ggplot2::aes(x=av,y=bv)) + ggplot2::theme_bw() +
    ggplot2::geom_line(ggplot2::aes(size=1.5)) + ggplot2::ylim(vec,vec1) +
    ggplot2::geom_point(data=df[c(1,nrow(df)),],ggplot2::aes(x=x,y=C,size=3,color=Year)) +
    ggplot2::geom_text(data=df[c(1,nrow(df)),],ggplot2::aes(x=x,y=C,size=2,color=Year,
                                                            label=Year,vjust=-1), show.legend = FALSE) +
    ggplot2::geom_point(data=df,ggplot2::aes(x=x,y=C,color=Year)) +
    ggplot2::geom_path(data=df,ggplot2::aes(x=x,y=C,color=Year)) +
    ggplot2::labs(title="SP curve and observed catch and SP",
                  x ="Spawning Biomass (SSB)", y = "Surplus Production (SP)") +
    ggplot2::guides(size="none",col=ggplot2::guide_legend(title="Observed catch")) +
    ggplot2::scale_color_gradient(breaks=c(Year[1],Year[length(Year)])) +
    ggplot2::theme(legend.position = c(.9,.85), legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5))

  if(!is.null(subtitle)){
    fitc_plot=fitc_plot+ggplot2::labs(subtitle=subtitle)
  }
  print(fitc_plot)

  if(plot_out==T){
    p <- grDevices::recordPlot()
    grDevices::jpeg("fit_catch.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()
  }

  # Compute reference points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # (Jacobson et al. 2002) and Winker et al. (2018)

  fit=RF(fit)

  # Extract

  Bmsy=fit$B_MSY
  Fmsy=fit$F_MSY

  # Plots RF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # F over years

  F_out=data$F_output
  if(is.na(data$F_input[1])) {F_inp=rep(NA,length(F_out))} else {F_inp=data$F_input}
  if(is.null(data$RP$F_MSY)) {Fmsy_inp=NA} else {Fmsy_inp=data$RP$F_MSY}
  Frel_out=F_out/Fmsy
  Frel_inp=F_inp/Fmsy_inp

  if(is.null(data$RP$B_MSY)) {Bmsy_inp=NA} else {Bmsy_inp=data$RP$B_MSY}
  Brel_out=B_aver/Bmsy
  Brel_inp=B_aver/Bmsy_inp

  max_f=max(c(F_inp,F_out,Fmsy,Fmsy_inp)*1.1,na.rm = TRUE)
  min_f=min(c(F_inp,F_out,Fmsy,Fmsy_inp),na.rm = TRUE)
  max_b=max(c(B_aver,Bmsy,Bmsy_inp)*1.1,na.rm = TRUE)
  min_b=min(c(B_aver,Bmsy,Bmsy_inp),na.rm = TRUE)
  max_fr=max(c(Frel_inp,Frel_out,1.1),na.rm = TRUE)
  min_fr=min(c(Frel_inp,Frel_out,0.9),na.rm = TRUE)
  max_br=max(c(Brel_inp,Brel_out,1.1),na.rm = TRUE)
  min_br=min(c(Brel_inp,Brel_out,0.9),na.rm = TRUE)

  if(control$method=="SSB"){bfac=c(rep("SSB from KBPM",length(Year)),rep("original SSB",length(Year)))
  } else {bfac=c(rep("B from KBPM",length(Year)),rep("original B",length(Year)))}
  if(control$method=="SSB"){
    btit="Spawning Biomass (SSB) over Years"
    baxis="Spawning biomass (SSB)"
    Brtit="Relative SSB over Years"} else {
      btit="Biomass over Years"
      baxis="Biomass"
      Brtit="Relative Biomass over Years"}


  plot_df=data.frame(Year=rep(Year,2),f=c(F_out,F_inp),fr=c(Frel_out,Frel_inp),
                     f_factor=c(rep("F from KBPM",length(Year)),rep("original F",length(Year))),
                     b=rep(B_aver,2),br=c(Brel_out,Brel_inp),b_factor=bfac)


  f_plot=ggplot2::ggplot(data=subset(plot_df,!is.na(f)),ggplot2::aes(x=Year,y=f,color=f_factor)) +
    ggplot2::theme_bw() + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::ylim(min_f,max_f) +
    ggplot2::labs(title="Fishing Mortality (F) over Years",y = "Fishing mortality (F)") +
    ggplot2::geom_hline(yintercept=c(Fmsy,Fmsy_inp),linetype="dashed",
                        color = c("#F8766D","#00BFC4"), na.rm=T) +
    ggplot2::annotate("text",x=Year[length(Year)]-1,y=Fmsy,label="Fmsy (KBPM)",
                      color = "#F8766D",size=3,vjust=-1) +
    ggplot2::guides(col=ggplot2::guide_legend(title="")) +
    ggplot2::theme(legend.position = c(.9,.95), legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5))
  if(!is.na(Fmsy_inp)){
    f_plot = f_plot +
      ggplot2::annotate("text",x=Year[1]+1,y=Fmsy_inp,label="original Fmsy",color = "#00BFC4",na.rm=T,size=3,vjust=-1)
  }
  print(f_plot)

  if(plot_out==T){
    p <- grDevices::recordPlot()
    grDevices::jpeg("F_absolute.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()
  }

  # F over years (relative)

  fr_plot=ggplot2::ggplot(data=subset(plot_df,!is.na(fr)),ggplot2::aes(x=Year,y=fr,color=f_factor)) +
    ggplot2::theme_bw() + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::ylim(min_fr,max_fr) + ggplot2::geom_hline(yintercept=1) +
    ggplot2::labs(title="Relative Fishing Mortality (F) over Years",y = "Fishing mortality (F)") +
    ggplot2::guides(col=ggplot2::guide_legend(title="")) +
    ggplot2::theme(legend.position = c(.9,.95), legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5))

  print(fr_plot)

  if(plot_out==T){
    p <- grDevices::recordPlot()
    grDevices::jpeg("F_relative.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()
  }

  # Biomass over years

  b_plot=ggplot2::ggplot(data=plot_df,ggplot2::aes(x=Year,y=b,color="#F8766D")) +
    ggplot2::theme_bw() + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::ylim(min_b,max_b) + ggplot2::labs(title=btit,y=baxis) +
    ggplot2::geom_hline(yintercept=c(Bmsy,Bmsy_inp),linetype="dashed",
                        color = c("#F8766D","#00BFC4"), na.rm=T) +
    ggplot2::annotate("text",x=Year[length(Year)]-1,y=Bmsy,label="Bmsy (KBPM)",
                      color = "#F8766D",size=3,vjust=-1) +
    ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5))
  if(!is.na(Bmsy_inp)){
    b_plot = b_plot +
      ggplot2::annotate("text",x=Year[1]+1,y=Bmsy_inp,label="original Bmsy",color = "#00BFC4",na.rm=T,size=3,vjust=-1)
  }
  print(b_plot)

  if(plot_out==T){
    p <- grDevices::recordPlot()
    grDevices::jpeg("B_absolute.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()
  }


  br_plot=ggplot2::ggplot(data=subset(plot_df,!is.na(br)),ggplot2::aes(x=Year,y=br,color=b_factor)) +
    ggplot2::theme_bw() + ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::ylim(min_br,max_br) +
    ggplot2::labs(title=Brtit ,y =baxis) + ggplot2::geom_hline(yintercept=1) +
    ggplot2::guides(col=ggplot2::guide_legend(title="")) +
    ggplot2::theme(legend.position = c(.9,.95), legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5))

  print(br_plot)

  if(plot_out==T){
    p <- grDevices::recordPlot()
    grDevices::jpeg("B_relative.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()
  }


  RP=list(K=fit$K,B_MSY=fit$B_MSY,F_MSY=fit$F_MSY,
          MSY=fit$MSY,MSYoverK=fit$MSYoverK)
  adjustment<-list(data=data,control=control,fit=list(
    Parameter_estimates=fit$Parameter_estimates,data=fit$data,RP=RP,optimr=fit$optimr))
  class(adjustment$fit)=class(fit)

  adjustment$fit$error=error(adjustment,plot_out=plot_out)

  if(plot_out==T){
    setwd(old_dir)}


  class(adjustment)="knobi"

  return(adjustment)
}
