#' @title KBPM retrospective analysis
#'
#' @description This function carries out the retrospective analysis evaluating the robustness of the KBPM fit.
#'
#' @param knobi_results A list containing the results of the KBPM fit. Object provided by \code{\link{knobi_fit}} function (main function).
#' @param yR Vector containing the years in which the catch time series ends in each of the retrospective analysis settings.
#' @param yR0 Optional argument. Vector containing the years in which the catch time series starts in each of the retrospective analysis settings. Equal length of 'yR' vector is required. By default it is assumed that the catch time series starts in the same year as in the original fit.
#' @param nR Number of retrospective patterns. Only required when 'yR' is not provided (if both arguments are included and error message is reported). See details.
#' @param plot_out Logical. TRUE means that a file with the plot of the retrospective fits is created. The default value is the input of this argument in the knobi_fit function.
#' @param plot_dir Optional directory for creating the folder and save the plots. Required when plot_out=TRUE. The default value is the input of this argument in the knobi_fit function.
#' @param plot_filename Optional name of the folder that will contain the plots. Required when plot_out=TRUE. The default value is the input of this argument in the knobi_fit function.
#'
#' @details If 'nR' is provided it specifies the number of fits to carry out. The first model considers the data deleting the last year and fits the surplus production curve, the next model deletes the two last years of the original data set and fits the SP curve, and the procedure continues until the last SP curve is fitted over the data excluding the last nR years. If the 'yR' argument is provided, the procedure is analogous but the number of years deleted in each of the retrospective fits is specified through the end years provided by this vector from the beginning of the time series or from the years specified in 'yR0'.
#'
#' @return A list containing the retrospective analysis results is provided. It includes the fits and the list of corresponding reference points.
#' The estimate surplus production curves of the retrospective analysis are plotted. The plot is shown in the plot window and saved (if plot_out=TRUE) on the provided directory or in the current directory.
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
#' library(knobi)
#'
#' # See knobi_fit example to obtain the knobi_results object
#' knobi_retrospectives<-knobi_retro(knobi_results,nR=5,plot_out=T)
#' knobi_retrospectives
#'
#' knobi_retro(knobi_results,yR=c(2010,2015),yR0=c(1995,2000))
#' }
#'
#' @export



knobi_retro<-function(knobi_results,yR=NULL,yR0=NULL,nR=NULL,plot_out=F,plot_filename=NULL,plot_dir=NULL){

  Year<-knobi_results$data$years
  lastyear<-max(Year)

  df <- data.frame(x = knobi_results$data$Average_Biomass,
                   y = knobi_results$data$SP, Year = Year)

  x<-knobi_results$data$Average_Biomass
  r<-knobi_results$fit$Parameter_estimates[[1]]
  K<-knobi_results$fit$Parameter_estimates[[2]]
  if (knobi_results$control$pella){p<-knobi_results$fit$Parameter_estimates[[3]]}

  cut<-K
  av <- seq(0, cut, length.out = 3*length(x))
  bv <- predict_model(knobi_results$fit); fit_base<-(knobi_results$fit)

  df_aux<-data.frame(av,bv)


  if(is.null(yR)==TRUE){

    if(is.null(nR)==TRUE){stop("You must provide values for nR or yR")}

    nR<-nR

    modelretro<-list()
    vec<-data.frame(matrix(0, ncol=5,nrow=nR))
    names_retro<-NULL

    for (i in 1:nR){
      newdf<-subset(df, Year<(lastyear-i+1))
      Data<-list(data=newdf, start_r=0.5, start_K=max(newdf$x))
      if (knobi_results$control$pella){
        Data$start_p<-1
        class(Data)<-"Pella"
      } else{class(Data)<-"Schaefer"}

      modelr <-fitting(Data)
      fit<-list(modelr$par)
      class(fit)<-class(Data)
      val<-RF(fit)

      vec[i,]<-c(val$K,
                 val$B_MSY,
                 val$F_MSY,
                 val$MSY,
                 val$MSYoverK)
      name_i<-paste(min(df$Year), "-", lastyear-i)
      names_retro<-c(names_retro,name_i)
      modelretro[[name_i]]<-fit

    }
  } else {

    if(is.null(nR)==FALSE){stop("You must provide values only for nR or yR")}

    nR<-length(yR)
    yR<-yR
    if(is.null(yR0)==TRUE){yR0<-rep(Year[1],length(yR))
    } else {
      if(length(yR)!=length(yR0)){stop("yR and yR0 must have the same length")}
      yR0<-yR0
    }

    Year<-knobi_results$data$years

    modelretro<-list()
    vec<-NULL; names_retro<-NULL

    for (i in 1:nR){
      ind<-which(Year %in% yR[i])
      newdf<-subset(df, Year<=Year[ind])
      ind2<-which(Year %in% yR0[i])
      newdf<-subset(newdf, Year>=Year[ind2])

      Data<-list(data=newdf, start_r=0.5, start_K=max(newdf$x))

      if (knobi_results$control$pella){
        Data$start_p<-1
        class(Data)<-"Pella"
      } else{class(Data)<-"Schaefer"}

      modelr <-fitting(Data)
      fit<-list(modelr$par)
      class(fit)<-class(Data)
      val<-RF(fit)

      name_i<-paste(Year[ind2], "-", Year[ind])
      names_retro<-c(names_retro,name_i)

      vec_i<-c(val$K,val$B_MSY,val$F_MSY,val$MSY,val$MSYoverK)
      vec<-rbind(vec,vec_i)

      modelretro[[name_i]]<-fit
    }

  }

  Retrospectives<-modelretro
  colnames(vec)<-c("K","B_MSY","F_MSY","MSY","MSYoverK")
  rownames(vec)<-names_retro
  Retrospectives$RP<-data.frame(vec)

  # Plots Retrospective  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(is.null(yR)==FALSE){nR<-length(yR)}
  df_plot<-data.frame(B=df_aux$av,SP=df_aux$bv,factor=rep(paste(min(df$Year), "-", lastyear),length(df_aux$av)))
  df_plot$factor<-as.character(df_plot$factor)

  for (i in 1:nR){
    fit<-modelretro[[i]]

    fit$data$SP<-df$y
    if(knobi_results$control$method=="Biomass"){fit$data$B<-df$x} else {fit$data$SSB<-df$x}
    x<-fit$data[[2]]
    r<-fit[[1]][1]
    K<-fit[[1]][2]
    if (knobi_results$control$pella){p<-fit[[1]][3]}

    cut<-K
    av <- seq(0, cut, length.out = 3*length(x))
    bretro <- predict_model(fit)
    aretro<-av
    df_plot<-data.frame(B=c(df_plot$B,aretro),SP=c(df_plot$SP,bretro),
                        factor=c(df_plot$factor,rep(names(modelretro[i]),3*length(x))))
    df_plot$factor<-as.character(df_plot$factor)
  }

  if(knobi_results$control$method=="SSB"){
    btit<-"SP curve and observed SSB and SP"
    baxis<-"Spawning biomass (SSB)"
    bleg<-"observed SSB"} else {
      btit<-"SP curve and observed Biomass and SP"
      baxis<-"Biomass"
      bleg<-"observed biomass"}


  max_y<-max(df_plot$SP,df$y)
  min_y<-min(df_plot$SP,df$y)

  retro_plot<-ggplot2::ggplot() + ggplot2::theme_bw() +
    ggplot2::geom_point(data=df[c(1,nrow(df)),],ggplot2::aes_(x=~x,y=~y,color=~bleg),size=3, show.legend=FALSE) +
    ggplot2::geom_text(data=df[c(1,nrow(df)),],ggplot2::aes_(x=~x,y=~y,color=~bleg,label=~Year,vjust=-1),
                       size=4, show.legend = FALSE) +
    ggplot2::geom_point(data=df,ggplot2::aes_(x=~x,y=~y,color=~bleg),show.legend=FALSE) +
    ggplot2::geom_path(data=df,ggplot2::aes_(x=~x,y=~y,color=~bleg)) +
    ggplot2::labs(title=btit, subtitle=knobi_results$data$Stock,
                  x =baxis, y = "Surplus Production (SP)") +
    ggplot2::geom_line(data=df_plot,ggplot2::aes_(x=~B,y=~SP,color=~factor)) + ggplot2::ylim(min_y,max_y) +
    ggplot2::guides(size="none",col=ggplot2::guide_legend(title="Retrospectives")) +
    ggplot2::theme(legend.position = c(.89,0.75), legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   axis.line=ggplot2::element_line())

  print(retro_plot)

  if(plot_out==TRUE){
    old_dir<-getwd()
    if (is.null(plot_dir)) {plot_dir<-knobi_results$control$plot_settings$plot_dir}
    setwd(plot_dir)
    if (is.null(plot_filename)){plot_filename<-knobi_results$control$plot_settings$plot_filename}
    if (plot_filename %in% list.dirs(full.names=FALSE)){
      setwd(paste0(plot_dir,"/",plot_filename))} else {
        dir.create(plot_filename)
        setwd(paste0(plot_dir,"/",plot_filename))}
    p <- grDevices::recordPlot()
    grDevices::jpeg("fits_retro.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()
    cat(paste0("Plot successfully saved in '",getwd(),"'"),"\n")
    setwd(old_dir)
  }

  return(Retrospectives)

}

