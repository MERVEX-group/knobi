#' @title Print a knobi object
#'
#' @description The default print method for a \code{\link{knobi_fit}} or a \code{\link{knobi_env}} object
#'
#' @param x,... Fitted model objects of class \code{knobi} produced by \code{knobi_fit()} or \code{knobi_env()}.
#'
#' @details Prints out the formula and the parameters estimates for the base KBPM fit or the environmental KBPM fit.
#'
#' @author
#' \itemize{
#' \item{Anxo Paz}
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' \item{M. Grazia Pennino}
#' }
#'
#' @seealso
#' \code{\link{knobi_fit}}, \code{\link{knobi_env}}
#'
#' @export
#'

print.knobi<-function(x, ...){

  if(is.null(x$model_env_Multiplicative)==T){
    if(x$control$pella==T){
      cat("\n Formula:\n","SP_t = (r/p)*B_t*(1-(B_t/K)^p) \n \n")
    } else {
      cat("\n Formula:\n","SP_t = r*B_t*(1-B_t/K) \n \n")}

    cat("Coefficients:\n")

    cat(names(x$fit$Parameter_estimates)[1]," ",x$fit$Parameter_estimates[1],"\n")
    cat(names(x$fit$Parameter_estimates)[2]," ",x$fit$Parameter_estimates[2],"\n")

    if(x$control$pella==T){
      cat(names(x$fit$Parameter_estimates)[3]," ",x$fit$Parameter_estimates[3],"\n \n")
    } else {cat("\n")}

  } else {

    if("p" %in% names(x$model_env_Multiplicative)){

      cat("\n Multiplicative model:\n","SP_t = (r/p)*B_t*(1-(B_t/K)^p)*exp(c%*%X_t) \n \n")
    } else {
      cat("\n Multiplicative model:\n","SP_t = r*B_t*(1-B_t/K)*exp(c%*%X_t) \n \n")
    }

    cat("Coefficients:\n")

    npms=length(x$model_env_Multiplicative)

    for(i in 1:npms){
      cat(names(x$model_env_Multiplicative)[i]," ",x$model_env_Multiplicative[i],"\n")
    }

    cat("\n")

    if("p" %in% names(x$model_env_Multiplicative)){

      cat("\n Additive model:\n","SP_t = (r/p)*B_t*(1-(B_t/K)^p)+c%*%X_tB_t \n \n")
    } else {
      cat("\n Additive model:\n","SP_t = r*B_t*(1-B_t/K)+c%*%X_tBt \n \n")
    }

    cat("Coefficients:\n")

    for(i in 1:npms){
      cat(names(x$model_env_Additive)[i]," ",x$model_env_Additive[i],"\n")
    }

    cat("\n")

  }
}
