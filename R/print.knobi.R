#' @title Print a knobi object
#'
#' @description The default print method for a \code{\link{knobi_fit}}, \code{\link{knobi_env}} or a \code{\link{knobi_proj}} object
#'
#' @param x,... Fitted model objects of class \code{knobi} produced by \code{knobi_fit()}, \code{knobi_env()} or \code{knobi_proj()}.
#'
#' @details Prints out the formula and the parameters estimates for the base KBPM fit or the environmental KBPM fit or the biomass and the surplus production estimated projections.
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
#' \code{\link{knobi_fit}}, \code{\link{knobi_env}}, \code{\link{knobi_proj}}
#'
#' @export
#'

print.knobi<-function(x, ...){

  if(is.null(x$fit)==FALSE){
    if(x$control$pella==TRUE){
      cat("\n Formula:\n","SP_t = (r/p)*B_t*(1-(B_t/K)^p) \n \n")
    } else {
      cat("\n Formula:\n","SP_t = r*B_t*(1-B_t/K) \n \n")}

    cat("Parameter estimates:\n")

    cat(names(x$fit$Parameter_estimates)[1]," ",x$fit$Parameter_estimates[1],"\n")
    cat(names(x$fit$Parameter_estimates)[2]," ",x$fit$Parameter_estimates[2],"\n")

    if(x$control$pella==TRUE){
      cat(names(x$fit$Parameter_estimates)[3]," ",x$fit$Parameter_estimates[3],"\n \n")
    } else {cat("\n")}

  }

  if(is.null(x$model_env_Multiplicative)==FALSE){

    if("p" %in% names(x$model_env_Multiplicative)){

      cat("\n Multiplicative model:\n","SP_t = (r/p)*B_t*(1-(B_t/K)^p)*exp(c*X_t) \n \n")
    } else {
      cat("\n Multiplicative model:\n","SP_t = r*B_t*(1-B_t/K)*exp(c*X_t) \n \n")
    }

    cat("Parameter estimates:\n")

    npms<-length(x$model_env_Multiplicative)

    for(i in 1:npms){
      cat(names(x$model_env_Multiplicative)[i]," ",x$model_env_Multiplicative[i],"\n")
    }

    cat("\n")

    if("p" %in% names(x$model_env_Multiplicative)){

      cat("\n Additive model:\n","SP_t = (r/p)*B_t*(1-(B_t/K)^p)+c*X_tB_t \n \n")
    } else {
      cat("\n Additive model:\n","SP_t = r*B_t*(1-B_t/K)+c*X_tBt \n \n")
    }

    cat("Parameter estimates:\n")

    for(i in 1:npms){
      cat(names(x$model_env_Additive)[i]," ",x$model_env_Additive[i],"\n")
    }

    cat("\n")

  }

  if(is.null(x$biomass)==FALSE){
    cat("\n Biomass projections: \n \n")
    print(subset(x$biomass[,-3],x$biomass$scenario!="input"))
    cat("\n")

    cat("\n Surplus Production projections: \n \n")
    print(subset(x$SP[,-3],x$SP$scenario!="input"))
    cat("\n \n")
  }

}
