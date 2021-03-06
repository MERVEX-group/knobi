% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knobi_proj.R
\name{knobi_proj}
\alias{knobi_proj}
\title{KBPM projections}
\usage{
knobi_proj(
  knobi_results,
  env_results = NULL,
  Ct = NULL,
  f = NULL,
  env = NULL,
  end_y = NULL,
  n_y = NULL,
  plot_out = F,
  plot_filename = NULL,
  plot_dir = NULL
)
}
\arguments{
\item{knobi_results}{The output object of \code{\link{knobi_fit}} function (main function).}

\item{env_results}{Optional. The output object of \code{\link{knobi_env}} function.}

\item{Ct}{Optional. Vector, data frame or matrix containing the values of the catch for the projected. Different catch scenarios are allowed. The values for each should be provided in each of the columns. The length vector or row number should be equal to the number of years. Projections can be based on selected catch or fishing mortality values, then only one of the arguments, Ct or f, must be introduced.}

\item{f}{Optional. Vector, data frame or matrix containing the values of the fishing mortality for the projected. Different catch scenarios are allowed. The values for each should be provided in each of the columns. The length vector or row number should be equal to the number of years. Projections can be based on selected catch or fishing mortality values, then only one of the arguments, Ct or f, must be introduced.}

\item{env}{Optional. If the multicovar argument of \code{\link{knobi_env}} is FALSE, a vector, data frame or matrix containing the values of the environmental covariates (unstandardized) for the projection years (rows) and the different catch or fishing mortality settings (columns).  On the other hand, if the multicovar argument of \code{\link{knobi_env}} is TRUE, the current argument must be a list, and each argument must be a data frame corresponding to each catch or fishing mortality setting containing the values of the environmental covariates for that scenario.}

\item{end_y}{Optional. The end year of our projections. If not provided, enter n_y in the argument above. If both are provided, the function uses n_y.}

\item{n_y}{Optional. The number of years for projections. If not provided, enter end_y in the next argument.}

\item{plot_out}{Logical. TRUE means that a file with the  environmental fit plots is created. By default this argument is FALSE.}

\item{plot_filename}{Optional name of the folder that will contain the plots. Required when plot_out=TRUE. The default value is the input of this argument in the knobi_fit function.}

\item{plot_dir}{Optional directory for creating the folder and save the plots. Required when plot_out=TRUE. The default value is the input of this argument in the knobi_fit function.}
}
\value{
A list containing the projection results is provided. \itemize{
\item base_model: Three-dimensional array containing projection information for each catch or fishing mortality scenario that is projected using the KBPM base model (see \code{\link{knobi_fit}} details). Rows correspond to years, columns to the derived quantities (biomass or SSB, surplus production, F and catch) and third dimension to the projection scenarios.
\item additive_model: Four-dimensional array containing projections information for each catch or fishing mortality scenario and for each environmental scenario considering the KBPM additive model (see \code{\link{knobi_env}} details). The rows correspond to the years, the columns to the derived quantities (biomass or SSB, surplus production, F and catch), the third dimension to the environmental projection scenarios and the fourth one to the catch scenarios. It is only returned if the arguments env_results and env are provided.
\item multiplicative_model: Four-dimensional array containing projections information for each catch or fishing mortality scenario and for each environmental scenario considering the KBPM multiplicative model (see \code{\link{knobi_env}} details). The rows correspond to the years, the columns to the derived quantities (biomass or SSB, surplus production, F and catch), the third dimension to the environmental projection scenarios and the fourth one to the catch scenarios. It is only returned if the arguments env_results and env are provided.
\item biomass: Data frame containing historical biomass and their projections for all scenarios.
\item catch: Data frame containing historical catch and their projection values for all scenarios.
\item f: Data frame containing the historical fishing mortality and their projection values for all scenarios.
\item SP: Data frame containing the historical surplus production and their projection values for all scenarios.}
The plots results are displayed in the plot window and are also saved (if plot_out="TRUE") in the  provided directory or in the same directory as \code{link{knobi_fit}}.
If environmental information is not provided, four plots are presented in a panel reporting the biomass, surplus production, catch and fishing mortality projections for each catch or fishing mortality scenario.
If environmental information is provided, the plots are presented in a panel for each catch or fishing mortality scenario reporting   biomass, surplus production, catch and fishing mortality projections in each of the environmental scenarios.
}
\description{
This function projects the stock biomass  (or stock spawning biomass) time series and then the surplus production  based on the selected catch or fishing mortality values for the projected years.
}
\examples{

\dontrun{

# First, run the example of knobi_fit and knobi_env function

# Then, create the data frame containing the selected catch for the projected years.
# In this illustration, catch values for all the projected years are constant and
# equal to the catch value in the last year multiplied by 1, 1.2 and 0.8.

catch=rep(knobi_results$data$Catch[length(knobi_results$data$Catch)],5)

Ct=data.frame(catch=catch,
              catch08=0.8*catch,
              catch12=1.2*catch)


# Then, create the data frame containing the environmental covariable values
# for the projected years.

env=data.frame(AMO1=c(0.2,0.3,0.4,0.5,0.6),
               AMO2=c(0.05,0.15,0.25,0.35,0.45),
               AMO3=c(-0.1,-0.2,-0.3,-0.4,-0.5))


# Based on the previous objects we can apply the projection function.

knobi_proj(knobi_results, knobi_environmental, Ct=Ct, n_y=5, env=env)

# An example without environmental information for the projections is given below.
# The number of projection years is also extended.

Ct=rbind(Ct,Ct[1:3,])
knobi_proj(knobi_results, Ct=Ct, end_y=2027)

 # Alternatively, projections can be based on fishing mortality.
 # The scenarios presented below have been created from the estimated F_msy of
 # knobi_fit analysis.

fmsy=knobi_results$fit$RP$F_MSY
ff=rep(fmsy,5)
f=data.frame(f=ff,f12=ff*1.2,f08=ff*0.8)

knobi_proj(knobi_results, f=f, n_y=5, env_results=env_results, env=env)

# The example without environmental information for the projections is given below.
# The number of projection years is also extended.

f=rbind(f,f[1:3,])
knobi_proj(knobi_results, f=f, end_y=2027)


# In case of multicovar=T in knobi_env, a list is required which each
# item is a data frame for each environmental scenario

env=list(climate_1=data.frame(AMO=c(0.2,0.2,0.3,0.3,0.4),
                              Tmax_Vigo=c(19,19,20,20,21)),
         climate_2=data.frame(AMO=c(0.2,0.3,0.4,0.5,0.6),
                              Tmax_Vigo=c(19,20,21,22,23)))

knobi_proj(knobi_results, knobi_environmental2, Ct=Ct[1:5,], n_y=5, env=env)
}

}
\author{
\itemize{
\item{Anxo Paz}
\item{Marta Cousido-Rocha}
\item{Santiago Cerviño López}
\item{M. Grazia Pennino}
}
}
