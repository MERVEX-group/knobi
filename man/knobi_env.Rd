% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knobi_env.R
\name{knobi_env}
\alias{knobi_env}
\title{KBPM environmental analysis}
\usage{
knobi_env(
  knobi_results,
  environmental,
  plot_out = F,
  plot_filename = NULL,
  plot_dir = NULL
)
}
\arguments{
\item{knobi_results}{The output object of \code{\link{knobi_fit}} function (main function).}

\item{environmental}{A list containing the following data and settings: \itemize{
\item data: data frame containing the values of each one of the environmental variable(s) in one column. Each row represents a year.
\item years: time series of years corresponding to the environmental variable(s).
\item lag: optional numerical vector providing the used lag value(s) in the relation among the surplus production and the environmental variable(s). This means that the SP_{t} is related to X_{t-lag} being X the environmental variable. The length of this argument must be equal to the number of environmental variables included.
\item nlag: if lag value is not provided, this argument is used to test all the lags smaller or equal to nlag (numerical vector) through cor.test function. The lag corresponding to the highest pearson correlation among the surplus production and the lagged environmental covariable values is considered in the environmental model.
\item start_c: optional numerical vector providing the start values of the environmental c parameter for the optimization of the additive and multiplicative models, respectively. By default, start_c=c(1,1). See details.
\item selected_var: optional character. By default, the fit is done using the environmental values according to the lag derived from the previous arguments. However, if this argument is equal to the name of the environmental variable no lag is applied to its values.
\item multicovar: optional logical. TRUE if you want to fit the environmental model including all the input environmental covariables, up to a maximum of 5. By default this argument is FALSE, which means that only the environmental covariable reporting the highest pearson correlation is included (after lagging it if corresponds).}}

\item{plot_out}{Logical. TRUE means that a file with the  environmental fit plots is created. By default this argument is FALSE.}

\item{plot_filename}{Optional name of the folder that will contain the plots. Required when plot_out=TRUE. The default value is the input of this argument in the knobi_fit function.}

\item{plot_dir}{Optional directory for creating the folder and save the plots. Required when plot_out=TRUE. The default value is the input of this argument in the knobi_fit function.}
}
\value{
A list containing the environmental analysis is provided. \itemize{
\item selected_lag: Estimated lag corresponding to the one reporting the highest correlation between the environmental variable and the surplus production. Derived if lag is not fixed.
\item fixed_lag: Input value of 'lag' argument.
\item lag_cor: Correlation between the environmental variable(s) value and the surplus production (after lagging the environmental one if corresponds).
\item selected_var: Environmental variable used in the fit, chosen by the user or the one derived from the highest pearson correlation procedure. In case that argument 'multicovar' is omitted, 'NULL' or equal to 'FALSE'.
\item model_env_Multiplicative: Estimates of the multiplicative model parameters.
\item model_env_Additive: Estimates of the additive model parameters.
\item ref_pts: Reference points (RPs) estimates for each model assuming X_t=0 (see vignettes).
\item scaled_environmental_var: Standardized variable used in the fit, with the 'scale' and 'center' attributes.
\item environmental_variables: Standardized covariables used in the fit (if 'multicovar=TRUE'), with the 'scale' and 'center' attributes.
\item plots3D: List with the 3D plots objects.
\item error: List of performance and accuracy: \itemize{
\item residuals: Pearson's residuals from the fit calculated as (observations-estimates)/sqrt(estimates) for each model (base KBPM, additive model and multiplicative model).
\item error_table: Array of performance and accuracy (observed vs. estimated) measures for each model: Standard error of the regression (SER), coefficient of determination (R-squared), adjusted coefficient of determination (adj-R-squared), Akaike information criterion (AIC), root-mean-squared error (RMSE), mean absolute percentage error (MAPE) and the value of the F statistic corresponding to the comparison of each environmental model respect to the base model (F-value) and its corresponding p-value (Pr(>F)).}}
Result plots are shown in the plot window and also saved (if plot_out="TRUE") on the provided directory or in the same directory as knobi_fit.
The first plot reports the correlation analysis between the environmental variable(s) and the KBPM data. The second one reports the fitted values of the base model (no environmental information) and of the environmental ones.
If multicovar=FALSE, 3D plots reporting the surplus production curve conditioned to a grid of environmental values are also reported.
}
\description{
Analyse and model the relationships between surplus production and environmental covariables to test whether productivity changes in response to environmental fluctuations.  Environmental effects are included as additive and multiplicative effects in the general KBPM formulation (see details).
}
\details{
Additive environmental model adds the following term on the right hand of equation (1) or (2) detailed in \code{\link{knobi_fit}} function: \eqn{cX_{t}B_{t}}, being \eqn{X_{t}}the environmental variable and \eqn{B_{t}} the biomass or SSB at time \eqn{t}.
Multiplicative environmental model multiplies the right hand of equation (1) or (2) by \eqn{exp(cX_{t})}.
The subscript t denotes the time (years).
}
\examples{

\dontrun{

# First, run the example of knobi_fit function

# Then, provide environmental data series

Env=data.frame(years=seq(1973,2020),AMO=c(-0.236,-0.441,-0.32,-0.385,-0.21,-0.201,-0.132,
-0.041,-0.098,-0.235,-0.093,-0.23,-0.29,-0.297,0.044,-0.028,-0.106,-0.061,-0.155,-0.242,
-0.234,-0.2,0.112,-0.082,0.028,0.349,0.094,0.004,0.095,0.041,0.207,0.182,0.268,0.242,
0.123,0.114,0.015,0.325,0.078,0.189,0.142,0.077,0.09,0.318,0.291,0.045,0.15,0.279),
TMax_Vigo=c(16.5,16.7,17.1,16.3,16.4,16.7,17.3,17.3,17.4,18,17.5,17.5,17.6,18,17,18.4,
17.8,19.6,19.1,18,17.9,17.7,17.7,19.4,18.2,19.7,19.2,18.6,18,18.3,18.5,18.8,18.7,18.7,
18.5,17.9,17.4,19.2,19,19.7,18,19.1,19.4,20,19.5,20.2,18.8,18.6))

# The environmental data series must start in the first year of the KBPM fit data
# minus the provided nlag or lag
years=knobi_results$data$years # See knobi_fit example to obtain the knobi_results object
ind=which(Env[,1]==years[1])
ind1=which(Env[,1]==years[length(years)])
nlag=5
Env=Env[(ind-nlag):ind1,]

# Now we create the environmental list
environmental=list()
environmental$data=data.frame(AMO=Env$AMO,Tmax=Env$TMax_Vigo)
environmental$years=Env$years
environmental$nlag=c(nlag,nlag)

knobi_environmental<-knobi_env(knobi_results,environmental)
knobi_environmental
knobi_environmental$plots3D$additive_plot

environmental$multicovar=T
knobi_env(knobi_results,environmental,plot_out=T)
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
