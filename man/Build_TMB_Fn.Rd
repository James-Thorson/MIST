% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Build_TMB_Fn.R
\name{Build_TMB_Fn}
\alias{Build_TMB_Fn}
\title{Build TMB object for geostatistical delta-GLMM}
\usage{
Build_TMB_Fn(TmbData, Version, estimate_phi = TRUE, Kappa = "constant",
  Parameters = "generate", Random = "generate", Map = "generate",
  use_REML = FALSE, loc_x = NULL, eigenbounds = c(Lower = -2, Upper =
  -0.01), TmbDir = system.file("executables", package = "MIST"),
  RunDir = getwd())
}
\arguments{
\item{Version}{a version number (see example for current default).}

\item{estimate_phi,}{Estimate the multiplicative offset between initial and equilibrium density}

\item{Kappa,}{controls kappa (decorrelation distance) of spatial variation in productivity among species, as well as the spatio-temporal "process error"
\describe{
  \item{"constant"}{Same value for all spatial plus the spatio-temporal component}
  \item{"spatial_vs_spatiotemporal"}{Same value for all spatial, but different for the spatio-temporal component}
  \item{"different"}{Different for all spatial components as well as the spatio-temporal component}
}}

\item{Parameters}{OPTIONAL, a tagged list of starting parameters}

\item{Random}{OPTIONAL, a character vector of random effects}

\item{Map}{OPTIONAL, a tagged list of parameters to either mirror or turn off}

\item{loc_x}{OPTIONAL, location for each sample used to generate plausible bounds for scale parameter}

\item{eigenbounds}{upper and lower bounds for the estimated eigenvalues of the community matrix}

\item{TmbDir}{OPTIONAL, a directory where the CPP file for the VAST model can be found locally}

\item{RunDir}{OPTIONAL, a directory where the CPP file is copied, copiled, and run (must have write privileges or else the function will crash)}

\item{TmbbData,}{a tagged list of data inputs generated by \code{Data_Fn}}

\item{Use_REML,}{OPTIONAL boolean whether to use maximum marginal likelihood or restricted maximum likelihood (termed "REML")}

\item{Use_REML,}{OPTIONAL boolean whether to use maximum marginal likelihood or restricted maximum likelihood (termed "REML")}

\item{silent}{Boolean, whether TMB should run silently}
}
\value{
Tagged list containing objects for running a VAST model
\describe{
  \item{Obj}{The built TMB object}
  \item{Upper}{A vector of upper bounds for parameters, optionally for use during optimization}
  \item{Lower}{A vector of lower bounds for parameters, optionally for use during optimization}
  \item{Parameters}{A tagged list of parameter starting values used when building Obj, which can be extracted, modified, and then put back into \code{Build_TMB_Fn} to define different starting values}
  \item{Map}{A taggged list of parameters to be turned off or mirrored, for similar use as Parameters}
  \item{Random}{A character vector of random effects, for similar use as Parameters}
}
}
\description{
\code{Build_TMB_Fn} builds a tagged list with everything necessary to run or interpret inputs for \code{MIST}
}

