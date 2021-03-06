% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Data_Fn.R
\name{Data_Fn}
\alias{Data_Fn}
\title{Build data input for MIST}
\usage{
Data_Fn(Version, obsmodel_p = NULL, b_i, s_i, t_i, p_i, a_x, MeshList,
  n_factors = 1, n_cointegrate = NULL, B_type = "Real_eigenvalue",
  startFromEquilibriumTF = FALSE, c_ktp = NULL, catchlogsd = 0.01,
  spatial_method = 0, CheckForErrors = TRUE)
}
\arguments{
\item{Version}{a version number (see example for current default).}

\item{obsmodel_p}{The observation model for each species p
\describe{
  \item{ObsModel=0}{Poisson}
  \item{ObsModel=1}{Lognormal}
  \item{ObsModel=2}{Zero-inflated lognormal}
  \item{ObsModel=3}{lognormal-Poisson}
  \item{ObsModel=4}{Normal}
}}

\item{b_i}{Sampled biomass per unit area for each observation i}

\item{s_i}{Spatial knot (e.g., grid cell) for each observation i}

\item{t_i}{Time interval (e.g., year) for each observation i}

\item{p_i}{Species for each observation i}

\item{a_x}{Area associated with each knot}

\item{MeshList,}{tagged list representing location information for the SPDE mesh hyperdistribution, i.e., from \code{SpatialDeltaGLMM::Spatial_Information_Fn}}

\item{n_factors}{Rank of covariance matrix for process error}

\item{n_cointegrate}{Number of regulatory relations in community matrix (default is number of species)}

\item{B_type}{\describe{
  \item{"Independent"}{Independent spatial-Gompertz dynamics for each species}
  \item{"Real_eigenvalue"}{Co-integration with eigenvalues restricted to real numbers}
  \item{"Complex_eigenvalue"}{Co-integration with eigenvalues including complex numbers}
}}

\item{startFromEquilibriumTF}{whether species start from equilibrium densities (i.e., turning of phi)}

\item{c_ktp}{harvest (biomass per area, in same units as \code{b_i} and \code{a_x}), where \code{c_ktp=0} involves no harvest (the default)}

\item{catchlogsd}{log-standard deviation of observed harvest specified via \code{c_ktp}}

\item{spatial_method}{DEPRECATED, always uses "Mesh" approximation}

\item{CheckForErrors}{Boolean, whether to check for errors in data inputs}
}
\value{
Tagged list containing inputs to function \code{MIST::Build_TMB_Fn()}
}
\description{
\code{Data_Fn} builds a tagged list of data inputs used by TMB for running the model
}

