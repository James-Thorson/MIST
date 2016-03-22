
#setwd( "C:/Users/James.Thorson/Desktop/Project_git/spatial_VAM/examples" )

# Install package
#devtools::install_github("james-thorson/spatial_VAM", auth_token=[ask Jim Thorson for code])
#devtools::install_github("james-thorson/utilities")
#devtools::install_github("james-thorson/spatial_DFA")

TmbDir = system.file("executables", package="SpatialVAM")
  #TmbDir = "C:/Users/James.Thorson/Desktop/Project_git/spatial_VAM/inst/executables/"

# load libraries
library( TMB )
library( INLA )
library( SpatialVAM )
#library( abind )
#library( RANN ) # nn2()
#library( ThorsonUtilities ) # nn2()

# directory
Date = Sys.Date()
  DateDir = paste0(getwd(),"/",Date)
  dir.create(DateDir)

##############
# Settings
##############

Version = "spatial_vam_v13"
# v2 -- added lognormal obs-model
# v3 -- exploring bug fixes for single-species model
# v4 -- added zero-inflated lognormal distribution, and option for diagonal Cov_pp (i.e., independence among species)
# v5 -- added lognormal-Poisson obs-model
# v6 -- added options for different ObsModel by species
# v7 -- added center of gravity computations, plus added diagnostic jnll_i
# v8 -- decompose B_pp = Alpha_pr * t(Beta_pr) a la co-integration
# v9 -- added penaly on eigenvalues, and also added eigendecomposition parameterization of B_pp
# v10 -- Switched logtauA_p to MargSigmaA_p so I can put bounds on it
# v11 -- switched logkappa to logkappa_z, where slots 0 through n_p-1 control range for spatial variation, and slot n_p controls range for spatio-temporal (SDFA) variation
# v12 -- changed eigen-cointegration parameterization slightly to hopefully improve convergence during simulation experiment
# v13 -- changed eigen-cointegration parameterization again to hopefully improve convergence during simulation experiment

# Other settings
Nfactors = 4        # Number of dynamic factors in process error
Ncointegrate = 4
Use_REML = FALSE
Kmeans_Config = list("n_x"=25, "nstart"=100, "iter.max"=1000)
Estimate_Phi = TRUE      # Phi is the offset of initial and equilibrium abundance
IndependentTF = FALSE # if TRUE, then Cov_pp and B_pp are diagonal
StartFromEquilibriumTF = FALSE    # Definition for expected density at initial year: TRUE=Expectation is equilibrium; FALSE=Expectation is productivity without accounting for interactions/density dependence
B_type = c("Cointegration", "Eigendecomposition", "Eigen-cointegration", "Eigen-cointegration_V2")[4]
Kappa_Independent = TRUE # Only matters if IndependentTF==TRUE
# Cointegration:  B_pp = Alpha_pr %*% t(Beta_pr)
# Eigendecomposition: B_pp = U_pr %*% diag(Lambda) %*% generalized_inverse(U_pr), where U_pr is the eigenvectors, and Lambda is the real eigenvalues

# Compile
file.copy( from=paste0(TmbDir,"/",Version,".cpp"), to=paste0(DateDir,"/",Version,".cpp") )
setwd( DateDir )
compile( paste0(Version,".cpp") )

# Derived
IndependentTF = FALSE # if TRUE, then Cov_pp and B_pp are diagonal
if( IndependentTF==TRUE ) B_type = "Cointegration"
if( IndependentTF==FALSE ) B_type = "Eigen-cointegration_V2"
ObsModel = c("Poisson", "LNP")[1]
Options_vec = c( "B_type"=switch(B_type, "Cointegration"=0, "Eigendecomposition"=1, "Eigen-cointegration"=2, "Eigen-cointegration_V2"=3), "IncludeAlpha"=1, "Cov_diagTF"=switch(as.character(IndependentTF), "TRUE"=1, "FALSE"=0), "StartVal"=switch(as.character(StartFromEquilibriumTF), "TRUE"=0, "FALSE"=1))
Nfactors_est = ifelse( IndependentTF==TRUE, 1, Nfactors)

# SimSettings
Nspecies = 4        # Only used in simulator
B_pp = cbind( c(0.7,0,-0.3,0), c(0,0.5,0,-0.5), c(-0.3,0,0.7,0), c(0,-0.5,0,0.5) )
SimSettings = list("n_species"=Nspecies, "n_years"=40, "n_years_burnin"=10, "start_from_equilibrium"=StartFromEquilibriumTF, "n_stations"=25, "n_samp_per_station"=8, "n_knots"=25, "Cov_pp"=0.05*diag(Nspecies), "B_pp"=B_pp, "phi_p"=rep(0,Nspecies), "SpatialScale"=0.1, "rho"=0.8, "SD_A"=0.01, "logMeanDens"=4, "SD_E"=0.01, "corr_E"=0, "ObsModel"="Poisson")

#########################
# Run model
#########################

### Load data
SimList = Sim_Fn( n_species=SimSettings[["n_species"]], n_years=SimSettings[["n_years"]], n_years_burnin=SimSettings[["n_years_burnin"]], n_stations=SimSettings[["n_stations"]], n_samp_per_station=SimSettings[["n_samp_per_station"]], n_knots=SimSettings[["n_knots"]], start_from_equilibrium=SimSettings[["start_from_equilibrium"]], Cov_pp=SimSettings[["Cov_pp"]], B_pp=SimSettings[["B_pp"]], phi_p=SimSettings[["phi_p"]], SpatialScale=SimSettings[["SpatialScale"]], rho=SimSettings[["rho"]], SD_A=SimSettings[["SD_A"]], logMeanDens=SimSettings[["logMeanDens"]], SD_E=SimSettings[["SD_E"]], corr_E=SimSettings[["corr_E"]], ObsModel=SimSettings[["ObsModel"]] )
DF = SimList$DF
Loc = SimList$Loc
a_x = rep(1,nrow(Loc))

# ObsModel
ObsModel_p = rep( switch(ObsModel,"Poisson"=0, "Lognormal"=1, "ZILN"=2, "LNP"=3), length(unique(DF[,'spp'])) )

# Make inputs
InputList = MakeInputs_Fn( Version=Version, options_vec=Options_vec, obsmodel_p=ObsModel_p, loc_x=Loc, a_x=a_x, data_frame=DF, n_factors=Nfactors_est, use_REML=Use_REML, independentTF=IndependentTF, estimate_phi=Estimate_Phi, n_cointegrate=Ncointegrate )

# Initialization
start_time = Sys.time()                                 # InputList[["TmbParams"]]  # SaveList$ParHat
dyn.load( dynlib(Version) )                                                         # log_tau=0.0,
obj = MakeADFun(data=InputList[["TmbData"]], DLL=Version, parameters=InputList[["TmbParams"]], random=InputList[["Random"]], map=InputList[["Map"]], hessian=FALSE, inner.control=list(maxit=1000) )

# First marginal likelihood and gradient
obj$fn( obj$par )
Gr = as.vector(obj$gr(obj$par))

# Bounds
EigenBounds = c("Lower"=-2, "Upper"=-0.001)
Upper = rep(Inf, length(obj$par) )
  if( B_type %in% c("Eigendecomposition") ) Upper[grep("Beta_pr",names(obj$par))] = EigenBounds['Upper']
  if( B_type %in% c("Eigen-cointegration") & Version%in%"spatial_vam_v11" ) Upper[grep("Alpha_pr",names(obj$par))[1+InputList$TmbData$n_p*0:InputList$TmbData$n_r]] = EigenBounds['Upper']
  if( B_type %in% c("Eigen-cointegration") & Version%in%"spatial_vam_v12" ) Upper[grep("Alpha_pr",names(obj$par))[1+(InputList$TmbData$n_p+1)*0:InputList$TmbData$n_r]] = EigenBounds['Upper']
  if( B_type %in% c("Eigen-cointegration_V2") ) Upper[grep("Beta_pr",names(obj$par))[length(grep("Beta_pr",names(obj$par)))-InputList$TmbData$n_r:1+1]] = EigenBounds['Upper']
  if( IndependentTF==TRUE ) Upper[grep("Alpha_pr",names(obj$par))] = EigenBounds['Upper']
Lower = rep(-Inf, length(obj$par) )
  if( B_type %in% c("Eigendecomposition") ) Lower[grep("Beta_pr",names(obj$par))] = EigenBounds['Lower']
  if( B_type%in%c("Eigen-cointegration") & Version%in%"spatial_vam_v11" ) Lower[grep("Alpha_pr",names(obj$par))[1+InputList$TmbData$n_p*0:InputList$TmbData$n_r]] = EigenBounds['Lower']
  if( B_type%in%c("Eigen-cointegration") & Version%in%"spatial_vam_v12" ) Lower[grep("Alpha_pr",names(obj$par))[1+(InputList$TmbData$n_p+1)*0:InputList$TmbData$n_r]] = EigenBounds['Lower']
  if( B_type %in% c("Eigen-cointegration_V2") ) Lower[grep("Beta_pr",names(obj$par))[length(grep("Beta_pr",names(obj$par)))-InputList$TmbData$n_r:1+1]] = EigenBounds['Lower']
  if( IndependentTF==TRUE ) Lower[grep("Alpha_pr",names(obj$par))] = EigenBounds['Lower']
  Lower[grep("logMargSigmaA_p",names(obj$par))] = log(0.001)

# Run model
for(i in 1:2) opt = nlminb(start=obj$env$last.par.best[-c(obj$env$random)], objective=obj$fn, gradient=obj$gr, upper=Upper, lower=Lower, control=list(eval.max=1e4, iter.max=1e4, trace=1) )
ParHat = obj$env$parList( opt$par )
opt[["final_gradient"]] = data.frame( "Name"=names(opt$par), "Lower"=Lower[1:length(opt$par)], "Est"=opt$par, "Upper"=Upper[1:length(opt$par)], "Gr"=obj$gr(opt$par))
opt[["AIC"]] = 2*length(opt$par) + 2*opt$objective

# Save stuff for now
Report = obj$report()
ParHat = obj$env$parList()
SD = sdreport(obj)

# Compute diagnostics
DerivedQuants = calc_derived_quants( Report=tryget("Report"), SD=tryget("SD"), InputList=tryget("InputList"), SimList=tryget("SimList") )

