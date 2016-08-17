
# Install package
if( !("SpatialVAM" %in% installed.packages()[,1])) devtools::install_github("james-thorson/SpatialVAM")  # , auth_token=[ask Jim Thorson for code]
if( !("FishData" %in% installed.packages()[,1])) devtools::install_github("james-thorson/FishData")
if( !("ThorsonUtilities" %in% installed.packages()[,1])) devtools::install_github("james-thorson/utilities")

# setwd( "C:/Users/James.Thorson/Desktop/Project_git/SpatialVAM/examples" )

# load libraries
library( TMB )
library( SpatialVAM )
TmbDir = system.file("executables", package="SpatialVAM")

# This is where all runs will be located
DateFile = paste(getwd(),'/',Sys.Date(),'_5species_EBS_Mesh/',sep='')
  dir.create(DateFile)

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
# v14 -- Added 2D AR1 option for spatial and spatio-temporal correlations

# Other settings
Nfactors_est = 5        # Number of dynamic factors in process error
Ncointegrate = 5
Use_REML = FALSE
Kmeans_Config = list("n_x"=25, "nstart"=100, "iter.max"=1000)
Estimate_Phi = TRUE      # Phi is the offset of initial and equilibrium abundance
IndependentTF = FALSE # if TRUE, then Cov_pp and B_pp are diagonal
StartFromEquilibriumTF = FALSE    # Definition for expected density at initial year: TRUE=Expectation is equilibrium; FALSE=Expectation is productivity without accounting for interactions/density dependence
B_type = c("Cointegration", "Eigendecomposition", "Eigen-cointegration", "Eigen-cointegration_V2")[4]
Kappa_Independent = TRUE # Only matters if IndependentTF==TRUE
ObsModel = c("Poisson", "LNP")[1]

# Derived
Options_vec = c( "B_type"=switch(B_type, "Cointegration"=0, "Eigendecomposition"=1, "Eigen-cointegration"=2, "Eigen-cointegration_V2"=3), "IncludeAlpha"=1, "Cov_diagTF"=switch(as.character(IndependentTF), "TRUE"=1, "FALSE"=0), "StartVal"=switch(as.character(StartFromEquilibriumTF), "TRUE"=0, "FALSE"=1))

# Determine region
Region = "Eastern_Bering_Sea"

# Decide on strata for use when calculating indices
strata.limits <- data.frame('STRATA'="All_areas")

# Decide on resolution
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )     # Samples: Do K-means on trawl locs; Domain: Do K-means on extrapolation grid
n_x = c(50, 100, 250, 500)[1] # Number of stations
grid_size_km = 50
Method = c("Grid", "Mesh")[2]

# Compile
file.copy( from=paste0(TmbDir,"/",Version,".cpp"), to=paste0(DateFile,"/",Version,".cpp") )
setwd( DateFile )
compile( paste0(Version,".cpp") )

#########################
# Run model
#########################

################
# Prepare data
# (THIS WILL VARY FOR DIFFERENT DATA SETS)
################

# Read or simulate trawl data
DF = FishData::download_catch_rates(survey="Eastern_Bering_Sea", species_set=5)
Data_Geostat = cbind( "spp"=DF[,"Sci"], "Year"=DF[,"Year"], "Catch_KG"=DF[,"Wt"], "AreaSwept_km2"=0.01, "Vessel"=0, "Lat"=DF[,"Lat"], "Lon"=DF[,"Long"] )

# Get extrapolation data
Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits )

# Calculate spatial information for SPDE mesh, strata areas, and AR1 process
Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn( grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile )
Data_Geostat = cbind( Data_Geostat, Spatial_List$loc_UTM, "knot_i"=Spatial_List$knot_i )

# Rename stuff
DF = ThorsonUtilities::rename_columns( Data_Geostat[,c('knot_i','spp','Year','Catch_KG')], newname=c('sitenum','spp','year','catch') )
Loc = Spatial_List[["loc_x"]]
a_x = Spatial_List[["a_xl"]][,1]

# ObsModel
ObsModel_p = rep( switch(ObsModel,"Poisson"=0, "Lognormal"=1, "ZILN"=2, "LNP"=3), length(unique(DF[,'spp'])) )

# Make inputs
InputList = MakeInputs_Fn( Version=Version, options_vec=Options_vec, obsmodel_p=ObsModel_p, loc_x=Loc, a_x=a_x, data_frame=DF, n_factors=Nfactors_est, use_REML=Use_REML, independentTF=IndependentTF, estimate_phi=Estimate_Phi, n_cointegrate=Ncointegrate )

# Initialization
dyn.load( dynlib(Version) )                                                         # log_tau=0.0,
obj = MakeADFun(data=InputList[["TmbData"]], DLL=Version, parameters=InputList[["TmbParams"]], random=InputList[["Random"]], map=InputList[["Map"]], hessian=FALSE, inner.control=list(maxit=1000) )

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
Opt = TMBhelper::Optimize( obj=obj, lower=Lower, upper=Upper, getsd=TRUE, savedir=DateDir )
Report = obj$report()
ParHat = obj$env$parList()

################
# Make diagnostic plots
################

# Plot index  # SpatialDeltaGLMM::
PlotIndex_Fn( DirName=DateDir, TmbData=InputList[["TmbData"]], Sdreport=Opt$SD, Year_Set=min(Data_Geostat[,'Year']):max(Data_Geostat[,'Year']), Years2Include=which(min(Data_Geostat[,'Year']):max(Data_Geostat[,'Year'])%in%sort(unique(Data_Geostat[,'Year']))), strata_names=names(strata.limits)[1], category_names=levels(Data_Geostat[,'spp']), use_biascorr=TRUE )

# Plot surface
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
Dim = c( "Nrow"=ceiling(sqrt(length(Years2Include))), "Ncol"=ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include)))) )
par( mfrow=Dim )
MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
PlotResultsOnMap_Fn(plot_set=3, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, category_names=levels(Data_Geostat[,'spp']), PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=paste0(DateDir,"Field_"), Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], mfrow=Dim, mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), Cex=MapDetails_List[["Cex"]], cex=1.8)

# Compute diagnostics
DerivedQuants = calc_derived_quants( Report=Report, SD=Opt$SD, InputList=InputList, species_names=levels(DF$spp) )

