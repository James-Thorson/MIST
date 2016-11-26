
# Install package
if( !("SpatialVAM" %in% installed.packages()[,1])) devtools::install_github("james-thorson/SpatialVAM")  # , auth_token=[ask Jim Thorson for code]
if( !("FishData" %in% installed.packages()[,1])) devtools::install_github("james-thorson/FishData")
if( !("ThorsonUtilities" %in% installed.packages()[,1])) devtools::install_github("james-thorson/utilities")

# setwd( "C:/Users/James.Thorson/Desktop/Project_git/SpatialVAM/examples" )

# load libraries
library( TMB )
library( SpatialVAM )
#TmbDir = system.file("executables", package="SpatialVAM")
TmbDir = "C:/Users/James.Thorson/Desktop/Project_git/SpatialVAM/inst/executables/"

# Source new files
FileSet = c("Build_TMB_Fn.R", "Data_Fn.R", "Summarize.R")
for( i in 1:length(FileSet)) source( paste0(TmbDir,"../../R/",FileSet[i]) )

# This is where all runs will be located
DateFile = paste(getwd(),'/',Sys.Date(),'_5species_EBS_Mesh/',sep='')
  dir.create(DateFile)

##############
# Settings
##############

Version = "spatial_vam_v14"
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
Estimate_Phi = TRUE      # Phi is the offset of initial and equilibrium abundance
StartFromEquilibriumTF = FALSE    # Definition for expected density at initial year: TRUE=Expectation is equilibrium; FALSE=Expectation is productivity without accounting for interactions/density dependence
B_type = c("Independent", "Real_eigenvalue", "Eigen-Complex_eigenvalue")[2]
Kappa = c("constant", "spatial_vs_spatiotemporal", "different")[1]  # Only works for new Data_Fn
EigenBounds = c("Lower"=-2, "Upper"=-0.001)
ObsModel = c("Poisson", "LNP", "ZILN")[3]

# Determine region
Region = "Eastern_Bering_Sea"

# Decide on strata for use when calculating indices
strata.limits <- data.frame('STRATA'="All_areas")

# Decide on resolution
Kmeans_Config = list("n_x"=25, "nstart"=100, "iter.max"=1000)
n_x = c(50, 100, 250, 500)[1] # Number of stations
grid_size_km = 50

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
Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn( grid_size_km=grid_size_km, n_x=n_x, Method="Mesh", Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile )
Data_Geostat = cbind( Data_Geostat, Spatial_List$loc_UTM, "knot_i"=Spatial_List$knot_i )

# ObsModel
ObsModel_p = rep( switch(ObsModel,"Poisson"=0, "Lognormal"=1, "ZILN"=2, "LNP"=3), length(unique(Data_Geostat[,'spp'])) )

# Make inputs
TmbData = Data_Fn( "Version"=Version, "obsmodel_p"=ObsModel_p, "n_cointegrate"=Ncointegrate, "b_i"=Data_Geostat[,'Catch_KG'], "s_i"=Data_Geostat[,'knot_i'], "t_i"=Data_Geostat[,'Year'], "p_i"=Data_Geostat[,'spp'], "a_x"=Spatial_List$a_xl[,1], "B_type"=B_type, "startFromEquilibriumTF"=FALSE, "spatial_method"=0, "MeshList"=Spatial_List$MeshList, "n_factors"=Nfactors_est )

# Initialization
TmbList = Build_TMB_Fn( "TmbData"=TmbData, "Version"=Version, "use_REML"=ifelse(is.na(Use_REML),TRUE,Use_REML), "loc_x"=Spatial_List$MeshList$loc_x, "estimate_phi"=Estimate_Phi, "Kappa"=Kappa, "eigenbounds"=EigenBounds, "TmbDir"=TmbDir, "RunDir"=DateFile )
obj = TmbList$Obj                            # "Parameters"=InputList$TmbParams,

# Run model
Opt = TMBhelper::Optimize( obj=obj, lower=TmbList$Lower, upper=TmbList$Upper, getsd=TRUE, savedir=DateFile )
Report = obj$report()
ParHat = obj$env$parList()

################
# Make diagnostic plots
################

# Settings
Year_Set = min(DF[,'year']):max(DF[,'year'])

# Plot index  # SpatialDeltaGLMM::
SpatialDeltaGLMM::PlotIndex_Fn( DirName=DateDir, TmbData=InputList[["TmbData"]], Sdreport=Opt$SD, Year_Set=min(Data_Geostat[,'Year']):max(Data_Geostat[,'Year']), Years2Include=which(min(Data_Geostat[,'Year']):max(Data_Geostat[,'Year'])%in%sort(unique(Data_Geostat[,'Year']))), strata_names=names(strata.limits)[1], category_names=levels(Data_Geostat[,'spp']), use_biascorr=TRUE )

# Get region-specific settings for plots
MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"="Gulf_of_St_Lawrence", "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
# Plot maps representing density or other variables
SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=ConfigDir, category_names=unique(DF$spp), Year_Set=Year_Set, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)

# Center of gravity
SpatialDeltaGLMM::Plot_range_shifts(Report=Report, TmbData=TmbData, Sdreport=SD, Znames=colnames(TmbData$Z_xm), PlotDir=DateDir, category_names=unique(DF$spp), Year_Set=Year_Set)

# Compute diagnostics
DerivedQuants = Summarize( Report=Report, SD=Opt$SD, InputList=InputList, species_names=levels(DF$spp) )

