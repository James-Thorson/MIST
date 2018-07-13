
#.onLoad <- function(libname, pkgname) {
#}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("###########################################################################################")
  packageStartupMessage("Loading package MIST, developed by James Thorson for the Northwest Fisheries Science Center")
  packageStartupMessage("###########################################################################################")
  if( !"INLA" %in% utils::installed.packages()[,1] ){
    packageStartupMessage("Installing INLA...")
    utils::install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")
  }
  #if( !"SpatialDeltaGLMM"%in%utils::installed.packages()[,1] ){
  #  packageStartupMessage("Installing package: SpatialDeltaGLMM...")
  #  devtools::install_github("nwfsc-assess/geostatistical_delta-GLMM")
  #}
  #if( !"TMB" %in% utils::installed.packages()[,1] ){
  #  packageStartupMessage("Installing TMB...")
  #  devtools::install_github("kaskr/adcomp/TMB")
  #}
  #if( !"TMBhelper" %in% utils::installed.packages()[,1] ){
  #  packageStartupMessage("Installing package: TMBhelper...")
  #  devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
  #}
}
