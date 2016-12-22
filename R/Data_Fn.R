
#' Build data input for SpatialVAM
#'
#' \code{Data_Fn} builds a tagged list of data inputs used by TMB for running the model
#'
#' @param Version a version number (see example for current default).
#' @param obsmodel_p The observation model for each species p
#' \describe{
#'   \item{ObsModel=0}{Poisson}
#'   \item{ObsModel=1}{Lognormal}
#'   \item{ObsModel=2}{Zero-inflated lognormal}
#'   \item{ObsModel=3}{lognormal-Poisson}
#'   \item{ObsModel=4}{Normal}
#' }
#' @param b_i Sampled biomass per unit area for each observation i
#' @param s_i Spatial knot (e.g., grid cell) for each observation i
#' @param t_i Time interval (e.g., year) for each observation i
#' @param p_i Species for each observation i
#' @param a_x Area associated with each knot
#' @param n_cointegrate Number of regulatory relations in community matrix (default is number of species)
#' @param n_factors Rank of covariance matrix for process error
#' @param B_type
#' \describe{
#'   \item{"Independent"}{Independent spatial-Gompertz dynamics for each species}
#'   \item{"Real_eigenvalue"}{Co-integration with eigenvalues restricted to real numbers}
#'   \item{"Complex_eigenvalue"}{Co-integration with eigenvalues including complex numbers}
#' }
#' @param startFromEquilibriumTF whether species start from equilibrium densities (i.e., turning of phi)
#' @param MeshList, tagged list representing location information for the SPDE mesh hyperdistribution, i.e., from \code{SpatialDeltaGLMM::Spatial_Information_Fn}
#' @param spatial_method DEPRECATED, always uses "Mesh" approximation
#' @param CheckForErrors Boolean, whether to check for errors in data inputs

#' @return Tagged list containing inputs to function \code{SpatialVAM::Build_TMB_Fn()}

#' @export
Data_Fn = function(Version, obsmodel_p=NULL,  b_i, s_i, t_i, p_i, a_x, n_cointegrate=NULL, n_factors=1,
  B_type="Real_eigenvalue", startFromEquilibriumTF=FALSE, MeshList, spatial_method=0, CheckForErrors=TRUE  ){

  # Assemble options vector
  options_vec = c( "B_type"=switch(B_type,"Independent"=0,"Real_eigenvalue"=3,"Complex_eigenvalue"=4), "IncludeAlpha"=1, "independentTF"=ifelse(B_type=="Independent",1,0), "StartVal"=1-as.numeric(startFromEquilibriumTF), "Spatial_Method"=spatial_method)

  # Expand a_x for auxiliary knots
  a_k = c(a_x, rep(0,MeshList$anisotropic_mesh$n-nrow(MeshList$loc_x)))

  # Check for bad data entry
  data_frame = cbind( "sitenum"=s_i, "year"=t_i, "spp"=p_i, "catch"=b_i )
  if( CheckForErrors==TRUE ){
    #if( !all(length(b_i)!=n_i | length(a_i)!=n_i | length(v_i)!=n_i | length(s_i)!=n_i | length(t_i)!=n_i ) stop("b_i, a_i, v_i, s_i, or t_i doesn't have length n_i")
  }

  # Fill in defaults
  if( is.null(n_cointegrate) || is.na(as.numeric(n_cointegrate)) ) n_cointegrate = length(unique(data_frame[,'spp']))
  if( n_cointegrate>length(unique(data_frame[,'spp'])) ) stop( "n_cointegrate can't be greater than the number of species")

  # Nonspatial
  if( Version%in%c("nonspatial_vam_v3","nonspatial_vam_v2","nonspatial_vam_v1")){
    # Pre-process to generate indices if using nonspatial
    Index_tp = tapply(DF$catch, INDEX=list(DF$year,DF$spp), FUN=mean)
    SE_Index_tp = tapply(DF$catch, INDEX=list(DF$year,DF$spp), FUN=function(vec){ sqrt(var(vec)/length(vec)) })
    SE_log_Index_tp = sqrt( log( (SE_Index_tp/Index_tp)^2 + 1))

    # Nonspatial data frame
    nonspatial_data_frame = cbind( "catch"=as.vector(Index_tp), "se_log_catch"=as.vector(SE_log_Index_tp), "year"=as.vector(row(Index_tp)), "spp"=as.vector(col(Index_tp)) )

    # Data
    if(Version%in%c("nonspatial_vam_v1")) Return = list("Options_vec"=options_vec, "n_i"=nrow(nonspatial_data_frame), "n_t"=length(unique(nonspatial_data_frame[,'year'])), "n_p"=length(unique(nonspatial_data_frame[,'spp'])), "n_j"=n_factors, "c_i"=nonspatial_data_frame[,'catch'], "se_log_c_i"=nonspatial_data_frame[,'se_log_catch'], "p_i"=as.numeric(nonspatial_data_frame[,'spp'])-1, "t_i"=nonspatial_data_frame[,'year']-1)
    if(Version%in%c("nonspatial_vam_v2")) Return = list("Options_vec"=options_vec, "ObsModel_p"=obsmodel_p, "n_i"=nrow(nonspatial_data_frame), "n_t"=length(unique(nonspatial_data_frame[,'year'])), "n_p"=length(unique(nonspatial_data_frame[,'spp'])), "n_j"=n_factors, "c_i"=nonspatial_data_frame[,'catch'], "se_log_c_i"=nonspatial_data_frame[,'se_log_catch'], "p_i"=as.numeric(nonspatial_data_frame[,'spp'])-1, "t_i"=nonspatial_data_frame[,'year']-1)
    if(Version%in%c("nonspatial_vam_v3")) Return = list("Options_vec"=options_vec, "ObsModel_p"=obsmodel_p, "n_i"=nrow(nonspatial_data_frame), "n_t"=length(unique(nonspatial_data_frame[,'year'])), "n_p"=length(unique(nonspatial_data_frame[,'spp'])), "n_r"=n_cointegrate, "n_j"=n_factors, "c_i"=nonspatial_data_frame[,'catch'], "se_log_c_i"=nonspatial_data_frame[,'se_log_catch'], "p_i"=as.numeric(nonspatial_data_frame[,'spp'])-1, "t_i"=nonspatial_data_frame[,'year']-1)
    #Data[['Options_vec']]['ObsModel'] = 1  # Switch to lognormal distribution
    if( "ObsModel_p" %in% names(Return)) Return[["ObsModel_p"]][] = 1
  } # End nonspatial

  # Spatial
  if( Version%in%c("spatial_vam_v14","spatial_vam_v13","spatial_vam_v12","spatial_vam_v11","spatial_vam_v10","spatial_vam_v9","spatial_vam_v8","spatial_vam_v7","spatial_vam_v6","spatial_vam_v5","spatial_vam_v4","spatial_vam_v3","spatial_vam_v2","spatial_vam_v1")){
    require( Matrix )
    # Sparse matrices for 2D AR1 process
    # Generate sparse matrices for precision matrix of 2D AR1 process
    M0 = Spatial_List$GridList[["M0"]]
    M1 = Spatial_List$GridList[["M1"]]
    M2 = Spatial_List$GridList[["M2"]]

    # Data
    # Necessary columns: sitenum, year, catch, spp
    if(Version%in%c("spatial_vam_v2","spatial_vam_v1")) Return = list("Options_vec"=options_vec, "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=diff(range(data_frame[,'year']))+1, "n_k"=MeshList$anisotropic_mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-1, "spde"=NULL)
    if(Version%in%c("spatial_vam_v5","spatial_vam_v4","spatial_vam_v3")) Return = list("Options_vec"=options_vec, "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=diff(range(data_frame[,'year']))+1, "n_k"=MeshList$anisotropic_mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-min(data_frame[,'year']), "spde"=NULL, "M0"=M0, "M1"=M1, "M2"=M2)
    if(Version%in%c("spatial_vam_v6")) Return = list("Options_vec"=options_vec, "ObsModel_p"=obsmodel_p, "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=diff(range(data_frame[,'year']))+1, "n_k"=MeshList$anisotropic_mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-min(data_frame[,'year']), "spde"=NULL, "M0"=M0, "M1"=M1, "M2"=M2)
    if(Version%in%c("spatial_vam_v7")) Return = list("Options_vec"=options_vec, "ObsModel_p"=obsmodel_p, "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=diff(range(data_frame[,'year']))+1, "n_k"=MeshList$anisotropic_mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-min(data_frame[,'year']), "a_k"=a_k, "Z_kl"=rbind(MeshList$loc_x,matrix(0,ncol=2,nrow=MeshList$anisotropic_mesh$n-nrow(MeshList$loc_x))), "spde"=NULL, "M0"=M0, "M1"=M1, "M2"=M2)
    if(Version%in%c("spatial_vam_v8")) Return = list("Options_vec"=options_vec, "ObsModel_p"=obsmodel_p, "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=diff(range(data_frame[,'year']))+1, "n_k"=MeshList$anisotropic_mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_r"=n_cointegrate, "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-min(data_frame[,'year']), "a_k"=a_k, "Z_kl"=rbind(MeshList$loc_x,matrix(0,ncol=2,nrow=MeshList$anisotropic_mesh$n-nrow(MeshList$loc_x))), "spde"=NULL, "M0"=M0, "M1"=M1, "M2"=M2)
    if(Version%in%c("spatial_vam_v13","spatial_vam_v12","spatial_vam_v11","spatial_vam_v10","spatial_vam_v9")) Return = list("Options_vec"=options_vec, "ObsModel_p"=obsmodel_p, "PenMult_z"=c(1000,10), "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=diff(range(data_frame[,'year']))+1, "n_k"=MeshList$anisotropic_mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_r"=n_cointegrate, "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-min(data_frame[,'year']), "a_k"=a_k, "Z_kl"=rbind(as.matrix(MeshList$loc_x),matrix(0,ncol=2,nrow=MeshList$anisotropic_mesh$n-nrow(MeshList$loc_x))), "spde"=NULL, "G0"=M0, "G1"=M1, "G2"=M2)
    if(Version%in%c("spatial_vam_v14")) Return = list("Options_vec"=options_vec, "ObsModel_p"=obsmodel_p, "PenMult_z"=c(1000,10), "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=diff(range(data_frame[,'year']))+1, "n_k"=MeshList$anisotropic_mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_r"=n_cointegrate, "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-min(data_frame[,'year']), "a_k"=a_k, "Z_kl"=rbind(as.matrix(MeshList$loc_x),matrix(0,ncol=2,nrow=MeshList$anisotropic_mesh$n-nrow(MeshList$loc_x))), "spde"=NULL, "M0"=M0, "M1"=M1, "M2"=M2)
    if("spde" %in% names(Return)) Return[['spde']] = list("n_s"=MeshList$anisotropic_spde$n.spde, "n_tri"=nrow(MeshList$anisotropic_mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$anisotropic_spde$param.inla$M0, "G0_inv"=INLA::inla.as.dgTMatrix(solve(MeshList$anisotropic_spde$param.inla$M0)) )

    # Changes necessary for 2D AR1 process
    if( length(options_vec)>4 && options_vec[5]==1 ){
      #Return[["n_k"]] = Return[["n_s"]]
    }
  } # End spatial

  return(Return)
}

