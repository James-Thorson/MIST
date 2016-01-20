
MakeInputs_Fn = function( Version, options_vec, obsmodel_p=NULL, loc_x, a_x, data_frame, n_cointegrate=NULL, n_factors=1, use_REML=FALSE, independentTF=FALSE, estimate_phi=TRUE ){

  # 
  rmatrix = function(nrow, ncol, mean=0, sd=1, diag=NA){
    Return = matrix(rnorm(nrow*ncol,mean=mean,sd=sd), nrow=nrow, ncol=ncol)
    if( !is.na(diag)) Return[cbind(1:min(nrow,ncol),1:min(nrow,ncol))] = diag
    return(Return)
  }
  seq_pos = function( n ) seq(from=1, to=n, length=n)
  
  # Build SPDE object using INLA
  MeshList = SpatialDeltaGLMM::Calc_Anisotropic_Mesh(loc_x=loc_x)
  a_k = c(a_x, rep(0,MeshList$mesh$n-nrow(loc_x)))

  # Check inputs in DF
  if( !all( c("sitenum","year","spp","catch") %in% colnames(data_frame)) ) stop( "data_frame must contain columns labeled sitenum, year, spp, and catch")

  # Fill in defaults
  if( is.null(n_cointegrate) | is.na(as.numeric(n_cointegrate)) ) n_cointegrate = length(unique(data_frame[,'spp']))
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
    if(Version%in%c("nonspatial_vam_v1")) Data = list("Options_vec"=options_vec, "n_i"=nrow(nonspatial_data_frame), "n_t"=length(unique(nonspatial_data_frame[,'year'])), "n_p"=length(unique(nonspatial_data_frame[,'spp'])), "n_j"=n_factors, "c_i"=nonspatial_data_frame[,'catch'], "se_log_c_i"=nonspatial_data_frame[,'se_log_catch'], "p_i"=as.numeric(nonspatial_data_frame[,'spp'])-1, "t_i"=nonspatial_data_frame[,'year']-1)
    if(Version%in%c("nonspatial_vam_v2")) Data = list("Options_vec"=options_vec, "ObsModel_p"=obsmodel_p, "n_i"=nrow(nonspatial_data_frame), "n_t"=length(unique(nonspatial_data_frame[,'year'])), "n_p"=length(unique(nonspatial_data_frame[,'spp'])), "n_j"=n_factors, "c_i"=nonspatial_data_frame[,'catch'], "se_log_c_i"=nonspatial_data_frame[,'se_log_catch'], "p_i"=as.numeric(nonspatial_data_frame[,'spp'])-1, "t_i"=nonspatial_data_frame[,'year']-1)
    if(Version%in%c("nonspatial_vam_v3")) Data = list("Options_vec"=options_vec, "ObsModel_p"=obsmodel_p, "n_i"=nrow(nonspatial_data_frame), "n_t"=length(unique(nonspatial_data_frame[,'year'])), "n_p"=length(unique(nonspatial_data_frame[,'spp'])), "n_r"=n_cointegrate, "n_j"=n_factors, "c_i"=nonspatial_data_frame[,'catch'], "se_log_c_i"=nonspatial_data_frame[,'se_log_catch'], "p_i"=as.numeric(nonspatial_data_frame[,'spp'])-1, "t_i"=nonspatial_data_frame[,'year']-1)
    #Data[['Options_vec']]['ObsModel'] = 1  # Switch to lognormal distribution
    if( "ObsModel_p" %in% names(Data)) Data[["ObsModel_p"]][] = 1

    # Parameters
    if(Version%in%c("nonspatial_vam_v2","nonspatial_vam_v1")) Params = list("alpha_p"=rep(0,Data$n_p), "phi_p"=rep(0,Data$n_p), "L_val"=rnorm(Data$n_j*Data$n_p-Data$n_j*(Data$n_j-1)/2), "B_pp"=matrix(rnorm(Data$n_p^2,sd=0.1),Data$n_p,Data$n_p), "logsigma_pz"=matrix(0,nrow=Data$n_p,2), "d_tp"=array(2,dim=unlist(Data[c('n_t','n_p')])))
    if(Version%in%c("nonspatial_vam_v3")) Params = list("alpha_p"=rep(0,Data$n_p), "phi_p"=rep(0,Data$n_p), "L_val"=rnorm(Data$n_j*Data$n_p-Data$n_j*(Data$n_j-1)/2), "Alpha_pr"=rbind(diag(-0.5,Data$n_r),rmatrix(nrow=Data$n_p-Data$n_r,ncol=Data$n_r,sd=0.01)), "Beta_pr"=rbind(diag(Data$n_r),rmatrix(nrow=Data$n_p-Data$n_r,ncol=Data$n_r,sd=0.2)), "logsigma_pz"=matrix(0,nrow=Data$n_p,2), "d_tp"=array(2,dim=unlist(Data[c('n_t','n_p')])))
    #

    # Random
    Random = c( "d_tp" )
    if(use_REML==TRUE) Random = c(Random, "alpha_p", "phi_p")  # , "B_pp"

    # Map
    Map = NULL
    # fix phi
    if( estimate_phi==FALSE ){
      Map[["phi_p"]] = factor( rep(NA,length(Params[["phi_p"]])) )
      Params[["phi_p"]][] = 0
    }
    # Fix Alpha_pr and Beta_pr for the eigendecomposition method
    if( Options_vec[["B_type"]]==1 ){
      # Alpha_pr, fix first row at one, so that the magnitude of each eigenvector is not colinear with the eigenvalues
      Map[["Alpha_pr"]] = array(1:prod(dim(Params[["Beta_pr"]])), dim=dim(Params[["Beta_pr"]]))
      Map[["Alpha_pr"]][1,] = NA
      Map[["Alpha_pr"]] = factor(Map[["Alpha_pr"]])
      Params[["Alpha_pr"]] = array( ifelse(!is.na(Map[["Alpha_pr"]]),Params[["Alpha_pr"]],1), dim=dim(Params[["Alpha_pr"]]))
      # Beta_pr, just estimate row column, which is interpreted as eigenvalues
      Map[["Beta_pr"]] = array(NA, dim=dim(Params[["Beta_pr"]]))
      Map[["Beta_pr"]][1,] = 1:ncol(Map[["Beta_pr"]])
      Map[["Beta_pr"]] = factor(Map[["Beta_pr"]])
      Params[["Beta_pr"]] = array( ifelse(!is.na(Map[["Beta_pr"]]),-0.5,0), dim=dim(Params[["Beta_pr"]]))  # B = U%*%L%*%solve(U) + I_pp, so -1<=eigenvalues<=0
    }
    # Identifiability restrictions on Beta_pr for co-integration method
    if( "Beta_pr"%in%names(Params) && Options_vec[["B_type"]]%in%c(0,2) ){
      Map[["Beta_pr"]] = factor(rbind( matrix(NA,nrow=Data$n_r,ncol=Data$n_r), matrix(seq_pos(Data$n_r*(Data$n_p-Data$n_r)),nrow=Data$n_p-Data$n_r,ncol=Data$n_r)))
    }
    # Better initial conditions for eigen-cointegration
    if( "Alpha_pr"%in%names(Params) && Options_vec[["B_type"]]%in%c(2) ){
      Params[["Alpha_pr"]][1,] = -0.5
    }
    # Make B_pp diagonal
    if( independentTF==TRUE ){
      if( "B_pp" %in% names(Params) ){
        Map[["B_pp"]] = matrix(NA, nrow=Data$n_p, ncol=Data$n_p)
        diag(Map[["B_pp"]]) = 1:Data$n_p
        Map[["B_pp"]] = factor( Map[["B_pp"]] )
        Params[["B_pp"]][] = 0
        Params[["L_val"]][] = 0
      }
      if( "Alpha_pr"%in%names(Params) && n_cointegrate==Data$n_p ){
        Map[["Alpha_pr"]] = matrix(NA, nrow=Data$n_p, ncol=Data$n_p)
        diag(Map[["Alpha_pr"]]) = 1:Data$n_p
        Map[["Alpha_pr"]] = factor( Map[["Alpha_pr"]] )
        Params[["Alpha_pr"]][] = 0
        Params[["L_val"]][] = 0
      }
    }
    # Observation model
    Map[["logsigma_pz"]] = matrix( 1:(2*Data$n_p), ncol=2, byrow=TRUE )
    for(p in 1:Data$n_p){
      if( Data$ObsModel_p[p]==0 ){
        Map[["logsigma_pz"]][p,] = c(NA,NA)
      }
      if( Data$ObsModel_p[p] %in% c(1,3,4) ){
        Map[["logsigma_pz"]][p,2] = NA
      }
    }
    Map[["logsigma_pz"]] = factor( Map[["logsigma_pz"]] )
  } # End nonspatial

  # Spatial
  if( Version%in%c("spatial_vam_v13","spatial_vam_v12","spatial_vam_v11","spatial_vam_v10","spatial_vam_v9","spatial_vam_v8","spatial_vam_v7","spatial_vam_v6","spatial_vam_v5","spatial_vam_v4","spatial_vam_v3","spatial_vam_v2","spatial_vam_v1")){
    # Data
    # Necessary columns: sitenum, year, catch, spp
    if(Version%in%c("spatial_vam_v2","spatial_vam_v1")) Data = list("Options_vec"=options_vec, "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=diff(range(data_frame[,'year']))+1, "n_k"=MeshList$mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-1, "spde"=NULL)
    if(Version%in%c("spatial_vam_v5","spatial_vam_v4","spatial_vam_v3")) Data = list("Options_vec"=options_vec, "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=diff(range(data_frame[,'year']))+1, "n_k"=MeshList$mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-min(data_frame[,'year']), "spde"=NULL, "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2)
    if(Version%in%c("spatial_vam_v6")) Data = list("Options_vec"=options_vec, "ObsModel_p"=obsmodel_p, "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=diff(range(data_frame[,'year']))+1, "n_k"=MeshList$mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-min(data_frame[,'year']), "spde"=NULL, "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2)
    if(Version%in%c("spatial_vam_v7")) Data = list("Options_vec"=options_vec, "ObsModel_p"=obsmodel_p, "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=diff(range(data_frame[,'year']))+1, "n_k"=MeshList$mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-min(data_frame[,'year']), "a_k"=a_k, "Z_kl"=rbind(loc_x,matrix(0,ncol=2,nrow=MeshList$mesh$n-nrow(loc_x))), "spde"=NULL, "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2)
    if(Version%in%c("spatial_vam_v8")) Data = list("Options_vec"=options_vec, "ObsModel_p"=obsmodel_p, "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=diff(range(data_frame[,'year']))+1, "n_k"=MeshList$mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_r"=n_cointegrate, "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-min(data_frame[,'year']), "a_k"=a_k, "Z_kl"=rbind(loc_x,matrix(0,ncol=2,nrow=MeshList$mesh$n-nrow(loc_x))), "spde"=NULL, "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2)
    if(Version%in%c("spatial_vam_v13","spatial_vam_v12","spatial_vam_v11","spatial_vam_v10","spatial_vam_v9")) Data = list("Options_vec"=options_vec, "ObsModel_p"=obsmodel_p, "PenMult_z"=c(0,10), "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=diff(range(data_frame[,'year']))+1, "n_k"=MeshList$mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_r"=n_cointegrate, "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-min(data_frame[,'year']), "a_k"=a_k, "Z_kl"=rbind(loc_x,matrix(0,ncol=2,nrow=MeshList$mesh$n-nrow(loc_x))), "spde"=NULL, "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2)
    if("spde" %in% names(Data)) Data[['spde']] = list("n_s"=MeshList$spde$n.spde, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$spde$param.inla$M0, "G0_inv"=inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)) )

    # Parameters
    if(Version=="spatial_vam_v1") Params = list("Hinput_z"=c(0,0), "logkappa"=log(1), "alpha_p"=rep(0,Data$n_p), "phi_p"=rep(0,Data$n_p), "logtauA_p"=rep(0,Data$n_p), "L_val"=rnorm(Data$n_j*Data$n_p-Data$n_j*(Data$n_j-1)/2), "B_pp"=matrix(rnorm(Data$n_p^2,sd=0.1),Data$n_p,Data$n_p), "d_ktp"=array(2,dim=unlist(Data[c('n_k','n_t','n_p')])), "Ainput_kp"=matrix(0,nrow=Data$n_k,ncol=Data$n_p))
    if(Version%in%c("spatial_vam_v3","spatial_vam_v2")) Params = list("Hinput_z"=c(0,0), "logkappa"=log(1), "alpha_p"=rep(0,Data$n_p), "phi_p"=rep(0,Data$n_p), "logtauA_p"=rep(0,Data$n_p), "L_val"=ifelse( is.na(fixdiag(Nrow=Data$n_p, Ncol=Data$n_j)), 1, 0), "B_pp"=diag(0.5,Data$n_p), "logsigma_pz"=matrix(0,nrow=Data$n_p,1), "d_ktp"=abind(SimList$d_stp,array(0,dim=c(Data$n_k-Data$n_s,Data$n_t,Data$n_p)),along=1), "Ainput_kp"=matrix(0,nrow=Data$n_k,ncol=Data$n_p))      # "d_ktp"=array(2,dim=unlist(Data[c('n_k','n_t','n_p')]))
    if(Version%in%c("spatial_vam_v4")) Params = list("Hinput_z"=c(0,0), "logkappa"=log(1), "alpha_p"=rep(0,Data$n_p), "phi_p"=rep(0,Data$n_p), "logtauA_p"=rep(0,Data$n_p), "L_val"=ifelse( is.na(fixdiag(Nrow=Data$n_p, Ncol=Data$n_j)), 1, 0), "B_pp"=diag(0.5,Data$n_p), "logsigma_pz"=matrix(0,nrow=Data$n_p,2), "d_ktp"=array(0,dim=c(Data$n_k,Data$n_t,Data$n_p)), "Ainput_kp"=matrix(0,nrow=Data$n_k,ncol=Data$n_p))      # "d_ktp"=array(2,dim=unlist(Data[c('n_k','n_t','n_p')]))
    if(Version%in%c("spatial_vam_v7","spatial_vam_v6","spatial_vam_v5")) Params = list("Hinput_z"=c(0,0), "logkappa"=log(1), "alpha_p"=rep(0,Data$n_p), "phi_p"=rep(0,Data$n_p), "logtauA_p"=rep(0,Data$n_p), "L_val"=ifelse( is.na(fixdiag(Nrow=Data$n_p, Ncol=Data$n_j)), 1, 0), "B_pp"=diag(0.5,Data$n_p), "logsigma_pz"=matrix(0,nrow=Data$n_p,2), "d_ktp"=array(0,dim=c(Data$n_k,Data$n_t,Data$n_p)), "Ainput_kp"=matrix(0,nrow=Data$n_k,ncol=Data$n_p), "delta_i"=rep(0,Data$n_i))      # "d_ktp"=array(2,dim=unlist(Data[c('n_k','n_t','n_p')]))
    if(Version%in%c("spatial_vam_v9","spatial_vam_v8")) Params = list("Hinput_z"=c(0,0), "logkappa"=log(1), "alpha_p"=rep(0,Data$n_p), "phi_p"=rep(0,Data$n_p), "logtauA_p"=rep(0,Data$n_p), "L_val"=ifelse( is.na(fixdiag(Nrow=Data$n_p, Ncol=Data$n_j)), 1, 0), "Alpha_pr"=rbind(diag(-0.5,Data$n_r),rmatrix(nrow=Data$n_p-Data$n_r,ncol=Data$n_r,sd=0.01)), "Beta_pr"=rbind(diag(Data$n_r),rmatrix(nrow=Data$n_p-Data$n_r,ncol=Data$n_r,sd=0.2)), "logsigma_pz"=matrix(0,nrow=Data$n_p,2), "d_ktp"=array(0,dim=c(Data$n_k,Data$n_t,Data$n_p)), "Ainput_kp"=matrix(0,nrow=Data$n_k,ncol=Data$n_p), "delta_i"=rep(0,Data$n_i))      # "d_ktp"=array(2,dim=unlist(Data[c('n_k','n_t','n_p')]))
    if(Version%in%c("spatial_vam_v10")) Params = list("Hinput_z"=c(0,0), "logkappa"=log(1), "alpha_p"=rep(0,Data$n_p), "phi_p"=rep(0,Data$n_p), "logMargSigmaA_p"=rep(0,Data$n_p), "L_val"=ifelse( is.na(fixdiag(Nrow=Data$n_p, Ncol=Data$n_j)), 1, 0), "Alpha_pr"=rbind(diag(-0.5,Data$n_r),rmatrix(nrow=Data$n_p-Data$n_r,ncol=Data$n_r,sd=0.01)), "Beta_pr"=rbind(diag(Data$n_r),rmatrix(nrow=Data$n_p-Data$n_r,ncol=Data$n_r,sd=0.2)), "logsigma_pz"=matrix(0,nrow=Data$n_p,2), "d_ktp"=array(0,dim=c(Data$n_k,Data$n_t,Data$n_p)), "Ainput_kp"=matrix(0,nrow=Data$n_k,ncol=Data$n_p), "delta_i"=rep(0,Data$n_i))      # "d_ktp"=array(2,dim=unlist(Data[c('n_k','n_t','n_p')]))
    if(Version%in%c("spatial_vam_v12","spatial_vam_v11")) Params = list("Hinput_z"=c(0,0), "logkappa_z"=rep(0,Data$n_p+1), "alpha_p"=rep(0,Data$n_p), "phi_p"=rep(0,Data$n_p), "logMargSigmaA_p"=rep(0,Data$n_p), "L_val"=ifelse( is.na(fixdiag(Nrow=Data$n_p, Ncol=Data$n_j)), 1, 0), "Alpha_pr"=rbind(rmatrix(nrow=Data$n_r,ncol=Data$n_r,sd=0.01,diag=-0.5),rmatrix(nrow=Data$n_p-Data$n_r,ncol=Data$n_r,sd=0.01)), "Beta_pr"=rbind(diag(Data$n_r),rmatrix(nrow=Data$n_p-Data$n_r,ncol=Data$n_r,sd=0.2)), "logsigma_pz"=matrix(0,nrow=Data$n_p,2), "d_ktp"=array(0,dim=c(Data$n_k,Data$n_t,Data$n_p)), "Ainput_kp"=matrix(0,nrow=Data$n_k,ncol=Data$n_p), "delta_i"=rep(0,Data$n_i))      # "d_ktp"=array(2,dim=unlist(Data[c('n_k','n_t','n_p')]))
    if(Version%in%c("spatial_vam_v13")) Params = list("Hinput_z"=c(0,0), "logkappa_z"=rep(0,Data$n_p+1), "alpha_p"=rep(0,Data$n_p), "phi_p"=rep(0,Data$n_p), "logMargSigmaA_p"=rep(0,Data$n_p), "L_val"=ifelse( is.na(fixdiag(Nrow=Data$n_p, Ncol=Data$n_j)), 1, 0), "Alpha_pr"=rbind(rmatrix(nrow=Data$n_r,ncol=Data$n_r,sd=0.01,diag=-0.5),rmatrix(nrow=Data$n_p-Data$n_r,ncol=Data$n_r,sd=0.01)), "Beta_pr"=rbind(diag(Data$n_r),rmatrix(nrow=Data$n_p-Data$n_r,ncol=Data$n_r,sd=0.2)), "logsigma_pz"=matrix(0,nrow=Data$n_p,2), "d_ktp"=array(0,dim=c(Data$n_k,Data$n_t,Data$n_p)), "Ainput_kp"=matrix(0,nrow=Data$n_k,ncol=Data$n_p), "delta_i"=rep(0,Data$n_i))      # "d_ktp"=array(2,dim=unlist(Data[c('n_k','n_t','n_p')]))

    # Random
    # Treating alpha_p, phi_p and B_pp as random (in REML) results in very slow inner optimization!  (100s of steps)
    Random = c( "Ainput_kp", "d_ktp" )
    if( Version%in%c("spatial_vam_v13","spatial_vam_v12","spatial_vam_v11","spatial_vam_v10","spatial_vam_v9","spatial_vam_v8","spatial_vam_v7","spatial_vam_v6","spatial_vam_v5")) Random = c(Random, "delta_i")
    if(use_REML==TRUE) Random = c(Random, "alpha_p", "phi_p")  # , "B_pp"

    # Map
    Map = NULL
    # Anisotropy
    Map[["Hinput_z"]] = factor( rep(NA,2) )
    # Observation model
    Map[["logsigma_pz"]] = matrix( 1:(2*Data$n_p), ncol=2, byrow=TRUE )
    Map[["delta_i"]] = 1:Data$n_i
    for(p in 1:Data$n_p){
      if( Data$ObsModel_p[p]==0 ){
        Map[["logsigma_pz"]][p,] = c(NA,NA)
        if( Version%in%c("spatial_vam_v13","spatial_vam_v12","spatial_vam_v11","spatial_vam_v10","spatial_vam_v9","spatial_vam_v8","spatial_vam_v7","spatial_vam_v6","spatial_vam_v5") ) Map[["delta_i"]][which((Data$p_i+1)==p)] = rep(NA,sum((Data$p_i+1)==p))
      }
      if( Data$ObsModel_p[p]==1 ){
        Map[["logsigma_pz"]][p,2] = NA
        if( Version%in%c("spatial_vam_v13","spatial_vam_v12","spatial_vam_v11","spatial_vam_v10","spatial_vam_v9","spatial_vam_v8","spatial_vam_v7","spatial_vam_v6","spatial_vam_v5") ) Map[["delta_i"]][which((Data$p_i+1)==p)] = rep(NA,sum((Data$p_i+1)==p))
      }
      if( Data$ObsModel_p[p]==2 ){
        if( Version%in%c("spatial_vam_v13","spatial_vam_v12","spatial_vam_v11","spatial_vam_v10","spatial_vam_v9","spatial_vam_v8","spatial_vam_v7","spatial_vam_v6","spatial_vam_v5") ) Map[["delta_i"]][which((Data$p_i+1)==p)] = rep(NA,sum((Data$p_i+1)==p))
        # Check number of zeros
        NumZero = tapply( Data$c_i, INDEX=Data$p_i, FUN=function(vec){sum(vec==0)})
        if( any(NumZero==0) ){
          Map[["logsigma_pz"]][,2] = ifelse( NumZero==0, NA, Map[["logsigma_pz"]][,2])
          Params[["logsigma_pz"]][,2] = ifelse( NumZero==0, 20, Params[["logsigma_pz"]][,2])
        }
      }
      if( Data$ObsModel_p[p]==3 ){
        Map[["logsigma_pz"]][p,2] = NA
      }
      if( Data$ObsModel_p[p]==4 ){
        Map[["logsigma_pz"]][p,2] = NA
        if( Version%in%c("spatial_vam_v13","spatial_vam_v12","spatial_vam_v11","spatial_vam_v10","spatial_vam_v9","spatial_vam_v8","spatial_vam_v7","spatial_vam_v6","spatial_vam_v5") ) Map[["delta_i"]][which((Data$p_i+1)==p)] = rep(NA,sum((Data$p_i+1)==p))
      }
    }
    Map[["logsigma_pz"]] = factor( Map[["logsigma_pz"]] )
    Map[["delta_i"]] = factor( Map[["delta_i"]] )
    # Fix Alpha_pr and Beta_pr for the eigendecomposition method
    if( Options_vec[["B_type"]]==1 ){
      # Alpha_pr, fix first row at one, so that the magnitude of each eigenvector is not colinear with the eigenvalues
      Map[["Alpha_pr"]] = array(1:prod(dim(Params[["Beta_pr"]])), dim=dim(Params[["Beta_pr"]]))
      Map[["Alpha_pr"]][1,] = NA
      Map[["Alpha_pr"]] = factor(Map[["Alpha_pr"]])
      Params[["Alpha_pr"]] = array( ifelse(!is.na(Map[["Alpha_pr"]]),Params[["Alpha_pr"]],1), dim=dim(Params[["Alpha_pr"]]))
      # Beta_pr, just estimate row column, which is interpreted as eigenvalues
      Map[["Beta_pr"]] = array(NA, dim=dim(Params[["Beta_pr"]]))
      Map[["Beta_pr"]][1,] = 1:ncol(Map[["Beta_pr"]])
      Map[["Beta_pr"]] = factor(Map[["Beta_pr"]])
      Params[["Beta_pr"]] = array( ifelse(!is.na(Map[["Beta_pr"]]),-0.5,0), dim=dim(Params[["Beta_pr"]]))  # B = U%*%L%*%solve(U) + I_pp, so -1<=eigenvalues<=0
    }
    # Identifiability restrictions on Beta_pr for co-integration method
    if( "Beta_pr"%in%names(Params) && Options_vec[["B_type"]]%in%c(0,2) ){
      Map[["Beta_pr"]] = factor(rbind( matrix(NA,nrow=Data$n_r,ncol=Data$n_r), matrix(seq_pos(Data$n_r*(Data$n_p-Data$n_r)),nrow=Data$n_p-Data$n_r,ncol=Data$n_r)))
    }
    if( "Beta_pr"%in%names(Params) && Options_vec[["B_type"]]%in%c(3) ){
      Map[["Beta_pr"]] = rbind( matrix(NA,nrow=Data$n_r,ncol=Data$n_r), matrix(seq_pos(Data$n_r*(Data$n_p-Data$n_r)),nrow=Data$n_p-Data$n_r,ncol=Data$n_r))
      Map[["Beta_pr"]][cbind(1:Data$n_r,1:Data$n_r)] = max(Map[["Beta_pr"]],0,na.rm=TRUE) + 1:Data$n_r
      Map[["Beta_pr"]] = factor(Map[["Beta_pr"]])
      Params[["Beta_pr"]][cbind(1:Data$n_r,1:Data$n_r)] = -0.5
    }
    # Better initial conditions for eigen-cointegration
    if( "Alpha_pr"%in%names(Params) && Options_vec[["B_type"]]%in%c(2) && Version%in%c("spatial_vam_v11") ){
      Params[["Alpha_pr"]][1,] = -0.5
    }
    if( "Alpha_pr"%in%names(Params) && Options_vec[["B_type"]]%in%c(3) ){
      Params[["Alpha_pr"]][] = 0.1 * rnorm(prod(dim(Params[["Alpha_pr"]])))
    }
    # Fix alpha if desired
    if( Data$Options_vec['IncludeAlpha']==0){
      Map[["Ainput_kp"]] = factor( array(NA,dim=dim(Params[["Ainput_kp"]])) )
      Map[["logtauA_p"]] = factor(NA)
    }
    # fix phi
    if( estimate_phi==FALSE ){
      Map[["phi_p"]] = factor( rep(NA,length(Params[["phi_p"]])) )
      Params[["phi_p"]][] = 0
    }
    # Make B_pp diagonal
    if( independentTF==TRUE ){
      if( "B_pp" %in% names(Params) ){
        Map[["B_pp"]] = matrix(NA, nrow=Data$n_p, ncol=Data$n_p)
        diag(Map[["B_pp"]]) = 1:Data$n_p
        Map[["B_pp"]] = factor( Map[["B_pp"]] )
        Params[["B_pp"]][] = 0
        Params[["L_val"]][] = 0
      }
      if( "Alpha_pr"%in%names(Params) && n_cointegrate==Data$n_p ){
        Map[["Alpha_pr"]] = matrix(NA, nrow=Data$n_p, ncol=Data$n_p)
        diag(Map[["Alpha_pr"]]) = 1:Data$n_p
        Map[["Alpha_pr"]] = factor( Map[["Alpha_pr"]] )
        Params[["Alpha_pr"]][] = 0
        Params[["L_val"]][] = 0
      }
    }
    # Fix logkappa_z at shared value by default
    if("logkappa_z" %in% names(Params)) Map[['logkappa_z']] = factor( rep(1,length(Params[["logkappa_z"]])) )
  } # End spatial
  
  Return = list("TmbData"=Data, "TmbParams"=Params, "Random"=Random, "Map"=Map, "MeshList"=MeshList)
  return(Return)
}
