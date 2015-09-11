
MakeInputs_Fn = function( options_vec, loc_x, data_frame, n_factors=1 ){

  # Build SPDE object using INLA
  MeshList = SpatialDeltaGLMM::Calc_Anisotropic_Mesh(loc_x=loc_x)

  if( !all( c("sitenum","year","spp","catch") %in% colnames(data_frame)) ) stop( "data_frame must contain columns labeled sitenum, year, spp, and catch")

  # Data
  # Necessary columns: sitenum, year, catch, spp
  if(Version%in%c("spatial_vam_v2","spatial_vam_v1")) Data = list("Options_vec"=options_vec, "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=length(unique(data_frame[,'year'])), "n_k"=MeshList$mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-1, "spde"=NULL)
  if(Version%in%c("spatial_vam_v5","spatial_vam_v4","spatial_vam_v3")) Data = list("Options_vec"=options_vec, "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=length(unique(data_frame[,'year'])), "n_k"=MeshList$mesh$n, "n_p"=length(unique(data_frame[,'spp'])), "n_j"=n_factors, "c_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'spp'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-min(data_frame[,'year']), "spde"=NULL, "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2)
  if("spde" %in% names(Data)) Data[['spde']] = list("n_s"=MeshList$spde$n.spde, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$spde$param.inla$M0, "G0_inv"=inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)) )

  # Parameters
  if(Version=="spatial_vam_v1") Params = list("Hinput_z"=c(0,0), "logkappa"=log(1), "alpha_p"=rep(0,Data$n_p), "phi_p"=rep(0,Data$n_p), "logtauA_p"=rep(0,Data$n_p), "L_val"=rnorm(Data$n_j*Data$n_p-Data$n_j*(Data$n_j-1)/2), "B_pp"=matrix(rnorm(Data$n_p^2,sd=0.1),Data$n_p,Data$n_p), "d_ktp"=array(2,dim=unlist(Data[c('n_k','n_t','n_p')])), "Ainput_kp"=matrix(0,nrow=Data$n_k,ncol=Data$n_p))
  if(Version%in%c("spatial_vam_v3","spatial_vam_v2")) Params = list("Hinput_z"=c(0,0), "logkappa"=log(1), "alpha_p"=rep(0,Data$n_p), "phi_p"=rep(0,Data$n_p), "logtauA_p"=rep(0,Data$n_p), "L_val"=ifelse( is.na(fixdiag(Nrow=Data$n_p, Ncol=Data$n_j)), 1, 0), "B_pp"=diag(0.5,Data$n_p), "logsigma_pz"=matrix(0,nrow=Data$n_p,1), "d_ktp"=abind(SimList$d_stp,array(0,dim=c(Data$n_k-Data$n_s,Data$n_t,Data$n_p)),along=1), "Ainput_kp"=matrix(0,nrow=Data$n_k,ncol=Data$n_p))      # "d_ktp"=array(2,dim=unlist(Data[c('n_k','n_t','n_p')]))
  if(Version%in%c("spatial_vam_v4")) Params = list("Hinput_z"=c(0,0), "logkappa"=log(1), "alpha_p"=rep(0,Data$n_p), "phi_p"=rep(0,Data$n_p), "logtauA_p"=rep(0,Data$n_p), "L_val"=ifelse( is.na(fixdiag(Nrow=Data$n_p, Ncol=Data$n_j)), 1, 0), "B_pp"=diag(0.5,Data$n_p), "logsigma_pz"=matrix(0,nrow=Data$n_p,2), "d_ktp"=array(0,dim=c(Data$n_k,Data$n_t,Data$n_p)), "Ainput_kp"=matrix(0,nrow=Data$n_k,ncol=Data$n_p))      # "d_ktp"=array(2,dim=unlist(Data[c('n_k','n_t','n_p')]))
  if(Version%in%c("spatial_vam_v5")) Params = list("Hinput_z"=c(0,0), "logkappa"=log(1), "alpha_p"=rep(0,Data$n_p), "phi_p"=rep(0,Data$n_p), "logtauA_p"=rep(0,Data$n_p), "L_val"=ifelse( is.na(fixdiag(Nrow=Data$n_p, Ncol=Data$n_j)), 1, 0), "B_pp"=diag(0.5,Data$n_p), "logsigma_pz"=matrix(0,nrow=Data$n_p,2), "d_ktp"=array(0,dim=c(Data$n_k,Data$n_t,Data$n_p)), "Ainput_kp"=matrix(0,nrow=Data$n_k,ncol=Data$n_p), "delta_i"=rep(0,Data$n_i))      # "d_ktp"=array(2,dim=unlist(Data[c('n_k','n_t','n_p')]))

  # Random
  # Treating alpha_p, phi_p and B_pp as random (in REML) results in very slow inner optimization!  (100s of steps)
  Random = c( "Ainput_kp", "d_ktp" )
  if( Version%in%c("spatial_vam_v5")) Random = c(Random, "delta_i")
  if(Use_REML==TRUE) Random = c(Random, "alpha_p", "phi_p")  # , "B_pp"

  # Map
  Map = NULL
  Map[["Hinput_z"]] = factor( rep(NA,2) )
  if(Data$Options_vec[1]==0){
    Map[["logsigma_pz"]] = factor( array(NA,dim=dim(Params[["logsigma_pz"]])) )
    if( Version%in%c("spatial_vam_v5")) Map[["delta_i"]] = factor( rep(NA,length(Params$delta_i)) )
  }
  if(Data$Options_vec[1]==1){
    Map[["logsigma_pz"]] = factor( cbind(1:Data$n_p,NA) )
    if( Version%in%c("spatial_vam_v5")) Map[["delta_i"]] = factor( rep(NA,length(Params$delta_i)) )
  }
  if(Data$Options_vec[1]==2){
    if( Version%in%c("spatial_vam_v5")) Map[["delta_i"]] = factor( rep(NA,length(Params$delta_i)) )
  }
  if(Data$Options_vec[1]==3){
    Map[["logsigma_pz"]] = factor( cbind(1:Data$n_p,NA) )
  }
  # Fix alpha if desired
  if( Data$Options_vec['IncludeAlpha']==0){
    Map[["Ainput_kp"]] = factor( array(NA,dim=dim(Params[["Ainput_kp"]])) )
    Map[["logtauA_p"]] = factor(NA)
  }
  # fix phi
  if( Estimate_Phi==FALSE ){
    Map[["phi_p"]] = factor( rep(NA,length(Params[["phi_p"]])) )
    Params[["phi_p"]][] = 0
  }
  # fix phi
  if( Data_set=="Swain_et_al" ){
    #Map[["phi_p"]] = factor( c(NA,NA,NA,1) )
    #Params[["phi_p"]][] = 0
  }
  # Make B_pp diagonal
  if( IndependentTF==TRUE ){
    Map[["B_pp"]] = matrix(NA, nrow=Data$n_p, ncol=Data$n_p)
    diag(Map[["B_pp"]]) = 1:Data$n_p
    Map[["B_pp"]] = factor( Map[["B_pp"]] )
    Params[["B_pp"]][] = 0
    Params[["L_val"]][] = 0
  }

  # Turn off additional stuff
  if(FALSE){
    Map[["logtauA_p"]] = factor( rep(NA,Data$n_p) )
    Map[["B_pp"]] = factor( array(NA,dim=dim(Params$B_pp)) )
    Params$B_pp[] = 0
  }

  Return = list("TmbData"=Data, "TmbParams"=Params, "Random"=Random, "Map"=Map)
  return(Return)
}
