Sim_Fn <-
function( n_species=4, n_years=20, n_stations=25, n_factors=2, B_pp=NULL, L_pj=NULL, phi_p=NULL, SpatialScale=0.1, SD_A=0.5, SD_E=0.2, rho=0.8, logMeanDens=1, RandomSeed=NA, Loc=NULL ){
  if( !is.na(RandomSeed) ) set.seed(RandomSeed) 
  require( RandomFields )
  
  # Parameters
  if( is.null(L_pj) ){
    L_pj = matrix( rnorm(n_factors*n_species), nrow=n_species, ncol=n_factors)
    for(i in 1:ncol(L_pj)){
      L_pj[seq(from=1,to=i-1,length=i-1),i] = 0
      if( L_pj[,i][which.max(abs(L_pj[,i]))]<0 ){
        L_pj[,i] = -1*L_pj[,i]
      }
    }
  }
  if( is.null(phi_p) ) phi_p = rnorm(n_species, mean=0, sd=1)
  alpha = rep(logMeanDens, n_species)

  if( is.null(B_pp) ){
    B_pp = matrix( rnorm(n_species*n_species,sd=0.1), nrow=n_species, ncol=n_species)
    diag(B_pp) = rho
  }

  # Spatial model
  if( is.null(Loc) ) Loc = cbind( "x"=runif(n_stations, min=0,max=1), "y"=runif(n_stations, min=0,max=1) )
  model_A <- RMgauss(var=SD_A^2, scale=SpatialScale)
  model_E <- RMgauss(var=SD_E^2, scale=SpatialScale)

  # Simulate Alpha
  A_sp = matrix(NA, nrow=n_stations, ncol=n_species)
  for(p in 1:n_species){
    A_sp[,p] = RFsimulate(model=model_A, x=Loc[,'x'], y=Loc[,'y'])@data[,1] + alpha[p]
  }
  
  # Simulate Epsilon
  E_stp = array(NA, dim=c(n_stations,n_years,n_species))
  D_stj = array(NA, dim=c(n_stations,n_years,n_factors))
  for(t in 1:n_years){
    for(j in 1:n_factors) D_stj[,t,j] = RFsimulate(model=model_E, x=Loc[,'x'], y=Loc[,'y'])@data[,1]
    E_stp[,t,] = D_stj[,t,] %*% t(L_pj)
  }

  # calculate log-density d_stp
  d_stp = dhat_stp = array(NA, dim=c(n_stations,n_years,n_species))
  # First year
  for(s in 1:n_stations){
  for(p in 1:n_species){
    dhat_stp[s,1,p] = phi_p[p] + A_sp[s,p]
    d_stp[s,1,p] = dhat_stp[s,1,p] + E_stp[s,1,p]
  }}
  # Project forward
  for(t in 2:n_years){
  for(s in 1:n_stations){
    dhat_stp[s,t,] = A_sp[s,p] + B_pp%*%dhat_stp[s,t-1,]
    d_stp[s,t,] = dhat_stp[s,t,] + E_stp[s,t,]
  }}
  
  # Simulate data
  DF = NULL
  for(t in 1:n_years){
  for(s in 1:n_stations){
  for(p in 1:n_species){
    Tmp = c("sitenum"=s, "spp"=p, "year"=t, "catch"=rpois(1,lambda=exp(d_stp[s,p,t])), 'waterTmpC'=0 )
    DF = rbind(DF, Tmp)
  }}}
  DF = data.frame(DF, row.names=NULL)
  DF[,'spp'] = factor( letters[DF[,'spp']] )
  if( n_species>26 ) stop( "problem with using letters")

  # Return stuff
  Sim_List = list("DF"=DF, "L_pj"=L_pj, "B_pp"=B_pp, "phi_p"=phi_p, "Loc"=Loc, "A_sp"=A_sp, "E_stp"=E_stp, "dhat_stp"=dhat_stp, "d_stp"=d_stp)
  return(Sim_List)
}
