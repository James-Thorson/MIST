#n_species=4; n_years=20; n_stations=25; B_pp=NULL; ObsModel="Poisson"; Cov_pp=NULL; phi_p=NULL; sdlog=0.1; SpatialScale=0.1; SD_A=0.5; SD_E=0.2; corr_E=0.5; rho=0.8; logMeanDens=1; RandomSeed=NA; Loc=NULL
#n_species=Nspecies; n_years=20; n_stations=30; phi_p=rep(0,Nspecies); SpatialScale=0.4; rho=0.5; SD_A=0.5; SD_E=0.2; corr_E=0.5; ObsModel=ObsModel 
Sim_Fn <-
function( n_species=4, n_years=20, n_stations=25, B_pp=NULL, ObsModel="Poisson", Cov_pp=NULL, B_params=c(0,0.2), phi_p=NULL, sdlog=0.1, SpatialScale=0.1, SD_A=0.5, SD_E=0.2, corr_E=0.5, rho=0.8, logMeanDens=1, RandomSeed=NA, Loc=NULL ){
  if( !is.na(RandomSeed) ) set.seed(RandomSeed) 
  require( RandomFields )
  
  # Covariance
  if( is.null(Cov_pp) ){
    Cov_pp = matrix(corr_E,n_species,n_species)
    diag(Cov_pp) = 1
    Cov_pp = (rep(SD_E,n_species)%o%rep(1,n_species)) * Cov_pp * (rep(1,n_species)%o%rep(SD_E,n_species))
  }
  L_pj = t(chol(Cov_pp))
  # Initial density relative to equilibrium
  if( is.null(phi_p) ) phi_p = rnorm(n_species, mean=0, sd=1)
  # Interaction matrix
  if( is.null(B_pp) ){
    B_pp = matrix( rnorm(n_species*n_species,mean=B_params[1],sd=B_params[2]), nrow=n_species, ncol=n_species)
    diag(B_pp) = rho
  }

  if(any( Mod(eigen(B_pp)$values)>1 )) stop( "B_pp is not stationary!")

  # Spatial model
  # Range = distance at which Correlation is approx. 0.1
  # Scale = Range/2
  if( is.null(Loc) ) Loc = cbind( "x"=runif(n_stations, min=0,max=1), "y"=runif(n_stations, min=0,max=1) )
  model_A <- RMgauss(var=SD_A^2, scale=SpatialScale)
  model_E <- RMgauss(var=1, scale=SpatialScale)     # Unit variance, so variance enters via Cov_pp (I confirmed that var input to RMgauss is the marginal variance)

  # Alpha
  alpha_p = (diag(Nspecies)-B_pp) %*% rep(logMeanDens,n_species)
  
  # Simulate Alpha
  A_sp = matrix(NA, nrow=n_stations, ncol=n_species)
  for(p in 1:n_species){
    A_sp[,p] = RFsimulate(model=model_A, x=Loc[,'x'], y=Loc[,'y'])@data[,1]
    A_sp[,p] = A_sp[,p] - mean(A_sp[,p]) + alpha_p[p]
  }
  
  # Stationary mean (Ives et al. 2003, Eq. 15)
  if( prod(eigen(B_pp)$values-1)==0 ){
    dinf_sp = NULL
    dzero_sp = A_sp
  }else{
    dzero_sp = dinf_sp = t( solve(diag(Nspecies)-B_pp) %*% t(A_sp))
  }

  # Simulate Epsilon
  D_stj = array(NA, dim=c(n_stations,n_years,n_species))
  E_stp = array(NA, dim=c(n_stations,n_years,n_species))
  for(t in 1:n_years){
    for(j in 1:n_species) D_stj[,t,j] = RFsimulate(model=model_E, x=Loc[,'x'], y=Loc[,'y'])@data[,1]
    E_stp[,t,] = D_stj[,t,] %*% t(L_pj)
  }
  # Sanity check on marginal process error variance: mean(apply(SimList$D_stj, MARGIN=2:3, FUN=function(vec){mean(vec^2)}))
  # Sanity check on transformed process error variance: mean(apply(SimList$E_stp, MARGIN=2:3, FUN=function(vec){sqrt(mean(vec^2))}))

  # calculate log-density d_stp
  d_stp = dhat_stp = array(NA, dim=c(n_stations,n_years,n_species))
  # First year
  for(s in 1:n_stations){
    dhat_stp[s,1,] = phi_p + dzero_sp[s,]
    d_stp[s,1,] = dhat_stp[s,1,] + E_stp[s,1,]
  }
  # Project forward
  for(t in 2:n_years){
  for(s in 1:n_stations){
    dhat_stp[s,t,] = A_sp[s,] + B_pp%*%d_stp[s,t-1,]
    d_stp[s,t,] = dhat_stp[s,t,] + E_stp[s,t,]
  }}
  
  # Simulate data
  DF = NULL
  for(t in 1:n_years){
  for(s in 1:n_stations){
  for(p in 1:n_species){
    Tmp = c("sitenum"=s, "spp"=p, "year"=t, "catch"=NA, 'waterTmpC'=0, 'lambda'=exp(d_stp[s,t,p]) )
    if(ObsModel=="Poisson") Tmp['catch'] = rpois(1,lambda=Tmp['lambda'])
    if(ObsModel=="Lognormal") Tmp['catch'] = rlnorm(1,meanlog=log(Tmp['lambda']),sdlog=sdlog)
    DF = rbind(DF, Tmp)
  }}}
  DF = data.frame(DF, row.names=NULL)
  DF[,'spp'] = factor( letters[DF[,'spp']] )
  if( n_species>26 ) stop( "problem with using letters")

  ### Explore simulated system
  # Stability (AR coefficient in direction of each eigenvector)
  Mod(eigen(B_pp)$values)
  
  # Stationary variance (Ives et al. 2003, Eq. 17)
  Vinf_pp = matrix( solve(diag(Nspecies^2) - kronecker(B_pp,B_pp)) %*% as.vector(Cov_pp), nrow=Nspecies, ncol=Nspecies)
  
  # Proportion of stationary covariance attributable to interactions B_pp (Ives et al. 2003, Eq. 24)
  approx_equal = function( n1,n2, tol=1e-6 ) ifelse( (abs(n1-n2)/mean(c(n1,n2)))<tol, TRUE, FALSE) 
  #if( !approx_equal(det(B_pp)^2,  det(Vinf_pp - Cov_pp)/det(Vinf_pp)) ) stop("Something wrong with MeanPropVar calculation")
  # Sanity check: det(Vinf - Cov_pp) / det(Vinf) # Should be equal, and is
  MeanPropVar = det(B_pp)^(2/Nspecies)    # proportion of variance due to interactions, calculated as (geometric) average per species
  
  # Reactivitiy (Reactivity<0, where 0 is non-stationary community and higher is more reactive, Ives et al. 2003 Eq. 25)
  trace_fn = function(mat) sum(diag(mat))
  Reactivity = -trace_fn(Cov_pp) / trace_fn(Vinf_pp)
  
  # Maximum reactivity (upper bound on reactivity, calculated from B without Vinf, Ives et al. 2003 Eq. 26)
  MaxReactivity = max(eigen( t(B_pp)%*%B_pp )$values) - 1
  
  # Plot total abundance
  par( mfrow=c(1,2), mar=c(3,3,0,0), mgp=c(1.75,0.25,0), tck=-0.02 )  
  matplot( apply( d_stp, MARGIN=2:3, FUN=sum), type="l", xlab="Year", ylab="Total log-abundance")
  plot( y=DF$catch, x=DF$lambda, col=rainbow(n_species)[as.numeric(DF$spp)], xlab="Expected count", ylab="Observed count")
    
  # Return stuff
  Sim_List = list("DF"=DF, "L_pj"=L_pj, "Cov_pp"=Cov_pp, "B_pp"=B_pp, "alpha_p"=alpha_p, "phi_p"=phi_p, "Loc"=Loc, "A_sp"=A_sp, "D_stj"=D_stj, "E_stp"=E_stp, "dhat_stp"=dhat_stp, "d_stp"=d_stp, "dinf_sp"=dinf_sp, "Vinf_pp"=Vinf_pp, "MeanPropVar"=MeanPropVar, "Reactivity"=Reactivity, "MaxReactivity"=MaxReactivity)
  return(Sim_List)
}
