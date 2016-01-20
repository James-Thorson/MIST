// Space time 
#include <TMB.hpp>

// trace of a matrix
template<class Type>
Type trace( matrix<Type> mat ){
  Type Return = 0;
  for(int i=0; i<mat.col(0).size(); i++) Return += mat(i,i); 
  return Return;
}

// 2nd power of a number
template<class Type>
Type square(Type x){ return x*x; }

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// function for logistic transform
template<class Type>
Type plogis(Type x){
  return 1.0 / (1.0 + exp(-x));
}

// dlognorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=false){
  Type Return;
  if(give_log==false) Return = dnorm( log(x), meanlog, sdlog, false) / x;
  if(give_log==true) Return = dnorm( log(x), meanlog, sdlog, true) - log(x);
  return Return;
}

// dzinflognorm
template<class Type>
Type dzinflognorm(Type x, Type meanlog, Type encounter_prob, Type log_notencounter_prob, Type sdlog, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true){
      if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
      if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
    }
  }else{
    if(give_log==false) Return = encounter_prob * dlognorm( x, meanlog, sdlog, false );
    if(give_log==true) Return = log(encounter_prob) + dlognorm( x, meanlog, sdlog, true );
  } 
  return Return;
}

// dzinfgamma, shape = 1/CV^2, scale = mean*CV^2
template<class Type>
Type dzinfgamma(Type x, Type posmean, Type encounter_prob, Type log_notencounter_prob, Type cv, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true){
      if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
      if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
    }
  }else{
    if(give_log==false) Return = encounter_prob * dgamma( x, pow(cv,-2), posmean*pow(cv,2), false );
    if(give_log==true) Return = log(encounter_prob) + dgamma( x, pow(cv,-2), posmean*pow(cv,2), true );
  } 
  return Return;
}

// dzinfnorm
template<class Type>
Type dzinfnorm(Type x, Type posmean, Type encounter_prob, Type log_notencounter_prob, Type cv, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true){
      if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
      if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
    }
  }else{
    if(give_log==false) Return = encounter_prob * dnorm( x, posmean, posmean*cv, false );
    if(give_log==true) Return = log(encounter_prob) + dnorm( x, posmean, posmean*cv, true );
  } 
  return Return;
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace R_inla;
  using namespace Eigen;
  using namespace density;

  // Settings
  DATA_FACTOR( Options_vec );
  // Slot 0 -- Method for assembling B_pp
  // Slot 1 -- include spatial variation in alpha
  // Slot 2 -- option for diagonal covariance (i.e., independence among species)
  // Slot 3 -- equilibrium distribution; 0=start from stationary mean; 1=start from variation but not density dependence
  DATA_FACTOR( ObsModel_p );
  // Slot 1-n_p -- distribution of data: 0=Poisson; 1=Lognormal; 2=Zero-inflated lognormal; 3=lognormal-Poisson; 4=Normal
  DATA_VECTOR( PenMult_z );  // 0: Penalty on eigenvalues; 1: Penalty on columns of Alpha_pr

  // Indices
  DATA_INTEGER(n_i);       // Total number of observations (i)
  DATA_INTEGER(n_s);	   // Number of stations (s)
  DATA_INTEGER(n_t);	   // Number of years (t)
  DATA_INTEGER(n_k);	   // Number of knots (k)
  DATA_INTEGER(n_p);   // Number of species (p)
  DATA_INTEGER(n_r);   // Rank of species interaction matrix B_pp (defining number of columns for alpha_pr and beta_pr)
  DATA_INTEGER(n_j);   // Number of dynamic factors in process error (j)

  // Data
  DATA_VECTOR( c_i );         // Count for observation
  DATA_FACTOR( p_i );       	// Species for observation
  DATA_FACTOR( s_i );       	// Site for observation
  DATA_FACTOR( t_i );       	// Year for observation

  // Derived quantities
  DATA_VECTOR( a_k );        // Area for each "real" stratum(km^2) in each stratum (zero for all knots k not associated with stations s)
  DATA_MATRIX( Z_kl );        // Derived quantity matrix

  // Aniso objects
  DATA_STRUCT(spde, spde_aniso_t);
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  // Fixed effects
  PARAMETER_VECTOR(Hinput_z); // Anisotropy parameters
  PARAMETER_VECTOR(logkappa_z);         // Controls range of spatial variation.  First n_p slots are independent for each spatial component (but can be fixed to be equal).  slot n_p is for spatio-temporal (SDFA) component
  PARAMETER_VECTOR(alpha_p);   // Mean of Gompertz-drift field
  PARAMETER_VECTOR(phi_p);              // Offset of beginning from equilibrium
  PARAMETER_VECTOR(logMargSigmaA_p);        // log-inverse SD of Alpha  // logtauA_p
  PARAMETER_VECTOR(L_val);    // Values in loadings matrix
  PARAMETER_MATRIX(Alpha_pr);   // error correction responses
  PARAMETER_MATRIX(Beta_pr);   // error correction loadings, B_pp = Alpha_pr %*% t(Beta_pr)
  PARAMETER_MATRIX(logsigma_pz);

  // Random effects
  PARAMETER_ARRAY(d_ktp);  // Spatial process variation
  PARAMETER_MATRIX(Ainput_kp);  // Spatial process variation
  PARAMETER_VECTOR(delta_i);

  // global stuff
  int n_l = Z_kl.row(0).size();
  Type jnll = 0;
  vector<Type> jnll_comp(4);
  jnll_comp.setZero();
  matrix<Type> Identity_pp(n_p, n_p);
  Identity_pp.setIdentity();
  
  // Anisotropy elements
  matrix<Type> H(2,2);
  H(0,0) = exp(Hinput_z(0));
  H(1,0) = Hinput_z(1);
  H(0,1) = Hinput_z(1);
  H(1,1) = (1+Hinput_z(1)*Hinput_z(1)) / exp(Hinput_z(0));

  // Calculate interactions matrix B_pp
  matrix<Type> Alpha_rp = Alpha_pr.transpose();
  matrix<Type> Beta_rp = Beta_pr.transpose();
  matrix<Type> B_pp(n_p, n_p);
  // Simple co-integration
  if( Options_vec(0)==0 ){
    B_pp = Alpha_pr * Beta_rp + Identity_pp;
  }
  // Simple eigen-decomposition
  if( Options_vec(0)==1 ){
    matrix<Type> L_rr(n_r, n_r);
    L_rr.setZero();
    for(int r=0; r<n_r; r++) L_rr(r,r) = Beta_pr(0,r);
    matrix<Type> Alpha_rr = Alpha_rp * Alpha_pr;
    matrix<Type> invAlpha_rp = atomic::matinv( Alpha_rr ) * Alpha_rp; // solve(t(mat)%*%mat) %*%t (mat)
    B_pp = Alpha_pr * L_rr;
    B_pp = B_pp * invAlpha_rp + Identity_pp;
    REPORT( invAlpha_rp );
    REPORT( L_rr );
  }
  // Hybrid eigen-cointegration
    // B_pp = Alpha_pp %*% t(Beta_pp) %*% L_pp %*% t(solve(Beta_pp)) %*% solve(Alpha_pp) + diag(n_p)
  if( Options_vec(0)==2 | Options_vec(0)==3 ){
    matrix<Type> Alpha_pp = Identity_pp;
    // Make Alpha_pp
    vector<Type> colnorm_r( n_r );
    colnorm_r.setZero();
    for(int p=0; p<n_p; p++){
    for(int r=0; r<n_r; r++){
      if( Options_vec(0)==2 ){
        if(p==r) Alpha_pp(p,r) = 1;
        if(p!=r) Alpha_pp(p,r) = Alpha_pr(p,r);
      }
      if( Options_vec(0)==3 ){
        Alpha_pp(p,r) = Alpha_pr(p,r);
      }
      colnorm_r(r) += square( Alpha_pp(p,r) );
    }}
    for(int p=0; p<n_p; p++){
    for(int r=0; r<n_r; r++){
      Alpha_pp(p,r) /= sqrt( colnorm_r(r) );
    }}
    // Make Beta_pp
    matrix<Type> Beta_pp = Identity_pp;
    for(int p=n_r; p<n_p; p++){
    for(int r=0; r<n_r; r++){
      Beta_pp(p,r) = Beta_pr(p,r);
    }}
    // Make L_pp
    matrix<Type> L_pp(n_p, n_p);
    L_pp.setZero();
    for(int r=0; r<n_r; r++){
      if( Options_vec(0)==2 ) L_pp(r,r) = Alpha_pr(r,r);
      if( Options_vec(0)==3 ) L_pp(r,r) = Beta_pr(r,r);
    }
    // Build B_pp
    matrix<Type> trans_Beta_pp = Beta_pp.transpose();
    matrix<Type> trans_invBeta_pp = atomic::matinv( Beta_pp ).transpose();
    matrix<Type> invAlpha_pp = atomic::matinv( Alpha_pp );
    B_pp = Alpha_pp * trans_Beta_pp;
    B_pp = B_pp * L_pp;
    B_pp = B_pp * trans_invBeta_pp;
    B_pp = B_pp * invAlpha_pp + Identity_pp;
    REPORT( Alpha_pp );
    REPORT( Beta_pp );
    REPORT( L_pp );
    // Penalize colnorm_r
    if( Options_vec(0)==3 ) jnll_comp(3) += PenMult_z(1) * ( log(colnorm_r)*log(colnorm_r) ).sum();
  }

  // Penalize unstable interaction matrices
  if( false ){
    vector< std::complex< Type > > eigenvalues_B_pp = B_pp.eigenvalues();
    vector< Type > real_eigenvalues_B_pp = eigenvalues_B_pp.real();
    vector< Type > imag_eigenvalues_B_pp = eigenvalues_B_pp.imag();
    if( PenMult_z(0) > 0 ){
      for( int p=0; p<n_p; p++){
        jnll_comp(3) += PenMult_z(0) * CppAD::CondExpGt( real_eigenvalues_B_pp(p), Type( 1.0), pow(real_eigenvalues_B_pp(p)-Type( 1.0),2), Type(0.0));
        jnll_comp(3) += PenMult_z(0) * CppAD::CondExpLt( real_eigenvalues_B_pp(p), Type(-1.0), pow(real_eigenvalues_B_pp(p)-Type(-1.0),2), Type(0.0));
      }
    }
    REPORT( real_eigenvalues_B_pp );
    REPORT( imag_eigenvalues_B_pp );
  }

  // Covariance via trimmed Cholesky
  matrix<Type> L_pj(n_p, n_j);
  matrix<Type> L_jp(n_j, n_p);
  L_pj.setZero();
  matrix<Type> Cov_pp(n_p, n_p);
  Cov_pp.setZero();
  int Count = 0;
  if( Options_vec(2)==0 ){
    // Assemble the loadings matrix (lower-triangular, loop through rows then columns)
    for(int j=0; j<n_j; j++){
    for(int p=0; p<n_p; p++){
      if(j<=p){
        L_pj(p,j) = L_val(Count);
        Count++;
      }
    }}
    // Calculate the covariance
    L_jp = L_pj.transpose();
    Cov_pp = L_pj*L_jp + Type(0.000001)*Identity_pp;  // additive constant to make Cov_pp invertible
  }
  // Diagonal covariance 
  if( Options_vec(2)==1 ){
    for(int p=0; p<n_p; p++){
      Cov_pp(p,p) = exp( L_val(p) );
    }
  }
  
  // Derived quantities related to spatial variation (logtauA=Spatial;  logtauE=Spatiotemporal)
  // calculation of logtauE_p and Range_pz depends upon whether we're treating species as independent or not
  matrix<Type> Range_pz(n_p, 2);
  vector<Type> logtauE_p(n_p);
  vector<Type> logtauA_p(n_p);
  vector<Type> MargSigmaA_p(n_p);
  for(int p=0; p<n_p; p++){
    logtauA_p(p) = 0.5*log(4*M_PI) - logMargSigmaA_p(p) - logkappa_z(p);
    Range_pz(p,0) = sqrt(8) / exp( logkappa_z(p) );
    if( Options_vec(2)==0 ){
      logtauE_p(p) = log(1/( exp(logkappa_z(n_p)) * sqrt(4*M_PI)) );
      Range_pz(p,1) = sqrt(8) / exp( logkappa_z(n_p) );
    }
    if( Options_vec(2)==1 ){
      logtauE_p(p) = log(1/( exp(logkappa_z(p)) * sqrt(4*M_PI)) );
      Range_pz(p,1) = sqrt(8) / exp( logkappa_z(p) );
    }
    MargSigmaA_p(p) = exp( logMargSigmaA_p(p) );
  }
  // Derived quantities related to GMRF
  Eigen::SparseMatrix<Type> Q_spatial;
  Eigen::SparseMatrix<Type> Q_spatiotemporal;
  MVNORM_t<Type> nll_mvnorm(Cov_pp);

  // Transform random fields
  matrix<Type> A_kp(n_k, n_p);
  for(int k=0; k<n_k; k++){
  for(int p=0; p<n_p; p++){
    A_kp(k,p) = Ainput_kp(k,p)/exp(logtauA_p(p)) + alpha_p(p);
  }}

  // Calculate mean of stationary distribution
  matrix<Type> TempMat_pp(n_p, n_p);
  matrix<Type> dinf_kp(n_k, n_p);
  TempMat_pp.setIdentity();
  TempMat_pp = TempMat_pp - B_pp;
  TempMat_pp = atomic::matinv( TempMat_pp );
  dinf_kp = A_kp * TempMat_pp.transpose();   // saves a transpose relative to t(Temp_pp * t(A_kp))
  
  // Calculate predicted density given stochastic process
  vector<Type> tmp1_p(n_p);
  vector<Type> tmp2_p(n_p);
  array<Type> dhat_ktp(n_k, n_t, n_p);
  // First year
  for(int k=0; k<n_k; k++){
  for(int p=0; p<n_p; p++){
    if(Options_vec(3)==0) dhat_ktp(k,0,p) = phi_p(p) + dinf_kp(k,p);
    if(Options_vec(3)==1) dhat_ktp(k,0,p) = phi_p(p) + A_kp(k,p);
  }}
  // Project forward
  for(int t=1; t<n_t; t++){
    for(int k=0; k<n_k; k++){
    for(int p=0; p<n_p; p++){
      dhat_ktp(k,t,p) = A_kp(k,p);
      for(int p1=0; p1<n_p; p1++) dhat_ktp(k,t,p) += B_pp(p,p1) * d_ktp(k,t-1,p1);
    }}
  }
  
  // Probability of random fields
  array<Type> Epsilon_kp(n_k, n_p);
  // Alpha (spatial variation in productivity)
  for(int p=0; p<n_p; p++){
    // SCALE function --  If doing it this way, A_kp=Alpha_kp (i.e., Alpha_kp not rescaled prior to being used in the likelihood)
    //if( Options_vec(1)==1 ) jnll_comp(0) += SCALE(nll_gmrf_spatial, exp(-logtauA_p(p)))(A_kp.col(p));
    Q_spatial = Q_spde(spde, exp(logkappa_z(p)), H);
    if( Options_vec(1)==1 ) jnll_comp(0) += GMRF(Q_spatial)(Ainput_kp.col(p));
  }
  // Epsilon (spatio-temporal process errors)
  for(int t=0; t<n_t; t++){
    for(int k=0; k<n_k; k++){
    for(int p=0; p<n_p; p++){
      Epsilon_kp(k,p) = d_ktp(k,t,p) - dhat_ktp(k,t,p);
    }}
    // Depends upon whether using independence or not
    if( Options_vec(2)==0 ){
      Q_spatiotemporal = Q_spde(spde, exp(logkappa_z(n_p)), H);
      jnll_comp(1) += SCALE(SEPARABLE(nll_mvnorm, GMRF(Q_spatiotemporal)), exp(-logtauE_p(0)))(Epsilon_kp);
    }
    if( Options_vec(2)==1 ){
      for(int p=0; p<n_p; p++){
        Q_spatiotemporal = Q_spde(spde, exp(logkappa_z(p)), H);
        jnll_comp(1) += SCALE( GMRF(Q_spatiotemporal), exp(L_val(p) - logtauE_p(p)))(Epsilon_kp.col(p));
      }
    }
  }

  // Probability of observations
  vector<Type> logchat_i(n_i);
  vector<Type> jnll_i(n_i);
  jnll_i.setZero();
  Type encounterprob;
  Type log_notencounterprob;  
  for(int i=0; i<n_i; i++){
    logchat_i(i) = d_ktp(s_i(i),t_i(i),p_i(i));
    if( !isNA(c_i(i)) ){                
      if( ObsModel_p(p_i(i))==0 ) jnll_i(i) = -1 * dpois( c_i(i), exp(logchat_i(i)), true );
      if( ObsModel_p(p_i(i))==1 ) jnll_i(i) = -1 * dlognorm( c_i(i), logchat_i(i), exp(logsigma_pz(p_i(i),0)), true );
      if( ObsModel_p(p_i(i))==2 ){
        encounterprob = ( 1.0 - exp(-1 * exp(logchat_i(i)) * exp(logsigma_pz(p_i(i),1))) );
        log_notencounterprob = -1 * exp(logchat_i(i)) * exp(logsigma_pz(p_i(i),1));
        jnll_i(i) = -1 * dzinflognorm( c_i(i), logchat_i(i)-log(encounterprob), encounterprob, log_notencounterprob, exp(logsigma_pz(p_i(i),0)), true);
      }
      if( ObsModel_p(p_i(i))==3 ) jnll_i(i) = -1 * dpois( c_i(i), exp(logchat_i(i)+delta_i(i)), true );
      if( ObsModel_p(p_i(i))==4 ) jnll_i(i) = -1 * dnorm( c_i(i), logchat_i(i), exp(logsigma_pz(p_i(i),0)), true );
    }
  }

  // Probability of overdispersion
  for(int i=0; i<n_i; i++){
    if( ObsModel_p(p_i(i))==3 ){
      jnll_i(i) -= dnorm( delta_i(i), Type(0.0), exp(logsigma_pz(p_i(i),0)), true);
    }
  }

  // Combine NLL
  jnll_comp(2) = jnll_i.sum();
  jnll = jnll_comp.sum();

  // Derived quantities -- Stationary variance
  int n2_p = n_p*n_p;
  matrix<Type> Kronecker_p2p2(n2_p,n2_p);
  Kronecker_p2p2 = kronecker( B_pp, B_pp );
  matrix<Type> InvDiff_p2p2(n2_p, n2_p);
  InvDiff_p2p2.setIdentity();
  InvDiff_p2p2 = InvDiff_p2p2 - Kronecker_p2p2;
  InvDiff_p2p2 = atomic::matinv( InvDiff_p2p2 );
  matrix<Type> Vinf_pp(n_p, n_p);
  Vinf_pp.setZero();
  for(int i=0; i<n_p; i++){
  for(int j=0; j<n_p; j++){
    int k = i + j*n_p;
    for(int i2=0; i2<n_p; i2++){
    for(int j2=0; j2<n_p; j2++){
      int k2 = i2 + j2*n_p;
      Vinf_pp(i,j) += InvDiff_p2p2(k,k2) * Cov_pp(i2,j2);
    }}
  }}
  REPORT( InvDiff_p2p2 );
  REPORT( Vinf_pp );
  REPORT( Kronecker_p2p2 );
  
  // Derived quantities -- reactivity = -trace(Cov_pp) / trace(Vinf_pp)
  Type Reactivity = -1 * trace(Cov_pp) / trace(Vinf_pp);
  REPORT( Reactivity );
  
  // Derived quantities -- geometric average proportion of stationary variance caused by interactions
  Type MeanPropVar = atomic::logdet(B_pp);
  MeanPropVar = exp( 2.0/n_p * MeanPropVar );
  REPORT( MeanPropVar );
  
  // Calculate indices
  array<Type> Index_ktp(n_k, n_t, n_p);
  matrix<Type> Index_tp(n_t, n_p);
  Index_tp.setZero();
  for(int k=0; k<n_k; k++){
  for(int t=0; t<n_t; t++){
  for(int p=0; p<n_p; p++){
    Index_ktp(k,t,p) = exp( d_ktp(k,t,p) ) * a_k(k) / 1000;  // Convert from kg to metric tonnes
    Index_tp(t,p) += Index_ktp(k,t,p); 
  }}}
  REPORT( Index_tp );
  ADREPORT( Index_tp );

  // Calculate other derived summaries
  // Each is the weighted-average X_xl over polygons (x) with weights equal to abundance in each polygon and time
  array<Type> mean_Z_tpl(n_t, n_p, n_l);
  mean_Z_tpl.setZero();
  int report_summary_TF = false;
  for(int t=0; t<n_t; t++){
  for(int p=0; p<n_p; p++){
  for(int l=0; l<n_l; l++){
    for(int k=0; k<n_k; k++){
      if( Z_kl(k,l)!=0 ){
        mean_Z_tpl(t,p,l) += Z_kl(k,l) * Index_ktp(k,t,p) / Index_tp(t,p);  
        report_summary_TF = true; 
      }
    }
  }}}
  if( report_summary_TF==true ){
    array<Type> cov_Z_tpll(n_t,n_p,n_l,n_l);
    cov_Z_tpll.setZero();
    for(int t=0; t<n_t; t++){
    for(int p=0; p<n_p; p++){
    for(int l1=0; l1<n_l; l1++){
    for(int l2=0; l2<n_l; l2++){
      for(int k=0; k<n_k; k++){
        cov_Z_tpll(t,p,l1,l2) += (Z_kl(k,l1)-mean_Z_tpl(t,p,l1))*(Z_kl(k,l2)-mean_Z_tpl(t,p,l2)) * Index_ktp(k,t,p)/Index_tp(t,p);  
      }
    }}}}
    REPORT( mean_Z_tpl );  
    ADREPORT( mean_Z_tpl );
    REPORT( cov_Z_tpll );  
  }
  
  // Parameters
  REPORT( alpha_p );
  REPORT( phi_p );
  REPORT( logkappa_z );
  REPORT( Hinput_z );
  REPORT( logsigma_pz );
  // Spatial field summaries
  REPORT( Range_pz );
  REPORT( Cov_pp );
  REPORT( B_pp );
  REPORT( Alpha_pr );
  REPORT( Beta_pr );
  REPORT( logtauE_p );
  REPORT( logtauA_p );
  REPORT( MargSigmaA_p );
  REPORT( L_pj );
  // Fields
  REPORT( dinf_kp );
  REPORT( A_kp );
  REPORT( d_ktp );
  REPORT( dhat_ktp );
  REPORT( Ainput_kp );
  // Sanity checks
  REPORT( L_jp );
  REPORT( M_PI );
  REPORT( logchat_i );
  // Derived summaries
  // Objective function components
  REPORT( jnll );
  REPORT( jnll_comp );
  REPORT( jnll_i );

  // Standard errors for derived quantities
  ADREPORT( B_pp );

  return jnll;
}
