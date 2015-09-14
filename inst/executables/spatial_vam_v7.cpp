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
  // Slot 0 -- DEFUNCT
  // Slot 1 -- include spatial variation in alpha
  // Slot 2 -- option for diagonal covariance (i.e., independence among species)
  // Slot 3 -- equilibrium distribution; 0=start from stationary mean; 1=start from variation but not density dependence
  DATA_FACTOR( ObsModel_p );
  // Slot 1-n_p -- distribution of data: 0=Poisson; 1=Lognormal; 2=Zero-inflated lognormal; 3=lognormal-Poisson; 4=Normal
  
  // Indices
  DATA_INTEGER(n_i);       // Total number of observations (i)
  DATA_INTEGER(n_s);	   // Number of stations (s)
  DATA_INTEGER(n_t);	   // Number of years (t)
  DATA_INTEGER(n_k);	   // Number of knots (k)
  DATA_INTEGER(n_p);   // Number of species (p)
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
  DATA_STRUCT(spde,spde_aniso_t);
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  // Fixed effects
  PARAMETER_VECTOR(Hinput_z); // Anisotropy parameters
  PARAMETER(logkappa);         // Controls range of spatial variation
  PARAMETER_VECTOR(alpha_p);   // Mean of Gompertz-drift field
  PARAMETER_VECTOR(phi_p);              // Offset of beginning from equilibrium
  PARAMETER_VECTOR(logtauA_p);        // log-inverse SD of Alpha
  PARAMETER_VECTOR(L_val);    // Values in loadings matrix
  PARAMETER_MATRIX(B_pp);   // Interaction matrix
  PARAMETER_MATRIX(logsigma_pz);

  // Random effects
  PARAMETER_ARRAY(d_ktp);  // Spatial process variation
  PARAMETER_MATRIX(Ainput_kp);  // Spatial process variation
  PARAMETER_VECTOR(delta_i);

  // global stuff
  int n_l = Z_kl.row(0).size();
  Type jnll = 0;
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();
  
  // Anisotropy elements
  matrix<Type> H(2,2);
  H(0,0) = exp(Hinput_z(0));
  H(1,0) = Hinput_z(1);
  H(0,1) = Hinput_z(1);
  H(1,1) = (1+Hinput_z(1)*Hinput_z(1)) / exp(Hinput_z(0));

  // Covariance via trimmed Cholesky
  matrix<Type> L_pj(n_p, n_j);
  matrix<Type> L_jp(n_j, n_p);
  L_pj.setZero();
  matrix<Type> Identity_pp(n_p, n_p);
  Identity_pp.setIdentity();
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
  
  // Derived quantities related to GMRF
  Type Range = sqrt(8) / exp( logkappa );
  Type logtauE = log(1/(exp(logkappa)*sqrt(4*M_PI)));
  vector<Type> MargSigmaA_p(n_p);
  for(int p=0; p<n_p; p++) MargSigmaA_p(p) = pow(4*M_PI,0.5) / exp(logtauA_p(p)) / exp(logkappa);    
  //Eigen::SparseMatrix<Type> Q = exp(4.0*logkappa)*G0 + Type(2.0)*exp(2.0*logkappa)*G1 + G2;
  Eigen::SparseMatrix<Type> Q = Q_spde(spde, exp(logkappa), H);
  GMRF_t<Type> nll_gmrf_spatial(Q);
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
  // Alpha
  for(int p=0; p<n_p; p++){
    // SCALE function --  If doing it this way, A_kp=Alpha_kp (i.e., Alpha_kp not rescaled prior to being used in the likelihood)
    //if( Options_vec(1)==1 ) jnll_comp(0) += SCALE(nll_gmrf_spatial, exp(-logtauA_p(p)))(A_kp.col(p));
    if( Options_vec(1)==1 ) jnll_comp(0) += nll_gmrf_spatial(Ainput_kp.col(p));
  }
  // Epsilon
  for(int t=0; t<n_t; t++){
    for(int k=0; k<n_k; k++){
    for(int p=0; p<n_p; p++){
      Epsilon_kp(k,p) = d_ktp(k,t,p) - dhat_ktp(k,t,p);
    }}
    jnll_comp(1) += SCALE(SEPARABLE(nll_mvnorm, nll_gmrf_spatial), exp(-logtauE))(Epsilon_kp); 
  }

  // Probability of observations
  vector<Type> logchat_i(n_i);
  Type encounterprob;
  Type log_notencounterprob;  
  for(int i=0; i<n_i; i++){
    logchat_i(i) = d_ktp(s_i(i),t_i(i),p_i(i));
    if( !isNA(c_i(i)) ){                
      if( ObsModel_p(p_i(i))==0 ) jnll_comp(2) -= dpois( c_i(i), exp(logchat_i(i)), true );
      if( ObsModel_p(p_i(i))==1 ) jnll_comp(2) -= dlognorm( c_i(i), logchat_i(i), exp(logsigma_pz(p_i(i),0)), true );
      if( ObsModel_p(p_i(i))==2 ){
        encounterprob = ( 1.0 - exp(-1 * exp(logchat_i(i)) * exp(logsigma_pz(p_i(i),1))) );
        log_notencounterprob = -1 * exp(logchat_i(i)) * exp(logsigma_pz(p_i(i),1));
        jnll_comp(2) -= dzinflognorm( c_i(i), logchat_i(i)-log(encounterprob), encounterprob, log_notencounterprob, exp(logsigma_pz(p_i(i),0)), true);
      }
      if( ObsModel_p(p_i(i))==3 ) jnll_comp(2) -= dpois( c_i(i), exp(logchat_i(i)+delta_i(i)), true );
      if( ObsModel_p(p_i(i))==4 ) jnll_comp(2) -= dnorm( c_i(i), logchat_i(i), exp(logsigma_pz(p_i(i),0)), true );
    }
  }

  // Probability of overdispersion
  if( Options_vec(0)==3 ){
    for(int i=0; i<n_i; i++){    
      jnll_comp(2) -= dnorm( delta_i(i), Type(0.0), exp(logsigma_pz(p_i(i),0)), true);
    }
  }

  // Combine NLL
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
  REPORT( logkappa );
  REPORT( Hinput_z );
  REPORT( logsigma_pz );
  // Spatial field summaries
  REPORT( Range );
  REPORT( Cov_pp );
  REPORT( B_pp );
  REPORT( jnll );
  REPORT( jnll_comp );
  REPORT( logtauE );
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
  
  return jnll;
}
