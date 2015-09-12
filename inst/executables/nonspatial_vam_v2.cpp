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
  DATA_INTEGER(n_t);	   // Number of years (t)
  DATA_INTEGER(n_p);   // Number of species (p)
  DATA_INTEGER(n_j);   // Number of dynamic factors in process error (j)

  // Data
  DATA_VECTOR( c_i );         // Count for observation
  DATA_VECTOR( se_log_c_i );         // Count for observation
  DATA_FACTOR( p_i );       	// Species for observation
  DATA_FACTOR( t_i );       	// Year for observation

  // Fixed effects
  PARAMETER_VECTOR(alpha_p);   // Mean of Gompertz-drift field
  PARAMETER_VECTOR(phi_p);              // Offset of beginning from equilibrium
  PARAMETER_VECTOR(L_val);    // Values in loadings matrix
  PARAMETER_MATRIX(B_pp);   // Interaction matrix
  PARAMETER_MATRIX(logsigma_pz);

  // Random effects
  PARAMETER_ARRAY(d_tp);  // Spatial process variation

  // global stuff
  Type jnll = 0;
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();
  
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
  
  // Derived quantities covariance
  MVNORM_t<Type> nll_mvnorm(Cov_pp);
  
  // Calculate mean of stationary distribution
  matrix<Type> TempMat_pp(n_p, n_p);
  vector<Type> dinf_p(n_p);
  TempMat_pp.setIdentity();
  TempMat_pp = TempMat_pp - B_pp;
  TempMat_pp = atomic::matinv( TempMat_pp );
  //dinf_p = ( alpha_p.matrix() * TempMat_pp.transpose() ).array();
  for(int p=0; p<n_p; p++) dinf_p(p) = ( alpha_p * TempMat_pp.transpose().col(p).array() ).sum();
  
  // Calculate predicted density given stochastic process
  vector<Type> tmp1_p(n_p);
  vector<Type> tmp2_p(n_p);
  array<Type> dhat_tp(n_t, n_p);
  // First year
  for(int p=0; p<n_p; p++){
    if(Options_vec(3)==0) dhat_tp(0,p) = phi_p(p) + dinf_p(p);
    if(Options_vec(3)==1) dhat_tp(0,p) = phi_p(p) + alpha_p(p);
  }
  // Project forward
  for(int t=1; t<n_t; t++){
    for(int p=0; p<n_p; p++){
      dhat_tp(t,p) = alpha_p(p);
      for(int p1=0; p1<n_p; p1++) dhat_tp(t,p) += B_pp(p,p1) * d_tp(t-1,p1);
    }
  }
  
  // Probability of random fields
  vector<Type> Epsilon_p(n_p);
  // Epsilon
  for(int t=0; t<n_t; t++){
    for(int p=0; p<n_p; p++){
      Epsilon_p(p) = d_tp(t,p) - dhat_tp(t,p);
    }
    jnll_comp(1) += nll_mvnorm(Epsilon_p); 
  }

  // Probability of observations
  vector<Type> logchat_i(n_i);
  Type encounterprob;
  Type log_notencounterprob;  
  for(int i=0; i<n_i; i++){
    logchat_i(i) = d_tp(t_i(i),p_i(i));
    if( !isNA(c_i(i)) ){                
      if(ObsModel_p(p_i(i))==0) jnll_comp(2) -= dpois( c_i(i), exp(logchat_i(i)), true );
      if(ObsModel_p(p_i(i))==1) jnll_comp(2) -= dlognorm( c_i(i), logchat_i(i), sqrt( exp(2*logsigma_pz(p_i(i),0)) + pow(se_log_c_i(i),2) ), true );
      if(ObsModel_p(p_i(i))==2){
        encounterprob = ( 1.0 - exp(-1 * exp(logchat_i(i)) * exp(logsigma_pz(p_i(i),1))) );
        log_notencounterprob = -1 * exp(logchat_i(i)) * exp(logsigma_pz(p_i(i),1));
        jnll_comp(2) -= dzinflognorm( c_i(i), logchat_i(i)-log(encounterprob), encounterprob, log_notencounterprob, exp(logsigma_pz(p_i(i),0)), true);
      }
      // Not implemented: if(ObsModel_p(p_i(i))==3) 
      if(ObsModel_p(p_i(i))==4) jnll_comp(2) -= dnorm( c_i(i), logchat_i(i), exp(logsigma_pz(p_i(i),0)), true );
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
  
  // Parameters
  REPORT( alpha_p );
  REPORT( phi_p );
  REPORT( logsigma_pz );
  // Spatial field summaries
  REPORT( Cov_pp );
  REPORT( B_pp );
  REPORT( jnll );
  REPORT( jnll_comp );
  REPORT( L_pj );
  // Fields
  REPORT( dinf_p );
  REPORT( d_tp );
  REPORT( dhat_tp );
  // Sanity checks
  REPORT( L_jp );
  REPORT( M_PI );
  REPORT( logchat_i );
  // Derived summaries
  
  return jnll;
}
