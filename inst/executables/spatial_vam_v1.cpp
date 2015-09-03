// Space time 
#include <TMB.hpp>

// 2nd power of a number
template<class Type>
Type square(Type x){ return pow(x,2.0); }

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
  // Slot 0 -- distribution of data

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

  // Aniso objects
  DATA_STRUCT(spde,spde_aniso_t);

  // Fixed effects
  PARAMETER_VECTOR(Hinput_z); // Anisotropy parameters
  PARAMETER(logkappa);         // Controls range of spatial variation
  PARAMETER_VECTOR(alpha_p);   // Mean of Gompertz-drift field
  PARAMETER_VECTOR(phi_p);              // Offset of beginning from equilibrium
  PARAMETER_VECTOR(logtauA_p);        // log-inverse SD of Alpha
  PARAMETER_VECTOR(L_val);    // Values in loadings matrix
  PARAMETER_MATRIX(B_pp);   // Interaction matrix
  PARAMETER_VECTOR(log_sigma_p);

  // Random effects
  PARAMETER_ARRAY(d_ktp);  // Spatial process variation
  PARAMETER_MATRIX(A_input);  // Spatial process variation
  PARAMETER_VECTOR(tilda_i);

  // global stuff
  Type jnll = 0;
  
  // Anisotropy elements
  matrix<Type> H(2,2);
  H(0,0) = exp(Hinput_z(0));
  H(1,0) = Hinput_z(1);
  H(0,1) = Hinput_z(1);
  H(1,1) = (1+square(Hinput_z(1))) / exp(Hinput_z(0));

  // Assemble the loadings matrix
  matrix<Type> L_pj(n_p, n_j);
  int Count = 0;
  for(int p=0; p<n_p; p++){
  for(int j=0; j<n_j; j++){
    if(p>=j){
      L_pj(p,j) = L_val(Count);
      Count++;
    }else{
      L_pj(p,j) = 0.0;
    }
  }}
  
  // Calculate the precision matrix among species
  matrix<Type> L_jp(n_j, n_p);
  for(int p=0; p<n_p; p++){
  for(int j=0; j<n_j; j++){
    L_jp(j,p) = L_pj(p,j);
  }}
  matrix<Type> Cov_pp(n_p, n_p);
  matrix<Type> Prec_pp(n_p, n_p);
  Cov_pp = L_jp * L_pj;
  //Prec_pp = atomic::matinv( Cov_pp );
  
  // Derived quantities related to GMRF
  Type Range = sqrt(8) / exp( logkappa );
  Type logtauE = log(1/(exp(logkappa)*sqrt(4*3.141592)));
  Eigen::SparseMatrix<Type> Q;
  Q = Q_spde(spde, exp(logkappa), H);
  GMRF_t<Type> nll_gmrf_spatial(Q);
  MVNORM_t<Type> nll_mvnorm(Cov_pp);
  
  // Transform random fields
  matrix<Type> A_kp(n_k, n_p);
  for(int k=0; k<n_k; k++){
  for(int p=0; p<n_p; p++){
    A_kp(k,p) = alpha_p(p) + A_input(k,p)/exp(logtauA_p(p));
  }}
  
  // Calculate expected density given stochastic process
  vector<Type> tmp1_p(n_p);
  vector<Type> tmp2_p(n_p);
  array<Type> dhat_ktp(n_k, n_t, n_p);
  // First year
  for(int s=0; s<n_s; s++){
  for(int p=0; p<n_p; p++){
    dhat_ktp(s,0,p) = phi_p(p) + A_kp(s,p);
  }}
  for(int k=n_s; k<n_k; k++){
  for(int p=0; p<n_p; p++){
    dhat_ktp(k,0,p) = 0.0;
  }}
  // Project forward
  for(int t=1; t<n_t; t++){
    for(int s=0; s<n_s; s++){
      for(int p=0; p<n_p; p++) tmp1_p(p) = dhat_ktp(s,t,p);
      for(int p=0; p<n_p; p++){
        tmp2_p(p) = ( B_pp.row(p).array() * tmp1_p.array() ).sum();
        tmp2_p(p) += A_kp(s,p);
      }
      for(int p=0; p<n_p; p++) dhat_ktp(s,t,p) = tmp2_p(p);
    }
    for(int k=n_s; k<n_k; k++){
    for(int p=0; p<n_p; p++){
      dhat_ktp(k,t,p) = 0.0;
    }}
  }
  
  // Probability of random fields
  array<Type> Epsilon_kp(n_k, n_p);
  for(int t=0; t<n_t; t++){
    // Epsilon
    for(int k=0; k<n_k; k++){
    for(int p=0; p<n_p; p++){
      Epsilon_kp(k,p) = d_ktp(k,t,p) - dhat_ktp(k,t,p);
    }}
    jnll += SCALE(SEPARABLE(nll_gmrf_spatial, nll_mvnorm), exp(-logtauE))(Epsilon_kp); 
    //SEPARABLE(AR1(rho_j(j)),GMRF(Q))(Epsilon_tmp);
    //jnll += SEPARABLE(nll_gmrf_spatial, nll_mvnorm)(Epsilon_kp); 
    // Alpha
    for(int p=0; p<n_p; p++){
      jnll += nll_gmrf_spatial(A_input.col(p));
    }
  }

  // Probability of observations
  vector<Type> logchat_i(n_i);
  for(int i=0; i<n_i; i++){
    logchat_i(i) = d_ktp(s_i(i),t_i(i),p_i(i));
    if( !isNA(c_i(i)) ){                
      jnll -= dpois( c_i(i), exp(logchat_i(i)), true );
    }
  }

  // Spatial field summaries
  REPORT( Range );
  // Fields
  REPORT( A_kp );
  REPORT( d_ktp );
  REPORT( dhat_ktp );
  REPORT( jnll );
  
  return jnll;
}
