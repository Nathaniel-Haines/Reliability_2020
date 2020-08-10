data {
  int N;
  int N_time; // # of timepoints
  int T_max;
  int T_subj[N,N_time];
  real delay_later[N,N_time,T_max];
  real amount_later[N,N_time,T_max];
  real delay_sooner[N,N_time,T_max];
  real amount_sooner[N,N_time,T_max];
  int choice[N,N_time,T_max]; // 0 for instant reward, 1 for delayed reward
}
parameters {
  // Group-level correlation matrix (cholesky factor for faster computation)
  cholesky_factor_corr[2] R_chol_k; 
  cholesky_factor_corr[2] R_chol_c;
  
  // Group-level parameter SDs
  vector<lower=0>[2] sigma_k;
  vector<lower=0>[2] sigma_c; 
  
  // Group-level parameter means
  vector[2] mu_k;
  vector[2] mu_c;
  
  // Individual-level parameters (before being transformed)
  matrix[2,N] k_pr; 
	matrix[2,N] c_pr;  
}
transformed parameters {
  // Individual-level parameter off-sets (for non-centered parameterization)
  matrix[2,N] k_tilde;
  matrix[2,N] c_tilde;
  
  // Individual-level parameters 
  matrix[N,2] k;
  matrix[N,2] c;
  
  // Construct inidividual offsets (for non-centered parameterization)
  k_tilde = diag_pre_multiply(sigma_k, R_chol_k) * k_pr;
  c_tilde = diag_pre_multiply(sigma_c, R_chol_c) * c_pr; 
  
  // Compute individual-level parameters from non-centered parameterization
  for (i in 1:N) {
    // Discounting rate (k) at time 1
    k[i,1] = exp(mu_k[1] + k_tilde[1,i]);
    // Discounting rate (k) at time 2
    k[i,2] = exp(mu_k[2] + k_tilde[2,i]);
    
    // Choice sensitivity/inverse temp. (c) at time 1
    c[i,1] = exp(mu_c[1] + c_tilde[1,i]);
    // Choice sensitivity/inverse temp. (c) at time 2
    c[i,2] = exp(mu_c[2] + c_tilde[2,i]);
  }
}
model {
  // Prior on cholesky factor of correlation (i.e. test-retest) matrix
  R_chol_k ~ lkj_corr_cholesky(1);
  R_chol_c ~ lkj_corr_cholesky(1); 
  
  // Priors on group-level means
  mu_k ~ normal(0, 1);
  mu_c ~ normal(0, 1);
  
  // Priors on group-level SDs
  sigma_k ~ normal(0, 0.2);
  sigma_c ~ normal(0, 0.2);
  
  // Priors on individual-level parameters
  to_vector(k_pr) ~ normal(0, 1);
  to_vector(c_pr) ~ normal(0, 1);
  
  // For each subject
  for (i in 1:N) {
    // Define values
    vector[T_subj[i,1]] ev_later_t1;
    vector[T_subj[i,1]] ev_sooner_t1;
    vector[T_subj[i,2]] ev_later_t2;
    vector[T_subj[i,2]] ev_sooner_t2;
  
    // vectorized likelihood for time 1
    ev_later_t1   = to_vector(amount_later[i,1,:])  ./ ( 1 + k[i,1] * to_vector(delay_later[i,1,:]));
    ev_sooner_t1  = to_vector(amount_sooner[i,1,:]) ./ ( 1 + k[i,1] * to_vector(delay_sooner[i,1,:]));
    choice[i,1,:] ~ bernoulli_logit( c[i,1] * (ev_later_t1 - ev_sooner_t1) );
    
    // vectorized likelihood for time 2
    ev_later_t2   = to_vector(amount_later[i,2,:])  ./ ( 1 + k[i,2] * to_vector(delay_later[i,2,:]));
    ev_sooner_t2  = to_vector(amount_sooner[i,2,:]) ./ ( 1 + k[i,2] * to_vector(delay_sooner[i,2,:]));
    choice[i,2,:] ~ bernoulli_logit( c[i,2] * (ev_later_t2 - ev_sooner_t2) );
  }
}
generated quantities { 
  // test-retest correlations
  corr_matrix[2] R_k;
  corr_matrix[2] R_c;

  // posterior predictions and log-likelihood
  real post_pred_t1[N, T_max];
  real post_pred_t2[N, T_max];
  // real log_lik[N, N_time, T_max];

	// Reconstruct correlation matrix from cholesky factor
  R_k = R_chol_k * R_chol_k';
  R_c = R_chol_c * R_chol_c';
  
  {
    // For each subject
    for (i in 1:N) {
      // Define values
      vector[T_subj[i,1]] ev_later_t1;
      vector[T_subj[i,1]] ev_sooner_t1;
      vector[T_subj[i,2]] ev_later_t2;
      vector[T_subj[i,2]] ev_sooner_t2;
    
      // vectorized likelihood for time 1
      ev_later_t1   = to_vector(amount_later[i,1,:])  ./ ( 1 + k[i,1] * to_vector(delay_later[i,1,:]));
      ev_sooner_t1  = to_vector(amount_sooner[i,1,:]) ./ ( 1 + k[i,1] * to_vector(delay_sooner[i,1,:]));
      post_pred_t1[i,:] = bernoulli_rng(inv_logit( c[i,1] * (ev_later_t1 - ev_sooner_t1)));
      
      // vectorized likelihood for time 2
      ev_later_t2   = to_vector(amount_later[i,2,:])  ./ ( 1 + k[i,2] * to_vector(delay_later[i,2,:]));
      ev_sooner_t2  = to_vector(amount_sooner[i,2,:]) ./ ( 1 + k[i,2] * to_vector(delay_sooner[i,2,:]));
      post_pred_t2[i,:] = bernoulli_rng(inv_logit( c[i,2] * (ev_later_t2 - ev_sooner_t2)));
    }
  }
} 
