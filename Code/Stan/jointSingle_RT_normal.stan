data {
	int N;      // # of subjects
	int N_cond; // # of conditions
	int N_time; // # of timepoints
	int T_max;  // max # of trials across subjects
	int T_subj[N, N_cond, N_time];  // # of trials within subjects, conditions, timepoints
	real RT[N, N_cond, N_time, T_max]; // Reaction times for each subject, condition, timepoint, and trial
}
parameters {
  // Group-level correlation matrix (cholesky factor for faster computation)
  cholesky_factor_corr[4] L_R_mu; 
  cholesky_factor_corr[4] L_R_sigma;
  
  // Group-level parameter means
  vector[4] mu_mean;     
  vector[4] sigma_mean;
  
  // Group-level parameter SDs
  vector<lower=0>[4] mu_sd;
  vector<lower=0>[4] sigma_sd;
  
  // Individual-level parameters (before being transformed)
  matrix[4,N] mu_i_pr;   
	matrix[4,N] sigma_i_pr;
}
transformed parameters {
  // Individual-level parameter off-sets (for non-centered parameterization)
  matrix[4,N] mu_i_tilde;
  matrix[4,N] sigma_i_tilde;
  
  // Individual-level parameters 
  matrix[N,4] mu_i;
  matrix[N,4] sigma_i;
  
  // Construct inidividual offsets (for non-centered parameterization)
  mu_i_tilde = diag_pre_multiply(mu_sd, L_R_mu) * mu_i_pr;
  sigma_i_tilde = diag_pre_multiply(sigma_sd, L_R_sigma) * sigma_i_pr; 
  
  // Compute individual-level parameters from non-centered parameterization
  for (i in 1:N) {
    // Congruent at time 1
    mu_i[i,1] = mu_mean[1] + mu_i_tilde[1,i];
    // Congruent at time 2
    mu_i[i,3] = mu_mean[3] + mu_i_tilde[3,i];
    
    // Congruent at time 1
    sigma_i[i,1] = exp(sigma_mean[1] + sigma_i_tilde[1,i]);
    // Congruent at time 2
    sigma_i[i,3] = exp(sigma_mean[3] + sigma_i_tilde[3,i]);
    
    // Inongruent at time 1
    mu_i[i,2] = mu_mean[2] + mu_i_tilde[2,i];
    // Inongruent at time 2
    mu_i[i,4] = mu_mean[4] + mu_i_tilde[4,i];
    
    // Inongruent at time 1
    sigma_i[i,2] = exp(sigma_mean[2] + sigma_i_tilde[2,i]);
    // Inongruent at time 2
    sigma_i[i,4] = exp(sigma_mean[4] + sigma_i_tilde[4,i]);
  }
}
model {
  // Prior on cholesky factor of correlation matrix
  L_R_mu    ~ lkj_corr_cholesky(1);
  L_R_sigma ~ lkj_corr_cholesky(1); 
  
  // Priors on group-level means 
  mu_mean     ~ normal(0, 1);
  sigma_mean  ~ normal(0, 1);
  
  // Priors on group-level SDs
  mu_sd     ~ normal(0, 1);
  sigma_sd  ~ normal(0, 1);
  
  // Priors on individual-level parameters
  to_vector(mu_i_pr)    ~ normal(0, 1);
  to_vector(sigma_i_pr) ~ normal(0, 1);
  
  // For each subject
  for (i in 1:N) {
    // Congruent at time 1
    RT[i,1,1,1:T_subj[i,1,1]] ~ normal(mu_i[i,1], sigma_i[i,1]);
    // Incongruent at time 1
    RT[i,2,1,1:T_subj[i,2,1]] ~ normal(mu_i[i,2], sigma_i[i,2]);
    // Congruent at time 2
    RT[i,1,2,1:T_subj[i,1,2]] ~ normal(mu_i[i,3], sigma_i[i,3]);
    // Incongruent at time 2
    RT[i,2,2,1:T_subj[i,2,2]] ~ normal(mu_i[i,4], sigma_i[i,4]);
  }
}
generated quantities { 
  // test-retest correlations
  corr_matrix[4] R_mu;
  corr_matrix[4] R_sigma;
  
  // posterior predictions and log-likelihood
  real post_pred_c1_t1[N, T_max];
  real post_pred_c1_t2[N, T_max];
  real post_pred_c2_t1[N, T_max];
  real post_pred_c2_t2[N, T_max];
  real log_lik[N, N_cond, N_time, T_max];
  
	// Reconstruct correlation matrix from cholesky factor
  R_mu = L_R_mu * L_R_mu';
  R_sigma = L_R_sigma * L_R_sigma';
  
  // initialize LL and post_pred arrays to -1
  post_pred_c1_t1 = rep_array(-1.0, N, T_max);
  post_pred_c1_t2 = rep_array(-1.0, N, T_max);  
  post_pred_c2_t1 = rep_array(-1.0, N, T_max);  
  post_pred_c2_t2 = rep_array(-1.0, N, T_max);
  for (i in 1:N) {
    log_lik[i,,,] = rep_array(-1.0, N_cond, N_time, T_max);
  }
  
  // Generate posterior predictions
  for (i in 1:N) {
    // Congruent at time 1
    for (t in 1:T_subj[i,1,1]) {
      post_pred_c1_t1[i,t] = normal_rng(mu_i[i,1], sigma_i[i,1]);
      log_lik[i,1,1,t] = normal_lpdf(RT[i,1,1,t] | mu_i[i,1], sigma_i[i,1]);
    }
    // Incongruent at time 1
    for (t in 1:T_subj[i,2,1]) {
      post_pred_c2_t1[i,t] = normal_rng(mu_i[i,2], sigma_i[i,2]);
      log_lik[i,2,1,t] = normal_lpdf(RT[i,2,1,t] | mu_i[i,2], sigma_i[i,2]);
    }
    // Congruent at time 2
    for (t in 1:T_subj[i,1,2]) {
      post_pred_c1_t2[i,t] = normal_rng(mu_i[i,3], sigma_i[i,3]);
      log_lik[i,1,2,t] = normal_lpdf(RT[i,1,2,t] | mu_i[i,3], sigma_i[i,3]);
    }
    // Incongruent at time 2
    for (t in 1:T_subj[i,2,2]) {
      post_pred_c2_t2[i,t] = normal_rng(mu_i[i,4], sigma_i[i,4]);
      log_lik[i,2,2,t] = normal_lpdf(RT[i,2,2,t] | mu_i[i,4], sigma_i[i,4]);
    }
  }
} 
