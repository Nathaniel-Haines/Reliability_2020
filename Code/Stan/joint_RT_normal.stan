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
  cholesky_factor_corr[2] L_R_mu; 
  cholesky_factor_corr[2] L_R_sigma;
  
  // Group-level parameter means
  vector[2] mu_mean_base;
  vector[2] mu_mean_delta;        
  vector[2] sigma_mean_base;
  vector[2] sigma_mean_delta;    
  
  // Group-level parameter SDs
  vector<lower=0>[2] mu_sd_base;
  vector<lower=0>[2] mu_sd_delta;
  vector<lower=0>[2] sigma_sd_base;
  vector<lower=0>[2] sigma_sd_delta; 
  
  // Individual-level parameters (before being transformed)
  matrix[N,2] mu_i_base_pr; 
  matrix[2,N] mu_i_delta_pr;  
  matrix[N,2] sigma_i_base_pr; 
	matrix[2,N] sigma_i_delta_pr;
}
transformed parameters {
  // Individual-level parameter off-sets (for non-centered parameterization)
  matrix[2,N] mu_i_delta_tilde;
  matrix[2,N] sigma_i_delta_tilde;
  
  // Individual-level parameters 
  matrix[N,2] mu_i_base;
  matrix[N,2] mu_i_delta;
  matrix[N,2] sigma_i_base;
  matrix[N,2] sigma_i_delta;
  
  // Construct inidividual offsets (for non-centered parameterization)
  mu_i_delta_tilde = diag_pre_multiply(mu_sd_delta, L_R_mu) * mu_i_delta_pr;
  sigma_i_delta_tilde = diag_pre_multiply(sigma_sd_delta, L_R_sigma) * sigma_i_delta_pr; 
  
  // Compute individual-level parameters from non-centered parameterization
  for (i in 1:N) {
    // Congruent at time 1
    mu_i_base[i,1] = mu_mean_base[1] + mu_sd_base[1] * mu_i_base_pr[i,1];
    // Congruent at time 2
    mu_i_base[i,2] = mu_mean_base[2] + mu_sd_base[2] * mu_i_base_pr[i,2];
    
    // Congruent at time 1
    sigma_i_base[i,1] = sigma_mean_base[1] + sigma_sd_base[1] * sigma_i_base_pr[i,1];
    // Congruent at time 2
    sigma_i_base[i,2] = sigma_mean_base[2] + sigma_sd_base[2] * sigma_i_base_pr[i,2];
    
    // Condition effect on mu at time 1
    mu_i_delta[i,1] = mu_mean_delta[1] + mu_i_delta_tilde[1,i];
    // Condition effect on mu at time 2
    mu_i_delta[i,2] = mu_mean_delta[2] + mu_i_delta_tilde[2,i];
    
    // Condition effect on SD at time 1
    sigma_i_delta[i,1] = sigma_mean_delta[1] + sigma_i_delta_tilde[1,i];
    // Condition effect on SD at time 2
    sigma_i_delta[i,2] = sigma_mean_delta[2] + sigma_i_delta_tilde[2,i];
  }
}
model {
  // Prior on cholesky factor of correlation matrix
  L_R_mu    ~ lkj_corr_cholesky(1);
  L_R_sigma ~ lkj_corr_cholesky(1); 
  
  // Priors on group-level means 
  mu_mean_base     ~ normal(0, 1);
  mu_mean_delta    ~ normal(0, 1);
  sigma_mean_base  ~ normal(0, 1);
  sigma_mean_delta ~ normal(0, 1); 
  
  // Priors on group-level SDs
  mu_sd_base     ~ normal(0, 1);
  mu_sd_delta    ~ normal(0, 1);
  sigma_sd_base  ~ normal(0, 1);
  sigma_sd_delta ~ normal(0, 1);
  
  // Priors on individual-level parameters
  to_vector(mu_i_base_pr)     ~ normal(0, 1);
  to_vector(mu_i_delta_pr)    ~ normal(0, 1);
  to_vector(sigma_i_base_pr)  ~ normal(0, 1); 
  to_vector(sigma_i_delta_pr) ~ normal(0, 1); 
  
  // For each subject
  for (i in 1:N) {
    // Congruent at time 1
    RT[i,1,1,1:T_subj[i,1,1]] ~ normal(mu_i_base[i,1], 
                                       exp(sigma_i_base[i,1]));
    // Incongruent at time 1
    RT[i,2,1,1:T_subj[i,2,1]] ~ normal(mu_i_base[i,1] + mu_i_delta[i,1], 
                                       exp(sigma_i_base[i,1] + sigma_i_delta[i,1]));
    // Congruent at time 2
    RT[i,1,2,1:T_subj[i,1,2]] ~ normal(mu_i_base[i,2], 
                                       exp(sigma_i_base[i,2]));
    // Incongruent at time 2
    RT[i,2,2,1:T_subj[i,2,2]] ~ normal(mu_i_base[i,2] + mu_i_delta[i,2], 
                                       exp(sigma_i_base[i,2] + sigma_i_delta[i,2]));
  }
}
generated quantities { 
  // test-retest correlations
  corr_matrix[2] R_mu;
  corr_matrix[2] R_sigma;
  
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
      post_pred_c1_t1[i,t] = normal_rng(mu_i_base[i,1], 
                                        exp(sigma_i_base[i,1]));
      log_lik[i,1,1,t] = normal_lpdf(RT[i,1,1,t] | 
                                     mu_i_base[i,1], 
                                     exp(sigma_i_base[i,1]));
    }
    // Incongruent at time 1
    for (t in 1:T_subj[i,2,1]) {
      post_pred_c2_t1[i,t] = normal_rng(mu_i_base[i,1] + mu_i_delta[i,1], 
                                        exp(sigma_i_base[i,1] + sigma_i_delta[i,1]));
      log_lik[i,2,1,t] = normal_lpdf(RT[i,2,1,t] | 
                                     mu_i_base[i,1] + mu_i_delta[i,1], 
                                     exp(sigma_i_base[i,1] + sigma_i_delta[i,1]));
    }
    // Congruent at time 2
    for (t in 1:T_subj[i,1,2]) {
      post_pred_c1_t2[i,t] = normal_rng(mu_i_base[i,2], 
                                        exp(sigma_i_base[i,2]));
      log_lik[i,1,2,t] = normal_lpdf(RT[i,1,2,t] | 
                                     mu_i_base[i,2], 
                                     exp(sigma_i_base[i,2]));
    }
    // Incongruent at time 2
    for (t in 1:T_subj[i,2,2]) {
      post_pred_c2_t2[i,t] = normal_rng(mu_i_base[i,2] + mu_i_delta[i,2], 
                                        exp(sigma_i_base[i,2] + sigma_i_delta[i,2]));
      log_lik[i,2,2,t] = normal_lpdf(RT[i,2,2,t] | 
                                     mu_i_base[i,2] + mu_i_delta[i,2], 
                                     exp(sigma_i_base[i,2] + sigma_i_delta[i,2]));
    }
  }
} 
