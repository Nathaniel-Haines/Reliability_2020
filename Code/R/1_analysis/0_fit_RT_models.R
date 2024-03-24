library(cmdstanr)
library(foreach)
library(posterior)

# Compile models
m0 <- cmdstan_model("Code/Stan/joint_RT_normal.stan")
m1 <- cmdstan_model("Code/Stan/joint_RT_lognormal.stan")
m2 <- cmdstan_model("Code/Stan/joint_RT_shiftlnorm.stan")
m3 <- cmdstan_model("Code/Stan/joint_RT_shiftlnorm_mix.stan")

# Load preprocessed data
stan_data <- readRDS("Data/1_Preprocessed/stan_ready_all.rds")

# Lists for parellel model fitting 
models <- list(m0, m1, m2, m3)
model_names <- c("normal", "lognormal", "shift_lognormal", "shift_lognormal_mix")
data_names <- c("Study1-Flanker", "Study2-Flanker", "Study1-Stroop",
                "Study2-Stroop", "Study3-Posner", "Study1a-IAT",
                "Study2b-IAT")

# Iterate through models and tasks
results <- foreach(m=seq_along(model_names), .combine = "c") %do% {
  foreach(d=data_names, .combine = "c") %do% {
    if (model_names[m]=="shift_lognormal" & d=="Study3-Posner") {
      tmp_seed  <- 1
      tmp_adapt <- .9
    } else {
      tmp_seed  <- 43201
      tmp_adapt <- .8
    }
    fit <- models[[m]]$sample( 
      stan_data[[d]], 
      iter_sampling = 750, 
      iter_warmup = 250, 
      chains = 8,
      parallel_chains = 8,
      seed = tmp_seed,
      adapt_delta = tmp_adapt,
      init = 0
    )
    rvars_pars <- as_draws_rvars(
      fit$draws(
        c(
          "R_mu", "R_sigma", 
          "mu_i_delta", "sigma_i_delta", "sigma_i_base",
          "post_pred_c1_t1", "post_pred_c1_t2", 
          "post_pred_c2_t1", "post_pred_c2_t2"
        )
      )
    )
    pars <- lapply(rvars_pars, draws_of)
    saveRDS(pars, file = paste0("Data/2_Fitted/fit_", d, "_", model_names[m],".rds"))
    rm(fit, pars); gc()
    model_names[m]
  }
}
