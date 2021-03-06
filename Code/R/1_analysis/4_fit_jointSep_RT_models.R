library(rstan)
library(foreach)

# Compile models
m0 <- stan_model("Code/Stan/jointSep_RT_normal.stan")
m1 <- stan_model("Code/Stan/jointSep_RT_lognormal.stan")
m2 <- stan_model("Code/Stan/jointSep_RT_shiftlnorm.stan")

# Load preprocessed data
stan_data <- readRDS("Data/1_Preprocessed/stan_ready_all.rds")

# Lists for parellel model fitting 
models <- list(m0, m1, m2)
model_names <- c("normal", "lognormal", "shift_lognormal")
data_names <- c("Study1-Flanker", "Study2-Flanker", "Study1-Stroop", 
                "Study2-Stroop", "Study3-Posner", "Study1a-IAT", 
                "Study2b-IAT")

# Random seed and adapt_delta
tmp_seed  <- 1
tmp_adapt <- .8

# Iterate through models and tasks
results <- foreach(m=seq_along(model_names), .combine = "c") %do% {
  foreach(d=data_names, .combine = "c") %do% {
    if (d=="Study2b-IAT" & m==1) {
      tmp_seed  <- 2
      tmp_adapt <- .8
    } else {
      tmp_seed  <- 1
      tmp_adapt <- .8
    }
    fit <- sampling(models[[m]], 
                    data   = stan_data[[d]], 
                    iter   = 3000, 
                    warmup = 1000, 
                    chains = 3, 
                    cores  = 3,
                    seed   = tmp_seed,
                    control = list(adapt_delta = tmp_adapt))
    saveRDS(fit, file = paste0("Data/2_Fitted/fit_", d, "_jointSep_", model_names[m],".rds"))
    rm(fit)
    gc()
    model_names[m]
  }
}
