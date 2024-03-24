library(cmdstanr)
library(dplyr)
library(foreach)
library(posterior)

## Generative model
stan_data <- readRDS("Data/1_Preprocessed/stan_ready_all.rds")
m0 <- cmdstan_model("Code/Stan/joint_DDT_hyperbolic.stan")

fit <- m0$sample( 
  stan_data$`Study1-DDT`, 
  iter_sampling = 750, 
  iter_warmup = 250, 
  chains = 8,
  parallel_chains = 8,
  seed = 43201
)

rvars_pars <- as_draws_rvars(
  fit$draws(
    c("R_k", "R_c", "k", "c")
  )
)
pars <- lapply(rvars_pars, draws_of)

# Save out fit
saveRDS(pars, file = "Data/2_Fitted/fit_Study1-DDT_hyperbolic.rds")

## Two-stage approach
# Softmax function
logsumexp <- function (x) {
  y <- max(x)
  y + log(sum(exp(x - y)))
}
softmax <- function (x) {
  exp(x - logsumexp(x))
}

# Define the log-likelihood function used for MLE
mle_ddt <- function(par, data, subj, time) {
  # number of trials for current run
  n_trials <- data$T_subj[subj, time]
  # loop through each trial and compute log-likelihood
  ll <- foreach(t=1:n_trials, .combine = "c") %do% {
    # discounting function
    ev_later  <- data$amount_later[subj,time,t] / ( 1 + par[1] * data$delay_later[subj,time,t])
    ev_sooner <- data$amount_sooner[subj,time,t] / ( 1 + par[1] * data$delay_sooner[subj,time,t])
    ev <- c(ev_sooner, ev_later)
    
    # Generate choice probability with softmax
    pr <- softmax(par[2] * ev);
    # Perturb for non-infinite logs
    pr[pr==1] <- .999999
    pr[pr==0] <- .000001
    
    # log probability of "true" simulated choice
    log(pr[(data$choice[subj,time,t]+1)])
  }
  # return the summed (minus) log-likelihood, because optim minimizes by default
  sum(-1*ll)
}

# Get data
ddt_dat <- stan_data$`Study1-DDT`

# Estimate MLE for k and c
mle_results <- foreach(subj=1:ddt_dat$N, .combine = "rbind") %do% {
  # set seed
  set.seed(43202)
  
  # Use optim to minimize the (minus) log-likelihood function
  # Time 1
  mle_t1 <- optim(
    par    = c(0.05, 0.5),    # Initial guess for k and c
    fn     = mle_ddt,         # Function we are minimizing
    method = "L-BFGS-B",      # Specific algorithm used
    lower  = c(0.0001, 0.0001),         # Lower bounds for k and c 
    upper  = c(1, 1.5),       # Upper bound for k and c 
    data   = ddt_dat,         # choice/task data
    subj   = subj,            # current subject 
    time   = 1                # timepoint
  )               
  # Time 2
  mle_t2 <- optim(
    par    = c(0.05, 0.5),    # Initial guess for k and c
    fn     = mle_ddt,         # Function we are minimizing
    method = "L-BFGS-B",      # Specific algorithm used
    lower  = c(0.0001, 0.0001),         # Lower bounds for k and c 
    upper  = c(1, 1.5),       # Upper bound for k and c 
    data   = ddt_dat,         # choice/task data
    subj   = subj,            # current subject 
    time   = 2                # timepoint
  )               
  
  # Return estimates
  data.frame(
    subj_num = rep(subj, 2),
    Time = 1:2,
    k = c(mle_t1$par[1], mle_t2$par[1]),
    c = c(mle_t1$par[2], mle_t2$par[2])
  )
}

# Save out MLE fit
saveRDS(mle_results, file = "Data/2_Fitted/fit_Study1-DDT_hyperbolic_MLE.rds")
