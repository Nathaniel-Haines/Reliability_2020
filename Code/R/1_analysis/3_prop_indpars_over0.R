library(rstan)
library(dplyr)
library(foreach)
library(hBayesDM)

setwd("~/Dropbox/Box/GitHub/Reliability_2020/")

# Task, model names, and parameters
models <- c("Study1-Stroop_normal", "Study1-Stroop_lognormal", 
           "Study1-Stroop_shift_lognormal")
pars  <- c("mu_i_delta", "sigma_i_delta")

# Function that determines if the 95% HDI is above 0
hdi_inc0 <- function(x) {
  hdi <- HDIofMCMC(x)
  return(hdi[1] > 0)
}

# Loop through, determine how many people show the "Stroop effect"
eff_dat <- foreach(i=seq_along(models), .combine = "rbind") %do% {
  # Extract generative model parameters
  post_pars <- rstan::extract(readRDS(paste0("Data/2_Fitted/fit_", 
                                             models[i], ".rds")), 
                              pars = pars)
  
  # Person-level parameters
  mu_delta  <- post_pars$mu_i_delta[,,1]
  sig_delta <- post_pars$sigma_i_delta[,,1]
  
  # Save out data
  data.frame(subj      = 1:dim(mu_delta)[2],
             model     = rep(models[i], dim(mu_delta)[2]),
             pr_mu_g0  = apply(mu_delta, 2, hdi_inc0),
             pr_sig_g0 = apply(sig_delta, 2, hdi_inc0)) 
}

# Determine number of individual-level difference parameters 
# with 95% HDIs above 0
sum(eff_dat$pr_mu_g0[eff_dat$model=="Study1-Stroop_normal"])
sum(eff_dat$pr_mu_g0[eff_dat$model=="Study1-Stroop_lognormal"])
sum(eff_dat$pr_mu_g0[eff_dat$model=="Study1-Stroop_shift_lognormal"])
sum(eff_dat$pr_sig_g0[eff_dat$model=="Study1-Stroop_normal"])
sum(eff_dat$pr_sig_g0[eff_dat$model=="Study1-Stroop_lognormal"])
sum(eff_dat$pr_sig_g0[eff_dat$model=="Study1-Stroop_shift_lognormal"])
