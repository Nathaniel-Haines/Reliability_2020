library(rstan)
library(dplyr)
library(bayesplot)
library(foreach)
library(ggplot2)
library(patchwork)

setwd("~/Dropbox/Box/GitHub/Reliability_2020/")

# Tasks and parameters within each model for plotting
tasks <- c("Study1-Stroop_normal", "Study1-Stroop_lognormal", 
           "Study1-Stroop_shift_lognormal")
pars  <- c("mu_i_delta", "sigma_mean_delta", "sigma_i_delta_tilde")

eff_dat <- foreach(i=seq_along(tasks), .combine = "rbind") %do% {
  # Extract generative model parameters
  post_pars <- rstan::extract(readRDS(paste0("Data/2_Fitted/fit_", 
                                             tasks[i], ".rds")), 
                              pars = pars)
  
  mu_delta  <- post_pars$mu_i_delta[,,1]
  sig_delta <- post_pars$sigma_mean_delta[,1] + post_pars$sigma_i_delta_tilde[,1,]
  
  data.frame(subj      = 1:dim(mu_delta)[2],
             model     = rep(tasks[i], dim(mu_delta)[2]),
             pr_mu_g0  = apply(mu_delta, 2, function(x) mean(x>0)),
             pr_sig_g0 = apply(sig_delta, 2, function(x) mean(x>0))) 
}

# Determine proportion of individual-level difference parameters 
# with 95% mass above 0
mean(eff_dat$pr_mu_g0[eff_dat$model=="Study1-Stroop_normal"]>=.95)
mean(eff_dat$pr_mu_g0[eff_dat$model=="Study1-Stroop_lognormal"]>=.95)
mean(eff_dat$pr_mu_g0[eff_dat$model=="Study1-Stroop_shift_lognormal"]>=.95)
mean(eff_dat$pr_sig_g0[eff_dat$model=="Study1-Stroop_normal"]>=.95)
mean(eff_dat$pr_sig_g0[eff_dat$model=="Study1-Stroop_lognormal"]>=.95)
mean(eff_dat$pr_sig_g0[eff_dat$model=="Study1-Stroop_shift_lognormal"]>=.95)
  
