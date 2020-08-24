library(mvtnorm)
library(rstan)
library(dplyr)
library(foreach)
library(patchwork)

setwd("~/Dropbox/Box/GitHub/Reliability_2020/")

# For transforming long data to stan-ready format
source("Code/R/0_Preprocessing/utils.R")

# Stan model for recovery
m1 <- stan_model("Code/Stan/joint_RT_lognormal.stan")

# Simulation parameters
set.seed(43202)
n_subj   <- c(10, 50, 50, 100, 100)
n_trials <- c(10, 50, 100, 50, 100)

## Group-level generative parameters
# test-retest
mu_cor <- seq(-1, 1, length.out = 15)
sd_cor <- seq(-1, 1, length.out = 15)
con_cor <- 0

# Means
mu_mean_base  <- c(-.5, -.55)
mu_mean_delta <- c(.1, .08)
mu_sd_base    <- c(.11, .12)
mu_sd_delta   <- c(.035, .029)
# Standard deviations
sigma_mean_base  <- c(-1.28, -1.32)
sigma_mean_delta <- c(.12, .11)
sigma_sd_base    <- c(.19, .2)
sigma_sd_delta   <- c(.11, .09)

results <- foreach(i=seq_along(n_subj), .combine = "rbind") %do% {
  foreach(r=seq_along(mu_cor), .combine = "rbind") %do% {
    ## Individual-level parameters
    # For mean when c = 1
    mu_c1 <- rmvnorm(n_subj[i], mu_mean_base, 
                     diag(mu_sd_base) %*%
                       matrix(c(1,con_cor,con_cor,1), nrow = 2) %*% 
                       diag(mu_sd_base))
    # For sd when c = 1
    sd_c1 <- rmvnorm(n_subj[i], sigma_mean_base, 
                     diag(sigma_sd_base) %*%
                       matrix(c(1,con_cor,con_cor,1), nrow = 2) %*% 
                       diag(sigma_sd_base)) %>%
      exp(.)
    
    # For mean when c = 2
    mu_delta <- rmvnorm(n_subj[i], mu_mean_delta, 
                        diag(mu_sd_delta) %*%
                          matrix(c(1,mu_cor[r],mu_cor[r],1), nrow = 2) %*% 
                          diag(mu_sd_delta))
    mu_c2 <- mu_c1 + mu_delta
    # For sd when c = 2
    sd_delta <- rmvnorm(n_subj[i], sigma_mean_delta, 
                        diag(sigma_sd_delta) %*%
                          matrix(c(1,sd_cor[r],sd_cor[r],1), nrow = 2) %*% 
                          diag(sigma_sd_delta))
    sd_c2 <- exp(log(sd_c1) + sd_delta)
    
    ## Loop through subjects and save results in long format
    sim_dat <- foreach(t=1:n_subj[i], .combine = "rbind") %do% {
      data.frame(subj_num  = rep(t, n_trials[i]*4),
                 Time      = c(rep(1, n_trials[i]*2),
                               rep(2, n_trials[i]*2)),
                 Condition = c(rep(0, n_trials[i]),
                               rep(2, n_trials[i]),
                               rep(0, n_trials[i]),
                               rep(2, n_trials[i])),
                 RT        = c(rlnorm(n_trials[i], mu_c1[t,1], sd_c1[t,1]),
                               rlnorm(n_trials[i], mu_c2[t,1], sd_c2[t,1]),
                               rlnorm(n_trials[i], mu_c1[t,2], sd_c1[t,2]),
                               rlnorm(n_trials[i], mu_c2[t,2], sd_c2[t,2])),
                 Correct   = rep(1, n_trials[i]*4))
    }
    
    # For storing results in stan-ready format
    stan_data <- pre_stp_flk_psr(sim_dat)
    
    # Fit model for parameter recovery
    fit <- sampling(m1, 
                    data   = stan_data, 
                    iter   = 1000, 
                    warmup = 250, 
                    chains = 2, 
                    cores  = 2, 
                    pars   = c("R_mu", "R_sigma"),
                    seed   = 43201)
    
    # Sample mean/sd cor
    wide_data <- sim_dat %>% 
      group_by(subj_num, Time) %>%
      summarize(mu_diff = mean(RT[Condition==2]) - mean(RT[Condition==0]),
                sd_diff = sd(RT[Condition==2]) - sd(RT[Condition==0])) %>%
      tidyr::pivot_wider(names_from = Time, values_from = mu_diff:sd_diff, names_prefix = "T")
    mu_cor_est <- cor.test(wide_data$mu_diff_T1, wide_data$mu_diff_T2)
    sd_cor_est <- cor.test(wide_data$sd_diff_T1, wide_data$sd_diff_T2)
    
    # Generative model estimates 
    pars <- rstan::extract(fit, pars = c("R_mu", "R_sigma"))
    mu_cor_mean <- mean(pars$R_mu[,1,2]) 
    mu_cor_hdi  <- hBayesDM::HDIofMCMC(pars$R_mu[,1,2])
    sd_cor_mean <- mean(pars$R_sigma[,1,2]) 
    sd_cor_hdi  <- hBayesDM::HDIofMCMC(pars$R_sigma[,1,2])
    
    # Return data in long format
    data.frame(n_subj    = rep(n_subj[i], 4),
               n_trials  = rep(n_trials[i], 4),
               samp_size = rep(paste0("Participants: ", n_subj[i], "\n",
                                      "Trials: ", n_trials[i]), 4),
               type      = c(rep("Two-Stage", 2),
                             rep("Generative", 2)),
               par       = rep(c("R_mu", "R_sigma"), 2),
               R_true    = rep(c(mu_cor[r], sd_cor[r]), 2),
               R_mu      = c(mu_cor_est$estimate, sd_cor_est$estimate,
                             mu_cor_mean, sd_cor_mean),
               R_lo      = c(mu_cor_est$conf.int[1], sd_cor_est$conf.int[1],
                             mu_cor_hdi[1], sd_cor_hdi[1]),
               R_hi      = c(mu_cor_est$conf.int[2], sd_cor_est$conf.int[2],
                             mu_cor_hdi[2], sd_cor_hdi[2]))
  }
}

saveRDS(results, file = "Data/2_Fitted/parameter_recovery.rds")
