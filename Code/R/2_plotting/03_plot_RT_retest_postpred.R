library(cowplot)
library(patchwork)
library(bayesplot)
library(ggplot2)
library(foreach)

source("Code/R/utils/plot_utils.R")

stan_data <- readRDS("Data/1_Preprocessed/stan_ready_all.rds")

fit_names <- c("normal", "lognormal", "shift_lognormal", "shift_lognormal_mix")
data_names <- c("Study1-Stroop", "Study2-Stroop", 
                "Study1-Flanker", "Study2-Flanker", 
                "Study3-Posner", 
                "Study1a-IAT", "Study2b-IAT")

xlim_map <- function(dataset) {
  return(
    switch(
      dataset,
      "Study1-Stroop" = c(0, 2.5),
      "Study2-Stroop" = c(0, 2.5),
      "Study1-Flanker" = c(0, 1.5),
      "Study2-Flanker" = c(0, 1.5),
      "Study3-Posner" = c(0, 1),
      "Study1a-IAT" = c(0, 4),
      "Study2b-IAT" = c(0, 4)
    )
  )
} 

# Results
results <- foreach(d=data_names) %do% {
  samp_data <- stan_data[[d]]
  n_subj <- samp_data$N
  
  set.seed(43202)  
  subjs <- sample(1:n_subj, 3, replace = F)
  
  # Read in all model fits for plotting
  fits <- list()
  for (i in fit_names) {
    # Read in data 
    fit <- readRDS(paste0("Data/2_Fitted/fit_", d, "_", i, ".rds"))
    
    # Only save random 3 subjects for plotting (to avoid RAM problems)
    fit$post_pred_c1_t1 <- fit$post_pred_c1_t1[,subjs,]
    fit$post_pred_c1_t2 <- fit$post_pred_c1_t2[,subjs,]
    fit$post_pred_c2_t1 <- fit$post_pred_c2_t1[,subjs,]
    fit$post_pred_c2_t2 <- fit$post_pred_c2_t2[,subjs,]
    
    # Return model fit
    fits[[i]] <- fit
    rm(fit)
  }
  
  # Test-retest plots
  mu_plot <- plot_retest(fits, parameter = "mu", 
                         task = paste0(d, " Strength"), 
                         samp_data = samp_data, legend = "top")
  sig_plot <- plot_retest(fits, parameter = "sig", 
                          task = paste0(d, " Variability"), 
                          samp_data = samp_data, legend = "none")
  
  tr_plot <- (mu_plot / sig_plot)
  
  # Posterior predictions for RTs
  subj <- subjs[1]
  p1 <- plot_RT(pars = fits[[1]], raw = samp_data, subj = subj, 
                subjs = subjs, xlim = xlim_map(d),
                legend = "none")
  p2 <- plot_RT(pars = fits[[2]], raw = samp_data, subj = subj, 
                subjs = subjs, xlim = xlim_map(d),
                legend = "right")
  p3 <- plot_RT(pars = fits[[3]], raw = samp_data, subj = subj, 
                subjs = subjs, xlim = xlim_map(d),
                legend = "none")
  p4 <- plot_RT(pars = fits[[4]], raw = samp_data, subj = subj, 
                subjs = subjs, xlim = xlim_map(d),
                legend = "none")
  pp_plot <- (p1 / p2 / p3 / p4)
  
  # Return list
  rm(fits); gc()
  list(test_retest = tr_plot, post_pred = pp_plot)
}

names(results) <- c("Study1-Stroop", "Study2-Stroop", 
                    "Study1-Flanker", "Study2-Flanker", 
                    "Study3-Posner", 
                    "Study1a-IAT", "Study2b-IAT")

stroop1 <- results$`Study1-Stroop`$test_retest | results$`Study1-Stroop`$post_pred
stroop2 <- results$`Study2-Stroop`$test_retest | results$`Study2-Stroop`$post_pred
flank1 <- results$`Study1-Flanker`$test_retest | results$`Study1-Flanker`$post_pred
flank2 <- results$`Study2-Flanker`$test_retest | results$`Study2-Flanker`$post_pred
posner3 <- results$`Study3-Posner`$test_retest | results$`Study3-Posner`$post_pred
iat1a <- results$`Study1a-IAT`$test_retest | results$`Study1a-IAT`$post_pred
iat2b <- results$`Study2b-IAT`$test_retest | results$`Study2b-IAT`$post_pred

ggsave(stroop1, filename = "Data/3_Plotted/stroop1.pdf", unit = "in",
       width = 12, height = 6)
ggsave(stroop2, filename = "Data/3_Plotted/stroop2.pdf", unit = "in",
       width = 12, height = 6)
ggsave(flank1, filename = "Data/3_Plotted/flank1.pdf", unit = "in",
       width = 12, height = 6)
ggsave(flank2, filename = "Data/3_Plotted/flank2.pdf", unit = "in",
       width = 12, height = 6)
ggsave(posner3, filename = "Data/3_Plotted/posner3.pdf", unit = "in",
       width = 12, height = 6)
ggsave(iat1a, filename = "Data/3_Plotted/iat1a.pdf", unit = "in",
       width = 12, height = 6)
ggsave(iat2b, filename = "Data/3_Plotted/iat2b.pdf", unit = "in",
       width = 12, height = 6)
