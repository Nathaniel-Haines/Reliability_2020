library(ggplot2)
library(dplyr)

rowSDs <- function(x, na.rm=T) {
  apply(x, 1, sd, na.rm = na.rm)
}

plot_retest <- function(pars, parameter, task, samp_data, legend) {
  # First, calculate sample test-retest Pearson's correlation (i.e. two-stage approach)
  par_trans <- switch(parameter, "mu" = rowMeans, "sig" = rowSDs)
  samp_data$RT[samp_data$RT==0] <- NA
  samp_t1 <- par_trans(samp_data$RT[,1,1,],na.rm=T) - par_trans(samp_data$RT[,2,1,],na.rm=T)
  samp_t2 <- par_trans(samp_data$RT[,1,2,],na.rm=T) - par_trans(samp_data$RT[,2,2,],na.rm=T)
  samp_cor <- cor.test(samp_t1, samp_t2)
  
  # Plot posterior correlation and sample test-retest
  est1 <- switch(parameter,
                 "mu" = fits[[1]]$R_mu[,1,2],
                 "sig" = fits[[1]]$R_sigma[,1,2])
  est2 <- switch(parameter,
                 "mu" = fits[[2]]$R_mu[,1,2],
                 "sig" = fits[[2]]$R_sigma[,1,2])
  est3 <- switch(parameter,
                 "mu" = fits[[3]]$R_mu[,1,2],
                 "sig" = fits[[3]]$R_sigma[,1,2])
  est4 <- switch(parameter,
                 "mu" = fits[[4]]$R_mu[,1,2],
                 "sig" = fits[[4]]$R_sigma[,1,2])
  ggplot() +
    geom_density(aes(x = est4, fill = "Shifted Lognormal Mix")) +
    geom_density(aes(x = est3, fill = "Shifted Lognormal "), alpha = .9) +
    geom_density(aes(x = est2, fill = "Lognormal "), alpha = .7) +
    geom_density(aes(x = est1, fill = " Normal "), alpha = .5) +
    scale_fill_manual(values = c("#edafaf", "#A25050", "#b5000c", "#700000"),
                      labels = c(" Normal ", "Lognormal ", "Shifted Lognormal ", "Shifted Lognormal Mix")) +
    geom_vline(xintercept = samp_cor$estimate, linetype = 2, size = 1) +
    geom_segment(aes(x = samp_cor$conf.int[1], xend = samp_cor$conf.int[2], y = 0, yend = 0), 
                 color = I("black"), size = 1.5) +
    xlab(bquote(.(task) ~ " Test-Retest")) +
    xlim(-1,1) +
    theme_minimal(base_size = 15) +
    theme(panel.grid = element_blank(),
          legend.position = legend, 
          legend.title = element_blank())
}

plot_RT <- function(pars, raw, subj, subjs, xlim, ylim, legend) {
  color_scheme_set("red")
  
  n_samples <- dim(pars$post_pred_c1_t1)[1]
  samp <- sample(1:n_samples, n_draws, F)
  
  pars1 <- pars$post_pred_c1_t1[samp,which(subj==subjs),]
  pars2 <- pars$post_pred_c2_t1[samp,which(subj==subjs),]
  
  raw1 <- raw$RT[subj,1,1,]
  raw2 <- raw$RT[subj,2,1,]
  
  rm1 <- which(raw1>0)
  rm2 <- which(raw2>0)
  
  raw1 <- raw1[rm1]
  raw2 <- raw2[rm2]
  pars1 <- pars1[,rm1]
  pars2 <- pars2[,rm2]
  
  p1 <- ppc_dens_overlay(y = raw1,
                         yrep = pars1) +
    geom_vline(xintercept = mean(raw1), linetype = 2, color = I("black")) +
    coord_cartesian(xlim = xlim) +
    xlab("Condition 1") +
    theme_minimal(base_size = 15) +
    theme(panel.grid = element_blank(),
          legend.position = "none")
  p2 <- ppc_dens_overlay(y = raw2,
                         yrep = pars2) +
    geom_vline(xintercept = mean(raw2), linetype = 2, color = I("black")) +
    coord_cartesian(xlim = xlim) +
    xlab("Condition 2") +
    theme_minimal(base_size = 15) +
    theme(panel.grid = element_blank(),
          legend.position = "none")
  p1 + p2
}
