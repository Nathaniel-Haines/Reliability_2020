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
  n_iters <- length(fits[[1]]$R_sigma[,1,2])
  est1 <- switch(parameter,
                 "mu" = foreach(i=1:n_iters, .combine = "c") %do% {
                   cor(fits[[1]]$mu_i[i,,2] - fits[[1]]$mu_i[i,,1], fits[[1]]$mu_i[i,,4] - fits[[1]]$mu_i[i,,3])
                   }, 
                 "sig" = foreach(i=1:n_iters, .combine = "c") %do% {
                   cor(fits[[1]]$sigma_i[i,,2] - fits[[1]]$sigma_i[i,,1], fits[[1]]$sigma_i[i,,4] - fits[[1]]$sigma_i[i,,3])
                 })
  est2 <- switch(parameter,
                 "mu" = foreach(i=1:n_iters, .combine = "c") %do% {
                   cor(fits[[2]]$mu_i[i,,2] - fits[[2]]$mu_i[i,,1], fits[[2]]$mu_i[i,,4] - fits[[2]]$mu_i[i,,3])
                 }, 
                 "sig" = foreach(i=1:n_iters, .combine = "c") %do% {
                   cor(fits[[2]]$sigma_i[i,,2] - fits[[2]]$sigma_i[i,,1], fits[[2]]$sigma_i[i,,4] - fits[[2]]$sigma_i[i,,3])
                 })
  est3 <- switch(parameter,
                 "mu" = foreach(i=1:n_iters, .combine = "c") %do% {
                   cor(fits[[3]]$mu_i[i,,2] - fits[[3]]$mu_i[i,,1], fits[[3]]$mu_i[i,,4] - fits[[3]]$mu_i[i,,3])
                 }, 
                 "sig" = foreach(i=1:n_iters, .combine = "c") %do% {
                   cor(fits[[3]]$sigma_i[i,,2] - fits[[3]]$sigma_i[i,,1], fits[[3]]$sigma_i[i,,4] - fits[[3]]$sigma_i[i,,3])
                 })
  qplot() +
    geom_density(aes(x = est3, fill = "Shifted Lognormal ")) +
    geom_density(aes(x = est2, fill = "Lognormal "), alpha = .9) +
    geom_density(aes(x = est1, fill = " Normal "), alpha = .7) +
    scale_fill_manual(values = c("#edafaf", "#b5000c", "#700000"),
                      labels = c(" Normal ", "Lognormal ", "Shifted Lognormal ")) +
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

plot_RT <- function(pars, raw, subj, subjs, n_draws, xlim, ylim, legend) {
  color_scheme_set("red")
  samp <- sample(1:900, n_draws, F)
  
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
