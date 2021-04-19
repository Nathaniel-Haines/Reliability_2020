library(rstan)
library(cowplot)
library(patchwork)
library(bayesplot)
library(gridExtra)
library(grid)
library(dplyr)
library(scales)

# Read in raw data
stan_data <- readRDS("Data/1_Preprocessed/stan_ready_all.rds")
dat <- stan_data$`Study1-Stroop`

# observed response times for participant 1 in congruent and inconguent conditions
con_obs_dat   <- dat$RT[1,1,1,]
incon_obs_dat <- dat$RT[1,2,1,]

# Posterior lognormal generative for congruent (person 1)
pt_p1_c1_t1_mus  <- pars$beta_con_mu[iters,subj[1],1]
pt_p1_c1_t1_sigs <- pars$beta_con_sig[iters,subj[1],1]
pt_p1_c1_t1 <- ggplot() +
  geom_histogram(aes(x = p1_c1_t1_dat, y = ..density..), 
                 bins = 30, fill = "gray") +
  coord_cartesian(xlim = c(0, 2.2), ylim = c(0, 3.2)) +
  ggtitle("Congruent") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5))
# Posterior lognormal generative for congruent (person 2)
pt_p2_c1_t1_mus  <- pars$beta_con_mu[iters,subj[2],1]
pt_p2_c1_t1_sigs <- pars$beta_con_sig[iters,subj[2],1]
pt_p2_c1_t1 <- ggplot() +
  geom_histogram(aes(x = p2_c1_t1_dat, y = ..density..), 
                 bins = 30, fill = "gray") +
  coord_cartesian(xlim = c(0, 2.2), ylim = c(0, 3.2)) +
  ggtitle("Congruent") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5))
# Posterior lognormal generative for incongruent (person 1) 
pt_p1_c2_t1_mus  <- pars$beta_con_mu[iters,subj[1],1] + pars$beta_delta_mu[iters,subj[1],1]
pt_p1_c2_t1_sigs <- pars$beta_delta_sig[iters,subj[1],1]
pt_p1_c2_t1 <- ggplot() +
  geom_histogram(aes(x = p1_c2_t1_dat, y = ..density..), 
                 bins = 30, fill = "gray") +
  coord_cartesian(xlim = c(0, 2.2), ylim = c(0, 3.2)) +
  ggtitle("Incongruent") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5))
# Posterior lognormal generative for incongruent (person 2) 
pt_p2_c2_t1_mus  <- pars$beta_con_mu[iters,subj[2],1] + pars$beta_delta_mu[iters,subj[2],1]
pt_p2_c2_t1_sigs <- pars$beta_delta_sig[iters,subj[2],1]
pt_p2_c2_t1 <- ggplot() +
  geom_histogram(aes(x = p2_c2_t1_dat, y = ..density..), 
                 bins = 30, fill = "gray") +
  coord_cartesian(xlim = c(0, 2.2), ylim = c(0, 3.2)) +
  ggtitle("Incongruent") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5))
# Posterior lognormal generative add lines
for (i in seq_along(iters)) {
  input1 <- paste0("geom_function(fun = dlnorm, 
                                 args = list(meanlog = pt_p1_c1_t1_mus[",i,"], 
                                             sdlog = pt_p1_c1_t1_sigs[",i,"]), 
                                 color = '#b5000c', alpha = .1, xlim = c(0,2.2))")
  input2 <- paste0("geom_function(fun = dlnorm, 
                                 args = list(meanlog = pt_p1_c2_t1_mus[",i,"], 
                                             sdlog = pt_p1_c2_t1_sigs[",i,"]), 
                                 color = '#b5000c', alpha = .1, xlim = c(0,2.2))")
  input3 <- paste0("geom_function(fun = dlnorm, 
                                 args = list(meanlog = pt_p2_c1_t1_mus[",i,"], 
                                             sdlog = pt_p2_c1_t1_sigs[",i,"]), 
                                 color = '#b5000c', alpha = .1, xlim = c(0,2.2))")
  input4 <- paste0("geom_function(fun = dlnorm, 
                                 args = list(meanlog = pt_p2_c2_t1_mus[",i,"], 
                                             sdlog = pt_p2_c2_t1_sigs[",i,"]), 
                                 color = '#b5000c', alpha = .1, xlim = c(0,2.2))")
  pt_p1_c1_t1 <- pt_p1_c1_t1 + eval(parse(text=input1))  
  pt_p1_c2_t1 <- pt_p1_c2_t1 + eval(parse(text=input2))  
  pt_p2_c1_t1 <- pt_p2_c1_t1 + eval(parse(text=input3))  
  pt_p2_c2_t1 <- pt_p2_c2_t1 + eval(parse(text=input4))  
}
# Save out plots
pr_p1_lnorm <- arrangeGrob(patchworkGrob(pr_p1_c1_t1 + pr_p1_c2_t1), 
                           bottom = textGrob(expression("Person 1 RTs"), 
                                             gp = gpar(fontsize = 18)), 
                           top = textGrob("Session 1", gp = gpar(fontsize = 20)))
pr_p2_lnorm <- arrangeGrob(patchworkGrob(pr_p2_c1_t1 + pr_p2_c2_t1), 
                           bottom = textGrob(expression("Person 47 RTs"), 
                                             gp = gpar(fontsize = 18)), 
                           top = textGrob("Session 1", gp = gpar(fontsize = 20)))
pt_p1_lnorm <- arrangeGrob(patchworkGrob(pt_p1_c1_t1 + pt_p1_c2_t1), 
                           bottom = textGrob(expression("Person 1 RTs"), 
                                             gp = gpar(fontsize = 18)), 
                           top = textGrob("Session 1", gp = gpar(fontsize = 20)))
pt_p2_lnorm <- arrangeGrob(patchworkGrob(pt_p2_c1_t1 + pt_p2_c2_t1), 
                           bottom = textGrob(expression("Person 47 RTs"), 
                                             gp = gpar(fontsize = 18)), 
                           top = textGrob("Session 1", gp = gpar(fontsize = 20)))

ggsave(pr_p1_lnorm, filename = "Data/3_Plotted/prior_p1_lnorm.pdf", unit = "in",
       width = 4, height = 3)
ggsave(pr_p2_lnorm, filename = "Data/3_Plotted/prior_p2_lnorm.pdf", unit = "in",
       width = 4, height = 3)
ggsave(pt_p1_lnorm, filename = "Data/3_Plotted/post_p1_lnorm.pdf", unit = "in",
       width = 4, height = 3)
ggsave(pt_p2_lnorm, filename = "Data/3_Plotted/post_p2_lnorm.pdf", unit = "in",
       width = 4, height = 3)

## Correlation plots
# Prior on correlation for mu
set.seed(7)
pr_corr_mu_dat <- data.frame(T1 = rnorm(50,0,1) + exp(rnorm(50,0,1)) * rnorm(50,0,1),
                             T2 = rnorm(50,0,1) + exp(rnorm(50,0,1)) * rnorm(50,0,1))
pr_corr_mu <- pr_corr_mu_dat %>%
  ggplot(aes(x = T1, y = T2)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "black", size = 1) +
  stat_ellipse(geom="polygon", type="norm", level=1/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=2/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=3/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=4/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=5/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=6/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=7/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=8/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=9/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=.99, size=0, alpha=1/10, fill="#b5000c") +
  geom_point(color = "#990000") +
  coord_cartesian(xlim = c(-5,5), ylim = c(-5,5)) +
  ggtitle(expression("Test-retest "~mu[Delta])) +
  xlab(expression("Session 1 "~mu[Delta])) +
  ylab(expression("Session 2 "~mu[Delta])) +
  theme_minimal(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        plot.title = element_text(hjust = .5))

# Prior on correlation for sigma
pr_corr_sig_dat <- data.frame(T1 = rnorm(50,0,1) + exp(rnorm(50,0,1)) * rnorm(50,0,1),
                              T2 = rnorm(50,0,1) + exp(rnorm(50,0,1)) * rnorm(50,0,1))
pr_corr_sig <- pr_corr_sig_dat %>%
  ggplot(aes(x = T1, y = T2)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "black", size = 1) +
  stat_ellipse(geom="polygon", type="norm", level=1/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=2/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=3/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=4/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=5/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=6/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=7/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=8/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=9/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=.99, size=0, alpha=1/10, fill="#b5000c") +
  geom_point(color = "#990000") +
  coord_cartesian(xlim = c(-5,5), ylim = c(-5,5)) +
  ggtitle(expression("Test-retest "~sigma[Delta])) +
  xlab(expression("Session 1 "~sigma[Delta])) +
  ylab(expression("Session 2 "~sigma[Delta])) +
  theme_minimal(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        plot.title = element_text(hjust = .5))

# Posterior on correlation for mu 
pt_corr_mu_dat <- data.frame(T1 = colMeans(pars$beta_delta_mu[,,1]),
                             T2 = colMeans(pars$beta_delta_mu[,,2]))
pt_corr_mu <- pt_corr_mu_dat %>%
  ggplot(aes(x = T1, y = T2)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "black", size = 1) +
  stat_ellipse(geom="polygon", type="norm", level=1/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=2/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=3/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=4/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=5/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=6/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=7/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=8/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=9/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=.99, size=0, alpha=1/10, fill="#b5000c") +
  geom_point(color = "#990000") +
  coord_cartesian(xlim = c(0,.2), ylim = c(0,.2)) +
  ggtitle(expression("Test-retest "~mu[Delta])) +
  xlab(expression("Session 1 "~mu[Delta])) +
  ylab(expression("Session 2 "~mu[Delta])) +
  theme_minimal(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        plot.title = element_text(hjust = .5))

# Posterior on correlation for sigma
pt_corr_sig_dat <- data.frame(T1 = colMeans(pars$beta_delta_sig[,,1]),
                              T2 = colMeans(pars$beta_delta_sig[,,2]))
pt_corr_sig <- pt_corr_sig_dat %>%
  ggplot(aes(x = T1, y = T2)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "black", size = 1) +
  stat_ellipse(geom="polygon", type="norm", level=1/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=2/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=3/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=4/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=5/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=6/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=7/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=8/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=9/10, size=0, alpha=1/10, fill="#b5000c") +
  stat_ellipse(geom="polygon", type="norm", level=.99, size=0, alpha=1/10, fill="#b5000c") +
  geom_point(color = "#990000") +
  coord_cartesian(xlim = c(0,.5), ylim = c(0,.5)) +
  ggtitle(expression("Test-retest "~sigma[Delta])) +
  xlab(expression("Session 1 "~sigma[Delta])) +
  ylab(expression("Session 2 "~sigma[Delta])) +
  theme_minimal(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        plot.title = element_text(hjust = .5))

# save multivariate normals
ggsave(pr_corr_mu, filename = "Data/3_Plotted/prior_corr_mu.pdf", unit = "in",
       width = 4, height = 4)
ggsave(pr_corr_sig, filename = "Data/3_Plotted/prior_corr_sig.pdf", unit = "in",
       width = 4, height = 4)
ggsave(pt_corr_mu, filename = "Data/3_Plotted/post_corr_mu.pdf", unit = "in",
       width = 4, height = 4)
ggsave(pt_corr_sig, filename = "Data/3_Plotted/post_corr_sig.pdf", unit = "in",
       width = 4, height = 4)
