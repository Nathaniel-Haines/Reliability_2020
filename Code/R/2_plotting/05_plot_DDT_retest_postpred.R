library(rstan)
library(dplyr)
library(patchwork)
library(bayesplot)

# Long format input data
dat <- readRDS("Data/1_Preprocessed/long_format_all.rds")

# Generative model results
fit <- readRDS("Data/2_Fitted/fit_Study1-DDT_hyperbolic.rds")
pars <- rstan::extract(fit)

# Two-stage approach results
mle_results <- readRDS("Data/2_Fitted/fit_Study1-DDT_hyperbolic_MLE.rds")

# test-retest correlations for k and c
k_cor <- cor.test(log(mle_results$k[mle_results$Time==1]), log(mle_results$k[mle_results$Time==2]))
c_cor <- cor.test(log(mle_results$c[mle_results$Time==1]), log(mle_results$c[mle_results$Time==2]))

# Plot
p1 <- qplot() +
  geom_density(aes(x = pars$R_k[,1,2], fill = "Hyperbolic"), alpha = .9) +
  scale_fill_manual(values = c("#b5000c"),
                    labels = c("Hyperbolic")) +
  geom_vline(xintercept = k_cor$estimate, linetype = 2, size = 1) +
  geom_segment(aes(x = k_cor$conf.int[1], xend = k_cor$conf.int[2], 
                   y = 0, yend = 0), color = I("black"), size = 1.5) +
  xlab(bquote("log(" ~ italic("k") ~ ") Test-Retest")) +
  xlim(-1,1) +
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.position = "none", 
        legend.title = element_blank())

p2 <- qplot() +
  geom_density(aes(x = pars$R_c[,1,2], fill = "Hyperbolic"), alpha = .9) +
  scale_fill_manual(values = c("#b5000c"),
                    labels = c("Hyperbolic")) +
  geom_vline(xintercept = c_cor$estimate, linetype = 2, size = 1) +
  geom_segment(aes(x = c_cor$conf.int[1], xend = c_cor$conf.int[2], 
                   y = 0, yend = 0), color = I("black"), size = 1.5) +
  xlab(bquote("log(" ~ italic("c") ~ ") Test-Retest")) +
  xlim(-1,1) +
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.position = "none", 
        legend.title = element_blank())


# Plotting parameters (choose subjects and # of posterior draws)
subjs <- c(17, 22, 44)
n_draws <- 100

# Create all plots
all_plots <- foreach(i=subjs) %do% {
  # Empirical indifference points
  indiff_dat <- dat$`Study1-DDT` %>%
    filter(subj_num==i & Time==1) %>%
    group_by(delay_later) %>%
    summarize(indiff = ifelse(all(choice==1), 800, 
                              min(amount_sooner[choice==0]))) %>%
    mutate(delay_norm  = delay_later / max(delay_later),
           indiff_norm = indiff / max(indiff))
  
  # Bayesian posterior draws of discounting curve
  k_draws <- foreach(d=1:n_draws, .combine = "rbind") %do% {
    subj_k <- pars$k[sample(1:3000, size = 1),i,1]
    subj_dis <- 1 / (1 + subj_k*(0:520))
    data.frame(k = rep(subj_k, 521),
               delay   = seq(0, 1, length.out = 521),
               subj_dis = subj_dis,
               Estimator = "HBA")
  }
  # Maximum likelihood discounting curve
  subj_k <- with(mle_results, k[subj_num==i & Time==1])
  subj_dis <- 1 / (1 + subj_k*(0:520))
  k_draws <- rbind(k_draws,
                   data.frame(k = rep(subj_k, 521),
                              delay   = seq(0, 1, length.out = 521),
                              subj_dis = subj_dis,
                              Estimator = "MLE"))
  # Plot them all
  ggplot() +
    geom_line(data = k_draws, 
              aes(x = delay, y = subj_dis, group = k, 
                  color = Estimator, linetype = Estimator,  
                  alpha = Estimator, size = Estimator)) +
    scale_color_manual(values = c("#C79999", "black")) +
    scale_linetype_manual(values = c(1, 2)) +
    scale_alpha_manual(values = c(.2, 1)) +
    scale_size_manual(values = c(.5, 1)) +
    geom_point(data = indiff_dat, 
               aes(x = delay_norm, y = indiff_norm),
               color = I("#8F2727"), size = 2) +
    xlab("Normalized Delay") +
    ylab("Normalized Amount") +
    ylim(0,1) +
    theme_cowplot(font_size = 15) +
    theme(strip.background = element_blank(),
          plot.title = element_text(hjust = .5),
          legend.position = "none")
}

# combine all plot panels and save
p3 <- (p1 / p2) | (all_plots[[1]] / all_plots[[2]] / all_plots[[3]])
p3
ggsave(p3, filename = "Data/3_Plotted/DDT1.pdf", unit = "in",
       width = 10, height = 6)
