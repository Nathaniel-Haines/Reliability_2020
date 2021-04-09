library(dplyr)
library(tidyr)
library(foreach)
library(ggplot2)
library(cowplot)
library(patchwork)

setwd("~/Dropbox/Box/GitHub/Reliability_2020/")

# Range of values for plotting changes in mean
time <- seq(0, 2.5, length.out = 1000)
mu_vals <- seq(-2, .5, length.out = 10)
plot_mu <- foreach(i=seq_along(mu_vals), .combine = "rbind") %do% {
  data.frame(mu          = rep(mu_vals[i], length(time)),
             Time        = time,
             Density     = dlnorm(time, mu_vals[i], .2))  
}
# Range of values for plotting changes in sd
sig_vals <- seq(.1, 1, length.out = 10)
plot_sig <- foreach(i=seq_along(sig_vals), .combine = "rbind") %do% {
  data.frame(sig         = rep(sig_vals[i], length(time)),
             Time        = time,
             Density     = dlnorm(time, -.5, sig_vals[i]))  
}
# Range of values for plotting changes in shift
time <- seq(0, 2.5, length.out = 1000)
delta_vals <- seq(.05, 1.5, length.out = 10)
plot_delta <- foreach(i=seq_along(delta_vals), .combine = "rbind") %do% {
  data.frame(delta       = rep(delta_vals[i], length(time)),
             Time        = time,
             Density     = brms::dshifted_lnorm(time, -1, .2, delta_vals[i]) ) 
}

# Plotting each parameter
p1 <- plot_mu %>%
  ggplot(aes(x = Time, y = Density, color = mu)) +
  geom_path() +
  scale_color_continuous(expression(mu), low = "#b5000c", high = "gray") +
  ggtitle(expression("Change in Difficulty (" ~ mu ~ ")")) +
  xlab("Response Time") +
  ylab("Density") +
  annotate(geom = "text", x = 1.5, y = 10, size = 5,
                label = expression(sigma ~ " = .2")) +
  theme_cowplot(font_size = 15) +
  theme(plot.title = element_text(hjust = .5))
p2 <- plot_sig %>%
  ggplot(aes(x = Time, y = Density, color = sig)) +
  geom_path() +
  scale_color_continuous(expression(sigma), low = "#b5000c", high = "gray") +
  ggtitle(expression("Change in Dispersion (" ~ sigma ~ ")")) +
  xlab("Response Time") +
  ylab("Density") +
  annotate(geom = "text", x = 1.5, y = 4, size = 5,
           label = expression(mu ~ " = -.5")) +
  theme_cowplot(font_size = 15) +
  theme(plot.title = element_text(hjust = .5))
p3 <- plot_delta %>%
  ggplot(aes(x = Time, y = Density, color = delta)) +
  geom_path() +
  scale_color_continuous(expression(delta), low = "#b5000c", high = "gray") +
  ggtitle(expression("Change in Shift (" ~ delta ~ ")")) +
  xlab("Response Time") +
  ylab("Density") +
  annotate(geom = "text", x = 2, y = 6.5, size = 5,
           label = expression(mu ~ " = -1, " ~ sigma ~ " = .2")) +
  theme_cowplot(font_size = 15) +
  theme(plot.title = element_text(hjust = .5))

p4 <- ((p1 | p2) / p3) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 18))

ggsave(p4, filename = "Data/3_Plotted/fig_4.pdf",
       dpi = 300, height = 6, width = 8, unit = "in")
