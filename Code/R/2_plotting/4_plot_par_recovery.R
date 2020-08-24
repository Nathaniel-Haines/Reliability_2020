library(dplyr)
library(ggplot2)

setwd("~/Dropbox/Box/GitHub/Reliability_2020/")

results <- readRDS("Data/2_Fitted/parameter_recovery.rds")

p1 <- results %>% 
  mutate(par = ifelse(par=="R_mu", "mu", "sigma")) %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype = 2, color = "gray") +
  geom_line(aes(x = R_true, y = R_mu, color = type)) +
  geom_ribbon(aes(x = R_true, ymin = R_lo, ymax = R_hi, fill = type), alpha = .3) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual("Model", values = c("#8F2727", "darkgray")) +
  scale_fill_manual("Model", values = c("#8F2727", "black")) +
  facet_grid(c("samp_size", "par")) +
  xlab("True Correlation") +
  ylab("Estimated Correlation") +
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank())

ggsave(p1, filename = "Data/3_Plotted/parameter_recovery.png",
       unit = "in", width = 7.5, height = 8, dpi = 300)

