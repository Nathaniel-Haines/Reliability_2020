library(ggplot2)

setwd("~/Dropbox/Box/GitHub/Reliability_2020/")

# Function for normal mixture density
dmixnorm <- function(x, mean1, mean2) {
  dnorm(x, mean = mean1)*.5 + dnorm(x, mean = mean2)*.5
}

# Color scheme  
colors <- c("lightgray", "#edafaf", "#b5000c", "#850000", "black")

# The "Kvam Quintet of Equal Mean Distributions" 
p1 <- ggplot() +
  geom_function(aes(color = "N(3, 1)"), fun = dnorm, 
                args = list(mean = 3), size = 1) +
  geom_function(aes(color = "LogN(.59, 1)"), fun = dlnorm, 
                args = list(meanlog = .59), size = 1) +
  geom_function(aes(color = "N(1, 1) + N(5, 1)"), fun = dmixnorm, 
                args = list(mean1 = 1, mean2 = 5), size = 1) +
  geom_function(aes(color = "Exponential(3)"), fun = dexp, 
                args = list(rate = 1/3), size = 1) +
  geom_function(aes(color = "Uniform(0, 6)"), fun = dunif, 
                args = list(min = 0, max = 6), size = 1) +
  geom_vline(xintercept = 3, linetype = 2, color = "darkgray") +
  scale_color_manual(values = colors[1:5]) +
  scale_x_continuous(breaks = seq(0,6,1), limits = c(0,6)) +
  theme_bw(base_size = 15) +
  labs(color = "Distribution") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(p1, filename = "Data/3_Plotted/Figure_2.pdf", unit = "in",
       width = 7, height = 4)
