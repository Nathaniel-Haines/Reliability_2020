library(patchwork)
library(EnvStats)

# Read in raw data
stan_data <- readRDS("Data/1_Preprocessed/stan_ready_all.rds")
dat <- stan_data$`Study1-Stroop`

# observed response times for participant 1 in congruent and inconguent conditions
con_obs_dat   <- dat$RT[8,1,1,]
incon_obs_dat <- dat$RT[8,2,1,]

# Two stage summary approach (the implicit generative model)
ts_mu_con    <- mean(con_obs_dat)
ts_mu_incon  <- mean(incon_obs_dat)
ts_sig_con   <- .001
ts_sig_incon <- .001

ts_con <- ggplot() +
  geom_histogram(aes(x = con_obs_dat, y = ..density..), 
                 bins = 30, fill = "gray") +
  geom_function(fun = dnorm, 
                args = list(mean = ts_mu_con, 
                            sd = ts_sig_con), 
                color = '#b5000c', n = 1000, size = 1, xlim = c(0,2)) +
  coord_cartesian(xlim = c(0, 2.2), ylim = c(0, 3.7)) +
  ggtitle("Congruent") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5))
ts_incon <- ggplot() +
  geom_histogram(aes(x = incon_obs_dat, y = ..density..), 
                 bins = 30, fill = "gray") +
  geom_function(fun = dnorm, 
                args = list(mean = ts_mu_incon, 
                            sd = ts_sig_incon), 
                color = '#b5000c', n = 1000, size = 1, xlim = c(0,2)) +
  coord_cartesian(xlim = c(0, 2.2), ylim = c(0, 3.7)) +
  ggtitle("Incongruent") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5))
tw_gen <- ts_con + ts_incon

# Normal generative model
norm_mu_con    <- mean(con_obs_dat)
norm_mu_incon  <- mean(incon_obs_dat)
norm_sig_con   <- sd(con_obs_dat)
norm_sig_incon <- sd(incon_obs_dat)

norm_con <- ggplot() +
  geom_histogram(aes(x = con_obs_dat, y = ..density..), 
                 bins = 30, fill = "gray") +
  geom_function(fun = dnorm, 
                args = list(mean = norm_mu_con, 
                            sd = norm_sig_con), 
                color = '#b5000c', n = 1000, size = 1, xlim = c(0,2)) +
  coord_cartesian(xlim = c(0, 2.2), ylim = c(0, 3.7)) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5))
norm_incon <- ggplot() +
  geom_histogram(aes(x = incon_obs_dat, y = ..density..), 
                 bins = 30, fill = "gray") +
  geom_function(fun = dnorm, 
                args = list(mean = norm_mu_incon, 
                            sd = norm_sig_incon), 
                color = '#b5000c', n = 1000, size = 1, xlim = c(0,2)) +
  coord_cartesian(xlim = c(0, 2.2), ylim = c(0, 3.7)) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5))
norm_gen <- norm_con + norm_incon

# Lognormal generative model
lnorm_mu_con    <- mean(log(con_obs_dat))
lnorm_mu_incon  <- mean(log(incon_obs_dat))
lnorm_sig_con   <- sd(log(con_obs_dat))
lnorm_sig_incon <- sd(log(incon_obs_dat))

lnorm_con <- ggplot() +
  geom_histogram(aes(x = con_obs_dat, y = ..density..), 
                 bins = 30, fill = "gray") +
  geom_function(fun = dlnorm, 
                args = list(meanlog = lnorm_mu_con, 
                            sdlog = lnorm_sig_con), 
                color = '#b5000c', n = 1000, size = 1, xlim = c(0,2)) +
  coord_cartesian(xlim = c(0, 2.2), ylim = c(0, 3.7)) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5))
lnorm_incon <- ggplot() +
  geom_histogram(aes(x = incon_obs_dat, y = ..density..), 
                 bins = 30, fill = "gray") +
  geom_function(fun = dlnorm, 
                args = list(meanlog = lnorm_mu_incon, 
                            sdlog = lnorm_sig_incon), 
                color = '#b5000c', n = 1000, size = 1, xlim = c(0,2)) +
  coord_cartesian(xlim = c(0, 2.2), ylim = c(0, 3.7)) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5))
lnorm_gen <- lnorm_con + lnorm_incon

# Shifted lognormal generative model
# First, estimate parameters using maximum likelihood
mle_shlnorm_con   <- elnorm3(con_obs_dat)
mle_shlnorm_incon <- elnorm3(incon_obs_dat) 

thresh <- min(mle_shlnorm_con$parameters[3], mle_shlnorm_incon$parameters[3])

slnorm_mu_con      <- mean(log(con_obs_dat-thresh))
slnorm_mu_incon    <- mean(log(incon_obs_dat-thresh))
slnorm_sig_con     <- sd(log(con_obs_dat-thresh))
slnorm_sig_incon   <- sd(log(incon_obs_dat-thresh))

slnorm_con <- ggplot() +
  geom_histogram(aes(x = con_obs_dat, y = ..density..), 
                 bins = 30, fill = "gray") +
  geom_function(fun = dlnorm3, 
                args = list(meanlog = slnorm_mu_con, 
                            sdlog = slnorm_sig_con,
                            threshold = thresh), 
                color = '#b5000c', n = 1000, size = 1, xlim = c(0,2)) +
  coord_cartesian(xlim = c(0, 2.2), ylim = c(0, 3.7)) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5))
slnorm_incon <- ggplot() +
  geom_histogram(aes(x = incon_obs_dat, y = ..density..), 
                 bins = 30, fill = "gray") +
  geom_function(fun = dlnorm3, 
                args = list(meanlog = slnorm_mu_incon, 
                            sdlog = slnorm_sig_incon,
                            threshold = thresh), 
                color = '#b5000c', n = 1000, size = 1, xlim = c(0,2)) +
  coord_cartesian(xlim = c(0, 2.2), ylim = c(0, 3.7)) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5))
slnorm_gen <- slnorm_con + slnorm_incon

# all the goodness
p1 <- tw_gen / norm_gen / lnorm_gen / slnorm_gen

ggsave(p1, filename = "Data/3_Plotted/fig_behavioral_models.pdf", unit = "in",
       width = 6, height = 7)
