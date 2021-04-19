library(rstan)
library(cowplot)
library(patchwork)
library(tidybayes)
library(hBayesDM)
library(dplyr)
library(foreach)

# function for later
rowSDs <- function(x, na.rm=T) {
  apply(x, 1, sd, na.rm = na.rm)
}

stan_data <- readRDS("Data/1_Preprocessed/stan_ready_all.rds")

fit_names <- c("normal", "lognormal", "shift_lognormal")
data_names <- c("Study1-Stroop", "Study2-Stroop", 
                "Study1-Flanker", "Study2-Flanker", 
                "Study3-Posner", 
                "Study1a-IAT", "Study2b-IAT")

# Results 
results <- foreach(d=data_names, .combine = "rbind") %do% {
  samp_data <- stan_data[[d]]
  tmp <- foreach(i=fit_names, .combine = "rbind") %do% {
    # Extract generative model estimates
    pars <- rstan::extract(readRDS(paste0("Data/2_Fitted/fit_", d, "_jointSep_", i, ".rds")),
                           pars = c("R_mu_delta", "R_sigma_delta"))
    
    # Compute posterior means and 95% HDIs for both mu and sigma
    data.frame(Model     = i,
               Task      = d,
               Parameter = c("mu", "sigma"),
               mean      = c(mean(pars$R_mu_delta[,1,2]), 
                             mean(pars$R_sigma_delta[,1,2])),
               lb        = c(HDIofMCMC(pars$R_mu_delta[,1,2])[1], 
                             HDIofMCMC(pars$R_sigma_delta[,1,2])[1]),
               ub        = c(HDIofMCMC(pars$R_mu_delta[,1,2])[2], 
                             HDIofMCMC(pars$R_sigma_delta[,1,2])[2]))
  }
  # Set 0 values to NA (they were originally set to 0 for Stan)
  samp_data$RT[samp_data$RT==0] <- NA
  
  # Mean contrasts
  means_t1  <- rowMeans(samp_data$RT[,1,1,],na.rm=T) - rowMeans(samp_data$RT[,2,1,],na.rm=T)
  means_t2  <- rowMeans(samp_data$RT[,1,2,],na.rm=T) - rowMeans(samp_data$RT[,2,2,],na.rm=T)
  means_cor <- cor.test(means_t1, means_t2)
  
  # SD contrasts
  sd_t1  <- rowSDs(samp_data$RT[,1,1,],na.rm=T) - rowSDs(samp_data$RT[,2,1,],na.rm=T)
  sd_t2  <- rowSDs(samp_data$RT[,1,2,],na.rm=T) - rowSDs(samp_data$RT[,2,2,],na.rm=T)
  sd_cor <- cor.test(sd_t1, sd_t2)
  
  # Combine two-stage with generative model estimates and return
  data.frame(Model     = "Two-Stage",
             Task      = d, 
             Parameter = c("mu", "sigma"),
             mean      = c(means_cor$estimate, 
                           sd_cor$estimate),
             lb        = c(means_cor$conf.int[1],
                           sd_cor$conf.int[1]),
             ub        = c(means_cor$conf.int[2],
                           sd_cor$conf.int[2])) %>%
    bind_rows(tmp)
}

# RT data plot
p1 <- results %>%
  mutate(Model = factor(Model, 
                        levels = c("Two-Stage","normal","lognormal","shift_lognormal"),
                        labels = c("Two-Stage","Normal","Lognormal","Shifted Lognormal")),
         Task = factor(Task,
                       levels = rev(c("Study1-Stroop", "Study2-Stroop", "Study1-Flanker", 
                                      "Study2-Flanker", "Study3-Posner", "Study1a-IAT", 
                                      "Study2b-IAT")),
                       labels = rev(c("Stroop Study 1", "Stroop Study 2", "Flanker Study 1", 
                                      "Flanker Study 2", "Posner Study 3", "IAT Study 1a", 
                                      "IAT Study 2b"))),
         Parameter = ifelse(Parameter=="mu", "mu[Delta]", "sigma[Delta]")) %>%
  ggplot(aes(y = Task, x = mean, xmin = lb, xmax = ub, color = Model)) +
  geom_pointinterval(position = position_dodge(width = .3), size = 2) +
  geom_vline(xintercept = -1, color = "darkgray") +
  geom_vline(xintercept = 1, color = "darkgray") +
  ggtitle("Joint Separate Model") +
  xlim(-1,1) + 
  xlab("Test-Retest Correlation") +
  ylab("Dataset") +
  scale_color_manual(labels = c("Two-Stage","Normal","Lognormal","Shifted Lognormal"),
                     values = c("black","#edafaf", "#b5000c", "#700000")) +
  facet_grid(cols = vars(Parameter),
             labeller = label_parsed) +
  theme_minimal(base_size=15) +
  theme(plot.title = element_text(hjust=.5), 
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank())

ggsave(p1, filename = "Data/3_Plotted/test_retest_all_jointSep.pdf", unit = "in",
       width = 8, height = 5)
