library(rstan)
library(dplyr)
library(tidyr)
library(foreach)
library(ggplot2)
library(patchwork)

rowSDs <- function(x, na.rm = T) {
  apply(x, 1, sd, na.rm = na.rm)
}

# Tasks and parameters within each model for plotting
tasks <- c("Study1-Stroop_normal", 
           "Study2b-IAT_normal", 
           "Study1-DDT_hyperbolic")
pars  <- list(c("mu_i_delta", "sigma_i_delta", "sigma_i_base"), 
              c("mu_i_delta", "sigma_i_delta", "sigma_i_base"),
              c("k", "c"))

# For response time mean/standard deviation contrasts
stan_data <- readRDS("Data/1_Preprocessed/stan_ready_all.rds")

# MLEs for hyperbolic model
mle_data <- readRDS("Data/2_Fitted/fit_Study1-DDT_hyperbolic_MLE.rds")

plot_data <- foreach(i=seq_along(tasks)) %do% {
  # Extract generative model parameters
  post_pars <- rstan::extract(readRDS(paste0("Data/2_Fitted/fit_", 
                                             tasks[i], ".rds")), 
                       pars = pars[[i]])
  
  # Loop through each parameter of interest
  plots <- foreach(p=pars[[i]][1:2]) %do% {
    # Compute parameter expectations (i.e. posterior means)
    if (p=="sigma_i_delta") {
      # This transformation ensure we look at the differences between
      # the sigma estimates on the scale of the SD for each condition, which
      # makes it comparable to the standard deviation contrast method
      sig_delta <- exp(post_pars$sigma_i_base + post_pars[[p]]) - exp(post_pars$sigma_i_base)
      par_pooled <- apply(sig_delta, c(2,3), mean) %>%
        as.data.frame() %>%
        mutate(subj_num = row_number(),
               pooled = "Yes")
    } else if (p %in% c("k", "c")) {
      par_pooled <- apply(log(post_pars[[p]]), c(2,3), mean) %>%
        as.data.frame() %>%
        mutate(subj_num = row_number(),
               pooled = "Yes")
    } else {
      par_pooled <- apply(post_pars[[p]], c(2,3), mean) %>%
        as.data.frame() %>%
        mutate(subj_num = row_number(),
               pooled = "Yes")
    }
    
    # Extract unpooled estimates
    if (tasks[i]=="Study1-DDT_hyperbolic") {
      # Unpooled MLEs
      par_unpooled <- mle_data %>%
        select(subj_num, Time, contains(p)) %>%
        mutate(Time = ifelse(Time == 1, "V1", "V2"),
               pooled = "No") %>% 
        spread(key = Time, value = p)
    } else {
      # Determine if computing unpooled mean or sd contrast
      diff_fun <- switch(p,
                         "mu_i_delta"    = rowMeans, 
                         "sigma_i_delta" = rowSDs)
      
      # Extract preprocessed data and compute contrasts
      samp_data <- stan_data[[strsplit(tasks[i], "_")[[1]][1]]]
      samp_data$RT[samp_data$RT==0] <- NA
      samp_t1 <- diff_fun(samp_data$RT[,2,1,], na.rm=T) - 
                 diff_fun(samp_data$RT[,1,1,], na.rm=T)
      samp_t2 <- diff_fun(samp_data$RT[,2,2,], na.rm=T) - 
                 diff_fun(samp_data$RT[,1,2,], na.rm=T)
      
      # Data.frame for unpooled estimates
      par_unpooled <- data.frame(subj_num = 1:(samp_data$N),
                                 pooled = "No",
                                 V1 = samp_t1,
                                 V2 = samp_t2)
    }
    # MLE parameter transformation function
    if (p %in% c("k", "c")) { 
      par_unpooled$V1 <- log(par_unpooled$V1)
      par_unpooled$V2 <- log(par_unpooled$V2)
    } 
    # Pooling figure
    par_pooling <- bind_rows(par_unpooled, par_pooled) %>%
      mutate(subj_num = as.factor(subj_num)) %>%
      ggplot(aes(x = V1, y = V2)) +
      geom_abline(intercept = 0, slope = 1, linetype = 2, color = "black", size = 1) +
      stat_ellipse(geom="polygon", type="norm", level=1/10, size=0, alpha=1/10, fill="gray") +
      stat_ellipse(geom="polygon", type="norm", level=2/10, size=0, alpha=1/10, fill="gray") +
      stat_ellipse(geom="polygon", type="norm", level=3/10, size=0, alpha=1/10, fill="gray") +
      stat_ellipse(geom="polygon", type="norm", level=4/10, size=0, alpha=1/10, fill="gray") +
      stat_ellipse(geom="polygon", type="norm", level=5/10, size=0, alpha=1/10, fill="gray") +
      stat_ellipse(geom="polygon", type="norm", level=6/10, size=0, alpha=1/10, fill="gray") +
      stat_ellipse(geom="polygon", type="norm", level=7/10, size=0, alpha=1/10, fill="gray") +
      stat_ellipse(geom="polygon", type="norm", level=8/10, size=0, alpha=1/10, fill="gray") +
      stat_ellipse(geom="polygon", type="norm", level=9/10, size=0, alpha=1/10, fill="gray") +
      stat_ellipse(geom="polygon", type="norm", level=.99, size=0, alpha=1/10, fill="gray") +
      geom_line(aes(group = subj_num), size = 1/4) +
      geom_point(aes(group = subj_num, color = pooled)) +
      scale_color_manual("Pooled?",
                         values = c("darkgray", "#990000")) +
      scale_fill_manual(values = c("darkgray", "#b5000c")) +
      scale_alpha_manual(values = c(1/10, 1/10)) +
      theme_minimal(base_size = 15) +
      theme(panel.grid = element_blank(),
            legend.position = "none")
    par_pooling
  }
  plots
}

p1 <- (plot_data[[1]][[1]] | plot_data[[1]][[2]]+xlim(-0.1, .19)+ylim(-.12, .18)) /
  (plot_data[[2]][[1]] | plot_data[[2]][[2]]) /
  (plot_data[[3]][[1]] | plot_data[[3]][[2]])

ggsave(p1, filename = "Data/3_Plotted/two-stage_vs_generative_pars.pdf", 
       unit = "in", width = 6, height = 8)
