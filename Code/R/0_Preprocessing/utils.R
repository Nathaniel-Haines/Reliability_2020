pre_stp_flk_psr <- function(x) {
  # Current study-task pair
  tmp_data <- na.omit(x) %>%
    filter(RT >= 0)
  # Number of subjects
  n_subj <- length(unique(tmp_data$subj_num))
  # Number of conditions 
  n_cond <- 2
  # Number of timepoints
  n_time <- 2
  # Determine number of trials within subjects, conditions, and timepoints
  T_subj_all <- tmp_data %>%
    filter(Condition != 1) %>%
    group_by(subj_num, Condition, Time) %>%
    summarize(n_trials = n())
  T_max <- T_subj_all %>%
  {max(.$n_trials)}
  
  # Create RT data array for stan; dims = (subject, condition, time, trial)
  RT <- correct <- array(0, dim = c(n_subj, n_cond, n_time, T_max))
  T_subj <- RT_min <- RT_max <- array(NA, dim = c(n_subj, n_cond, n_time))
  for (i in 1:n_subj) {
    # Number of trials per condition/timpoint
    c0_t1 <- T_subj[i, 1, 1] <- with(T_subj_all, n_trials[subj_num==i & Condition==0 & Time==1])
    c2_t1 <- T_subj[i, 2, 1] <- with(T_subj_all, n_trials[subj_num==i & Condition==2 & Time==1])
    c0_t2 <- T_subj[i, 1, 2] <- with(T_subj_all, n_trials[subj_num==i & Condition==0 & Time==2])
    c2_t2 <- T_subj[i, 2, 2] <- with(T_subj_all, n_trials[subj_num==i & Condition==2 & Time==2])
    
    # Choice and RTs for congruent condition at time 1
    RT[i, 1, 1, 1:c0_t1]      <- with(tmp_data, RT[subj_num==i & Condition==0 & Time==1])
    RT_min[i, 1, 1]           <- min(RT[i, 1, 1, 1:c0_t1])
    RT_max[i, 1, 1]           <- max(RT[i, 1, 1, 1:c0_t1])
    correct[i, 1, 1, 1:c0_t1] <- with(tmp_data, Correct[subj_num==i & Condition==0 & Time==1])
    # Choice and RTs for incongruent condition at time 1
    RT[i, 2, 1, 1:c2_t1]      <- with(tmp_data, RT[subj_num==i & Condition==2 & Time==1])
    RT_min[i, 2, 1]           <- min(RT[i, 2, 1, 1:c2_t1])
    RT_max[i, 2, 1]           <- max(RT[i, 2, 1, 1:c2_t1])
    correct[i, 2, 1, 1:c2_t1] <- with(tmp_data, Correct[subj_num==i & Condition==2 & Time==1])
    # Choice and RTs for congruent condition at time 2
    RT[i, 1, 2, 1:c0_t2]      <- with(tmp_data, RT[subj_num==i & Condition==0 & Time==2])
    RT_min[i, 1, 2]           <- min(RT[i, 1, 2, 1:c0_t2])
    RT_max[i, 1, 2]           <- max(RT[i, 1, 2, 1:c0_t2])
    correct[i, 1, 2, 1:c0_t2] <- with(tmp_data, Correct[subj_num==i & Condition==0 & Time==2])
    # Choice and RTs for incongruent condition at time 2
    RT[i, 2, 2, 1:c2_t2]      <- with(tmp_data, RT[subj_num==i & Condition==2 & Time==2])
    RT_min[i, 2, 2]           <- min(RT[i, 2, 2, 1:c2_t2])
    RT_max[i, 2, 2]           <- max(RT[i, 2, 2, 1:c2_t2])
    correct[i, 2, 2, 1:c2_t2] <- with(tmp_data, Correct[subj_num==i & Condition==2 & Time==2])
  }
  
  # Stan-ready data list
  list(N       = n_subj,
       N_cond  = n_cond,
       N_time  = n_time,
       T_max   = T_max,
       T_subj  = T_subj,
       RT      = RT,
       RT_min  = RT_min,
       RT_max  = RT_max,
       correct = correct)
}

pre_iat <- function(x) {
  # Current study-task pair
  tmp_data <- x
  # Number of subjects
  n_subj <- length(unique(tmp_data$subj_num))
  subj_ids <- unique(tmp_data$subj_num)
  # Number of conditions 
  n_cond <- 2
  # Number of timepoints
  n_time <- 2
  # Determine number of trials within subjects, conditions, and timepoints
  T_subj_all <- tmp_data %>%
    group_by(subj_num, Condition, Time) %>%
    summarize(n_trials = n())
  T_max <- T_subj_all %>%
  {max(.$n_trials)}
  
  # Create RT data array for stan; dims = (subject, condition, time, trial)
  RT <- correct <- array(0, dim = c(n_subj, n_cond, n_time, T_max))
  T_subj <- RT_min <- RT_max <-array(NA, dim = c(n_subj, n_cond, n_time))
  for (i in 1:n_subj) {
    # Number of trials per condition/timpoint
    c0_t1 <- T_subj[i, 1, 1] <- with(T_subj_all, n_trials[subj_num==subj_ids[i] & Condition==1 & Time==1])
    c2_t1 <- T_subj[i, 2, 1] <- with(T_subj_all, n_trials[subj_num==subj_ids[i] & Condition==2 & Time==1])
    c0_t2 <- T_subj[i, 1, 2] <- with(T_subj_all, n_trials[subj_num==subj_ids[i] & Condition==1 & Time==2])
    c2_t2 <- T_subj[i, 2, 2] <- with(T_subj_all, n_trials[subj_num==subj_ids[i] & Condition==2 & Time==2])
    
    # Choice and RTs for congruent condition at time 1
    RT[i, 1, 1, 1:c0_t1]      <- with(tmp_data, RT[subj_num==subj_ids[i] & Condition==1 & Time==1])
    RT_min[i, 1, 1]           <- min(RT[i, 1, 1, 1:c0_t1])
    RT_max[i, 1, 1]           <- max(RT[i, 1, 1, 1:c0_t1])
    correct[i, 1, 1, 1:c0_t1] <- with(tmp_data, Correct[subj_num==subj_ids[i] & Condition==1 & Time==1])
    # Choice and RTs for incongruent condition at time 1
    RT[i, 2, 1, 1:c2_t1]      <- with(tmp_data, RT[subj_num==subj_ids[i] & Condition==2 & Time==1])
    RT_min[i, 2, 1]           <- min(RT[i, 2, 1, 1:c2_t1])
    RT_max[i, 2, 1]           <- max(RT[i, 2, 1, 1:c2_t1])
    correct[i, 2, 1, 1:c2_t1] <- with(tmp_data, Correct[subj_num==subj_ids[i] & Condition==2 & Time==1])
    # Choice and RTs for congruent condition at time 2
    RT[i, 1, 2, 1:c0_t2]      <- with(tmp_data, RT[subj_num==subj_ids[i] & Condition==1 & Time==2])
    RT_min[i, 1, 2]           <- min(RT[i, 1, 2, 1:c0_t2])
    RT_max[i, 1, 2]           <- max(RT[i, 1, 2, 1:c0_t2])
    correct[i, 1, 2, 1:c0_t2] <- with(tmp_data, Correct[subj_num==subj_ids[i] & Condition==1 & Time==2])
    # Choice and RTs for incongruent condition at time 2
    RT[i, 2, 2, 1:c2_t2]      <- with(tmp_data, RT[subj_num==subj_ids[i] & Condition==2 & Time==2])
    RT_min[i, 2, 2]           <- min(RT[i, 2, 2, 1:c2_t2])
    RT_max[i, 2, 2]           <- max(RT[i, 2, 2, 1:c2_t2])
    correct[i, 2, 2, 1:c2_t2] <- with(tmp_data, Correct[subj_num==subj_ids[i] & Condition==2 & Time==2])
  }
  
  # Stan-ready data list
  list(N       = n_subj,
       N_cond  = n_cond,
       N_time  = n_time,
       T_max   = T_max,
       T_subj  = T_subj,
       RT      = RT,
       RT_min  = RT_min,
       RT_max  = RT_max,
       correct = correct)
}

pre_ddt <- function(x) {
  # Current study-task pair
  tmp_data <- x
  # Number of subjects
  n_subj <- length(unique(tmp_data$subj_num))
  # Number of timepoints
  n_time <- 2
  # Determine number of trials within subjects and timepoints
  T_subj_all <- tmp_data %>%
    group_by(subj_num, Time) %>%
    summarize(n_trials = n())
  T_max <- T_subj_all %>%
  {max(.$n_trials)}
  
  # Create data arrays for stan; dims = (subject, time, trial)
  choice <- amount_later <- amount_sooner <- delay_later <- delay_sooner <- 
    array(0, dim = c(n_subj, n_time, T_max))
  T_subj <- array(T_max, dim = c(n_subj, n_time))
  for (i in 1:n_subj) {
    # Amount Later at time 1
    amount_later[i, 1,] <- with(tmp_data, amount_later[subj_num==i & Time==1])
    amount_later[i, 2,] <- with(tmp_data, amount_later[subj_num==i & Time==2])
    # Amount Sooner at time 1
    amount_sooner[i, 1,] <- with(tmp_data, amount_sooner[subj_num==i & Time==1])
    amount_sooner[i, 2,] <- with(tmp_data, amount_sooner[subj_num==i & Time==2])
    # Delay Later at time 1
    delay_later[i, 1,] <- with(tmp_data, delay_later[subj_num==i & Time==1])
    delay_later[i, 2,] <- with(tmp_data, delay_later[subj_num==i & Time==2])
    # Delay Sooner at time 1
    delay_sooner[i, 1,] <- with(tmp_data, delay_sooner[subj_num==i & Time==1])
    delay_sooner[i, 2,] <- with(tmp_data, delay_sooner[subj_num==i & Time==2])
    # Choices at time 1
    choice[i, 1,] <- with(tmp_data, choice[subj_num==i & Time==1])
    choice[i, 2,] <- with(tmp_data, choice[subj_num==i & Time==2])
  }
  
  # Stan-ready data list
  list(N       = n_subj,
       N_time  = n_time,
       T_max   = T_max,
       T_subj  = T_subj,
       amount_later  = amount_later,
       amount_sooner = amount_sooner,
       delay_later  = delay_later,
       delay_sooner = delay_sooner,
       choice  = choice)
}
