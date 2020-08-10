library(dplyr)
library(foreach)
library(haven)

setwd("~/Dropbox/Box/GitHub/Reliability_2020/")

# task-specific preprocessing functions
source("Code/R/0_Preprocessing/utils.R")

# Data path and individual subject file names
data_path <- "Data/0_Raw"

# Original source 
hedge     <- c("Study1-Flanker", "Study2-Flanker", "Study1-Stroop", 
               "Study2-Stroop", "Study3-Posner")
gawronski <- c("Study1a-IAT", "Study2b-IAT")
ahn       <- c("Study1-DDT")

# Names of study files
studies <- c(hedge, gawronski, ahn)

# List of headers for individual studies from Hedge et al
headers <- list(Flanker = c("Block", "Trial", "Arrow", "Condition", "Correct", "RT"),
                Stroop  = c("Block", "Trial", "Unused", "Condition", "Correct", "RT"),
                Posner  = c("Block", "Trial", "SOA", "Condition", "Correct", "RT"))

# Loop through each study-task pair
long_data <- foreach(s=seq_along(studies)) %do% {
  # For Hedge data (cognitive tasks)
  if (studies[s] %in% hedge) {
    # List of files for time 1 and 2
    files_t1 <- list.files(file.path(data_path, studies[s]), pattern = "*1.csv")
    files_t2 <- list.files(file.path(data_path, studies[s]), pattern = "*2.csv")
    
    # Get task name
    task_name <- strsplit(studies[s], split = "-")[[1]][2]
    
    # Loop through each individual subject file
    foreach(i=seq_along(files_t1), .combine = "rbind") %do% {
      # For time 1
      tmp_t1 <- read.csv(file.path(data_path, studies[s], files_t1[i]), 
                         header = F, col.names = headers[[task_name]]) %>%
        mutate(subj_num = i,
               Time = 1,
               Task = task_name)
      # For time 2 (about 3 weeks apart)
      tmp_t2 <- read.csv(file.path(data_path, studies[s], files_t2[i]), 
                         header = F, col.names = headers[[task_name]]) %>%
        mutate(subj_num = i,
               Time = 2,
               Task = task_name)
      tmp <- rbind(tmp_t1, tmp_t2)
      if (task_name=="Stopsignal") {
        tmp %>%
          mutate(Condition = ifelse(Condition == 1, 2, Condition))
      } else {
        tmp
      }
    }
    
  # For Gawronski data (self-concept and race IAT)
  } else if (studies[s] %in% gawronski) {
    # Get task name and type
    task_name <- strsplit(studies[s], split = "-")[[1]][2]
    task_type <- strsplit(studies[s], split = "-")[[1]][1]
    
    # For self-concept IAT
    if (task_type == "Study1a") {
      tmp_t1 <- read_sav(file.path(data_path, studies[s], "iat1.sav")) %>% 
        # Must join with "start" .sav files to find unique subject ID
        full_join(read_sav(file.path(data_path, studies[s], "start1.sav")), 
                  by = "Subjt1") %>%
        mutate(subj_num = subjectnum) %>%
        filter(!is.na(subj_num)) %>%
        mutate(Time = 1,
               Task = task_name)
      tmp_t2 <- read_sav(file.path(data_path, studies[s], "iat2.sav")) %>%
        full_join(read_sav(file.path(data_path, studies[s], "start2.sav")), 
                  by = "Subjt2") %>%
        mutate(subj_num = subjectnum) %>%
        filter(!is.na(subj_num)) %>%
        mutate(Time = 2,
               Task = task_name)
    # For race IAT
    } else if (task_type == "Study2b") {
      tmp_t1 <- read_sav(file.path(data_path, studies[s], "iat1.sav")) %>% 
        # Must join with "start" .sav files to find unique subject ID
        full_join(read_sav(file.path(data_path, studies[s], "start1.sav")), 
                  by = "Subj1") %>%
        mutate(subj_num = code) %>%
        filter(!is.na(subj_num)) %>%
        mutate(Time = 1,
               Task = task_name)
      tmp_t2 <- read_sav(file.path(data_path, studies[s], "iat2.sav")) %>%
        full_join(read_sav(file.path(data_path, studies[s], "start2.sav")), 
                  by = "Subj2") %>%
        mutate(subj_num = code) %>%
        filter(!is.na(subj_num)) %>%
        mutate(Time = 2,
               Task = task_name)
    }
    # Combine Time 1 and 2 into single data.frame
    tmp_t1 %>%
      bind_rows(tmp_t2) %>%
      select(subj_num, Time, RT, Block, Correct, Task) %>%
      filter(Block %in% c(6,7,11,12) & RT > 0) %>%
      # Treat "practice" and "test" trials as single blocks
      mutate(Condition = case_when(Block == 6 ~ 1,
                                   Block == 7 ~ 1,
                                   Block == 11 ~ 2,
                                   Block == 12 ~ 2)) %>%
      select(-Block) %>%
      group_by(subj_num, Time, Condition, Task) %>%
      mutate(RT = RT/1000, 
             Correct = case_when(Correct=="True" ~ 1,
                                 Correct=="False" ~ 0)) %>%
      select(subj_num, Correct, RT, Condition, Time, Task) %>%
      na.omit() %>%
      group_by(subj_num) %>%
      filter(any(Time==1) & any(Time==2))
  } else if (studies[s] %in% ahn) {
    # Get file directories for each subject at time 1 and 2
    files_t1 <- substr(list.files(file.path(data_path, studies[s]), 
                                  pattern = "_S1"), 1, 6)
    files_t2 <- substr(list.files(file.path(data_path, studies[s]), 
                                  pattern = "_S2"), 1, 6)
    # Include subjects with data at both times
    subj_ids <- files_t1[files_t1 %in% files_t2]
    
    # Extract information from each timepoint
    foreach(i=seq_along(subj_ids), .combine = "rbind") %do% {
      ### For time 1
      # Subject file path
      path_t1 <- file.path(data_path, studies[s], paste0(subj_ids[i], "_S1"), 
                           "Staircase1")
      # Staircase design
      design_t1 <- read.table(file.path(path_t1, "design_exp.txt"), 
                              sep = ",", header = T) %>%
        filter(trialNum > 4)
      # Choices
      choice_t1 <- read.table(file.path(path_t1, "choices.txt"), 
                              sep = ",", header = T) %>%
        filter(trialNum > 4)
      # Combine data and assigm data columns for modeling
      tmp_t1 <- full_join(design_t1, choice_t1, by = "trialNum") %>%
        mutate(subj_num = i,
               Time = 1,
               Task = "DDT",
               trial_num     = trialNum - 4, # non-practice trials start at 5
               amount_later  = ifelse(leftMoney == 800, leftMoney, rightMoney),
               amount_sooner = ifelse(leftMoney == 800, rightMoney, leftMoney),
               delay_later   = ifelse(leftMoney == 800, leftTime, rightTime),
               delay_sooner  = ifelse(leftMoney == 800, rightTime, leftTime),
               choice        = ifelse(leftMoney == 800, 1 - choice, choice)) %>%
        select(subj_num, Time, Task, trial_num, 
               contains("amount"), contains("delay"), choice)
      
      ### For time 2
      # Subject file path
      path_t2 <- file.path(data_path, studies[s], paste0(subj_ids[i], "_S2"), 
                           "Staircase1")
      # Staircase design
      design_t2 <- read.table(file.path(path_t2, "design_exp.txt"), 
                              sep = ",", header = T) %>%
        filter(trialNum > 4)
      # Choices
      choice_t2 <- read.table(file.path(path_t2, "choices.txt"), 
                              sep = ",", header = T) %>%
        filter(trialNum > 4)
      # Combine data and assigm data columns for modeling
      tmp_t2 <- full_join(design_t2, choice_t2, by = "trialNum") %>%
        mutate(subj_num = i,
               Time = 2,
               Task = "DDT",
               trial_num     = trialNum - 4, # non-practice trials start at 5
               amount_later  = ifelse(leftMoney == 800, leftMoney, rightMoney),
               amount_sooner = ifelse(leftMoney == 800, rightMoney, leftMoney),
               delay_later   = ifelse(leftMoney == 800, leftTime, rightTime),
               delay_sooner  = ifelse(leftMoney == 800, rightTime, leftTime),
               choice        = ifelse(leftMoney == 800, 1 - choice, choice)) %>%
        select(subj_num, Time, Task, trial_num, 
               contains("amount"), contains("delay"), choice)
      
      # Combine timepoints (about 1 month apart)
      rbind(tmp_t1, tmp_t2)
    }
  }
}
names(long_data) <- studies

# Save out long-format data
saveRDS(long_data, file = "Data/1_Preprocessed/long_format_all.rds")

# Create Stan-ready data format
stan_data <- foreach(s=studies) %do% {
  # Current study-task pair
  tmp_data <- long_data[[s]]
  
  # Determine task preprocessing function
  pre_function <- switch(tmp_data$Task[1],
                         "Flanker"    = pre_stp_flk_psr,
                         "Stroop"     = pre_stp_flk_psr,
                         "Posner"     = pre_stp_flk_psr,
                         "IAT"        = pre_iat,
                         "DDT"        = pre_ddt)
  # Run preprocessing
  pre_function(tmp_data)
}
names(stan_data) <- studies

# Save out stan-ready data
saveRDS(stan_data, file = "Data/1_Preprocessed/stan_ready_all.rds")
