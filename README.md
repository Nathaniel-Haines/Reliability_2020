# Instructions for Use

## General Information

Data, R, and Stan code to reproduce analyses and figures from "Theoretically Informed Generative Models Can Advance the Psychological and Brain Sciences: Lessons from the Reliability Paradox" by Nathaniel Haines, Peter D. Kvam, Louis Irving, Colin Tucker Smith, Theodore P. Beauchaine, Mark A. Pitt, Woo-Young Ahn, and Brandon M. Turner.

All analyses codes are in the Code/R/1_analysis/ directory. In there, individual R scripts are numbered in the order they should be run. The Stan model codes are all located in Code/Stan/. Note that fitted models are not included in the Data/2_Fitted/ directory, as they are large in size and could not be uploaded to the repository. Once all analysis scripts are run, plotting codes to reproduce the figures can be found in the Code/R/2_plotting/ directory.

The preprocessed data (located in Data/1_Preprocessed/) are all included and ready to be fit in Stan as is. However, we also include the trial-level raw data from each of the studies we used for re-analyzing the tasks in Data/0_Raw/. There, you will find each dataset in its own folder, along with a README file for each of the studies. All R scripts used to preprocess these raw data are included in Code/R/0_Preprocessing/, and they are numbered in the order in which they should be run to reproduce the results (there is only a single master script that needs to be run). See details below for more specific instructions.

Please reach out if you have any questions!

# Specific Instructions

## Necessary Software Installations

First, you will need to ensure that Stan and rstan are properly downloaded and installed. You can find a guide on installing these here: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started. 

After checking that rstan is installed properly, ensure that the following packages are installed (all versions numbers indicate those used on my own machine, but the most up-to-date versions will likely also work fine). All analyses were conducted with R version 3.6.0: 

- EnvStats_2.4.0
- tidybayes_2.3.1
- gridExtra_2.3
- scales_1.1.1
- bayesplot_1.7.0
- tidyr_1.1.3
- cowplot_0.9.4
- hBayesDM_0.7.2 
- Rcpp_1.0.6
- mvtnorm_1.1-1
- haven_2.3.1
- doParallel_1.0.14
- iterators_1.0.10
- patchwork_1.0.0
- foreach_1.4.7
- dplyr_1.0.5
- rstan_2.21.2
- StanHeaders_2.21.0-7
- ggplot2_3.3.3

## Running scripts

After installing necessary packages, codes should be run in sequence as described above. First, run pre-processing scripts, and then run analysis and plotting scripts. After starting your R environment, set the working directory to the main project directory:

`setwd("path_to_my_version/Reliability_2020")`

Next, you can begin to source the files in sequence, beginning with the pre-processing code: 

`source("Code/R/0_Preprocessing/0_make_stan_ready_all.R")`

Then, you can run the analysis scripts in order, starting with the response time models:

`source("Code/R/1_analysis/0_fit_RT_models.R")`

Run the delay discounting models and other analysis scripts in this same way. After running all analysis scripts, you can then run the plotting scripts. Some of these rely on results from the analyses in the previous stage, and others do not. You can run them in the same way as the scripts above. e.g., for the plot in the main text showing different distributions with the same mean:

`source("Code/R/2_plotting/00_Kvam_quintet.R")`

After running all the plotting scripts, the figures can be found in the "Data/3_Plotted/" directory. All fitted model objects can similarly be found in "Data/2_Fitted/".



