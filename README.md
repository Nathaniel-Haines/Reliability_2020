# Reliability_2020
Data, R, and Stan code to reproduce analyses and figures from "Learning from the Reliability Paradox: How Theoretically Informed Generative Models Can Advance the Social, Behavioral, and Brain Sciences" by Nathaniel Haines, Peter D. Kvam, Louis Irving, Colin Tucker Smith, Theodore P. Beauchaine, Mark A. Pitt, Woo-Young Ahn, and Brandon M. Turner.

All analyses codes are in the Code/R/1_analysis/ directory. In there, individual R scripts are numbered in the order they should be run. The Stan model codes are all located in Code/Stan/. Note that fitted models are not included in the Data/2_Fitted/ directory, as they are large in size and could not be uploaded to the repository. Once all analysis scripts are run, plotting codes to reproduce the figures can be found in the Code/R/2_plotting/ directory.

The preprocessed data (located in Data/1_Preprocessed/) are all included and ready to be fit in Stan as is. However, we also include the trial-level raw data from each of the studies we used for re-analyzing the tasks in Data/0_Raw/. There, you will find each dataset in its own folder, along with a README file for each of the studies. All R scripts used to preprocess these raw data are included in Code/R/0_Preprocessing/, and they are numbered in the order in which they should be run to reproduce the results (there is only a single master script that needs to be run).

Please reach out if you have any questions!
