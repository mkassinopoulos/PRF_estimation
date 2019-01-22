# PRF_estimation

This repository provides scripts for our work presented in:

Kassinopoulos, M., Mitsis, G.D., 2019. Identification of Physiological Response Functions to Correct for Fluctuations in Resting-State fMRI related to Heart Rate and Respiration. bioRxiv 512855. https://doi.org/10.1101/512855

In this work, we present a novel framework for estimating physiological response functions (PRFs) in resting-state fMRI and we demonstrate the feasibility of this framework on data from the Human Connectome Project (HCP). Specifically, we use physiological recordings acquired during the scan as well as the global signal (i.e. average timeseries across all voxels in the brain) to extract the cardiac and respiration response functions. These PRFs are subsequently used in the analysis to correct for BOLD fluctuations induced by changes in heart rate and respiration. In addition, in this work we present probabilistic maps of brain areas affected by cardiac pulsatility, heart rate and breathing pattern derived from 41 healthy young subjects of the HCP. These maps can be found in the following SharePoint link:

https://mcgill-my.sharepoint.com/:f:/g/personal/michalis_kassinopoulos_mail_mcgill_ca/EnY1Q8ny18ZKjEB-I_eqXJ0BT8zm4TsfGC6IIx06H9GswA?e=X7FTBg


## Extraction of PRF curves and associated physiological regressors
*PRF_sc_optimize_parameters.m* is a function that receives the heart rate, respiratory flow (variable extracted from respiratory signal), the global signal, and the parameters of the PRFs and returns the corresponding PRF curves as well as the regressors that can be incorporated later in the general linear model as nuissance regressors.

*PRF_sc_HCP.m* is a script that demonstrates how the PRF_sc_optimize_parameters.m is used to find the optimal parameters using numerical optimization techniques. This script uses data from HCP that can be found in the SharePoint link shown earlier. *PRF_sc_HCP.m* consists of three sections. In the first section, the user selects a scan out of 164 scans from 41 subjects (these subjects were chosen based on the quality of their physiological recordings). The second section finds the optimal parameters using a genetic algorithm as well as a gradient-based algorithm. Finally, the last section presents the PRF curves and the associated physiological regressors of the model which can be exported by the user in a textfile and used in fMRI software toolboxes (e.g. FSL, SPM).


## Author

Michalis Kassinopoulos, PhD candidate
Graduate Program in Biological and Biomedical Engineering, McGill University, Montreal, Canada
mkassinopoulos@gmail.com

Date: 15-Jan-2019

Please do not hesitate to contact me if you have any questions related to the use of these scripts.

Michalis
