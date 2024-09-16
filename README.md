# VERGE

This repository contains the code for the VERGE (Varying Effects Regression with Graph Estimation) algorithm as described in the manuscript accepted by Biometrics:

**"Bayesian network-guided sparse regression with flexible varying effects."**


Authors: Yangfan Ren, Christine B. Peterson, and Marina Vannucci

## Overview

The VERGE algorithm is implemented in Matlab and is a Bayesian hierarchical regression model that enables the selection of both network-linked predictor variables and covariates that modify the predictor effects.

## Contents

- `run_all.m`: the main file for running the VERGE algorithm
- `MCMC_all.m`: the MCMC algorithm
- `data_gen.m`: the file for generating simulation data
- `pred.m`: the file for prediction on test data
- `functions`: all the helper functions that are used in MCMC (need to be in the same folder with other files to run MCMC)
- `Combo`: files for Combo case study. The original source of this data was the COMBO study described in Wu et al. (2011). Here, we share the version of the data as processed by Zhang et al. (2021).
    - `run_all.m`: the main file for running the Combo data (needs to be used with MCMC_all.m and helper functions)
    - `check_combo.R`: results and plots for the case study
    - `X.csv, Y.csv, Z.csv`: predictors, response, and covariates for the Combo data
    - `results_combo.m`: saved results for the case study

## References

- Wu, G. D., Chen, J., Hoffmann, C., Bittinger, K., Chen, Y.-Y., et al. (2011). Linking long-term dietary patterns with gut microbial enterotypes. *Science* 334, 105–108.  
- Zhang, L., Shi, Y., Jenq, R. R., Do, K.-A., et al. (2021). Bayesian compositional regression with structured priors for microbiome feature selection. *Biometrics* 77, 824–838.

## Copyright

Yangfan Ren
