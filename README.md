# LAS: Log-adjusted shrinkage prior for signal estimation

This repository provides R code implementing spatially clustered regression for spatial data analysis, as proposed by the following paper.


[Hamura, H., Irie, K. and Sugasawa, S. (2020). Shrinkage with robustness: log-adjusted heavy-tailed priors. *arXiv:2001.08465*](https://arxiv.org/abs/2001.08465)

The repository includes the following files.

* `LAS-function.R` : Script implementing the proposed method
* `Simulation.R` : Script applying the proposed prios and some existing ones to simulated datasets 
* `Data_Analysis.R` : Script applying shrinkage priors to summary statistics of prostate cancer data
* `prostz.RData` : Prostate cancer data (used in `Data_Analysis.R`)


