# occPlus 0.1.0

* Initial CRAN release.
* Multi-species occupancy model for environmental DNA (eDNA) read-count data,
  following the statistical framework of Griffin et al. (2020).
* Two model-fitting backends:
    - `runOccPlus()`: Stan-based Hamiltonian Monte Carlo via **rstan**.
    - `runMCMCOccPlus()`: Custom Gibbs sampler implemented in R/C++ (Rcpp).
* Accounts for false-positive and false-negative observation errors at both the
  field-sampling and laboratory-amplification stages of eDNA surveys.
* Latent factor model captures species co-occurrence patterns and enables
  concurrent ordination alongside occupancy estimation.
* Suite of visualisation functions for occupancy rates, detection rates,
  covariate effects, read-intensity distributions, and inter-species correlation
  matrices.
* Two vignettes: a quickstart guide (`occPlus`) and a model-interpretation
  walkthrough (`model-interpretation`).
