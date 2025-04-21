# Overview:

This project is designed to provide a pipeline for fitting mixture models using Stan. 
The project features the the fit_model() function in gtwrapper.R, which handles the preprocessing, model fitting, and postprocessing of the mixture models.

There are four .stan files corresponding to the logistic, poisson, linear and gamma regressions. 
The mixture models implemented in this package currently include linear-linear, poisson-poisson, gamma-gamma, and logistic-logistic models.

# Installation
To use this package, you need to have the rstan package installed.

# Functionality

This framework allows you to:
- Simulate random data or load existing data
- Define a mixture model
- Fit the model via Stan
- Postprocess and summarize results

# Parameters

fit_model(
  formula,           # Model formula (like y ~ X1 + X2)
  p_family,          # Distributional family (e.g., "linear", "logistic", "poisson", "gamma")
  data,              # "random", data.frame, or CSV/RDS filepath (use "random" for synthetic data)
  components,        # Optional: specify mixture components for a mixture model (e.g., c("linear", "poisson"))
  result_type,       # 0 = raw z-samples (0/1 matrix); 1 = postprocessed summary (posterior stats)
  iterations,        # Total MCMC samples per chain (default: 10000)
  burning_iterations,# Warmup iterations per chain (default: 1000)
  chains,            # Number of MCMC chains (default: 2)
  seed               # Random seed (can be "random" to generate a new seed, or an integer)
)

# Output

The function fit_model returns the posterior samples or a matrix of results based on the specified result_type. 

# Dependencies

rstan
parallel
gtpreprocessing.R
gtpostprocessing.R