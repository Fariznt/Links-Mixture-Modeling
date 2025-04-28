# Overview:
This project provides a comprehensive pipeline for fitting Bayesian mixture models using Stan, including both standard regression and survival analysis. 

The framework consists of two main functions:

- fit_model() in gtwrapper.R - Handles preprocessing, fitting, and postprocessing for:
    Linear, logistic, Poisson, and Gamma regression
    Mixtures of these models (linear-linear, Poisson-Poisson, etc.)

- fit_survival_model() in gtsurvival.R - Specialized for survival analysis with:
    Weibull and Gamma survival models
    Support for censored data and time truncation

The package includes corresponding .stan files for each model type, implementing efficient MCMC sampling.

# Installation
To use this package, you need to have the rstan package installed. Furthermore for gtsurvival, you need to have the survival package installed.
To do this run the following commands in your terminal:
install.packages(c("rstan", "survival"))

# Functionality
This framework enables you to:

Simulate random data or load existing datasets
- Define mixture components (including survival models)
- Fit models via Stan with customizable MCMC settings
- Postprocess results including:
- Parameter estimates and credible intervals
- Component assignment probabilities
- Model diagnostics

For survival models specifically:
- Handle right-censored survival data
- Account for time truncation when needed
- Analyze mixture survival distributions


# Function Summary

[gtwrapper.R]
fit_model <- function(
  formula,                    # Model formula (e.g., y ~ X1 + X2)
  p_family,                   # Distribution family: "linear", "logistic", "poisson", or "gamma"
  data,                       # "random", data.frame, or filepath to CSV/RDS
  components = NULL,          # Optional: vector of components for mixture models (e.g., c("linear", "poisson"))
  result_type = 0,            # 0 = raw samples, 1 = postprocessed summary
  iterations = 10000,         # Total MCMC iterations per chain
  burning_iterations = 1000,  # Warmup iterations per chain
  chains = 2,                 # Number of MCMC chains
  seed = 123                  # Random seed ("random" or integer)
) 
  """
  Fit a Bayesian regression model with optional mixture components
  
  Fits various GLMs using Stan with support for:
  - Linear, logistic, Poisson, and Gamma regression
  - Mixture models when components are specified
  - Synthetic data generation
  
  Arguments:
    formula         Model formula specifying response and predictors.
                    Example: y ~ X1 + X2
                    
    p_family        Base distribution family:
                    - "linear": Normal/Gaussian regression
                    - "logistic": Binomial logistic regression  
                    - "poisson": Poisson regression
                    - "gamma": Gamma regression
                    
    data            Input data specification:
                    - "random": generates synthetic data
                    - data.frame: uses provided data
                    - character: path to data file (CSV/RDS)
                    
    components      Optional vector specifying mixture components.
                    Example: c("linear", "poisson") creates a 2-component
                    mixture model. When NULL, fits standard GLM.
                    
    result_type     Output format:
                    - 0: Returns raw posterior samples matrix
                    - 1: Returns processed summary statistics
                    
    iterations      Total MCMC iterations per chain (default: 10000)
    
    burning_iterations  Warmup iterations per chain (default: 1000)
    
    chains          Number of MCMC chains (default: 2)
    
    seed            Random seed specification:
                    - "random": generates new random seed
                    - integer: uses fixed seed for reproducibility
  
  Returns:
    For result_type=0: Matrix of posterior samples
    For result_type=1: List containing:
      - Parameter estimates (posterior means and CIs)
      - Model fit statistics
      - For mixtures: component assignments
  
  Example Usage:
  fit_result <- fit_model(formula = y ~ X1 + X2,
                        p_family = "linear",
                        data = "random",
                        components = c("linear", "linear"),
                        result_type = 1,
                        iterations = 3000,
                        burning_iterations = 1000,
                        chains = 50,
                        seed = 123)
    


[gtsurvival.R]
fit_survival_model <- function(
  formula,                    # Survival formula (like Surv(time, status) ~ X1 + X2)
  p_family,                   # Distribution family: "weibull" or "gamma"
  data,                       # "random", data.frame, or filepath (CSV/RDS)
  result_type,                # 0 = raw samples matrix, 1 = postprocessed summary
  iterations = 10000,         # Total MCMC iterations per chain
  burning_iterations = 1000,  # Warmup iterations per chain
  chains = 2,                 # Number of MCMC chains
  seed = "random",            # Random seed ("random" or integer)
  truncation = 0,             # 1 = enable truncation, 0 = disable (default)
  status_column = NULL        # Name of status column if not in formula
) 
  
  Function: Fit a Bayesian survival model with optional truncation
  
  Fits either a Weibull or Gamma survival model using Stan, with support for:
  - Right-censored data
  - Time truncation (when truncation=1)
  - Mixture components (via p_family specification)
  
  Arguments:
    formula         Survival formula using Surv() syntax. Example: 
                    Surv(time, status) ~ age + treatment
                    
    p_family        Distribution family for survival times:
                    - "weibull": Weibull proportional hazards model
                    - "gamma": Gamma accelerated failure time model
                    
    data            Input data: 
                    - "random": generates synthetic data
                    - data.frame: existing dataset  
                    - character: path to CSV/RDS file
                    
    result_type     Output format:
                    - 0: Returns raw posterior samples (z matrix)
                    - 1: Returns postprocessed summary statistics
                    
    iterations      Total MCMC iterations per chain (default: 10000)
    
    burning_iterations  Warmup iterations per chain (default: 1000)
    
    chains          Number of MCMC chains (default: 2)
    
    seed            Random seed:
                    - "random": generates new seed
                    - integer: fixed seed for reproducibility
                    
    truncation      Time truncation:
                    - 0: No truncation (default)
                    - 1: Apply truncation between (0, max(y))
                    
    status_column   Optional name of status column if not specified in formula
  
  Returns:
    For result_type=0: Matrix of posterior samples
    For result_type=1: List containing:
      - Parameter estimates (means and credible intervals)
      - Cluster assignments
      - Model fit statistics

  Example Usage:
    result <- fit_survival_model(
      formula = Surv(y, status) ~ X2, 
      p_family = "weibull",
      data = "random",
      result_type = 1,
      iterations = 500,
      burning_iterations = 250,
      chains = 2,
      seed = "random",
      truncation = 0,
      status_column = "status")