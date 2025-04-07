library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# source of helper functions
source("gtpreprocessing.R")
source("gtpostprocessing.R")


generate_stan <- function(components, formula, data) {
  
  # checks that inputs are as expected
  if (identical(components, c("linear", "linear"))) {
    model_frame <- model.frame(formula, data)
    y <- model.response(model_frame)
    X <- model.matrix(formula, model_frame)[, -1, drop = FALSE]
    K <- ncol(X)
    N <- nrow(X)
    
    # Fixed Stan code for linear-linear mixture
    stan_code <- paste(
      "data {",
      "  int<lower=1> N;             // Number of data points",
      "  int<lower=1> K;             // Number of predictors",
      "  matrix[N, K] X;             // Predictor matrix",
      "  vector[N] y;                // Response vector",
      "}",
      "parameters {",
      "  real<lower=0, upper=1> theta; // Mixture weight for the first component",
      "  real mu1;                     // Mean of the first component",
      "  real mu2;                     // Mean of the second component",
      "  real<lower=0> sigma1;         // Standard deviation of the first component",
      "  real<lower=0> sigma2;         // Standard deviation of the second component",
      "  vector[K] beta1;              // Regression coefficients for the first component",
      "  vector[K] beta2;              // Regression coefficients for the second component",
      "}",
      "model {",
      "  // Priors",
      "  mu1 ~ normal(0, 5);",
      "  mu2 ~ normal(0, 5);",
      "  sigma1 ~ cauchy(0, 2.5);",
      "  sigma2 ~ cauchy(0, 2.5);",
      "  beta1 ~ normal(0, 5);",
      "  beta2 ~ normal(0, 5);",
      "  theta ~ beta(1, 1);          ",
      "",
      "  // Mixture model likelihood",
      "  for (n in 1:N) {",
      "    target += log_sum_exp(",
      "      log(theta) + normal_lpdf(y[n] | mu1 + X[n] * beta1, sigma1),",
      "      log1m(theta) + normal_lpdf(y[n] | mu2 + X[n] * beta2, sigma2)",
      "    );",
      "  }",
      "}",
      "generated quantities {",
      "  int<lower=1, upper=2> z[N];    // Mixture membership",
      "  for (n in 1:N) {",
      "    // Calculate unnormalized log probabilities for each component",
      "    real log_prob1 = log(theta) + normal_lpdf(y[n] | mu1 + X[n] * beta1, sigma1);",
      "    real log_prob2 = log1m(theta) + normal_lpdf(y[n] | mu2 + X[n] * beta2, sigma2);",
      "    ",
      "    // Normalize probabilities using softmax",
      "    vector[2] prob = softmax([log_prob1, log_prob2]');",
      "    ",
      "    // Sample z[n] based on the posterior probabilities",
      "    z[n] = categorical_rng(prob);",
      "  }",
      "}",
      sep = "\n"
    )
    
    # Create filename
    stan_file <- "gtmix_linear.stan"
    writeLines(stan_code, stan_file)
    return(stan_file)
    
  } else if (identical(components, c("poisson", "poisson"))) {
    model_frame <- model.frame(formula, data)
    y <- model.response(model_frame)
    X <- model.matrix(formula, model_frame)[, -1, drop = FALSE]
    K <- ncol(X)
    N <- nrow(X)
    
    # Poisson-Poisson mixture code
    stan_code <- paste(
      "data {",
      "  int<lower=1> N;             // Number of observations",
      "  int<lower=0> y[N];          // Poisson response variable (counts)",
      "  int<lower=1> K;             // Number of predictors",
      "  matrix[N, K] X;             // Predictor matrix",
      "}",
      "",
      "parameters {",
      "  real<lower=0, upper=1> theta;           // Mixing proportions (constrained to sum to 1)",
      "  vector[K] beta1;            // Regression coefficients for component 1",
      "  vector[K] beta2;            // Regression coefficients for component 2",
      "}",
      "",
      "model {",
      "  vector[N] log_lik1;         // Log-likelihood for component 1",
      "  vector[N] log_lik2;         // Log-likelihood for component 2",
      "  ",
      "  // Linear predictors for each component",
      "  vector[N] eta1 = X * beta1; // Linear predictor for component 1",
      "  vector[N] eta2 = X * beta2; // Linear predictor for component 2",
      "",
      "  // Calculate log-likelihoods for each component",
      "  for (n in 1:N) {",
      "    log_lik1[n] = poisson_log_lpmf(y[n] | eta1[n]); // Component 1",
      "    log_lik2[n] = poisson_log_lpmf(y[n] | eta2[n]); // Component 2",
      "    target += log_mix(theta, log_lik1[n], log_lik2[n]);",
      "  }",
      "",
      "  // Priors for regression coefficients",
      "  beta1 ~ normal(0, 5);",
      "  beta2 ~ normal(0, 5);",
      "  theta ~ beta(1,1);",
      "}",
      "",
      "generated quantities {",
      "  int<lower=1, upper=2> z[N];    // Mixture membership",
      "  for (n in 1:N) {",
      "    // Calculate unnormalized log probabilities for each component",
      "    real log_prob1 = log(theta) + poisson_log_lpmf(y[n] | X[n] * beta1);",
      "    real log_prob2 = log1m(theta) + poisson_log_lpmf(y[n] | X[n] * beta2);",
      "    ",
      "    // Normalize probabilities using softmax",
      "    vector[2] prob = softmax([log_prob1, log_prob2]');",
      "    ",
      "    // Sample z[n] based on the posterior probabilities",
      "    z[n] = categorical_rng(prob);",
      "  }",
      "}",
      sep = "\n"
    )
    stan_file <- "gtmix_poisson.stan"
    writeLines(stan_code, stan_file)
    return(stan_file)
    
  } else if (identical(components, c("gamma", "gamma"))) {
    stan_code <- paste(
      "data {",
      "  int<lower=1> N;               // Number of observations",
      "  int<lower=1> K;               // Number of predictors",
      "  vector<lower=0>[N] y;         // Response variable (positive values)",
      "  matrix[N, K] X;               // Predictor matrix",
      "}",
      "",
      "parameters {",
      "  real<lower=0, upper=1> theta; // Mixing proportions (must sum to 1)",
      "  vector[K] beta1;              // Regression coefficients for component 1",
      "  vector[K] beta2;              // Regression coefficients for component 2",
      "  real<lower=0> phi1;           // Shape parameter for component 1",
      "  real<lower=0> phi2;           // Shape parameter for component 2",
      "}",
      "",
      "model {",
      "  vector[N] mu1 = exp(X * beta1);  // Mean of Gamma for component 1",
      "  vector[N] mu2 = exp(X * beta2);  // Mean of Gamma for component 2",
      "  vector[N] log_lik1;",
      "  vector[N] log_lik2;",
      "",
      "  // Calculate log-likelihoods for each component",
      "  for (n in 1:N) {",
      "    log_lik1[n] = gamma_lpdf(y[n] | phi1, phi1 / mu1[n]);",
      "    log_lik2[n] = gamma_lpdf(y[n] | phi2, phi2 / mu2[n]);",
      "    target += log_mix(theta, log_lik1[n], log_lik2[n]);",
      "  }",
      "",
      "  // Priors for regression coefficients and mix proportion",
      "  beta1 ~ normal(0, 5);",
      "  beta2 ~ normal(0, 5);",
      "  theta ~ beta(1,1);",
      "  // Priors for shape parameters",
      "  phi1 ~ exponential(1);",
      "  phi2 ~ exponential(1);",
      "}",
      "",
      "generated quantities {",
      "  int<lower=1, upper=2> z[N];    // Mixture membership",
      "  vector[N] mu1 = exp(X * beta1);  // Recompute mu1 for generated quantities",
      "  vector[N] mu2 = exp(X * beta2);  // Recompute mu2 for generated quantities",
      "  for (n in 1:N) {",
      "    // Calculate unnormalized log probabilities for each component",
      "    real log_prob1 = log(theta) + gamma_lpdf(y[n] | phi1, phi1 / mu1[n]);",
      "    real log_prob2 = log1m(theta) + gamma_lpdf(y[n] | phi2, phi2 / mu2[n]);",
      "    ",
      "    // Normalize probabilities using softmax",
      "    vector[2] prob = softmax([log_prob1, log_prob2]');",
      "    ",
      "    // Sample z[n] based on the posterior probabilities",
      "    z[n] = categorical_rng(prob);",
      "  }",
      "}",
      sep = "\n"
    )
    stan_file <- "gtmix_gamma.stan"
    writeLines(stan_code, stan_file)
    return(stan_file)
    
  } else if (identical(components, c("logistic", "logistic"))) {
    model_frame <- model.frame(formula, data)
    y <- model.response(model_frame)
    X <- model.matrix(formula, model_frame)[, -1, drop = FALSE]  # Remove intercept
    K <- ncol(X)
    N <- nrow(X)
    
    # Logistic-logistic mixture Stan code
    stan_code <- paste(
      "data {",
      "  int<lower=1> N;           // Number of observations",
      "  int<lower=1> K;           // Number of predictors",
      "  matrix[N, K] X;           // Predictor matrix",
      "  int<lower=0, upper=1> y[N]; // Binary outcome",
      "}",
      "",
      "parameters {",
      "  real<lower=0, upper=1> theta;     // Mixing proportion (for component 1)",
      "  vector[K] beta1;                 // Regression coefficients for component 1",
      "  vector[K] beta2;                 // Regression coefficients for component 2",
      "}",
      "",
      "model {",
      "  vector[N] log_lik1;  // Log-likelihood contributions from component 1",
      "  vector[N] log_lik2;  // Log-likelihood contributions from component 2",
      "",
      "  // Priors",
      "  theta ~ beta(1, 1);           // Uniform prior on mixing proportion",
      "  beta1 ~ normal(0, 5);         // Priors for regression coefficients (component 1)",
      "  beta2 ~ normal(0, 5);         // Priors for regression coefficients (component 2)",
      "",
      "  // Mixture model likelihood",
      "  for (n in 1:N) {",
      "    log_lik1[n] = bernoulli_logit_lpmf(y[n] | dot_product(X[n], beta1));",
      "    log_lik2[n] = bernoulli_logit_lpmf(y[n] | dot_product(X[n], beta2));",
      "    target += log_mix(theta, log_lik1[n], log_lik2[n]);",
      "  }",
      "}",
      "",
      "generated quantities {",
      "  int<lower=1, upper=2> z[N];      // Mixture membership for each observation",
      "  for (n in 1:N) {",
      "    vector[2] log_weights;",
      "    log_weights[1] = log(theta) + bernoulli_logit_lpmf(y[n] | dot_product(X[n], beta1));",
      "    log_weights[2] = log1m(theta) + bernoulli_logit_lpmf(y[n] | dot_product(X[n], beta2));",
      "    z[n] = categorical_rng(softmax(log_weights)); // Sample membership",
      "  }",
      "}",
      sep = "\n"
    )
    
    stan_file <- "gtmix_logistic.stan"
    writeLines(stan_code, stan_file)
    return(stan_file)
    
  } else {
    stop("Invalid mixture inputs. Must be: linear, poisson, gamma, or logistic")
  }
}

# main wrapper function
fit_model <- function(formula, p_family, data, components = NULL, result_type = 0, iterations = 10000, burning_iterations = 1000, chains = 2, seed = 123) {
  
  # Parse seed and generate random if "random" passed in
  if (seed == "random") {
    seed <- sample.int(1e6, 1)
  } else {
    seed <- as.integer(seed)
  }
  
  # Process the results based on the result_type
  if (!(result_type %in% c(0, 1))) {
    stop("Invalid result_type! Choose 0 for matrix or 1 for posterior samples.")
  }

  # Check what type of data should be loaded & handle NA values
  if (identical(data, "random")) {
    print("Generating synthetic data...")
    data <- generate_synthetic_data(p_family, seed)
    print(head(data))  # first few rows of generated data
    print(dim(data))   # dimensions of generated data
  } else {
    if (!file.exists(data)) {
      stop("Error: The specified CSV file does not exist.")
    } else {
      data <- read.csv(data)
      if (anyNA(data)) {
        na_locations <- which(is.na(data), arr.ind = TRUE)  # Get row and column locations of NA values
        stop("Error: NA values found in the data. Locations of NA values:\n", 
             paste("Row:", na_locations[, 1], "Column:", na_locations[, 2], collapse = "\n"))
      }
    }
  }

  # Check for missing data
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("The input data is not a valid data frame or is empty.")
  }

  # Check if need to generate stan file
  if (!is.null(components)) {
    stan_file <- generate_stan(components, formula, data)
  } else {
    stan_file <- switch(p_family,
                        "linear" = "gtlinear.stan",
                        "logistic" = "gtlogistic.stan",
                        "poisson" = "gtpoisson.stan",
                        "gamma" = "gtgamma.stan",
                        stop("Unknown family! Choose from: 'linear', 'logistic', 'poisson', or 'gamma'"))
  }

  # Prepare the data
  if (ncol(data) <= 1) {
    stop("Data should have at least one predictor and one response variable.")
  }
  
  #Handle the formula & prepare the data
  model_frame <- model.frame(formula, data)
  y <- model.response(model_frame)
  X <- model.matrix(formula, model_frame)[, -1]
  
  stan_data <- list(
    N = nrow(data),
    K = ncol(data) - 1,
    X = X,
    y = y
  )

  # Load the Stan model
  stan_model <- rstan::stan_model(file = stan_file)

  # Fit the model using sampling
  fit <- sampling(stan_model, data = stan_data, iter = iterations, warmup = burning_iterations, 
                  chains = chains, seed = seed)

  # Process results based on result_type
  result <- get_results(p_family, fit, result_type)
  return(result)
}


##########################################################
# Linear mixture fit test with random data

formula <- y ~ X1 + X2
fit_result <- fit_model(formula = formula,
                        p_family = "linear",
                        data = "random",
                        components = c("linear", "linear"),
                        result_type = 1,
                        iterations = 3000,
                        burning_iterations = 1000,
                        chains = 2,
                        seed = 123)