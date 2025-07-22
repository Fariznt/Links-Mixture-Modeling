#' Fits Bayesian mixture models using Stan with various distribution families (package loader)
#' @keywords internal 
#' "_PACKAGE"
#' @name LinksMixtureModeling
#' @import rstan
#' @importFrom stats model.frame model.response model.matrix
#' @importFrom utils read.csv
NULL

.onLoad <- function(libname, pkgname) {
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
}

#' Creates synthetic data for testing mixture models.
#'
#' @param family Distribution family ("linear", "logistic", "poisson", "gamma")
#' @param seed Random seed
#' @return A data frame with synthetic data
#' @keywords internal
generate_synthetic_mixture_data <- function(family, seed, priors) {
  set.seed(seed)
  
  N <- 200
  K <- 2
  
  if (family == "linear") {
    # LINEAR data generation
    beta1 <- c(1, 2)
    beta2 <- c(3, 4)
    sigma1 <- 1
    sigma2 <- 2
    
    X <- cbind(1, runif(N, -2, 2))
    z <- rbinom(N, size = 1, prob = 0.8) + 1  # true link (z=1) or false link (z=2)
    z <- ifelse(z == 1, 0, 1)
    y <- numeric(N)
    
    for (i in 1:N) {
      if (z[i] == 1) {
        y[i] <- rnorm(1, mean = sum(X[i, ] * beta1), sd = sigma1)
      } else {
        y[i] <- rnorm(1, mean = sum(X[i, ] * beta2), sd = sigma2)
      }
    }
    
  } else if (family == "logistic") {
    # LOGISTIC data generation
    beta1 <- c(0.5, -1)
    beta2 <- c(0.8, 1)
    X <- matrix(rnorm(N * K), nrow = N, ncol = K)
    z <- rbinom(N, size = 1, prob = 0.6) + 1
    logit <- function(x) exp(x) / (1 + exp(x))
    y <- numeric(N)
    
    for (i in 1:N) {
      if (z[i] == 1) {
        y[i] <- rbinom(1, 1, prob = logit(sum(X[i, ] * beta1)))
      } else {
        y[i] <- rbinom(1, 1, prob = logit(sum(X[i, ] * beta2)))
      }
    }
    
  } else if (family == "poisson") {
    # POISSON data generation
    beta1 <- c(1, 2)
    beta2 <- c(2, 3)
    X <- cbind(1, runif(N, -2, 2))
    lambda1 <- exp(X %*% beta1)
    lambda2 <- exp(X %*% beta2)
    z <- rbinom(N, 1, 0.6)
    y <- ifelse(z == 1, rpois(N, lambda1), rpois(N, lambda2))
    
  } else if (family == "gamma") {
    # GAMMA data generation
    beta1 <- c(0.5, 1.2)
    phi1 <- 5
    beta2 <- c(2.0, -0.8)
    phi2 <- 2
    X <- cbind(1, runif(N, -2, 2))
    z <- rbinom(N, 1, 0.6)
    y <- numeric(N)
    
    for (i in 1:N) {
      if (z[i] == 1) {
        mu <- exp(X[i, ] %*% beta1)
        y[i] <- rgamma(1, shape = phi1, rate = phi1 / mu)
      } else {
        mu <- exp(X[i, ] %*% beta2)
        y[i] <- rgamma(1, shape = phi2, rate = phi2 / mu)
      }
    }
    
  } else {
    stop("Unknown family type. Choose 'linear', 'logistic', 'poisson', or 'gamma'.")
  }
  
  return(data.frame(X, y))
}


#' Generate stan code for mixture models
#'
#' @param components Vector of component types
#' @param formula Model formula
#' @param data Input data
#' @return Path to generated Stan file
#' @keywords internal
generate_stan <- function(components, formula, data, priors) {
  
  # Process prior list to get variable and function definitions for concatenation
  # during stan generation
  definitions = get_stan_definitions(priors)
  variable_definitions = definitions[["variable_defs"]]
  function_definitions = definitions[["function_defs"]]
  
  # checks that inputs are as expected
  if (identical(components, c("linear", "linear"))) {
    model_frame <- stats::model.frame(formula, data)
    y <- stats::model.response(model_frame)
    X <- stats::model.matrix(formula, model_frame)[, -1, drop = FALSE]
    K <- ncol(X)
    N <- nrow(X)
  
    # Fixed Stan code for linear-linear mixture
    stan_code <- paste(
      "functions {",
      function_definitions,
      "}",
      "data {",
      "  int<lower=1> N;                 // Number of data points",
      "  int<lower=1> M;                 // Number of response variables",
      "  array[M] int<lower=1> K_per_m;  // Array of predictor counts for each response",
      "  int<lower=1> K_sum;             // Sum of all values in K; total predictor count",
      "  matrix[N, M] y_all;             // Matrix of M response variables",
      "  matrix[N, K_sum] X_all;        // All predictor matrices combined column-wise",
      "}",
      "",
      "transformed data {",
      variable_definitions, # defines any hyperparameter variables used in prior string
      "}",
      "",
      "parameters {",
      "  real<lower=0, upper=1> theta;  // Mixture weight for the first component",
      "  vector[M] mu1;                 // Means of the first component",
      "  vector[M] mu2;                 // Means of the second component",
      "  vector<lower=0>[M] sigma1;     // Standard deviations of the first component",
      "  vector<lower=0>[M] sigma2;     // Standard deviations of the second component",
      "  array[M] vector[K[m]] beta1;   // Regression coefficients for the first component",
      "  array[M] vector[K[m]] beta2;   // Regression coefficients for the second component",
      "}",
      "model {",
      "  // Priors",
      "  mu1 ~ ", priors$mu1, ";",
      "  mu2 ~ ", priors$mu2, ";",
      "  sigma1 ~ ", priors$sigma1, ";",
      "  sigma2 ~ ", priors$sigma2, ";",
      "  beta1 ~ ", priors$beta1, ";",
      "  beta2 ~ ", priors$beta2, ";",
      "  theta ~ ", priors$theta, ";",
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
      "transformed data {",
      variable_definitions, # defines any hyperparameter variables used in prior string
      "}",
      "",
      "parameters {",
      "  real<lower=0, upper=1> theta;           // Mixing proportions (constrained to sum to 1)",
      "  vector[K] beta1;            // Regression coefficients for component 1",
      "  vector[K] beta2;            // Regression coefficients for component 2",
      "}",
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
      "  beta1 ~ ", priors$beta1, ";",
      "  beta2 ~ ", priors$beta2, ";",
      "  theta ~ ", priors$theta, ";",
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
      "transformed data {",
      variable_definitions, # defines any hyperparameter variables used in prior string
      "}",
      "",
      "parameters {",
      "  real<lower=0, upper=1> theta; // Mixing proportions (must sum to 1)",
      "  vector[K] beta1;              // Regression coefficients for component 1",
      "  vector[K] beta2;              // Regression coefficients for component 2",
      "  real<lower=0> phi1;           // Shape parameter for component 1",
      "  real<lower=0> phi2;           // Shape parameter for component 2",
      "}",
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
      "  beta1 ~ ", priors$beta1, ";",
      "  beta2 ~ ", priors$beta2, ";",
      "  theta ~ ", priors$theta, ";",
      "  // Priors for shape parameters",
      "  phi1 ~ ", priors$phi1, ";",
      "  phi2 ~ ", priors$phi2, ";",
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
      "transformed data {",
      variable_definitions, # defines any hyperparameter variables used in prior string
      "}",
      "",
      "parameters {",
      "  real<lower=0, upper=1> theta;     // Mixing proportion (for component 1)",
      "  vector[K] beta1;                 // Regression coefficients for component 1",
      "  vector[K] beta2;                 // Regression coefficients for component 2",
      "}",
      "model {",
      "  vector[N] log_lik1;  // Log-likelihood contributions from component 1",
      "  vector[N] log_lik2;  // Log-likelihood contributions from component 2",
      "",
      "  // Priors",
      "  // Uniform prior on mixing proportion:",
      "  theta ~ ", priors$theta, ";",
      "  // Prior for regression coefficients (component 1):",
      "  beta1 ~ ", priors$beta1, ";",
      "  // Prior for regression coefficients (component 2):",
      "  beta2 ~ ", priors$beta2, ";",
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

#' Process model results from fit object
#'
#' @param family Distribution family
#' @param fit Stan fit object
#' @param return_type 0 for matrix output, 1 for posterior samples
#' @return Processed results
#' @keywords internal
get_mixture_results <- function(family, fit, return_type) {
  posterior <- extract(fit)
  z_samples <- posterior$z            # get z-samples from fit
  beta1.p <- posterior$beta1          # get beta1 from fit
  beta2.p <- posterior$beta2          # get bata2 from fit
  
  if (family == "gamma") {
    phi1.p <- posterior$phi1
    phi2.p <- posterior$phi2
  } else {
    phi1.p <- NULL
    phi2.p <- NULL
  }

  if (is.null(posterior$z) || is.null(posterior$beta1) || is.null(posterior$beta2)) {
  stop("Error: Missing expected parameters in posterior")
  }

  # helper function for flipping z_samples and processing params
  process_parameters <- function(z_samples, beta1.p, beta2.p, phi1.p = NULL, phi2.p = NULL) {
    # Flip z_samples (0 -> 1, 1 -> 0)
    for (i in seq_len(nrow(z_samples))) {
      z_samples[i, ] <- ifelse(z_samples[i, ] == 1, 0, 1)
      if (mean(z_samples[i, ]) < 0.5) {
        z_samples[i, ] <- 1 - z_samples[i, ]
        temp <- beta1.p[i, ]
        beta1.p[i, ] <- beta2.p[i, ]
        beta2.p[i, ] <- temp

        if (!is.null(phi1.p) && !is.null(phi2.p)) {  # Check if family is gamma
          temp_phi <- phi1.p[i]
          phi1.p[i] <- phi2.p[i]
          phi2.p[i] <- temp_phi
        }
      }
    }
    
      # Calculate statistics
      mean_beta1 <- apply(beta1.p, 2, mean)
      mean_beta2 <- apply(beta2.p, 2, mean)
      var_beta1 <- apply(beta1.p, 2, var)
      var_beta2 <- apply(beta2.p, 2, var)
      ci_beta1 <- apply(beta1.p, 2, function(x) quantile(x, probs = c(0.25, 0.95)))
      ci_beta2 <- apply(beta2.p, 2, function(x) quantile(x, probs = c(0.25, 0.95)))

    if (return_type == 0) {
      print("Printing z_samples matrix:")
      print(z_samples)  # Print the matrix z_samples
    } else if (return_type == 1) {
      print("Printing betas and related statistics:")
      # Print each beta with a descriptive label
      cat("Mean beta1: ", mean_beta1, "\n")
      cat("Mean beta2: ", mean_beta2, "\n")
      cat("Variance beta1: ", var_beta1, "\n")
      cat("Variance beta2: ", var_beta2, "\n")
      cat("95% CI beta1: ", ci_beta1, "\n")
      cat("95% CI beta2: ", ci_beta2, "\n")
    }
    
    return(list(
        z_samples = z_samples,
        mean_beta1 = mean_beta1,
        mean_beta2 = mean_beta2,
        var_beta1 = var_beta1,
        var_beta2 = var_beta2,
        ci_beta1 = ci_beta1,
        ci_beta2 = ci_beta2
      ))
  }
  # Call process_parameters helper function
  result_list <- process_parameters(z_samples, beta1.p, beta2.p, phi1.p, phi2.p)
  return(result_list)
}
