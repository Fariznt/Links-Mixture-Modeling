#' Dynamically generate a two‐component survival‐mixture Stan program
#'
#' @param p_family     One of "weibull" or "gamma"
#' @param formula      A Surv() formula, e.g. Surv(time, status) ~ X1 + X2
#' @param data         A data.frame or the string "random"
#' @param truncation   0 or 1: whether to apply truncation
#' @param status_column Name of the column holding 0/1 censoring flags
#' @return             A single string containing the full Stan program
#' @keywords           internal
generate_survival_stan <- function(p_family,
                                   formula,
                                   data,
                                   priors,
                                   truncation    = 0,
                                   status_column = "status") {
  
  # Process prior list to get variable and function definitions for concatenation
  # during stan generation
  definitions = get_stan_definitions(priors)
  variable_definitions = definitions[["variable_defs"]]
  function_definitions = definitions[["function_defs"]]
  
  model_frame <- stats::model.frame(formula, data)
  surv_obj    <- stats::model.response(model_frame)
  time        <- surv_obj[,"time"]
  status      <- surv_obj[, status_column]
  X           <- stats::model.matrix(formula, model_frame)[, -1, drop=FALSE]
  K           <- ncol(X)
  N           <- nrow(X)
  
  if (p_family == "weibull") {
    stan_code <- paste(
      "functions {",
      function_definitions,
      "}",
      "data {",
      "  int<lower=1> N;",
      "  int<lower=1> K;",
      "  matrix[N, K] X;",
      "  vector<lower=0>[N] y;",
      "  array[N] int<lower=0,upper=1> status;",
      "  int<lower=0,upper=1> use_truncation;",
      "  real<lower=0> truncation_lower;",
      "  real<lower=0> truncation_upper;",
      "}",
      "",
      "transformed data {",
      variable_definitions, # defines any hyperparameter variables used in prior string
      "}",
      "",
      "parameters {",
      "  vector[K] beta1;",
      "  vector[K] beta2;",
      "  real<lower=0> shape1;",
      "  real<lower=0> shape2;",
      "  real<lower=0> scale1;",
      "  real<lower=0> scale2;",
      "  real<lower=0,upper=1> theta;",
      "}",
      "model {",
      "  beta1 ~ ", priors$beta1,";",
      "  beta2 ~ ", priors$beta2,";",
      "  shape1 ~ ", priors$shape1,";",
      "  shape2 ~ ", priors$shape2,";",
      "  scale1 ~ ", priors$scale1,";",
      "  scale2 ~ ", priors$scale2,";",
      "  theta ~ ", priors$theta,";",
      "  ",
      "  for (n in 1:N) {",
      "    real log_p1 = log(theta);",
      "    real log_p2 = log(1 - theta);",
      "    real linpred1 = X[n] * beta1;",
      "    real linpred2 = X[n] * beta2;",
      "    real scale_factor1 = scale1 * exp(-linpred1);",
      "    real scale_factor2 = scale2 * exp(-linpred2);",
      "    ",
      "    if (status[n] == 1) {",
      "      log_p1 += weibull_lpdf(y[n] | shape1, scale_factor1);",
      "      log_p2 += weibull_lpdf(y[n] | shape2, scale_factor2);",
      "    } else {",
      "      log_p1 += weibull_lccdf(y[n] | shape1, scale_factor1);",
      "      log_p2 += weibull_lccdf(y[n] | shape2, scale_factor2);",
      "    }",
      "    target += log_sum_exp(log_p1, log_p2);",
      "  }",
      "}",
      "generated quantities {",
      "  array[N] int<lower=1,upper=2> z;",
      "  vector[N] log_lik;",
      "  ",
      "  for (n in 1:N) {",
      "    vector[2] log_p;",
      "    real linpred1 = X[n] * beta1;",
      "    real linpred2 = X[n] * beta2;",
      "    real scale_factor1 = scale1 * exp(-linpred1);",
      "    real scale_factor2 = scale2 * exp(-linpred2);",
      "    ",
      "    log_p[1] = log(theta) + weibull_lpdf(y[n] | shape1, scale_factor1);",
      "    log_p[2] = log(1 - theta) + weibull_lpdf(y[n] | shape2, scale_factor2);",
      "    ",
      "    z[n] = categorical_rng(softmax(log_p));",
      "    log_lik[n] = log_sum_exp(log_p[1], log_p[2]);",
      "  }",
      "}",
      sep = "\n"
    )
    stan_file <- "gtmix_weibull.stan"
    
  } else {
    stan_code <- paste(
      "functions {",
      function_definitions,
      "}",
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
      "  beta1 ~ ", priors$beta1,";",
      "  beta2 ~ ", priors$beta2,";",
      "  theta ~ ", priors$theta,";",
      "  // Priors for shape parameters",
      "  phi1 ~ ", priors$phi1,";",
      "  phi2 ~ ", priors$phi2,";",
      "}",
      "generated quantities {",
      "  array[N] int<lower=1, upper=2> z; // Mixture membership",
      "  vector[N] mu1 = exp(X * beta1);  // Recompute mu1 for generated quantities",
      "  vector[N] mu2 = exp(X * beta2);  // Recompute mu2 for generated quantities",
      "  for (n in 1:N) {",
      "    real log_prob1 = log(theta) + gamma_lpdf(y[n] | phi1, phi1 / mu1[n]);",
      "    real log_prob2 = log1m(theta) + gamma_lpdf(y[n] | phi2, phi2 / mu2[n]);",
      "    vector[2] prob = softmax([log_prob1, log_prob2]');",
      "    z[n] = categorical_rng(prob);",
      "  }",
      "}",
      sep = "\n"
    )
    stan_file <- "gtmix_gamma.stan"
  }
  
  writeLines(c(stan_code, ""), stan_file)
  stan_file
}

#' Creates synthetic data for testing survival models.
#'
#' @param family Distribution family ("weibull", "gamma")
#' @param seed Random seed
#' @return A data frame with synthetic data
#' @keywords internal
#' 
generate_synthetic_survival_data <- function(family, seed) {
  set.seed(seed)
  N <- 200
  X <- cbind(1, runif(N, -2, 2)) # Intercept and 1 covariate
  colnames(X) <- c("X1", "X2")
  
  # True parameters (shared for both distributions)
  beta1 <- c(0.5, 1.2) # Intercept and slope for group 1
  beta2 <- c(1.5, -0.8) # Intercept and slope for group 2
  z_true <- rbinom(N, 1, 0.6) # 1 = group 1, 0 = group 2
  y <- numeric(N)
  
  if (family == "gamma") {
    # Gamma parameters (shape = phi)
    phi1 <- 3.0
    phi2 <- 2.0
    
    for (i in 1:N) {
      if (z_true[i] == 1) {
        mu <- exp(X[i, ] %*% beta1)
        y[i] <- rgamma(1, shape = phi1, rate = phi1 / mu)
      } else {
        mu <- exp(X[i, ] %*% beta2)
        y[i] <- rgamma(1, shape = phi2, rate = phi2 / mu)
      }
    }
  } else if (family == "weibull") {
    # Weibull parameters
    shape1 <- 2.0 # alpha for group 1
    shape2 <- 1.5 # alpha for group 2
    
    for (i in 1:N) {
      if (z_true[i] == 1) {
        scale <- exp(-X[i, ] %*% beta1)
        y[i] <- rweibull(1, shape = shape1, scale = scale)
      } else {
        scale <- exp(-X[i, ] %*% beta2)
        y[i] <- rweibull(1, shape = shape2, scale = scale)
      }
    }
  } else {
    stop("Unknown family type. Use 'gamma' or 'weibull'.")
  }
  
  # Realistic censoring mechanism (time-dependent)
  # Changed quantile range to be within [0,1]
  censoring_time <- runif(N,
                          min = quantile(y, 0.4),
                          max = quantile(y, 0.9)
  ) # Now properly between 0 and 1
  status <- as.numeric(y <= censoring_time) # 1 = event, 0 = censored
  y <- pmin(y, censoring_time) # Observed time is min(event, censoring)
  
  return(data.frame(X, y = y, status = status, z_true = z_true))
}

#' Process model results
#'
#' @param family Distribution family
#' @param fit Stan fit object
#' @param return_type 0 for matrix output, 1 for posterior samples
#' @param status_vector vector for status column
#' @return Processed results
#' @keywords internal
get_survival_results <- function(family, fit, return_type = 0, status_vector = NULL) {
  posterior <- extract(fit)
  
  # Core parameters
  z_samples <- posterior$z
  beta1.p <- posterior$beta1
  beta2.p <- posterior$beta2
  
  # Family-specific
  if (family == "gamma") {
    phi1.p <- posterior$phi1
    phi2.p <- posterior$phi2
  } else if (family == "weibull") {
    shape1.p <- posterior$shape1
    shape2.p <- posterior$shape2
  } else {
    stop("Unsupported family. Use 'gamma' or 'weibull'.")
  }
  
  if (is.null(z_samples) || is.null(beta1.p) || is.null(beta2.p)) {
    stop("Missing required parameters (z, beta1, beta2) in Stan output.")
  }
  
  process_parameters <- function() {
    n_iter <- nrow(z_samples)
    
    for (i in seq_len(n_iter)) {
      z_samples[i, ] <- 1 - z_samples[i, ]
      if (mean(z_samples[i, ]) < 0.5) {
        z_samples[i, ] <- 1 - z_samples[i, ]
        tmp_beta <- beta1.p[i, ]
        beta1.p[i, ] <- beta2.p[i, ]
        beta2.p[i, ] <- tmp_beta
        
        if (family == "gamma") {
          tmp_phi <- phi1.p[i]
          phi1.p[i] <- phi2.p[i]
          phi2.p[i] <- tmp_phi
        }
        
        if (family == "weibull") {
          tmp_shape <- shape1.p[i]
          shape1.p[i] <- shape2.p[i]
          shape2.p[i] <- tmp_shape
        }
      }
    }
    
    stats <- list(
      mean_beta1 = apply(beta1.p, 2, mean),
      mean_beta2 = apply(beta2.p, 2, mean),
      var_beta1 = apply(beta1.p, 2, var),
      var_beta2 = apply(beta2.p, 2, var),
      ci_beta1 = apply(beta1.p, 2, function(x) quantile(x, probs = c(0.25, 0.95))),
      ci_beta2 = apply(beta2.p, 2, function(x) quantile(x, probs = c(0.25, 0.95)))
    )
    
    if (family == "gamma") {
      stats$mean_phi1 <- mean(phi1.p)
      stats$mean_phi2 <- mean(phi2.p)
      stats$ci_phi1 <- quantile(phi1.p, probs = c(0.25, 0.95))
      stats$ci_phi2 <- quantile(phi2.p, probs = c(0.25, 0.95))
    }
    
    if (family == "weibull") {
      stats$mean_shape1 <- mean(shape1.p)
      stats$mean_shape2 <- mean(shape2.p)
      stats$ci_shape1 <- quantile(shape1.p, probs = c(0.25, 0.95))
      stats$ci_shape2 <- quantile(shape2.p, probs = c(0.25, 0.95))
    }
    
    stats$z_samples <- z_samples
    
    # Include posterior mode assignment per observation
    if (!is.null(status_vector)) {
      z_mode <- apply(z_samples, 2, function(x) {
        as.integer(round(mean(x)))  # modal class approx
      })
      stats$z_assignment <- z_mode
      stats$status <- status_vector
      
      stats$cross_tab <- table(Assignment = z_mode, Status = status_vector)
    }
    
    return(stats)
  }
  
  stats <- process_parameters()
  
  # Output
  if (return_type == 0) {
    cat("Posterior z_samples:\n")
    print(stats$z_samples)
  } else if (return_type == 1) {
    cat("Post-processed beta statistics:\n")
    print(stats[c("mean_beta1", "mean_beta2", "ci_beta1", "ci_beta2")])
    
    if (family == "gamma") {
      cat("Gamma phi stats:\n")
      print(stats[c("mean_phi1", "mean_phi2", "ci_phi1", "ci_phi2")])
    }
    if (family == "weibull") {
      cat("Weibull shape stats:\n")
      print(stats[c("mean_shape1", "mean_shape2", "ci_shape1", "ci_shape2")])
    }
    
    if (!is.null(status_vector)) {
      cat("Assignment vs. Status Table:\n")
      print(stats$cross_tab)
    }
  } else {
    stop("Invalid return_type: Use 0 (matrix) or 1 (summary).")
  }
  
  return(stats)
}