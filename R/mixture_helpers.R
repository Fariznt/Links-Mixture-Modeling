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
    
    # Linear-linear mixture code
    stan_code <- paste(
      "functions {",
      function_definitions,
      "}",
      "data {",
      "  int<lower=1> N;                 // Number of data points",
      "  int<lower=1> M;                 // Number of response variables",
      "  matrix[N, M] y_all;             // Matrix of M response variables",
      "  array[M] int<lower=1> K_per_m;  // Array of predictor counts for each response",
      "  int<lower=1> K_sum;             // Sum of all values in K; total predictor count",
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
      "  vector[K_sum] beta1_flat;      // Regression coefficients for the first component",
      "  vector[K_sum] beta2_flat;      // Regression coefficients for the second component",
      "}",
      "model {",
      "  // Priors",
      "  mu1 ~ ", priors$mu1, ";",
      "  mu2 ~ ", priors$mu2, ";",
      "  sigma1 ~ ", priors$sigma1, ";",
      "  sigma2 ~ ", priors$sigma2, ";",
      "  beta1_flat ~", priors$beta1, ";",
      "  beta2_flat ~", priors$beta2, ";",
      "  theta ~ ", priors$theta, ";",
      "",
      "  // Mixture model likelihood",
      "  for (n in 1:N) {",
      "    real log_prob_comp1 = log(theta);",
      "    real log_prob_comp2 = log1m(theta);",
      "    int current_pos = 1;",
      "    for (m in 1:M) {",
      "      row_vector[K_per_m[m]] X_n_m = X_all[n, current_pos:(current_pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] beta1_m = beta1_flat[current_pos:(current_pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] beta2_m = beta2_flat[current_pos:(current_pos + K_per_m[m] - 1)];",
      "      log_prob_comp1 += normal_lpdf(y_all[n, m] | mu1[m] + X_n_m * beta1_m, sigma1[m]);",
      "      log_prob_comp2 += normal_lpdf(y_all[n, m] | mu2[m] + X_n_m * beta2_m, sigma2[m]);",
      "      current_pos += K_per_m[m];",
      "    }",
      "    target += log_sum_exp(log_prob_comp1, log_prob_comp2);",
      "  }",
      "}",
      "generated quantities {",
      "  int<lower=1, upper=2> z[N];",
      "  for (n in 1:N) {",
      "    real log_prob1 = log(theta);",
      "    real log_prob2 = log1m(theta);",
      "    int current_pos = 1;",
      "    for (m in 1:M) {",
      "      row_vector[K_per_m[m]] X_n_m = X_all[n, current_pos:(current_pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] beta1_m = beta1_flat[current_pos:(current_pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] beta2_m = beta2_flat[current_pos:(current_pos + K_per_m[m] - 1)];",
      "      log_prob1 += normal_lpdf(y_all[n, m] | mu1[m] + X_n_m * beta1_m, sigma1[m]);",
      "      log_prob2 += normal_lpdf(y_all[n, m] | mu2[m] + X_n_m * beta2_m, sigma2[m]);",
      "      current_pos += K_per_m[m];",
      "    }",
      "    vector[2] log_probs;",
      "    log_probs[1] = log_prob1;",
      "    log_probs[2] = log_prob2;",
      "    vector[2] prob = softmax(log_probs);",
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
    # Poisson-Poisson mixture code
    stan_code <- paste(
      "functions {",
      function_definitions,
      "}",
      "data {",
      "  int<lower=1> N;                 // Number of data points",
      "  int<lower=1> M;                 // Number of response variables",
      "  int<lower=0> y_all[N,M];        // Poisson counts (N obs × M responses)",
      "  array[M] int<lower=1> K_per_m;  // Array of predictor counts for each response",
      "  int<lower=1> K_sum;             // Sum of all values in K; total predictor count",
      "  matrix[N, K_sum] X_all;         // All predictor matrices combined column-wise",
      "}",
      "",
      "transformed data {",
      variable_definitions,
      "}",
      "",
      "parameters {",
      "  real<lower=0, upper=1> theta;  // Mixing proportions (constrained to sum to 1)",
      "  vector[K_sum] beta1_flat;      // concatenated betas for comp 1",
      "  vector[K_sum] beta2_flat;      // concatenated betas for comp 2",
      "}",
      "model {",
      "  // Priors",
      "  beta1_flat ~ ", priors$beta1, ";",
      "  beta2_flat ~ ", priors$beta2, ";",
      "  theta ~ ", priors$theta, ";",
      "",
      "  // Mixture model likelihood",
      "  for (n in 1:N) {",
      "    real log_prob_comp1 = log(theta);",
      "    real log_prob_comp2 = log1m(theta);",
      "    int current_pos = 1;",
      "    for (m in 1:M) {",
      "      row_vector[K_per_m[m]] X_n_m = X_all[n, current_pos:(current_pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] beta1_m = beta1_flat[current_pos:(current_pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] beta2_m = beta2_flat[current_pos:(current_pos + K_per_m[m] - 1)];",
      "",
      "      // log-rate (eta) for each component",
      "      real eta1 = dot_product(X_n_m, beta1_m);",
      "      real eta2 = dot_product(X_n_m, beta2_m);",
      "",
      "      //add each equation's contribution",
      "      log_prob_comp1 += poisson_log_lpmf(y_all[n, m] | eta1);",
      "      log_prob_comp2 += poisson_log_lpmf(y_all[n, m] | eta2);",
      "",
      "      current_pos += K_per_m[m];",
      "    }",
      "    target += log_sum_exp(log_prob_comp1, log_prob_comp2);",
      "  }",
      "}",
      "",
      "generated quantities {",
      "  int<lower=1, upper=2> z[N];    // Mixture membership",
      "  for (n in 1:N) {",
      "  // initialize log‑probs for the two components",
      "    real log_prob1 = log(theta);",
      "    real log_prob2 = log1m(theta);",
      "    int current_pos = 1;",
      "",
      "    // accumulate each formula’s contribution",
      "    for (m in 1:M) {",
      "      row_vector[K_per_m[m]] X_n_m = X_all[n, current_pos:(current_pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] beta1_m = beta1_flat[current_pos:(current_pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] beta2_m = beta2_flat[current_pos:(current_pos + K_per_m[m] - 1)];",
      "",
      "      real eta1 = dot_product(X_n_m, beta1_m);",
      "      real eta2 = dot_product(X_n_m, beta2_m);",
      "",
      "      log_prob1 += poisson_log_lpmf(y_all[n, m] | eta1);",
      "      log_prob2 += poisson_log_lpmf(y_all[n, m] | eta2);",
      "",
      "      current_pos += K_per_m[m];",
      "    }",
      "    vector[2] log_probs = [log_prob1, log_prob2]';",
      "",
      "    // Normalize probabilities using softmax",
      "    vector[2] prob = softmax(log_probs);",
      "",
      "    // Sample z[n] based on the posterior probabilities",
      "    z[n] = categorical_rng(prob);",
      "  }",
      "}",
      sep = "\n"
    )
    stan_file <- "gtmix_poisson.stan"
    writeLines(c(stan_code, ""), stan_file)
    return(stan_file)
    
  } else if (identical(components, c("gamma", "gamma"))) {
    stan_code <- paste(
      "functions {",
      function_definitions,
      "}",
      "data {",
      "  int<lower=1> N;                 // Number of data points",
      "  int<lower=1> M;                 // Number of response variables",
      "  matrix<lower=0>[N, M] y_all;     // Positive responses for each obs × formula",
      "  array[M] int<lower=1> K_per_m;   // Array of predictor counts per response",
      "  int<lower=1> K_sum;             // Total number of predictors",
      "  matrix[N, K_sum] X_all;         // All predictor matrices combined",
      "}",
      "",
      "transformed data {",
      variable_definitions,          # any hyper‐params you defined
      "}",
      "",
      "parameters {",
      "  real<lower=0,upper=1> theta;      // Mixture weight",
      "  vector[K_sum] beta1_flat;        // Flattened betas for comp1",
      "  vector[K_sum] beta2_flat;        // Flattened betas for comp2",
      "  vector<lower=0>[M] phi1;         // Shape params for comp1",
      "  vector<lower=0>[M] phi2;         // Shape params for comp2",
      "}",
      "",
      "model {",
      "  // Priors",
      "  beta1_flat ~ ", priors$beta1, ";",
      "  beta2_flat ~ ", priors$beta2, ";",
      "  theta ~ ", priors$theta, ";",
      "  phi1 ~ ", priors$phi1, ";",
      "  phi2 ~ ", priors$phi2, ";",
      "",
      "  // Mixture‐model likelihood",
      "  for (n in 1:N) {",
      "    real log_prob_comp1 = log(theta);",
      "    real log_prob_comp2 = log1m(theta);",
      "    int current_pos = 1;",
      "    for (m in 1:M) {",
      "      row_vector[K_per_m[m]] X_n_m = ",
      "        X_all[n, current_pos:(current_pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] b1 = ",
      "        beta1_flat[current_pos:(current_pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] b2 = ",
      "        beta2_flat[current_pos:(current_pos + K_per_m[m] - 1)];",
      "      real mu1 = exp(dot_product(X_n_m, b1));",
      "      real mu2 = exp(dot_product(X_n_m, b2));",
      "      log_prob_comp1 += gamma_lpdf(y_all[n, m] | phi1[m], phi1[m] / mu1);",
      "      log_prob_comp2 += gamma_lpdf(y_all[n, m] | phi2[m], phi2[m] / mu2);",
      "      current_pos += K_per_m[m];",
      "    }",
      "    target += log_sum_exp(log_prob_comp1, log_prob_comp2);",
      "  }",
      "}",
      "",
      "generated quantities {",
      "  int<lower=1,upper=2> z[N];        // Mixture membership",
      "  for (n in 1:N) {",
      "    real log_prob1 = log(theta);",
      "    real log_prob2 = log1m(theta);",
      "    int current_pos = 1;",
      "    for (m in 1:M) {",
      "      row_vector[K_per_m[m]] X_n_m = ",
      "        X_all[n, current_pos:(current_pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] b1 = ",
      "        beta1_flat[current_pos:(current_pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] b2 = ",
      "        beta2_flat[current_pos:(current_pos + K_per_m[m] - 1)];",
      "      real mu1 = exp(dot_product(X_n_m, b1));",
      "      real mu2 = exp(dot_product(X_n_m, b2));",
      "      log_prob1 += gamma_lpdf(y_all[n, m] | phi1[m], phi1[m] / mu1);",
      "      log_prob2 += gamma_lpdf(y_all[n, m] | phi2[m], phi2[m] / mu2);",
      "      current_pos += K_per_m[m];",
      "    }",
      "    vector[2] log_probs = [log_prob1, log_prob2]';",
      "    z[n] = categorical_rng(softmax(log_probs));",
      "  }",
      "}",
      sep = "\n"
    )
    stan_file <- "gtmix_gamma.stan"
    writeLines(stan_code, stan_file)
    return(stan_file)
    
  } else if (identical(components, c("logistic", "logistic"))) {
    # Logistic–Logistic mixture code (multi‐formula)
    stan_code <- paste(
      "functions {",
      function_definitions,
      "}",
      "data {",
      "  int<lower=1> N;                 // Number of observations",
      "  int<lower=1> M;                 // Number of response variables",
      "  int<lower=0,upper=1> y_all[N,M]; // Binary outcomes for each obs×response",
      "  array[M] int<lower=1> K_per_m;  // Predictors per response",
      "  int<lower=1> K_sum;             // Total predictors across all formulas",
      "  matrix[N, K_sum] X_all;         // Combined predictor matrix (col‑wise)",
      "}",
      "",
      "transformed data {",
      variable_definitions,
      "}",
      "",
      "parameters {",
      "  real<lower=0,upper=1> theta;    // Mixing weight",
      "  vector[K_sum] beta1_flat;       // All comp‑1 betas concatenated",
      "  vector[K_sum] beta2_flat;       // All comp‑2 betas concatenated",
      "}",
      "",
      "model {",
      "  // Priors",
      "  beta1_flat ~ ", priors$beta1, ";",
      "  beta2_flat ~ ", priors$beta2, ";",
      "  theta      ~ ", priors$theta,  ";",
      "",
      "  // Mixture likelihood",
      "  for (n in 1:N) {",
      "    real log_p1 = log(theta);",
      "    real log_p2 = log1m(theta);",
      "    int pos = 1;",
      "    for (m in 1:M) {",
      "      row_vector[K_per_m[m]] X_nm = ",
      "        X_all[n, pos:(pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] b1 = ",
      "        beta1_flat[pos:(pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] b2 = ",
      "        beta2_flat[pos:(pos + K_per_m[m] - 1)];",
      "      log_p1 += bernoulli_logit_lpmf(y_all[n,m] | ",
      "                          dot_product(X_nm, b1));",
      "      log_p2 += bernoulli_logit_lpmf(y_all[n,m] | ",
      "                          dot_product(X_nm, b2));",
      "      pos += K_per_m[m];",
      "    }",
      "    target += log_sum_exp(log_p1, log_p2);",
      "  }",
      "}",
      "",
      "generated quantities {",
      "  int<lower=1,upper=2> z[N];",
      "  for (n in 1:N) {",
      "    real lp1 = log(theta);",
      "    real lp2 = log1m(theta);",
      "    int pos = 1;",
      "    for (m in 1:M) {",
      "      row_vector[K_per_m[m]] X_nm = ",
      "        X_all[n, pos:(pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] b1 = ",
      "        beta1_flat[pos:(pos + K_per_m[m] - 1)];",
      "      vector[K_per_m[m]] b2 = ",
      "        beta2_flat[pos:(pos + K_per_m[m] - 1)];",
      "      lp1 += bernoulli_logit_lpmf(y_all[n,m] | ",
      "                         dot_product(X_nm, b1));",
      "      lp2 += bernoulli_logit_lpmf(y_all[n,m] | ",
      "                         dot_product(X_nm, b2));",
      "      pos += K_per_m[m];",
      "    }",
      "    vector[2] lps = [lp1, lp2]';",
      "    vector[2] probs = softmax(lps);",
      "    z[n] = categorical_rng(probs);",
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
  beta1.p <- posterior$beta1_flat          # get beta1 from fit
  beta2.p <- posterior$beta2_flat          # get bata2 from fit
  
  if (family == "gamma") {
    phi1.p <- posterior$phi1
    phi2.p <- posterior$phi2
  } else {
    phi1.p <- NULL
    phi2.p <- NULL
  }
  
  if (is.null(posterior$z) || is.null(posterior$beta1) || is.null(posterior$beta2)) {
    stop("Error: Missing expected parameters in posterior", posterior$z, posterior$beta1, posterior$beta2, "end") #TODO fix this line
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
