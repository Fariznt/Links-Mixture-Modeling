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

#' Creates synthetic data for testing mixture models with multiple formulas.
#'
#' @param family Distribution family ("linear", "logistic", "poisson", or "gamma")
#' @param formulas List of model formulas (e.g., list(y1 ~ X1 + X2, y2 ~ X1 + X3))
#' @param N The number observations to generate, default 200
#' @param mix_prop Mixture proportion; the relative frequency of the more frequent component
#' in the two-component mixture. In data linkage, the true link component proportion.
#' @param seed Seed used to initialize pseudorandom number generators
#' @param beta1 Value of all true coefficients of first component (default: 1)
#' @param beta2 Value of all true coefficients of second component (default: 2)
#' @param phi1 For gamma distribution, value of all true shapes of first component. (default: 5)
#' @param phi2 For gamma distribution, value of all true shapes of second component. (default: 2)
#' @details It's especially recommended to define beta1, beta2, phi1, and phi2 beyond
#' default values when generating data for multiple formulas to mimic real-world
#' heterogeneity of response variables.
#' @return A list with:
#'   - data: data frame of predictors and responses
#'   - latent_values: list containing z and component parameters
#' @keywords internal
generate_synthetic_glm_data <- function(family, formulas, 
                                            N = 200, 
                                            mix_prop = 0.8, 
                                            seed = NULL,
                                            beta1 = NULL,
                                            beta2 = NULL,
                                            phi1 = NULL,
                                            phi2 = NULL) {
  # Argument validity checks
  if (!is.null(beta1) && length(beta1) != length(formulas)) {
    stop("beta1, if defined, must be a vector of the same length as formulas")
  }
  if (!is.null(beta2) && length(beta2) != length(formulas)) {
    stop("beta2, if defined, must be a vector of the same length as formulas")
  }
  
  if (!is.null(seed)) { # if no seed is passed, R will make one
    set.seed(seed)
  }
  
  # Correct if formulas is a single formula, not a list
  if (!is.list(formulas)) {
    formulas <- list(formulas)
  }
  
  # pull out right-side of formulas as list of character vectors
  all_rhs = lapply(formulas, function(f) all.vars(f[[3]])) 
  # combine and save unique predictors
  predictor_names <- unique(unlist(all_rhs)) 
  
  # Simulate predictors uniformly between -1 and 1
  data <- as.data.frame(
    lapply(predictor_names, function(var) runif(N, -1, 1))
  )
  names(data) <- predictor_names # assign col names
  
  # Latent mixture assignments (0 or 1)
  z <- rbinom(N, 1, mix_prop)
  
  # Prepare expected results for output
  latent_values <- list(
    z = z, # vector of latent mixture membership
    true_params = list() # list of true parameter values (coefficients, shapes) 
  )

  # Loop over each formula to generate responses
  for (f in formulas) {
    response_name <- as.character(f[[2]]) # LHS of formula
    preds <- all.vars(f[[3]]) # RHS of formula as vector
    
    # Define true coefficients for this formula
    i = which(formulas == f) # index of this formula
    
    # betas
    if (is.null(beta1)) {
      beta1.f <- rep(1, length(preds) + 1)
    } else {
      beta1.f <- rep(beta1[i], length(preds) + 1)
    }
    if (is.null(beta2)) {
      beta2.f <- rep(2, length(preds) + 1)
    } else {
      beta2.f <- rep(beta2[i], length(preds) + 1)
    }
    
    if (is.null(phi1)) {
      phi1.f <- 5
    } else {
      phi1.f <- phi1[i]
    }
    if (is.null(phi2)) {
      phi2.f <- 2
    } else {
      phi2.f <- phi2[i]
    }
    
    # Build design matrix with intercept
    X <- cbind(1, as.matrix(data[preds]))
    y <- numeric(N)
    
    # Generate outcomes per observation
    for (i in seq_len(N)) {
      eta_raw <- sum(X[i, ] * if (z[i] == 0) beta1.f else beta2.f)
      # clamp eta into [-10, 10] to keep mu and therefore rate positive finite
      eta <- pmin(pmax(eta_raw, -10), 10)
      
      mu_raw   <- exp(eta)
      # also clamp mu
      mu <- pmin(pmax(mu_raw, 1e-6), 1e6)
      
      y[i] <- switch(
        family,
        linear   = rnorm(1, mean = eta, sd = 1),
        logistic = rbinom(1, 1, prob = 1 / (1 + exp(-eta))), # TODO check
        poisson  = rpois(1, lambda = exp(eta)), # TODO check
        gamma    = rgamma( # TODO check
          1,
          shape = if (z[i] == 0) phi1.f else phi1.f,
          rate  = if (z[i] == 0) phi1.f / mu else phi2.f / mu
        ),
        stop("Invalid family: ", family)
      )
    }
    
    # Append response to data
    data[[response_name]] <- y
    
    # Record parameters used for this response
    params <- list(beta1 = beta1.f, beta2 = beta2.f)
    if (family == "gamma") {
      params$phi1 <- phi1.f
      params$phi2 <- phi2.f
    }
    latent_values$true_params[[response_name]] <- params
  }
  
  # Return both simulated dataset and ground-truth info
  return(list(
    data = data,
    latent_values = latent_values
  ))
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

#' Process model results from fit object and generate summaries
#'
#' @param family Distribution family
#' @param posterior Posterior extracted from fit object returned by stampler
#' @param formulas List of formulas used on the stan model
#' @return Processed results depending on return_type
#' @keywords internal
process_results <- function(family, posterior, formulas) {
  results <- list() # Holds all processed values
  
  # Extract results from fit
  z_samples <- posterior$z
  beta1_flat <- posterior$beta1_flat
  beta2_flat <- posterior$beta2_flat
  if (family == "gamma") {
    phi1 <- posterior$phi1
    phi2 <- posterior$phi2
  } else {
    phi1 <- NULL
    phi2 <- NULL
  }
  
  results$raw_z_samples <- z_samples
  
  # Process formulas into usable format for summary generation
  processed_formulas <- setNames(
    # list values (list of predictors)
    lapply(formulas, function(f) all.vars(f[[3]])), # f pulls out all RHS vars (predictors) in formula list
    # list keys (response vars)
    vapply(formulas, function(g) as.character(g[[2]]), "") # g pulls out all LHS vars (response vars)
  )
  
  # Process flat variables
  beta1_list <- list() # list index corresponds to response variable...
  beta2_list <- list() # ...w/ each element's col = predictor, row=iteration
  
  predictor_counts <- lengths(processed_formulas)
  response_count <- length(predictor_counts)
  pos = 1
  for (response_var in names(processed_formulas)) {
    predictor_count <- length(processed_formulas[[response_var]])
    beta1_list[[response_var]] <- beta1_flat[, pos:(pos + predictor_count - 1)]
    beta2_list[[response_var]] <- beta2_flat[, pos:(pos + predictor_count - 1)]
    pos <- pos + predictor_count
  }
  
  # Label switching (correcting switching of observations during sampling)
  for (i in seq_len(nrow(z_samples))) {
    z_samples[i, ] <- ifelse(z_samples[i, ] == 1, 0, 1)
    if (mean(z_samples[i, ]) < 0.5) {
      z_samples[i, ] <- 1 - z_samples[i, ]
      
      for (response_var in names(processed_formulas)) {
        temp_beta <- beta1_list[[response_var]][i, ]
        beta1_list[[response_var]][i, ] <- beta2_list[[response_var]][i, ]
        beta2_list[[response_var]][i, ] <- temp_beta
        if (!is.null(phi1) && !is.null(phi2)) {  # Check if family is gamma
          j = which(names(processed_formulas) == response_var) # index assoc. with this response variable
          temp_phi <- phi1[i,j]
          phi1[i,j] <- phi2[i,j]
          phi2[i,j] <- temp_phi
        }
      }
    }
  }
  
  # <-- Posterior component membership probability table -->

  # Probabilities of membership in the more frequent cluster, i.e.
  # component with the larger mixture weight. Same as probability of true link
  # in data linkage application
  component_membership_probabilities <- data.frame(
    membership_probability = colMeans(z_samples)
  )
  results$component_membership_probabilities <- component_membership_probabilities

  
  # <-- Process parameter statistics -->
  param_stats <- list()
  for (response_var in names(processed_formulas)) {
    param_stats[[response_var]] = list()
    
    # <- Process betas ->
    preds <- processed_formulas[[response_var]]
    num_of_preds <- length(preds)
    
    # Make empty data frames to load
    comp1_beta_summary <- summary_template(preds, num_of_preds) 
    comp2_beta_summary <- summary_template(preds, num_of_preds)
    
    # populate the data frames with beta summary stats by row
    for (r in seq_along(preds)) {
      # add stats as row to 1st component
      comp1_beta_summary[r, ] <- get_stats(beta1_list[[response_var]][, r])
      # add stats as row to 2nd component
      comp2_beta_summary[r, ] <- get_stats(beta2_list[[response_var]][, r])
    }
    
    # Add to results
    param_stats[[response_var]][["component1_beta_summary"]] <- comp1_beta_summary
    param_stats[[response_var]][["component2_beta_summary"]] <- comp2_beta_summary
    
    # <- Process phis ->
    if (!is.null(phi1) && !is.null(phi2)) {
      index = which(names(processed_formulas) == response_var)
      
      comp1_phi_summary <- summary_template("shape", 1)
      comp1_phi_summary[1, ] <- get_stats(phi1[[index]]) 
      
      comp2_phi_summary <- summary_template("shape", 1)
      comp2_phi_summary[1, ] <- get_stats(phi2[[index]]) 
      
      param_stats[[response_var]][["component1_phi_summary"]] <- comp1_phi_summary
      param_stats[[response_var]][["component2_phi_summary"]] <- comp2_phi_summary
    }
  }
  results <- c(results, param_stats)
  
  return(results)
}

#' process_results helper that creates data frame template for summary statistics
#' @param num_of_preds Total number of rows for which summary being made
#' @param preds Vector of row names
#' @return Empty data frame with named columns
#' @keywords internal
summary_template <- function(row_names, row_count) {
  # create empty data.frame with one row per predictor for beta params
  data.frame(
    mean = numeric(row_count),
    variance = numeric(row_count),
    lower_ci_95 = numeric(row_count),
    upper_ci_95 = numeric(row_count),
    lower_ci_50 = numeric(row_count),
    upper_ci_50 = numeric(row_count),
    row.names = row_names,
    stringsAsFactors = FALSE
  )
}


#' process_results helper that takes a numeric vector of draws for a parameter 
#' component (a real parameter itself, or an element of a array parameter, etc.)
#' and returns summary statistics as a vector values 
#' 
#' @param v A numeric vector of draws. This is part of a parameter associated with a 
#'    single response variable (plus a single predictor in case of betas)
#' @return A vector of statistics formatted as
#'        c(<mean>, <variance>, <lower 95>, <upper 95>, <lower 50>, <upper 50>)
#'        where lower 95 is the lower bound of the 95% credible interval,
#'        upper 50 is the upper bound of the 50% credible interval, etc.
#' @keywords internal
get_stats <- function(v) {
  mean <- mean(v)
  variance <- var(v)
  q95 <- quantile(v, c(0.025, 0.975))
  q50 <- quantile(v, c(0.25,  0.75))
  
  c(mean, variance, q95[1], q95[2], q50[1], q50[2])
}
