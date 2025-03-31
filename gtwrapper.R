library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# source of helper functions
source("gtpreprocessing.R")
source("gtpostprocessing.R")


generate_mixture_stan <- function(components, formula, data) {
  
  # checks that all models are valid
  valid_components <- c("linear", "poisson", "gamma", "logistic")
  if (!all(components %in% valid_components)) {
    stop("Invalid component type. Must be: linear, poisson, gamma, or logistic")
  }
  
  # parse al the necessary info
  model_frame <- model.frame(formula, data)
  y <- model.response(model_frame)
  X <- model.matrix(formula, model_frame)[, -1, drop = FALSE]  # Remove intercept
  K <- ncol(X)
  N <- nrow(X)
  
  response_type <- if ("logistic" %in% components) "int<lower=0, upper=1>" else "real"
  if (response_type == "int<lower=0, upper=1>" && !all(y %in% c(0, 1))) {
    stop("Logistic component requires binary response (0/1)")
  }
  
  get_component_likelihood <- function(type, i) {
    switch(type,
           "linear" = sprintf("normal_lpdf(y[n] | mu%d + X[n] * beta%d, sigma%d)", i, i, i),
           "poisson" = sprintf("poisson_log_lpmf(y[n] | dot_product(X[n], beta%d))", i),
           "gamma" = sprintf("gamma_lpdf(y[n] | shape%d, shape%d/exp(dot_product(X[n], beta%d)))", i, i, i),
           "logistic" = sprintf("bernoulli_logit_lpmf(y[n] | dot_product(X[n], beta%d))", i)
    )
  }
  
  param_blocks <- unlist(lapply(seq_along(components), function(i) {
    type <- components[i]
    switch(type,
           "linear" = c(
             sprintf("real mu%d;", i),
             sprintf("real<lower=0> sigma%d;", i),
             sprintf("vector[K] beta%d;", i)
           ),
           "poisson" = sprintf("vector[K] beta%d;", i),
           "gamma" = c(
             sprintf("real<lower=0> shape%d;", i),
             sprintf("vector[K] beta%d;", i)
           ),
           "logistic" = sprintf("vector[K] beta%d;", i)
    )
  }))
  
  prior_blocks <- unlist(lapply(seq_along(components), function(i) {
    type <- components[i]
    switch(type,
           "linear" = c(
             sprintf("mu%d ~ normal(0, 5);", i),
             sprintf("sigma%d ~ cauchy(0, 2.5);", i),
             sprintf("beta%d ~ normal(0, 5);", i)
           ),
           "poisson" = sprintf("beta%d ~ normal(0, 5);", i),
           "gamma" = c(
             sprintf("shape%d ~ gamma(0.1, 0.1);", i),
             sprintf("beta%d ~ normal(0, 5);", i)
           ),
           "logistic" = sprintf("beta%d ~ normal(0, 5);", i)
    )
  }))
  
  if (length(components) == 2) {
    likelihood_code <- sprintf(
      "for (n in 1:N) {\n    target += log_mix(theta,\n      %s,\n      %s);\n  }",
      get_component_likelihood(components[1], 1),
      get_component_likelihood(components[2], 2)
    )
  } else {
    stop("Currently only supporting 2-component mixtures")
  }
  
  gen_quant_code <- sprintf(
    "array[N] int<lower=1, upper=2> z;\n  for (n in 1:N) {\n    real log_prob1 = log(theta) + %s;\n    real log_prob2 = log1m(theta) + %s;\n    z[n] = categorical_rng(softmax([log_prob1, log_prob2]'));\n  }",
    get_component_likelihood(components[1], 1),
    get_component_likelihood(components[2], 2)
  )
  
  # final stan code
  stan_code <- paste(
    "data {",
    "  int<lower=1> N;",
    "  int<lower=1> K;",
    "  matrix[N, K] X;",
    "  vector[N] y;",
    "}",
    "",
    "parameters {",
    if (length(components) > 2) "  simplex[length(components)] theta;" else "  real<lower=0, upper=1> theta;",
    paste(" ", param_blocks, collapse = "\n"),
    "}",
    "",
    "model {",
    if (length(components) > 2) "  theta ~ dirichlet(rep_vector(1.0, length(components)));" else "  theta ~ beta(1, 1);",
    paste(" ", prior_blocks, collapse = "\n"),
    likelihood_code,
    "}",
    "",
    "generated quantities {",
    gen_quant_code,
    "}",
    sep = "\n"
  )
  
  # Create filename
  stan_file <- paste0("gtmix_", paste0(components, collapse = "_"), ".stan")
  writeLines(stan_code, stan_file)
  return(stan_file)
}

# main wrapper function
fit_model <- function(formula, p_family, data, mixture_components = NULL, result_type = 0, iterations = 10000, burning_iterations = 1000, chains = 2, seed = 123) {
  
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

  # Choose the stan file based on family
  # stan_file <- switch(p_family,
                      # "linear" = "gtlinear.stan",
                      # "logistic" = "gtlogistic.stan",
                      # "poisson" = "gtpoisson.stan",
                      # "gamma" = "gtgamma.stan",
                      # stop("Unknown family! Choose from: 'linear', 'logistic', 'poisson', or 'gamma'"))
  
  if (!is.null(mixture_components)) {
    if (!is.null(p_family)) {
      warning("Both p_family and mixture_components specified. Using mixture_components.")
    }
    stan_file <- generate_mixture_stan(mixture_components, formula, data)
    p_family <- "mixture" # can be anything really
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

fit_result <- fit_model(formula = y ~ X1 + X2,
                        "linear",
                        data = "random",
                        mixture_components = c("linear", "linear"))

##########################################################
# Linear Fit Test with Random data (temporary)

formula <- y ~ X1 + X2
fit_result <- fit_model(formula = formula,
                        p_family = "linear",
                        data = "random",
                        result_type = 1,
                        iterations = 1000,
                        burning_iterations = 1000,
                        chains = 2,
                        seed = 123)




load_data <- function(data_input) {
  if (is.data.frame(data_input)) { # Deals with data already in data.frame format
    return (data_input)
    
  } else if (data_input == "random") {
    return (data_input)
    
  } else if (file.exists(data_input)) { # Deals with csv or the rds data objects Prof Gutman mentioned
    
    if (grepl("\\.csv$", data_input)) { # csv
      data <- read.csv(data_input)
      
    } else if (grepl("\\.rds$", data_input)) { #rds
      data <- readRDS(data_input)
      if (!is.data.frame(data)) { # checks that it is a data.frame object stored inside
        stop("Error: .rds file does not contain a data.frame")
      }
    } else {
      stop("Error: Unsupported file format. Use .csv or .rds")
    }
    return(data)
  } else {
    stop("Error: Data input must be a data.frame or a valid file path")
  }
}