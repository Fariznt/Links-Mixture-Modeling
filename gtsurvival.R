library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# source of helper functions
source("gtpreprocessing.R")
source("gtpostprocessing.R")

# main wrapper function for survival models with truncation param
fit_survival_model <- function(formula, p_family, data, result_type, iterations, burning_iterations, chains, seed, truncation = 0) {
  
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
  stan_file <- switch(p_family,
                      "weibull" = "gtweibull.stan",
                      "gamma" = "gtgamma.stan",
                      stop("Unknown family! Choose from: 'weibull' or 'gamma'"))

  # Prepare the data
  if (ncol(data) <= 1) {
    stop("Data should have at least one predictor and one response variable.")
  }
  
  # Handle the formula & prepare the data
  model_frame <- model.frame(formula, data)
  y <- model.response(model_frame)
  X <- model.matrix(formula, model_frame)[, -1]
  
  # Add truncation if truncation = 1
  if (truncation == 1) {
    truncation_lower <- 0       # lower truncation bound
    truncation_upper <- max(y)  # upper truncation bound, could adjust this based on the data or user input
    stan_data <- list(
      N = nrow(data),
      K = ncol(data) - 1,
      X = X,
      y = y,
      truncation_lower = truncation_lower,
      truncation_upper = truncation_upper
    )
  } else {
    stan_data <- list(
      N = nrow(data),
      K = ncol(data) - 1,
      X = X,
      y = y
    )
  }

  # Load the Stan model
  stan_model <- rstan::stan_model(file = stan_file)

  # Fit the model using sampling
  fit <- sampling(stan_model, data = stan_data, iter = iterations, warmup = burning_iterations, 
                  chains = chains, seed = seed)

  # Process results based on result_type
  result <- get_results(p_family, fit, result_type)
  return(result)
}
