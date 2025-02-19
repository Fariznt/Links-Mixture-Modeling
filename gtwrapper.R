library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# main wrapper function
fit_model <- function(formula, family, data, iterations, burning_iterations, chains, seed) {
  
  # Parses seed and generates random if "random" passed in
  if (seed == "random") {
    seed <- sample.int(1e6, 1)
  } else {
    seed <- as.integer(seed)
  }
  
  # Checks what type of data should be loaded & NA values
  if (identical(data, "random")) {
    source("gtrun.R")
    data <- generate_synthetic_data(family, seed)
    
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
  
  
  # Choose the post-processing function based on family
  if (family == "linear") {
    stan_file <- file.path("gtlinear.stan")   # relative paths for now
  } else if (family == "logistic") {
    stan_file <- file.path("gtlogistic.stan")
  } else if (family == "poisson") {
    stan_file <- file.path("gtpoisson.stan")
  } else if (family == "gamma") {
    stan_file <- file.path("gtgamma.stan")
  } else {
    stop("Unknown family! Choose from: 'linear', 'logistic', 'poission', or 'gamma'")
  }

  # Check for missing or NA values in the dataset:
  # returns an array of all the missing values with the same size as the data set
  missing_data <- which(is.na(data), arr.ind = TRUE)

  if (length(missing_data) > 0) {
    stop(
      "Error: Missing values found in the following rows and columns:\n",
      paste(apply(missing_data, 1, function(x) paste("Row", x[1], "Column", x[2])), collapse = "\n")
    )
  }

  # Prepare the data
  stan_data <- list(
    N = nrow(data), # Number of observations
    K = ncol(data) - 1, # Number of predictors
    X = as.matrix(data[, -ncol(data)]), # Predictor matrix
    y = data[, ncol(data)] # Response vector
  )

  # Load the Stan model
  stan_model <- rstan::stan_model(file = stan_file)

  # Fit the model using sampling
  fit <- sampling(stan_model, data = stan_data, iter = iterations, warmup = burning_iterations, chains = chains, seed = seed)
  
  # TODO: figure out what needs to be returned from postprocessing function

  return(list(fit = fit))
}


args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript gtwrapper.R formula family data (optional: iterations, burning_iterations, chains, seed)")
}

# assign default values
default_iterations <- 10000
default_burning_iterations <- 1000
default_chains <- 2
default_seed <- 123

# assign argument inputs 
formula <- as.formula(args[1])
family <- args[2]
data <- args[3]
iterations <- ifelse(length(args) < 4, default_iterations, as.integer(args[4]))
burning_iterations <- ifelse(length(args) < 5, default_burning_iterations, as.integer(args[5]))
chains <- ifelse(length(args) < 6, default_chains, as.integer(args[6]))
seed <- ifelse(length(args) < 7, default_seed, ifelse(args[7] == "random", sample.int(1e6, 1), as.integer(args[7])))


# Run the model
fit_result <- fit_model(formula, family, data, iterations, burning_iterations, chains, seed)

print(fit_result$fit)
