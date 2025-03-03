library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("gtpostprocessing.R")

# generates synthetic data (encapsulates previous preprocessing gtrun.R functionality in a function)
generate_synthetic_data <- function(family, seed) {
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



# main wrapper function
fit_model <- function(formula, p_family, data, iterations, burning_iterations, chains, seed, result_type) {
  
  # Parses seed and generates random if "random" passed in
  if (seed == "random") {
    seed <- sample.int(1e6, 1)
  } else {
    seed <- as.integer(seed)
  }
  
  # Checks what type of data should be loaded & NA values
  if (identical(data, "random")) {
    source("gtrun.R")
    data <- generate_synthetic_data(p_family, seed)
    
  } else {
    if (!is.data.frame(data)) {
      stop("Error: Data must be a data frame.")
    }
    if (anyNA(data)) {
      na_locations <- which(is.na(data), arr.ind = TRUE)  # Get row and column locations of NA values
      stop("Error: NA values found in the data. Locations of NA values:\n", 
            paste("Row:", na_locations[, 1], "Column:", na_locations[, 2], collapse = "\n"))
    }
  }


  # Choose the post-processing function based on family
  if (p_family == "linear") {
    stan_file <- "gtlinear.stan"
  } else if (p_family == "logistic") {
    stan_file <- "gtlogistic.stan"
  } else if (p_family == "poisson") {
    stan_file <- "gtpoisson.stan"
  } else if (p_family == "gamma") {
    stan_file <- "gtgamma.stan"
  } else {
    stop("Unknown family! Choose from: 'linear', 'logistic', 'poisson', or 'gamma'")
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
  
  # Process the results based on the result_type
  if (result_type == 0) {
    # return 0/1 matrix (link)
    result_type <- get_results(p_family, fit, return_type = 0)
  } else if (result_type == 1) {
    # Return posterior samples for the betas
    result_type <- get_results(p_family, fit, return_type = 1)
  } else {
    stop("Invalid result_type! Choose 0 for matrix or 1 for posterior samples.")
  }

  return(list(fit = fit, result_type = result_type))
}


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

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript gtwrapper.R formula family data (optional: iterations, burning_iterations, chains, seed, result_type)")
}

# assign argument inputs 
formula <- as.formula(args[1])
p_family <- as.character(args[2])
print(p_family)
data <- load_data(args[3])

# assign argument inputs or default values if necessary
iterations <- ifelse(length(args) < 4, 10000, as.integer(args[4]))
burning_iterations <- ifelse(length(args) < 5, 1000, as.integer(args[5]))
chains <- ifelse(length(args) < 6, 2, as.integer(args[6]))
seed <- ifelse(length(args) < 7, 123, ifelse(args[7] == "random", sample.int(1e6, 1), as.integer(args[7])))
result_type <- ifelse(length(args) < 8, 1, as.integer(args[8]))

# Run the model
fit_result <- fit_model(formula, p_family, data, iterations, burning_iterations, chains, seed, result_type)

print(fit_result$fit)
print(fit_result$result_type)
