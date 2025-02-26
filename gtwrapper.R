library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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

# function to either return 0,1 matrix or post process results.
# return_type: 0 for matrix, 1 for post process, default is matrix
get_results <- function(family, fit, return_type) {
  posterior <- extract(fit)
  
  # helper function to process the common steps for each family
  process_parameters <- function(z_samples, beta1.p, beta2.p, phi1.p = NULL, phi2.p = NULL) {
    # Flip z_samples (0 -> 1, 1 -> 0)
    for (i in 1:nrow(z_samples)) {
      z_samples[i, ] <- ifelse(z_samples[i, ] == 1, 0, 1)
      if (mean(z_samples[i, ]) < 0.5) {
        z_samples[i, ] <- 1 - z_samples[i, ]
        temp <- beta1.p[i, ]
        beta1.p[i, ] <- beta2.p[i, ]
        beta2.p[i, ] <- temp
        if (family == "gamma") {
          temp_phi <- phi1.p[i]
          phi1.p[i] <- phi2.p[i]
          phi2.p[i] <- temp_phi
        }
      }
    }
    
    # Compute posterior statistics for beta1 and beta2
    mean_beta1 <- apply(beta1.p, 2, mean)
    mean_beta2 <- apply(beta2.p, 2, mean)
    var_beta1 <- apply(beta1.p, 2, var)
    var_beta2 <- apply(beta2.p, 2, var)
    ci_beta1 <- apply(beta1.p, 2, function(x) quantile(x, probs = c(0.25, 0.95)))
    ci_beta2 <- apply(beta2.p, 2, function(x) quantile(x, probs = c(0.25, 0.95)))
    
    # If the family is gamma, include the phi parameters
    if (family == "gamma") {
      mean_phi1 <- mean(phi1.p)
      mean_phi2 <- mean(phi2.p)
      var_phi1 <- var(phi1.p)
      var_phi2 <- var(phi2.p)
      ci_phi1 <- quantile(phi1.p, probs = c(0.25, 0.95))
      ci_phi2 <- quantile(phi2.p, probs = c(0.25, 0.95))
      return(list(mean_beta1 = mean_beta1, mean_beta2 = mean_beta2,
                  var_beta1 = var_beta1, var_beta2 = var_beta2,
                  ci_beta1 = ci_beta1, ci_beta2 = ci_beta2,
                  mean_phi1 = mean_phi1, mean_phi2 = mean_phi2,
                  var_phi1 = var_phi1, var_phi2 = var_phi2,
                  ci_phi1 = ci_phi1, ci_phi2 = ci_phi2))
    } else {
      return(list(mean_beta1 = mean_beta1, mean_beta2 = mean_beta2,
                  var_beta1 = var_beta1, var_beta2 = var_beta2,
                  ci_beta1 = ci_beta1, ci_beta2 = ci_beta2))
    }
  }
  
  # Common return type for 0-1 matrix
  if (return_type == 0) {
    if (family == "gamma") {
      return(matrix(c(1, 1, 1, 1), nrow = 1, ncol = 4))  # For gamma family
    } else {
      return(matrix(c(1, 1), nrow = 1, ncol = 2))  # For other families
    }
  }
  
  # Only perform post-processing if return_type == 1
  if (return_type == 1) {
    z_samples <- posterior$z
    beta1.p <- posterior$beta1
    beta2.p <- posterior$beta2
    phi1.p <- ifelse(family == "gamma", posterior$phi1, NULL)
    phi2.p <- ifelse(family == "gamma", posterior$phi2, NULL)
    
    return(process_parameters(z_samples, beta1.p, beta2.p, phi1.p, phi2.p))
  }
}



# main wrapper function
fit_model <- function(formula, family, data, iterations, burning_iterations, chains, seed, result_type) {
  
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
  
  # Process the results based on the result_type
  if (result_type == 0) {
    # return 0/1 matrix (link)
    result_type <- get_results(family, fit, return_type = 0)
  } else if (result_type == 1) {
    # Return posterior samples for the betas
    result_type <- get_results(family, fit, return_type = 1)
  } else {
    stop("Invalid result_type! Choose 0 for matrix or 1 for posterior samples.")
  }

  return(list(fit = fit, result_type = result_type))
}


args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript gtwrapper.R formula family data (optional: iterations, burning_iterations, chains, seed, result_type)")
}

# assign default values
default_iterations <- 10000
default_burning_iterations <- 1000
default_chains <- 2
default_seed <- 123
default_result_type <- 1  # Default to 1 (posterior samples)

# assign argument inputs 
formula <- as.formula(args[1])
family <- args[2]
data <- args[3]
iterations <- ifelse(length(args) < 4, default_iterations, as.integer(args[4]))
burning_iterations <- ifelse(length(args) < 5, default_burning_iterations, as.integer(args[5]))
chains <- ifelse(length(args) < 6, default_chains, as.integer(args[6]))
seed <- ifelse(length(args) < 7, default_seed, ifelse(args[7] == "random", sample.int(1e6, 1), as.integer(args[7])))
result_type <- ifelse(length(args) < 8, default_result_type, as.integer(args[8]))

# Run the model
fit_result <- fit_model(formula, family, data, iterations, burning_iterations, chains, seed, result_type)

print(fit_result$fit)
print(fit_result$result_type)
