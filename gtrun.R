## Generates synthetic data for linear, logistic, poisson and gamma models, then fits stan model for each and gives results.
## Note that results directly from stan model are not usable as they may have labels switched. instead we post-process them first and compute posterior means, variances and confidence limits manually

library(rstan)

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
