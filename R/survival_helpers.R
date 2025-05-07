#' Creates synthetic data for testing survival models.
#'
#' @param family Distribution family ("weibull", "gamma")
#' @param seed Random seed
#' @return A data frame with synthetic data
#' @keywords internal
#' 
generate_synthetic_data <- function(family, seed) {
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
get_results <- function(family, fit, return_type = 0, status_vector = NULL) {
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
