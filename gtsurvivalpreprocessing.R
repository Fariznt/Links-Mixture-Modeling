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
