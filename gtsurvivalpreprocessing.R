generate_synthetic_data <- function(family, seed) {
  set.seed(seed)

  N <- 200
  X <- cbind(1, runif(N, -2, 2))  # Intercept and 1 covariate
  beta1 <- c(0.5, 1.2)
  beta2 <- c(1.5, -0.8)
  z <- rbinom(N, 1, 0.6)  # 1 = group 1, 0 = group 2
  y <- numeric(N)
  status <- rbinom(N, 1, 0.85)  # ~15% censoring

  if (family == "gamma") {
    # GAMMA survival times
    phi1 <- 3.0
    phi2 <- 2.0

    for (i in 1:N) {
      if (z[i] == 1) {
        mu <- exp(X[i, ] %*% beta1)
        y[i] <- rgamma(1, shape = phi1, rate = phi1 / mu)
      } else {
        mu <- exp(X[i, ] %*% beta2)
        y[i] <- rgamma(1, shape = phi2, rate = phi2 / mu)
      }
    }

  } else if (family == "weibull") {
    # WEIBULL survival times
    shape1 <- 2.0
    shape2 <- 1.5

    for (i in 1:N) {
      if (z[i] == 1) {
        scale <- exp(X[i, ] %*% beta1)
        y[i] <- rweibull(1, shape = shape1, scale = scale)
      } else {
        scale <- exp(X[i, ] %*% beta2)
        y[i] <- rweibull(1, shape = shape2, scale = scale)
      }
    }

  } else {
    stop("Unknown family type. Use 'gamma' or 'weibull'.")
  }

  return(data.frame(X, y = y, status = status))
}