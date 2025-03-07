# # function to either return 0,1 matrix or post process results.
# # return_type: 0 for matrix, 1 for post process, default is matrix
library(rstan)

get_results <- function(family, fit, return_type) {
  posterior <- extract(fit)
  z_samples <- posterior$z            # get z-samples from fit
  beta1.p <- posterior$beta1          # get beta1 from fit
  beta2.p <- posterior$beta2          # get bata2 from fit
  
  if (family == "gamma") {
    phi1.p <- posterior$phi1
    phi2.p <- posterior$phi2
  } else {
    phi1.p <- NULL
    phi2.p <- NULL
  }

  if (is.null(posterior$z) || is.null(posterior$beta1) || is.null(posterior$beta2)) {
  stop("Error: Missing expected parameters in posterior")
  }

  # helper function for flipping z_samples and processing params
  process_parameters <- function(z_samples, beta1.p, beta2.p, phi1.p = NULL, phi2.p = NULL) {
    # Flip z_samples (0 -> 1, 1 -> 0)
    for (i in seq_len(nrow(z_samples))) {
      z_samples[i, ] <- ifelse(z_samples[i, ] == 1, 0, 1)
      if (mean(z_samples[i, ]) < 0.5) {
        z_samples[i, ] <- 1 - z_samples[i, ]
        temp <- beta1.p[i, ]
        beta1.p[i, ] <- beta2.p[i, ]
        beta2.p[i, ] <- temp

        if (!is.null(phi1.p) && !is.null(phi2.p)) {  # Check if family is gamma
          temp_phi <- phi1.p[i]
          phi1.p[i] <- phi2.p[i]
          phi2.p[i] <- temp_phi
        }
      }
    }

    if (return_type == 0) {
      return (z_samples)
    } else {
      mean_beta1 <- apply(beta1.p, 2, mean)
      mean_beta2 <- apply(beta2.p, 2, mean)
      var_beta1 <- apply(beta1.p, 2, var)
      var_beta2 <- apply(beta2.p, 2, var)
      ci_beta1 <- apply(beta1.p, 2, function(x) quantile(x, probs = c(0.25, 0.95)))
      ci_beta2 <- apply(beta2.p, 2, function(x) quantile(x, probs = c(0.25, 0.95)))
      return(list(mean_beta1 = mean_beta1, mean_beta2 = mean_beta2, var_beta1 = var_beta1, var_beta2 = var_beta2,
                  ci_beta1 = ci_beta1, ci_beta2 = ci_beta2))
    }
  }

  # Call process_parameters helper function
  result <- process_parameters(z_samples, beta1.p, beta2.p, phi1.p, phi2.p)
  return(result)
}