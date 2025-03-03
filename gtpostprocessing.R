# function to either return 0,1 matrix or post process results.
# return_type: 0 for matrix, 1 for post process, default is matrix
get_results <- function(family, fit, return_type) {
  posterior <- extract(fit)
  
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