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
