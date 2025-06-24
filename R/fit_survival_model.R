#' Fits a survival model using Stan with various distribution families.
#'
#' @param formula Model formula (e.g., y ~ x1 + x2)
#' @param p_family Distribution family ("weibull", "gamma")
#' @param data Input data frame or "random" for synthetic data
#' @param hyperparameters Named list (or NULL). Expected keys:  
#'   `beta  = list(mean, sd)`,  
#'   `theta = list(alpha, beta)`,  
#'   `phi   = list(rate)` (gamma survival only),  
#'   `shape = list(alpha, beta)`, `scale = list(alpha, beta)` (weibull survival only).  
#'   mean, sd, alpha, beta, rate are length-2 vectors (one value per mixture component);  
#'   if a scalar is supplied it is used for both components.  
#'   NULL or missing elements trigger weakly-informative default priors.
#' @param result_type 0 for matrix output, 1 for posterior samples,
#' @param iterations Total number of MCMC iterations, default is
#' @param burning_iterations Number of burn-in iterations
#' @param chains Number of MCMC chains
#' @param seed Random seed (integer or "random")
#' @param truncation 0 or 1 indicating whether to use truncated distributions
#' @param status_column indicates which column is the status one
#' @return Results depending on result_type
#' @export
#' @importFrom rstan stan_model sampling
#' @importFrom stats model.frame model.response model.matrix
#' @importFrom utils read.csv
#' @importFrom survival Surv
#' @examples
#' \dontrun{
#' # Weibull mixture example
#' surv_fit <- LinksMixtureModeling::fit_survival_model(
#' formula = Surv(y, status) ~ X2, 
#' p_family = "weibull",
#' data = "random",
#' result_type = 1,
#' iterations = 500,
#' burning_iterations = 250,
#' chains = 2,
#' seed = "random",
#' truncation = 0,
#' status_column = "status") }

fit_survival_model <- 
  function(formula, 
           p_family, 
           data, 
           hyperparameters = NULL,
           result_type, 
           iterations, 
           burning_iterations, 
           chains, 
           seed, 
           truncation = 0, 
           status_column = NULL) {
    
    # Define default hyperparameter values
    defaults <- list(
      beta  = list(mean  = c(0, 0), sd = c(5, 5)),
      theta = list(alpha = 1, beta  = 1),
      phi   = list(rate  = c(1, 1)),
      shape = list(alpha = c(2, 2), beta = c(1, 1)), # alpha is shape, beta is rate
      scale = list(alpha = c(2, 2), beta = c(1, 1)) # alpha is shape, beta is rate
    )
    
    # Set default hyperparameter values if not passed in
    if (is.null(hyperparameters)) { # if nothing passed in, use all defaults
      hyperparameters <- defaults
    } else if (!is.list(hyperparameters)) { # very simple validation of prior format
      stop("Invalid argument: 'priors' must be a named list of lists.")
    } else { # use hyperparameters passed in, with default values for missing ones
      hyperparameters <- utils::modifyList(defaults, hyperparameters)
    }
    
  # Parse seed and generate random if "random" passed in
  if (seed == "random") {
    seed <- sample.int(1e6, 1)
  } else {
    seed <- as.integer(seed)
  }
  
  # Process the results based on the result_type
  if (!(result_type %in% c(0, 1))) {
    stop("Invalid result_type! Choose 0 for matrix or 1 for posterior samples.")
  }
  
  # Check what type of data should be loaded & handle NA values
  if (identical(data, "random")) {
    print("Generating synthetic data...")
    data <- generate_synthetic_survival_data(p_family, seed)
    print(head(data)) # first few rows of generated data
    print(dim(data)) # dimensions of generated data
  } else {
    if (!file.exists(data)) {
      stop("Error: The specified CSV file does not exist.")
    } else {
      data <- read.csv(data)
      if (anyNA(data)) {
        na_locations <- which(is.na(data), arr.ind = TRUE) # Get row and column locations of NA values
        stop(
          "Error: NA values found in the data. Locations of NA values:\n",
          paste("Row:", na_locations[, 1], "Column:", na_locations[, 2], collapse = "\n")
        )
      }
    }
  }
  
  # Check for missing data
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("The input data is not a valid data frame or is empty.")
  }
  
  # Ensure the status column is handled
  if (!is.null(status_column)) {
    if (status_column %in% colnames(data)) {
      data$status <- ifelse(data[[status_column]] == 0, 0, 1) # 0 if censored, 1 if event timed out
    } else {
      stop("The specified status column does not exist in the data.")
    }
  } else {
    stop("Please specify the 'status_column' parameter.")
  }
  
  pkg_path <- system.file(package = "LinksMixtureModeling")
  
  # Choose the stan file based on family
  stan_file <- switch(p_family,
                      "weibull" = file.path(pkg_path, "stan", "gtweibull.stan"),
                      "gamma" = file.path(pkg_path, "stan", "gtgamma.stan"),
                      stop("Unknown family! Choose from: 'weibull' or 'gamma'")
  )
  
  # Prepare the data for Stan model
  if (ncol(data) <= 1) {
    stop("Data should have at least one predictor and one response variable.")
  }
  
  # Handle the formula & prepare the data
  surv_obj <- eval(formula[[2]], envir = data)
  y <- surv_obj[, "time"]
  status <- surv_obj[, "status"]
  
  # Create design matrix (remove intercept if present)
  X <- model.matrix(formula, data)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, -1, drop = FALSE]  # Remove intercept
  }
  
  # Prepare hyperparameter data
  # -- beta
  if (length(hyperparameters$beta$mean) == 2) {
    beta1_loc   <- hyperparameters$beta$mean[1]
    beta2_loc   <- hyperparameters$beta$mean[2]
  } else {
    beta1_loc   <- beta2_loc <- hyperparameters$beta$mean
  }
  if (length(hyperparameters$beta$sd) == 2) {
    beta1_scale <- hyperparameters$beta$sd[1]
    beta2_scale <- hyperparameters$beta$sd[2]
  } else {
    beta1_scale <- beta2_scale <- hyperparameters$beta$sd
  }
  # -- theta
  theta_alpha <- hyperparameters$theta$alpha
  theta_beta  <- hyperparameters$theta$beta
  # -- phi
  if (is.list(hyperparameters$phi$rate)) {
    phi1_rate <- hyperparameters$phi$rate[1]
    phi2_rate <- hyperparameters$phi$rate[2]
  } else {
    phi1_rate <- phi2_rate <- hyperparameters$phi$rate
  }
  # -- shape
  if (length(hyperparameters$shape$alpha) == 2) {
    shape1_alpha   <- hyperparameters$shape$alpha[1]
    shape2_alpha   <- hyperparameters$shape$alpha[2]
  } else {
    shape1_alpha   <- shape2_alpha <- hyperparameters$shape$alpha
  }
  if (length(hyperparameters$shape$beta) == 2) {
    shape1_beta   <- hyperparameters$shape$beta[1]
    shape2_beta   <- hyperparameters$shape$beta[2]
  } else {
    shape1_beta <- shape2_beta <- hyperparameters$shape$beta
  }
  # -- scale
  if (length(hyperparameters$scale$alpha) == 2) {
    scale1_alpha   <- hyperparameters$scale$alpha[1]
    scale2_alpha   <- hyperparameters$scale$alpha[2]
  } else {
    scale1_alpha   <- scale2_alpha <- hyperparameters$scale$alpha
  }
  if (length(hyperparameters$scale$beta) == 2) {
    scale1_beta   <- hyperparameters$scale$beta[1]
    scale2_beta   <- hyperparameters$scale$beta[2]
  } else {
    scale1_beta   <- scale2_beta <- hyperparameters$scale$beta
  }

  # Prepare core Stan data
  stan_data <- list(
    N = nrow(data),
    K = ncol(X),
    X = X,
    y = y,
    status = status,
    use_truncation = ifelse(truncation == 1, 1, 0),
    truncation_lower = if(truncation == 1) 0 else 0.1,
    truncation_upper = if(truncation == 1) max(y) else max(y)*2,
    
    # prior hyperparameters
    beta1_loc    = beta1_loc,
    beta2_loc    = beta2_loc,
    beta1_scale  = beta1_scale,
    beta2_scale  = beta2_scale,
    
    theta_alpha  = theta_alpha,
    theta_beta   = theta_beta
  )
  
  # Append to stan_data the p_family-specific prior hyperparameters
  if (p_family == "weibull") {
    stan_data <- c(stan_data, list(
      shape1_alpha = shape1_alpha,
      shape2_alpha = shape2_alpha,
      shape1_beta = shape1_beta,
      shape2_beta = shape2_beta,
      scale1_alpha = scale1_alpha,
      scale2_alpha = scale2_alpha,
      scale1_beta = scale1_beta,
      scale2_beta = scale2_beta
    ))
  } else if (p_family == "gamma") {
    stan_data <- c(stan_data, list(
      phi1_rate = phi1_rate,
      phi2_rate = phi2_rate
    ))
  } else {
    stop("Unsupported p_family: ", p_family)
  }
  
  # Load the Stan model
  stan_model <- rstan::stan_model(file = stan_file)
  
  # Fit the model using sampling
  fit <- sampling(stan_model,
                  data = stan_data, iter = iterations, warmup = burning_iterations,
                  chains = chains, seed = seed
  )
  
  posterior <- extract(fit)
  print(names(posterior)) # Check the parameters available
  
  # Process results based on result_type
  result <- get_survival_results(p_family, fit, result_type, status_vector = data$status)
  return(result)
}