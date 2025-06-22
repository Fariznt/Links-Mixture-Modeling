#' Fits a Bayesian mixture model using Stan with various distribution families.
#'
#' @param formula Model formula (e.g., y ~ x1 + x2)
#' @param p_family Distribution family ("linear", "logistic", "poisson", "gamma")
#' @param data Input data frame or "random" for synthetic data
#' @param components Vector specifying mixture components (e.g., c("linear", "linear"))
#' @param prior Named list (or NULL).  Expected keys:
#'   `mu = list(mean, sd)`, `beta = list(mean, sd)`,
#'   `sigma = list(scale)`, `theta = list(alpha, beta)`,
#'   `phi = list(rate)` (gamma only).
#'   mean, sd, scale, rate are length-2 vectors, where each value corresponds to
#'   one of the two components in the mixture. If a scalar is passed, it
#'   is used for both components. NULL or missing
#'   values triggers weakly-informative defaults.
#' @param result_type 0 for matrix output, 1 for posterior samples,
#' @param iterations Total number of MCMC iterations, default is
#' @param burning_iterations Number of burn-in iterations
#' @param chains Number of MCMC chains
#' @param seed Random seed (integer or "random")
#' @return Results depending on result_type
#' @export
#' @importFrom rstan stan_model sampling
#' @importFrom stats model.frame model.response model.matrix
#' @importFrom utils read.csv
#' @examples
#' \dontrun{
#' # Linear mixture example
#' fit <- fit_model(y ~ x1 + x2, "linear", data = df, 
#'                 components = c("linear", "linear"))
#' }

fit_model <- 
  function(formula, 
           p_family, 
           data, 
           components, 
           hyperparameters = NULL,
           result_type = 0, 
           iterations = 10000, 
           burning_iterations = 1000, 
           chains = 2, 
           seed = 123) {

  # Define default hyperparameter values
  defaults <- list(
    mu    = list(mean  = c(0, 0), sd = c(5, 5)),
    beta  = list(mean  = c(0, 0), sd = c(5, 5)),
    sigma = list(scale = c(2.5, 2.5)),
    theta = list(alpha = 1, beta  = 1),
    phi   = list(rate  = c(1, 1))
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
    data <- generate_synthetic_mixture_data(p_family, seed)
    print(head(data))  # first few rows of generated data
    print(dim(data))   # dimensions of generated data
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

  # Check for missing data
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("The input data is not a valid data frame or is empty.")
  }

  # Generate stan file
  stan_file <- generate_stan(components, formula, data)

  # Prepare the data~
  if (ncol(data) <= 1) {
    stop("Data should have at least one predictor and one response variable.")
  }
  
  #Handle the formula & prepare the data
  model_frame <- model.frame(formula, data)
  y <- model.response(model_frame)
  X <- model.matrix(formula, model_frame)[, -1]

  # Prepare hyperparameter data
  # -- mu
  if (length(hyperparameters$mu$mean) == 2) {
    mu1_loc <- hyperparameters$mu$mean[1]
    mu2_loc <- hyperparameters$mu$mean[2]
  } else {
    mu1_loc <- mu2_loc <- hyperparameters$mu$mean
  }
  if (length(hyperparameters$mu$sd) == 2) {
    mu1_scale <- hyperparameters$mu$sd[1]
    mu2_scale <- hyperparameters$mu$sd[2]
  } else {
    mu1_scale <- mu2_scale <- hyperparameters$mu$sd
  }
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
  # -- sigma (scale-only)
  if (length(hyperparameters$sigma$scale) == 2) {
    sigma1_scale <- hyperparameters$sigma$scale[1]
    sigma2_scale <- hyperparameters$sigma$scale[2]
  } else {
    sigma1_scale <- sigma2_scale <- hyperparameters$sigma$scale
  }
  
  # -- theta
  theta_alpha <- hyperparameters$theta$alpha
  theta_beta  <- hyperparameters$theta$beta
  
  # -- phi (Gamma rate parameters)
  if (is.list(hyperparameters$phi$rate)) {
    phi1_rate <- hyperparameters$phi$rate[1]
    phi2_rate <- hyperparameters$phi$rate[2]
  } else {
    phi1_rate <- phi2_rate <- hyperparameters$phi$rate
  }
  
  stan_data <- list(
    N = nrow(data),
    K = ncol(data) - 1,
    X = X,
    y = y,
    
    mu1_loc = mu1_loc,
    mu2_loc = mu2_loc,
    mu1_scale = mu1_scale,
    mu2_scale = mu2_scale,
    
    beta1_loc    = beta1_loc,
    beta2_loc    = beta2_loc,
    beta1_scale  = beta1_scale,
    beta2_scale  = beta2_scale,
    
    sigma1_scale = sigma1_scale,
    sigma2_scale = sigma2_scale,
    
    theta_alpha  = theta_alpha,
    theta_beta   = theta_beta,
    
    phi1_rate    = phi1_rate,
    phi2_rate    = phi2_rate
  )

  # Load the Stan model
  stan_model <- rstan::stan_model(file = stan_file)

  # Fit the model using sampling
  fit <- sampling(stan_model, data = stan_data, iter = iterations, warmup = burning_iterations, 
                  chains = chains, seed = seed)

  # Process results based on result_type
  result <- get_mixture_results(p_family, fit, result_type)
  return(result)
}


##########################################################
# Linear mixture fit test with random data

# formula <- y ~ X1 + X2
# fit_result <- fit_model(formula = formula,
                        # p_family = "linear",
                        # data = "random",
                        # components = c("linear", "linear"),
                        # result_type = 1,
                        # iterations = 3000,
                        # burning_iterations = 1000,
                        # chains = 2,
                        # seed = 123)