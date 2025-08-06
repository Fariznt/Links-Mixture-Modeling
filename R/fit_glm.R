#' Fits a Bayesian mixture model using Stan with various distribution families.
#'
#' @param formula Model formula (e.g., y ~ x1 + x2)
#' @param p_family Distribution family ("linear", "logistic", "poisson", "gamma")
#' @param data Input data frame or "random" for synthetic data
#' @param priors Named list (or NULL) of strings defining priors and optionally
#' hyperparameter definitions used in the prior definition. NULL priors or missing
#' list elements defining needed priors triggers weakly-informative defaults.
#' @param iterations Total number of MCMC iterations, default is
#' @param warmup_iterations Number of burn-in iterations
#' @param chains Number of MCMC chains
#' @param seed Random seed (integer or "random")
#' @param diagnostics If TRUE, returns results in a list that also includes performance
#' metrics and (when data == "random") the latent values used for data generation. 
#' FALSE by default.
#' @return Results depending on result_type
#' @export
#' @importFrom rstan stan_model sampling
#' @importFrom stats model.frame model.response model.matrix
#' @importFrom utils read.csv
#' @examples
#' \dontrun{
#' # Linear mixture example, assuming cov_matrix is a defined 2x2 R matrix
#' lin_fit <- LinksMixtureModeling::fit_glm(
#'  formula             = y ~ X1 + X2,        # y, X1, X2 auto-generated
#'  p_family            = "linear",
#'  data                = "random",           # triggers built-in generator
#'  priors              = list(beta1_sigma = cov_matrix, 
#'                             beta1_loc = c(7,8), 
#'                             beta1 = "multi_normal(beta1_loc, beta1_sigma)",
#'                             mu1 = "normal(0,4)",
#'                             mu2_loc = 0,
#'                             mu2_scale = 6,
#'                             mu2 = "normal(mu2_loc, mu2_scale)"),
#'  iterations          = 2000,               # 1000 warm-up + 1000 sampling
#'  warmup_iterations   = 1000,
#'  chains              = 2,
#'  seed                = 123
#')
#' }
fit_glm <- 
  function(formulas,
           p_family, 
           data, 
           priors = NULL,
           iterations = 10000, 
           warmup_iterations = 1000, 
           chains = 2, 
           seed = 123,
           diagnostics = FALSE) {

  validate_args(priors, p_family)
    
  # Set default priors if some or all not passed in
  priors <- fill_defaults(priors, p_family, 'glm')
  
  # Correct/prepare user input when formulas is just 1 formula
  if (!is.list(formulas)) { 
    formulas <- list(formulas)
  }

  # Parse seed and generate random if "random" passed in
  if (seed == "random") {
    seed <- sample.int(1e6, 1)
  } else {
    seed <- as.integer(seed)
  }

  # Check what type of data should be loaded & handle NA values
  if (identical(data, "random")) {
    print("Generating synthetic data...")
    generated <- generate_synthetic_glm_data(p_family, formulas, seed)
    data <- generated$data
  } else if (is.character(data) && length(data) == 1) {
    # data is a fie path
    if (!file.exists(data)) {
      stop("Error: The specified CSV file does not exist.")
    }
    data <- read.csv(data)
    if (anyNA(data)) {
      na_locations <- which(is.na(data), arr.ind = TRUE)  # Get row and column locations of NA values
      stop("Error: NA values found in the data. Locations of NA values:\n", 
           paste("Row:", na_locations[, 1], "Column:", na_locations[, 2], collapse = "\n"))
    }
  } else {
    if (!is.data.frame(data) || nrow(data) == 0) {
      stop("`data` must be either:\n",
           "\"random\"\n",
           "A single string path to a CSV file\n",
           "A nonempty data.frame")
    }
  }

  # Generate stan file
  stan_file <- generate_stan(c(p_family, p_family), formula, data, priors)

  # Prepare the data~
  if (ncol(data) <= 1) {
    stop("Data should have at least one predictor and one response variable.")
  }
  
  # Prepare data/formulas to pass into stan model
  
  # Initial declaration of stan data variables
  N <- nrow(data)
  M <- length(formulas)
  
  K_per_m <- integer(0) 
  y_all <- matrix(nrow = N, ncol = 0) 
  X_all <- matrix(nrow = N, ncol = 0)
  K_sum <- 0

  # Continue constructing values of M, K_per_m, X_all, and y_all
  for (i in 1:length(formulas)) {
    model_frame <- model.frame(formulas[[i]], data)
    model_matrix <- model.matrix(formulas[[i]], model_frame)
    
    # append predictor count for this formula onto vector of all predictor counters
    predictor_count <- ncol(model_matrix) - 1 # -1 to disclude the intercept column
    K_per_m <- c(K_per_m, predictor_count) 
    
    # Add to sum of total predictors
    K_sum <- K_sum + predictor_count
    
    # append response vector for this formula (column-wise) to matrix of response vectors
    y <- model.response(model_frame) # response vector for this formula
    y_all <- cbind(y_all, y)
    
    # append predictor matrix for this formula (column-wise) to matrix of all predictors
    X <- model_matrix[, -1, drop = FALSE] # form matrix for this formula by dropping intercept column
    X_all <- cbind(X_all, X) # append to merged matrix
  }

  stan_data <- list(
    N = N, # num of observations
    M = M, # num of response vars
    y_all = y_all, # matrix of response vectors
    K_per_m = array(K_per_m), # vector of predictors per response var
    X_all = X_all, # col-wise merged matrix of all predictor matrices
    K_sum = K_sum # total num of predictors
  )

  # Load the Stan model
  compile_time <- system.time({
    stan_model <- rstan::stan_model(file = stan_file)
  })[["elapsed"]]
  

  # Fit the model using sampling
  sampling_time <- system.time({
    fit <- sampling(stan_model, data = stan_data, iter = iterations, warmup = warmup_iterations, 
                    chains = chains, seed = seed, verbose = TRUE)
  })[["elapsed"]]

  
  posterior <- extract(fit)
  
  # Stop if sampler failed
  if (length(posterior$z) == 0) {
    stop("Model compiled but sampling failed to produce any draws. \n",
         "Internal cause: All chains resulted in zero saved iterations,",
         " either due to invalid input (e.g. when iterations ≤ burn-in iterations)",
         " or due to crashing/termination during sampling.")
  }

  if (diagnostics) {
    output <- list(
      results = process_results(p_family, posterior, formulas),
      compile_time = compile_time,
      sampling_time = sampling_time
    )
    if (exists("generated", inherits = FALSE)) {
      output$latent_values <- generated$latent_values
      output$generated_data <- generated$data
    }
  } else {
    output <- process_results(p_family, posterior, formulas)
  }
  
  output
}


##########################################################
# Linear mixture fit test with random data

# formula <- y ~ X1 + X2
# fit_result <- fit_model(formula = formula,
                        # p_family = "linear",
                        # data = "random",
                        # result_type = 1,
                        # iterations = 3000,
                        # warmup_iterations = 1000,
                        # chains = 2,
                        # seed = 123)