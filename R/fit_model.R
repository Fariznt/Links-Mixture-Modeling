#' Fits a Bayesian mixture model using Stan with various distribution families.
#'
#' @param formula Model formula (e.g., y ~ x1 + x2)
#' @param p_family Distribution family ("linear", "logistic", "poisson", "gamma")
#' @param data Input data frame or "random" for synthetic data
#' @param components Vector specifying mixture components (e.g., c("linear", "linear"))
#' @param priors Named list (or NULL) of strings defining priors and optionally
#' hyperparameter definitions used in the prior definition. NULL priors or missing
#' list elements defining needed priors triggers weakly-informative defaults.
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
#' # Linear mixture example, assuming cov_matrix is a defined R matrix
#' lin_fit <- LinksMixtureModeling::fit_model(
#'  formula             = y ~ X1 + X2,        # y, X1, X2 auto-generated
#'  p_family            = "linear",
#'  data                = "random",           # triggers built-in generator
#'  components          = c("linear", "linear"),
#'  priors              = list(beta1_sigma = cov_matrix, 
#'                             beta1_loc = c(7,8), 
#'                             beta1 = "multi_normal(beta1_loc, beta1_sigma)",
#'                             mu1 = "normal(0,4)",
#'                             mu2_loc = 0,
#'                             mu2_scale = 6,
#'                             mu2 = "normal(mu2_loc, mu2_scale)"),
#'  result_type         = 1,                  # post-processed summary
#'  iterations          = 2000,               # 1000 warm-up + 1000 sampling
#'  burning_iterations  = 1000,
#'  chains              = 2,
#'  seed                = 123
#')
#' }
fit_model <- 
  function(formulas, 
           p_family, 
           data, 
           components, 
           priors = NULL,
           result_type = 0, 
           iterations = 10000, 
           burning_iterations = 1000, 
           chains = 2, 
           seed = 123) {
  
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
  stan_file <- generate_stan(components, formula, data, priors)

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
    X <- as.matrix(model_matrix[, -1, drop = FALSE]) # form matrix for this formula by dropping intercept column
    X_all <- cbind(X_all, X) # append to merged matrix
  }

  stan_data <- list(
    N = N, # num of observations
    M = M, # num of response vars
    y_all = y_all, # matrix of response vectors
    K_per_m = K_per_m, # vector of predictors per response var
    X_all = X_all, # col-wise merged matrix of all predictor matrices
    K_sum = K_sum # total num of predictors
    # TODO: probably possible to remove K_sum, and calculate in stan
  )

  # Load the Stan model
  stan_model <- rstan::stan_model(file = stan_file)

  # Fit the model using sampling
  fit <- sampling(stan_model, data = stan_data, iter = iterations, warmup = burning_iterations, 
                  chains = chains, seed = seed, verbose = TRUE)

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