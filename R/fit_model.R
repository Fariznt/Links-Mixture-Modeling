#' Fits a Bayesian mixture model using Stan with various distribution families.
#'
#' @param formula Model formula (e.g., y ~ x1 + x2)
#' @param p_family Distribution family ("linear", "logistic", "poisson", "gamma")
#' @param data Input data frame or "random" for synthetic data
#' @param components Vector specifying mixture components (e.g., c("linear", "linear"))
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

fit_model <- function(formula, p_family, data, components, result_type = 0, iterations = 10000, burning_iterations = 1000, chains = 2, seed = 123) {
  
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
    data <- generate_synthetic_data(p_family, seed)
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

  # Prepare the data
  if (ncol(data) <= 1) {
    stop("Data should have at least one predictor and one response variable.")
  }
  
  #Handle the formula & prepare the data
  model_frame <- model.frame(formula, data)
  y <- model.response(model_frame)
  X <- model.matrix(formula, model_frame)[, -1]
  
  stan_data <- list(
    N = nrow(data),
    K = ncol(data) - 1,
    X = X,
    y = y
  )

  # Load the Stan model
  stan_model <- rstan::stan_model(file = stan_file)

  # Fit the model using sampling
  fit <- sampling(stan_model, data = stan_data, iter = iterations, warmup = burning_iterations, 
                  chains = chains, seed = seed)

  # Process results based on result_type
  result <- get_results(p_family, fit, result_type)
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