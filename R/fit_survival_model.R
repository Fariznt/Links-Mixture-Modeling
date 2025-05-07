#' Fits a survival model using Stan with various distribution families.
#'
#' @param formula Model formula (e.g., y ~ x1 + x2)
#' @param p_family Distribution family ("weibull", "gamma")
#' @param data Input data frame or "random" for synthetic data
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

fit_survival_model <- function(formula, p_family, data, result_type, iterations, burning_iterations, chains, seed, truncation = 0, status_column = NULL) {
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
  
  # Choose the stan file based on family
  stan_file <- switch(p_family,
                      "weibull" = system.file("stan", "gtweibull.stan", package = "LinksMixtureModeling"),
                      "gamma" = system.file("stan", "gtgamma.stan", package = "LinksMixtureModeling"),
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
  
  # Prepare Stan data
  stan_data <- list(
    N = nrow(data),
    K = ncol(X),
    X = X,
    y = y,
    status = status,
    use_truncation = ifelse(truncation == 1, 1, 0),
    truncation_lower = if(truncation == 1) 0 else 0.1,
    truncation_upper = if(truncation == 1) max(y) else max(y)*2
  )
  
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