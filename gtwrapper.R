library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# main wrapper function
fit_model <- function(formula, family, data, iterations, burning_iterations, chains, seed) {
  # Choose the stan file based on family
  if (family == "linear") {
    stan_file <- "gtlinear.stan" # relative paths for now
  } else if (family == "logistic") {
    stan_file <- "gtlogistic.stan"
  } else if (family == "poisson") {
    stan_file <- "gtpoisson.stan"
  } else if (family == "gamma") {
    stan_file <- "gtgamma.stan"
  } else {
    stop("Unknown family! Choose from: 'linear', 'logistic', 'poission', or 'gamma'")
  }

  # Check for missing or NA values in the dataset:
  # returns an array of all the missing values with the same size as the data set
  missing_data <- which(is.na(data), arr.ind = TRUE)

  if (length(missing_data) > 0) {
    stop(
      "Error: Missing values found in the following rows and columns:\n",
      paste(apply(missing_data, 1, function(x) paste("Row", x[1], "Column", x[2])), collapse = "\n")
    )
  }

  # Prepare the data
  stan_data <- list(
    N = nrow(data), # Number of observations
    K = ncol(data) - 1, # Number of predictors
    X = as.matrix(data[, -ncol(data)]), # Predictor matrix
    y = data[, ncol(data)] # Response vector
  )

  # Load the Stan model
  stan_model <- rstan::stan_model(file = stan_file)

  # Fit the model using sampling
  fit <- sampling(stan_model, data = stan_data, iter = iterations, warmup = burning_iterations, chains = chains, seed = seed)

  # TO DO: the post processing stuff, based on each family type





  return(list(fit = fit))
}

# User prompts
formula <- as.formula(readline("Enter the formula (e.g., y ~ x1 + x2): "))
family <- readline("Enter the family (linear, logistic, poisson, gamma): ")
data_path <- readline("Enter the name of your dataset (CSV format, must be in working directory): ")
iterations <- as.integer(readline("Enter the number of iterations (default 1000): "))
burning_iterations <- as.integer(readline("Enter the number of burning iterations (default 1000): "))
chains <- as.integer(readline("Enter the number of chains (default 2): "))
seed_input <- readline("Enter the seed, or 'RAN' for a random seed (default 123): ")

# Set seed value: If 'RAN' is input, generate a random seed, otherwise use input or default to 123
seed <- ifelse(seed_input == "RAN", sample.int(1e6, 1),
  ifelse(seed_input == "", 123, as.integer(seed_input))
)

# Reading in the dataset
data <- read.csv(data_path)

# Run the model
fit_result <- fit_model(formula, family, data,
  iterations = ifelse(is.null(iterations) || iterations == "", 1000, iterations),
  burning_iterations = ifelse(is.null(burning_iterations) || burning_iterations == "", 1000, burning_iterations),
  chains = ifelse(is.null(chains) || chains == "", 2, chains),
  seed = seed
)

print(fit_result$fit)