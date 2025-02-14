library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# main wrapper function
fit_model <- function(formula, family, data, iterations, burning_iterations, chains, seed) {
  
  # generate a random seed if one is not given
  seed_input <- as.integer(seed)
  if (is.null(seed)) {
    seed <- sample.int(1e6, 1) 
  }
  
  # Choose the post-processing function based on family
  if (family == "linear") {
    stan_file <- "gtlinear.stan"   # relative paths for now
  } else if (family == "logistic") {
    stan_file <- "gtlogistic.stan"
  } else if (family == "poisson") {
    stan_file <- "gtpoisson.stan"
  } else if (family == "gamma") {
    stan_file <- "gtgamma.stan"
  } else {
    stop("Unknown family! Choose from: 'linear', 'logistic', 'poission', or 'gamma'")
  }
  
  # Prepare the data 
  stan_data <- list(
    N = nrow(data),                      # Number of observations
    K = ncol(data) - 1,                  # Number of predictors
    X = as.matrix(data[, -ncol(data)]),  # Predictor matrix
    y = data[, ncol(data)]               # Response vector
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
seed <- readline("Enter the seed, or random for random seed (default 123): ")

# Reading in the dataset
data <- read.csv(data_path)

# Run the model
fit_result <- fit_model(formula, family, data, 
                        iterations = ifelse(iterations == "", 1000, iterations), 
                        burning_iterations = ifelse(burning_iterations == "", 1000, burning_iterations), 
                        chains = ifelse(chains == "", 2, chains), 
                        seed = ifelse())

print(fit_result$fit)
