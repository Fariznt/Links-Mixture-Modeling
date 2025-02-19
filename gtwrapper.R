library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# main wrapper function
fit_model <- function(formula, family, data, iterations = 1000, burning_iterations = 500, chains = 1, seed = NULL) {
  
  # generate a random seed if one is not given
  if (is.null(seed)) {
    seed <- sample.int(1e6, 1) 
  }
  
  # Choose the stan file based on family
  if (family == "linear") {
    stan_file <- "gtlinear.stan"         # relative paths for now
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
