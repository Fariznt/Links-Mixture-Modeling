library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# source of helper functions
source("gtpreprocessing.R")
source("gtpostprocessing.R")

# main wrapper function
fit_model <- function(formula, p_family, data, result_type, iterations, burning_iterations, chains, seed) {
  
  # Parse seed and generate random if "random" passed in
  if (seed == "random") {
    seed <- sample.int(1e6, 1)
  } else {
    seed <- as.integer(seed)
  }
  
  print("RESULT TYPE IS")
  print(result_type)
  
  print("FAMILy IS")
  print(p_family)
  
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

  # Choose the stan file based on family
  stan_file <- switch(p_family,
                      "linear" = "gtlinear.stan",
                      "logistic" = "gtlogistic.stan",
                      "poisson" = "gtpoisson.stan",
                      "gamma" = "gtgamma.stan",
                      stop("Unknown family! Choose from: 'linear', 'logistic', 'poisson', or 'gamma'"))

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
}


##########################################################
# Linear Fit Test with Random data (temporary)

formula <- y ~ X1 + X2
fit_result <- fit_model(formula = formula,
                        p_family = "linear",
                        data = "random",
                        result_type = 1,
                        iterations = 10000,
                        burning_iterations = 1000,
                        chains = 2,
                        seed = 123)




load_data <- function(data_input) {
  if (is.data.frame(data_input)) { # Deals with data already in data.frame format
    return (data_input)
    
  } else if (data_input == "random") {
    return (data_input)
    
  } else if (file.exists(data_input)) { # Deals with csv or the rds data objects Prof Gutman mentioned
    
    if (grepl("\\.csv$", data_input)) { # csv
      data <- read.csv(data_input)
      
    } else if (grepl("\\.rds$", data_input)) { #rds
      data <- readRDS(data_input)
      if (!is.data.frame(data)) { # checks that it is a data.frame object stored inside
        stop("Error: .rds file does not contain a data.frame")
      }
    } else {
      stop("Error: Unsupported file format. Use .csv or .rds")
    }
    return(data)
  } else {
    stop("Error: Data input must be a data.frame or a valid file path")
  }
}


##########################################################
# command line prompting
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript gtwrapper.R formula family data result_type (optional: iterations, burning_iterations, chains, seed)")
}

# assign default values
default_iterations <- 1000
default_burning_iterations <- 1000
default_chains <- 2
default_seed <- 123
default_result_type <- 0  # Default to 0 (matrix)

cat("args[4]:", args[4], "\n")

# assign argument inputs 
formula <- as.formula(args[1])
p_family <- args[2]
data <- args[3]
result_type <- ifelse(length(args) < 4, default_result_type, as.integer(args[4]))
iterations <- ifelse(length(args) < 5, default_iterations, as.integer(args[5]))
burning_iterations <- ifelse(length(args) < 6, default_burning_iterations, as.integer(args[6]))
chains <- ifelse(length(args) < 7, default_chains, as.integer(args[7]))
seed <- ifelse(length(args) < 8, default_seed, ifelse(args[8] == "random", sample.int(1e6, 1), as.integer(args[8])))

# Run the model
fit_result <- fit_model(formula, p_family, data, result_type, iterations, burning_iterations, chains, seed)

print(fit_result$fit)
print(fit_result$result_type)