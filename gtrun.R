## Generates synthetic data for linear, logistic, poisson and gamma models, then fits stan model for each and gives results. 
## Note that results directly from stan model are not usable as they may have labels switched. instead we post-process them first and compute posterior means, variances and confidence limits manually


rm(list=ls())
# For reproducibility
set.seed(123) 


# LINEAR
# Primitives
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='~/Desktop/utra/Links-Mixture-Modeling/gtlinear.stan')

# Parameters
N <- 200                # Number of observations
K <- 2                  # Number of predictors including intercept
theta <- 0.8           # Proportion of observations that are true links
beta1 <- c(1, 2)       # Regression coefficients for true links
beta2 <- c(3, 4)       # Regression coefficients for false links
sigma1 <- 1             # Standard deviation for the true links
sigma2 <- 2             # Standard deviation for the false links

# Simulate predictor matrix
X <- cbind(1, runif(N, -2, 2))  # Intercept and a single predictor

# Simulate mixture memberships
z <- rbinom(N, size = 1, prob = theta) + 1 #true link (z=1) or false link (z=2)
z <- ifelse(z == 1, 0, 1)

# Simulate outcomes
y <- numeric(N)
for (i in 1:N) {
  if (z[i] == 1) {
    y[i] <- rnorm(1, mean = sum(X[i, ] * beta1), sd = sigma1)
  } else {
    y[i] <- rnorm(1, mean = sum(X[i, ] * beta2), sd = sigma2)
  }
}


# Prepare data for Stan
stan_data <- list(
  N = N,       # Number of observations
  K = K,                # Number of predictors
  X = X, # Predictor matrix
  y = y            # Outcome vector
)


# Fit the model
fit <- sampling(stanDso, data = stan_data, iter = 1e4,chains = 1,seed = 123)

# Print summary of results from stan fit; but cant use as may have switched labels. instead post-process first and compute results manually
#print(fit, pars = c("theta", "mu1", "mu2", "sigma1", "sigma2", "beta1", "beta2"))

# Extract posterior
posterior <- extract(fit)

# Extract membership samples
z_samples = posterior$z

# Extract beta samples
beta1.p = posterior$beta1
beta2.p = posterior$beta2

# Post process beta samples
for (i in 1:nrow(z_samples)) {
  
  # Relabel
  z_samples[i, ] <- ifelse(z_samples[i,] == 1, 0, 1)
  
  if (mean(z_samples[i,]) < 0.5) {  # Check if label 1 is less dominant
    # Swap labels in z is not
    z_samples[i,] <- 1 - z_samples[i, ]
    
    # Swap coefficients for beta
    temp <- beta1.p[i,]
    beta1.p[i, ] <- beta2.p[i, ]
    beta2.p[i, ] <- temp 
  }
}

# Compute posterior mean, variance and confidence limits for coefficent in each mixture component
mean_beta1 <- apply(beta1.p,2,mean)
mean_beta2 <- apply(beta2.p,2,mean)
var_beta1 <- apply(beta1.p,2,var)
var_beta2 <- apply(beta2.p,2,var)
ci_beta1 <- apply(beta1.p,2,function(x) quantile(x, probs = c(0.25, 0.95)))
ci_beta2 <- apply(beta2.p,2,function(x) quantile(x, probs = c(0.25, 0.95)))




# LOGISTIC
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='~/Desktop/utra/Links-Mixture-Modeling/gtlogistic.stan')


# Data
# Simulate a dataset
set.seed(123)
N <- 200                  # Number of observations
K <- 2                    # Number of predictors
theta <- 0.6              # Proportion of component 1
beta1 <- c(0.5, -1)         # Coefficients for component 1
beta2 <- c(0.8, 1)         # Coefficients for component 2

# Generate predictors
X <- matrix(rnorm(N * K), nrow = N, ncol = K)

# Generate mixture memberships
z <- rbinom(N, size = 1, prob = theta) + 1 # 1 for component 1, 2 for component 2

# Generate binary outcomes
logit <- function(x) exp(x) / (1 + exp(x))
y <- numeric(N)
for (i in 1:N) {
  if (z[i] == 1) {
    y[i] <- rbinom(1, 1, prob = logit(sum(X[i, ] * beta1)))
  } else {
    y[i] <- rbinom(1, 1, prob = logit(sum(X[i, ] * beta2)))
  }
}

# Prepare Stan data
stan_data <- list(
  N = N,
  K = K,
  X = X,
  y = y
)

# Fit the model
fit <- sampling(stanDso, data = stan_data, iter = 20000,chains = 1,seed = 123)


# Print summary of results
print(fit, pars = c("theta", "beta1", "beta2"))
z_samples <- extract(fit, "z")$z  # Mixture membership samples





## POISSON
rm(list=ls())
# primitives
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='~/Desktop/utra/Links-Mixture-Modeling/gtpoisson.stan')


# Number of observations
N <- 200
K <- 2
X <- cbind(1, runif(N, -2, 2))  # Intercept and a single predictor

# True parameters for two components
beta1 <- c(1, 2)
beta2 <- c(2, 3)

# Mixing proportion
theta <- 0.6

# Generate data from the mixture
lambda1 <- exp(X %*% beta1)
lambda2 <- exp(X %*% beta2)
z <- rbinom(N, 1, theta)  # Latent component memberships
y <- ifelse(z == 1, rpois(N, lambda1), rpois(N, lambda2))

# Save data for Stan
stan_data <- list(N = N, K = K, X = X, y = y)

# Fit the model
fit <- sampling(stanDso, data = stan_data, iter = 50000,chains = 1,seed = 123)


# Print summary of results from stan fit; but cant use as may have switched labels. instead post-process first and compute results manually
#print(fit, pars = c("theta", "beta1", "beta2"))

#### Now to post-process for possible label switching
# Extract posterior
posterior <- extract(fit)

# Extract membership samples
#z_samples <- extract(fit, "z")$z  
z_samples = posterior$z

# Extract beta samples
beta1.p = posterior$beta1
beta2.p = posterior$beta2

# Post process beta samples
for (i in 1:nrow(z_samples)) {
  
  # Relabel
  z_samples[i, ] <- ifelse(z_samples[i,] == 1, 0, 1)
  
  if (mean(z_samples[i,]) < 0.5) {  # Check if label 1 is dominant
    # Swap labels in z
    z_samples[i,] <- 1 - z_samples[i, ]
    
    # Swap coefficients for beta
    temp <- beta1.p[i,]
    beta1.p[i, ] <- beta2.p[i, ]
    beta2.p[i, ] <- temp 
  }
}

# Compute posterior mean, variance and confidence limits for coefficent in each mixture component
mean_beta1 <- apply(beta1.p,2,mean)
mean_beta2 <- apply(beta2.p,2,mean)
var_beta1 <- apply(beta1.p,2,var)
var_beta2 <- apply(beta2.p,2,var)
ci_beta1 <- apply(beta1.p,2,function(x) quantile(x, probs = c(0.25, 0.95)))
ci_beta2 <- apply(beta2.p,2,function(x) quantile(x, probs = c(0.25, 0.95)))

print(mean_beta1)
print(mean_beta2)
print(var_beta1)
print(var_beta2)
print(ci_beta1)
print(ci_beta1)

print(z_samples)




## GAMMA
rm(list=ls())
# primitives
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='~/Desktop/utra/Links-Mixture-Modeling/gtgamma.stan')

# For reproducibility
set.seed(123)

# Parameters for links
beta1 <- c(0.5, 1.2)   # Intercept and slope for component 1
phi1 <- 5              # Shape parameter (1/dispersion)

# Parameters for nonlinks
beta2 <- c(2.0, -0.8)  # Intercept and slope for component 2
phi2 <- 2              # Shape parameter (1/dispersion)

# Mixing proportion
theta <- 0.6  # proportion of links

# Generate data
N <- 200  
K <- 2
X <- cbind(1, runif(N, -2, 2))  # Intercept and a single predictor

# Generate component memberships
z <- rbinom(N, 1, theta)  # Latent component memberships

# Generate response (Y) from Gamma distribution
y <- numeric(N)
for (i in 1:N) {
  if (z[i] == 1) {
    mu <- exp(X[i, ] %*% beta1)
    y[i] <- rgamma(1, shape = phi1, rate = phi1 / mu)
  } else {
    mu <- exp(X[i, ] %*% beta2)
    y[i] <- rgamma(1, shape = phi2, rate = phi2 / mu)
  }
}

# Prepare data for Stan
stan_data <- list(
  N = N,
  K = K,                
  y = y,
  X = X
)


# Fit the model
fit <- sampling(stanDso, data = stan_data, iter = 50000,chains = 1,seed = 123)

# Print summary of results from stan fit; but cant use as may have switched labels. instead post-process first and compute results manually
#print(fit, pars = c("theta", "beta1", "beta2"))

#### Now to post-process for possible label switching
# Extract posterior
posterior <- extract(fit)

# Extract membership samples
#z_samples <- extract(fit, "z")$z  
z_samples = posterior$z

# Extract beta samples
beta1.p = posterior$beta1
beta2.p = posterior$beta2
phi1 = posterior$phi1
phi2 = posterior$phi2

# Post process beta samples
for (i in 1:nrow(z_samples)) {
  
  # Relabel
  z_samples[i, ] <- ifelse(z_samples[i,] == 1, 0, 1)
  
  #swap
  if (mean(z_samples[i,]) < 0.5) {  # Check if label 1 is dominant
    # Swap labels in z
    z_samples[i,] <- 1 - z_samples[i, ]
    
    # Swap coefficients for beta
    temp <- beta1.p[i,]
    beta1.p[i, ] <- beta2.p[i, ]
    beta2.p[i, ] <- temp 
    
    # Swap coefficients for shape (dispersion)
    temp <- phi1[i]
    phi1[i] <- phi2[i]
    phi2[i] <- temp 
  }
}

# Compute posterior mean, variance and confidence limits for coefficent in each mixture component
mean_beta1 <- apply(beta1.p,2,mean)
mean_beta2 <- apply(beta2.p,2,mean)
var_beta1 <- apply(beta1.p,2,var)
var_beta2 <- apply(beta2.p,2,var)
ci_beta1 <- apply(beta1.p,2,function(x) quantile(x, probs = c(0.25, 0.95)))
ci_beta2 <- apply(beta2.p,2,function(x) quantile(x, probs = c(0.25, 0.95)))
