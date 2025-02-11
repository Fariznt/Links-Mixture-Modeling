data {
  int<lower=1> N;             // Number of observations
  array[N] int<lower=0> y;         // Poisson response variable (counts)
  int<lower=1> K;             // Number of predictors
  matrix[N, K] X;             // Predictor matrix
}

parameters {
  real<lower=0, upper=1> theta;           // Mixing proportions (constrained to sum to 1)
  vector[K] beta1;            // Regression coefficients for component 1
  vector[K] beta2;            // Regression coefficients for component 2
}

model {
  vector[N] log_lik1;         // Log-likelihood for component 1
  vector[N] log_lik2;         // Log-likelihood for component 2
  
  // Linear predictors for each component
  vector[N] eta1 = X * beta1; // Linear predictor for component 1
  vector[N] eta2 = X * beta2; // Linear predictor for component 2

  // Calculate log-likelihoods for each component
  for (n in 1:N) {
    log_lik1[n] = poisson_log_lpmf(y[n] | eta1[n]); // Component 1
    log_lik2[n] = poisson_log_lpmf(y[n] | eta2[n]); // Component 2
    target += log_mix(theta, log_lik1[n], log_lik2[n]);
  }

  // Priors for regression coefficients
  beta1 ~ normal(0, 5);
  beta2 ~ normal(0, 5);
  theta ~ beta(1,1);
}

generated quantities {
  array[N] int<lower=1, upper=2> z;    // Mixture membership
  for (n in 1:N) {
    // Calculate unnormalized log probabilities for each component
    real log_prob1 = log(theta) + poisson_log_lpmf(y[n] | X[n] * beta1);
    real log_prob2 = log1m(theta) + poisson_log_lpmf(y[n] | X[n] * beta2);
    
    // Normalize probabilities using softmax
    vector[2] prob = softmax([log_prob1, log_prob2]');
    
    // Sample z[n] based on the posterior probabilities
    z[n] = categorical_rng(prob);
  }
}

