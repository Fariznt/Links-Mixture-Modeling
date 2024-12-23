data {
  int<lower=1> N;           // Number of observations
  int<lower=1> K;           // Number of predictors
  matrix[N, K] X;           // Predictor matrix
  int<lower=0, upper=1> y[N]; // Binary outcome
}

parameters {
  real<lower=0, upper=1> theta;     // Mixing proportion (for component 1)
  vector[K] beta1;                 // Regression coefficients for component 1
  vector[K] beta2;                 // Regression coefficients for component 2
}

model {
  vector[N] log_lik1;  // Log-likelihood contributions from component 1
  vector[N] log_lik2;  // Log-likelihood contributions from component 2

  // Priors
  theta ~ beta(1, 1);           // Uniform prior on mixing proportion
  beta1 ~ normal(0, 5);         // Priors for regression coefficients (component 1)
  beta2 ~ normal(0, 5);         // Priors for regression coefficients (component 2)

  // Mixture model likelihood
  for (n in 1:N) {
    log_lik1[n] = bernoulli_logit_lpmf(y[n] | dot_product(X[n], beta1));
    log_lik2[n] = bernoulli_logit_lpmf(y[n] | dot_product(X[n], beta2));
    target += log_mix(theta, log_lik1[n], log_lik2[n]);
  }
  
}

generated quantities {
  int<lower=1, upper=2> z[N];      // Mixture membership for each observation
  for (n in 1:N) {
    vector[2] log_weights;
    log_weights[1] = log(theta) + bernoulli_logit_lpmf(y[n] | dot_product(X[n], beta1));
    log_weights[2] = log1m(theta) + bernoulli_logit_lpmf(y[n] | dot_product(X[n], beta2));
    z[n] = categorical_rng(softmax(log_weights)); // Sample membership
  }
}
